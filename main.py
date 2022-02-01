# coding=utf-8
# !python3

"""This file ... is the main file :) """
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os
import glob
import time
import query_creation
import slurmarry_creation
import output_processing
import assembly_healing


def create_queries(project_dir, cov_file_path, assembly_file, low_cov_start, low_cov_end, min_query_len,
                   queries_per_file):
    """
    Function for the query creation part, including low coverage region detection, merging of low cov. regions that are
    close together and the query file creation.
    Creates a new storage_dir if necessary, creates a input and output dir (overwrite if already present)
    :param project_dir: Should include / at the end
    :param cov_file_path: path to the coverage file
    :param assembly_file: path to the assembly file
    :param low_cov_start: threshold for introducing a low cov region
    :param low_cov_end: threshold for ending a low cov regions
    :param min_query_len: minimum query length in bp
    :param queries_per_file: amount of queries placed per diamond input file
    :return: None
    """

    # Stores the relevant data of the current run
    run_information = ["CLCR create_queries run \t\t" + time.ctime(time.time())]      # Initialising

    start_time = time.time()

    # Detect the low coverage regions
    low_cov_regions = query_creation.detect_regions(cov_file_path, low_cov_start, low_cov_end)

    # Count the amount of low cov. regions before merging
    region_count = 0
    for scaffold in low_cov_regions:
        region_count += len(scaffold[1])
    run_information.append("Detected low cov. regions:\t\t" + str(region_count))

    storage_files_dir_path = project_dir + "storage_files/"
    # Create a storage dir, if it not exists
    os_command = " if ! [ -d " + storage_files_dir_path + " ] ; then mkdir " + storage_files_dir_path + "; fi"
    os.system(os_command)

    # Create the "original" low cov regions file
    low_cov_storage_tsv_path = storage_files_dir_path + "original_low_cov_regions.tsv"
    query_creation.create_low_cov_tsv_file(low_cov_regions, low_cov_storage_tsv_path)

    # Merge low cov regions, that are closer together than the min_query_length/2
    low_cov_regions = query_creation.merge_close_regions(low_cov_regions, int(min_query_len/2))

    # Count the regions after merging
    region_count = 0
    for scaffold in low_cov_regions:
        region_count += len(scaffold[1])
    run_information.append("Low cov. regions after merging:\t" + str(region_count))

    # Delete the query_files directory with all files in it, if it already exists
    query_files_dir_path = project_dir + "query_files/"
    os_command = " if [ -d " + query_files_dir_path + " ] ; then rm -r " + query_files_dir_path + "; fi"
    os.system(os_command)

    # Create the new query_files directory
    os_command = "mkdir " + query_files_dir_path
    os.system(os_command)

    # Delete the diamond_output directory with all files in it, if it already exists
    diamond_output_dir_path = project_dir + "diamond_output/"
    os_command = " if [ -d " + diamond_output_dir_path + " ] ; then rm -r " + diamond_output_dir_path + "; fi"
    os.system(os_command)

    # Create the new diamond_output directory
    os_command = "mkdir " + diamond_output_dir_path
    os.system(os_command)

    # Create the query files
    fasta_count, created_queries = query_creation.query_files_creation(low_cov_regions, assembly_file, min_query_len,
                                                                       query_files_dir_path, queries_per_file)

    run_information.append("Created query .fasta files:\t\t" + str(fasta_count))
    run_information.append("Created queries:\t\t\t" + str(created_queries))

    # Calculate the runtime
    run_time = time.strftime("%Hh%Mm%Ss", time.gmtime((time.time() - start_time)))
    # Append the runtime info at the second position
    run_information = [run_information[0], str("Runtime:\t\t\t\t" + run_time)] + run_information[1:]

    # Create run information file
    run_info_file = open((storage_files_dir_path + "query_creation_" + time.strftime("%Y%m%d-%H%M%S") + ".txt"), "w")
    for line in run_information:
        run_info_file.write((line + "\n"))

    run_info_file.close()

    return None


def prepare_slurm_run(project_dir, protein_database, auto_run):
    """
    Creates the slurmarray file for the current run, if auto_run ist set True, the job will be automatically submitted
    to the computer cluster
    :param project_dir: path to the project directory
    :param protein_database: path to the protein database
    :param auto_run: automatic job submitting if True
    :return: None
    """

    # Delete the slurm directory with all files in it, if it already exists
    slurm_dir_path = project_dir + "slurm_files/"
    os_command = " if [ -d " + slurm_dir_path + " ] ; then rm -r " + slurm_dir_path + "; fi"
    os.system(os_command)

    # Create the new slurm directory
    os_command = "mkdir " + slurm_dir_path
    os.system(os_command)

    # Relevant dirs for diamond
    input_dir = project_dir + "query_files/"
    output_dir = project_dir + "diamond_output/"

    # Create slurmarray file
    slurmarry_creation.create_slurmarry(protein_database, input_dir, output_dir, slurm_dir_path)

    # Start slurm job array
    if auto_run:
        # Count the input files
        input_file_list = glob.glob(input_dir + "/temp_in_*.fasta")
        file_count = len(input_file_list)

        # Submit job
        os_command = "sbatch --array=1-" + str(file_count) + " " + slurm_dir_path + "CLCR_slurmarray.slurm"
        os.system(os_command)

    return None


def create_healed_assembly(project_dir, unhealed_assembly, dynamic_threshold_dist):
    """
    This function creates a healed version of the handed over assembly. Log file is saved in the storage_files directory
    :param project_dir: path of the project dir
    :param unhealed_assembly: path to the original assembly
    :param dynamic_threshold_dist: The max_detect_distance defines the distance from
                            a detected frameshift position to the original low cov. region, where a frameshift is still
                            considered and not excluded in the further analysis. Values around 10 might be reasonable.
    :return: None
    """

    # Stores the relevant data of the current run
    run_information = ["CLCR create_healed_assembly run \t\t" + time.ctime(time.time())]   # Initialising

    start_time = time.time()

    # Read in the diamond results
    output_dir = project_dir + "diamond_output/"
    all_diamond_results = output_processing.read_in_diamond_output(output_dir)

    # Read in the original low cov. regions
    low_cov_storage_tsv = project_dir + "storage_files/original_low_cov_regions.tsv"
    stored_low_cov_regions = query_creation.read_in_low_cov_tsv_file(low_cov_storage_tsv)

    # Filter out relevant frameshift positions
    temp_list_1 = output_processing.filter_out_relevant_results(all_diamond_results, dynamic_threshold_dist,
                                                                stored_low_cov_regions)
    considered_diamond_hits_list, healing_region_list, considered_frameshifts_count, frameshifts_detected, \
    insertion_count, deletion_count = temp_list_1

    # Create healing data file
    healing_data_path = project_dir + "storage_files/healing_data.tsv"
    assembly_healing.create_detailed_healing_information_file(considered_diamond_hits_list, healing_data_path)

    # Create healed assembly file
    new_assembly_dir = project_dir + "healed_assembly/"
    temp_list_2 = assembly_healing.heal_assembly_file(healing_region_list, unhealed_assembly, new_assembly_dir)
    new_fna_file_path, simplified_distance_distribution = temp_list_2

    # Calculate the runtime
    run_time = time.strftime("%Hh%Mm%Ss", time.gmtime((time.time() - start_time)))
    # Append relevant information
    run_information.append("Runtime:\t\t\t\t" + run_time)
    run_information.append("Queries with at least one frameshift: \t" + str(len(healing_region_list)))
    run_information.append("Queries with at least one Diamond hit: \t" + str(len(all_diamond_results)))
    run_information.append("Considered frameshifts in the healing process: \t" + str(considered_frameshifts_count))
    # Regardless if they are in the original low cov or not, putative intron transition frameshift excluded
    # run_information.append("Considered frameshifts by the overlapping heuristic: \t" + str(frameshifts_detected))
    run_information.append("Healed insertions: \t" + str(insertion_count))
    run_information.append("Healed deletion: \t" + str(deletion_count))

    # Create run information file
    storage_files_dir_path = project_dir + "storage_files/"
    run_info_file = open((storage_files_dir_path + "create_healed_assembly" + time.strftime("%Y%m%d-%H%M%S") + ".txt"),
                         "w")
    for line in run_information:
        run_info_file.write((line + "\n"))

    run_info_file.close()

    return None


def main():
    print("CLCR MAIN CALLED")

    # New program structure version

    """# Query creation
    project_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_01.02.2022/"
    cov_file_path = "/share/gluster/NOTLOESUNG/freya/T_draco/t_draco.pbc.wgs_short.txt"
    assembly_file = "/share/gluster/assemblies/Tdraco/" \
                    "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"
    low_cov_start = 15
    low_cov_end = 18
    min_query_len = 500
    queries_per_file = 5000

    create_queries(project_dir, cov_file_path, assembly_file, low_cov_start, low_cov_end, min_query_len,
                   queries_per_file)"""

    """# Create and submit slurm array job
    project_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_01.02.2022/"
    protein_database = "/home/johannes/Desktop/trachinus_draco/protein_db/protein_db.dmnd"
    auto_run = False
    prepare_slurm_run(project_dir, protein_database, auto_run)"""

    # Create healed assembly file
    project_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_01.02.2022/"
    unhealed_assembly = "/share/gluster/assemblies/Tdraco/" \
                   "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"
    create_healed_assembly(project_dir, unhealed_assembly)


    # Read in the diamond results

    output_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/output_files/"

    all_diamond_results = output_processing.read_in_diamond_output(output_dir)

    low_cov_storage_tsv = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                          "storage_files/original_low_cov_regions.tsv"
    stored_low_cov_regions = query_creation.read_in_low_cov_tsv_file(low_cov_storage_tsv)

    consid_diamond_hits_list, healing_region_list = output_processing.filter_out_relevant_results(all_diamond_results,
                                                                                                  10,
                                                                                                  stored_low_cov_regions
                                                                                                  )

    healing_data_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                        "storage_files/healing_data.tsv"

    assembly_healing.create_detailed_healing_information_file(consid_diamond_hits_list, healing_data_path)

    print(len(healing_region_list), " queries with at least one frameshift found")
    print(len(all_diamond_results), " queries had at least one Diamond hit")

    old_assembly = "/share/gluster/assemblies/Tdraco/" \
                   "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    new_assembly_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_27.01.2022/" \
                       "healed_assembly/"

    new_fna_file_path, simplified_distance_distribution = assembly_healing.heal_assembly_file(healing_region_list, old_assembly,
                                                                                              new_assembly_dir)
    
    #"""


if __name__ == '__main__':
    main()

