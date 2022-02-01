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
    run_information = [run_information[0], str("Runtime:\t\t\t" + run_time)] + run_information[1:]

    # Create run information file
    run_info_file = open((storage_files_dir_path + "query_creation_" + time.strftime("%Y%m%d-%H%M%S") + ".txt"), "w")
    for line in run_information:
        run_info_file.write((line + "\n"))

    run_info_file.close()

    return None


def prepare_slurm_run(project_dir, protein_database, auto_run):
    """

    :param project_dir:
    :param protein_database:
    :param auto_run:
    :return:
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
    output_dir = project_dir + "diamond_output"
    slurm_file = slurm_dir_path + "CLCR_slurmarray.slurm"

    # Create slurmarray file
    slurmarry_creation.create_slurmarry(protein_database, input_dir, output_dir, slurm_file)

    # Start slurm job array
    if auto_run:
        # Count the input files
        input_file_list = glob.glob(input_dir + "/temp_in_*.txt")
        file_count = len(input_file_list)

        # Submit job
        os_command = "sbatch --array=1-" + str(file_count) + slurm_file
        os.system(os_command)

    return None


def main():
    print("CLCR MAIN CALLED")

    # New program structure version

    # Query creation
    project_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_01.02.2022/"
    cov_file_path = "/share/gluster/NOTLOESUNG/freya/T_draco/t_draco.pbc.wgs_short.txt"
    assembly_file = "/share/gluster/assemblies/Tdraco/" \
                    "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"
    low_cov_start = 15
    low_cov_end = 18
    min_query_len = 500
    queries_per_file = 5000

    create_queries(project_dir, cov_file_path, assembly_file, low_cov_start, low_cov_end, min_query_len,
                   queries_per_file)

    # Create and submit slurm array job




    """# Sending the slurm jobarray to the cluster
    protein_database = "/home/johannes/Desktop/trachinus_draco/protein_db/protein_db.dmnd"
    input_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/query_files/"
    output_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/output_files/"

    slurmarry_creation.create_slurmarry(protein_database, input_dir, output_dir)

    # sbatch --array=1-51 CLCR_slurmarray.slurm
    # """

    """# Read in the diamond results

    output_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/output_files/"

    all_diamond_results = output_processing.read_in_diamond_output(output_dir)

    low_cov_storage_tsv = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                          "storage_files/original_low_cov_regions.tsv"
    stored_low_cov_regions = query_creation.read_in_low_cov_tsv_file(low_cov_storage_tsv)

    considered_diamond_hits_list, healing_region_list = output_processing.filter_out_relevant_results(all_diamond_results, 10,
                                                                                                      stored_low_cov_regions)

    healing_data_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                        "storage_files/healing_data.tsv"

    assembly_healing.create_detailed_healing_information_file(considered_diamond_hits_list, healing_data_path)

    #output_path = "/home/johannes/Desktop/considered_hits_len_distribution.png"
    #merged_len_dist = [[5, 71975], [25, 29666], [50, 24082], [100, 32263], [250, 35836], [500, 31516], [1000, 17666],
    #                    [5000, 9168], [float("inf"), 695]]
    #length_distribution = outputP.considered_diamond_hit_length_distribution_plot(considered_diamond_hits_list,
    #                                                                              output_path, merged_len_dist)
    #print(length_distribution)

    print(len(healing_region_list), " queries with at least one frameshift found")
    print(len(all_diamond_results), " queries had at least one Diamond hit")

    # Assembly file
    #assembly_file = "/share/gluster/assemblies/Tdraco/" \
    #                "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    # Create queries who might introduce putative wrong healing:
    #wrong_healing_query_dir = "/home/johannes/Desktop/trachinus_draco/putative_wrong_healed_queries/"
    #compD.custom_query_files_creation(diamond_query_list, assembly_file, 0, wrong_healing_query_dir, 1)

    # sbatch --array=1-53 CLCR_slurmarray.slurm

    old_assembly = "/share/gluster/assemblies/Tdraco/" \
                   "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    new_assembly_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_27.01.2022/" \
                       "healed_assembly/"

    new_fna_file_path, simplified_distance_distribution = assembly_healing.heal_assembly_file(healing_region_list, old_assembly,
                                                                                              new_assembly_dir)

    #print(simplified_distance_distribution)
    
    #"""


if __name__ == '__main__':
    main()

