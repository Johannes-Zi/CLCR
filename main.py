# coding=utf-8

"""This file ... is the main file :) """
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import detectR
import compD
import outputP
import outputC


def main():
    print("MAIN CALLED")

    """# Short read low cov. region detection:
    short_read_coverage_file = "/share/gluster/NOTLOESUNG/freya/T_draco/t_draco.pbc.wgs_short.txt"
    short_low_cov_regions = detectR.detect_regions(short_read_coverage_file, 15, 18)
    region_count = 0
    for scaffold in short_low_cov_regions:
        region_count += len(scaffold[1])
    print("Total short read low cov. regions before merging: ", region_count)

    # Create length distribution plot of the
    # len_distribution = low_cov_length_distribution_plot(short_low_cov_regions,
    #                                                    "/home/johannes/Desktop/low_cov_len_distribution_raw.png")

    short_low_cov_regions = detectR.merge_close_regions(short_low_cov_regions, 250)

    # Create length distribution plot of the
    # len_distribution = low_cov_length_distribution_plot(short_low_cov_regions,
    #                                                    "/home/johannes/Desktop/low_cov_len_distribution_merged.png")

    region_count = 0
    for scaffold in short_low_cov_regions:
        region_count += len(scaffold[1])
    print("Total short read low cov. regions after merging: ", region_count)

    # Create queries
    assembly_file = "/share/gluster/assemblies/Tdraco/" \
                    "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    # Create short read low cov. queries:
    short_read_query_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries/"
    compD.query_files_creation(short_low_cov_regions, assembly_file, 500, short_read_query_dir, 5000)
    # """

    """# Sending the slurm jobarray to the cluster
    protein_database = "/home/johannes/Desktop/trachinus_draco/protein_db/protein_db.dmnd"
    input_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries/"
    output_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries_output/"

    compD.create_slurmarry(protein_database, input_dir, output_dir)

    # sbatch --array=1-51 CLCR_slurmarray.slurm 
    # sbatch --array=1-5 CLCR_slurmarray.slurm
    # """

    # Read in the diamond resutls

    output_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries_output/"

    all_diamond_results = outputP.read_in_diamond_output(output_dir)

    considered_diamond_hits_list, healing_region_list = outputP.filter_out_relevant_results(all_diamond_results, 10)

    # output_path = "/home/johannes/Desktop/considered_hits_len_distribution.png"
    # length_distribution = outputP.considered_diamond_hit_length_distribution_plot(considered_diamond_hits_list,
    #                                                                              output_path)
    # print(length_distribution)

    print(len(healing_region_list), " queries with at least one frameshift found")
    print(len(all_diamond_results), " queries had at least one Diamond hit")

    old_assembly = "/share/gluster/assemblies/Tdraco/" \
                   "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    new_assembly_dir = "/home/johannes/Desktop/trachinus_draco/healed_assembly/"

    outputC.heal_assembly_file(healing_region_list, old_assembly, new_assembly_dir)
    
    #"""


if __name__ == '__main__':
    main()

