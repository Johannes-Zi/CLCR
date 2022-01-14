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
    coverage_file_path = "/share/gluster/NOTLOESUNG/freya/T_draco/t_draco.pbc.wgs_short.txt"
    low_cov_regions = detectR.detect_regions(coverage_file_path, 15, 18)
    region_count = 0
    for scaffold in low_cov_regions:
        region_count += len(scaffold[1])
    print("Total short read low cov. regions before merging: ", region_count)

    low_cov_storage_tsv = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/storage_files/original_low_cov_regions.tsv"
    detectR.create_low_cov_tsv_file(low_cov_regions, low_cov_storage_tsv)

    # Create length distribution plot of the
    # len_distribution = low_cov_length_distribution_plot(low_cov_regions,
    #                                                    "/home/johannes/Desktop/low_cov_len_distribution_raw.png")

    low_cov_regions = detectR.merge_close_regions(low_cov_regions, 250)
    #"""

    """# Create length distribution plot of the
    print("Calculated length distribution")
    len_distribution = detectR.low_cov_length_distribution_plot(low_cov_regions,
                                                       "/home/johannes/Desktop/low_cov_len_distribution_merged.png")
    print(len_distribution)
    #"""

    """region_count = 0
    for scaffold in low_cov_regions:
        region_count += len(scaffold[1])
    print("Total short read low cov. regions after merging: ", region_count)

    # Create queries
    assembly_file = "/share/gluster/assemblies/Tdraco/" \
                    "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    # Create short read low cov. queries:
    query_dir_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/query_files/"
    compD.query_files_creation(low_cov_regions, assembly_file, 500, query_dir_path, 5000)
    # """

    """# Sending the slurm jobarray to the cluster
    protein_database = "/home/johannes/Desktop/trachinus_draco/protein_db/protein_db.dmnd"
    input_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/query_files/"
    output_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/output_files/"

    compD.create_slurmarry(protein_database, input_dir, output_dir)

    # sbatch --array=1-51 CLCR_slurmarray.slurm
    # """

    # Read in the diamond results

    """output_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/output_files/"

    all_diamond_results = outputP.read_in_diamond_output(output_dir)

    low_cov_storage_tsv = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/storage_files/original_low_cov_regions.tsv"
    stored_low_cov_regions = detectR.read_in_low_cov_tsv_file(low_cov_storage_tsv)

    considered_diamond_hits_list, healing_region_list = outputP.filter_out_relevant_results(all_diamond_results, 10,
                                                                                            stored_low_cov_regions)

    healing_data_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/storage_files/healing_data.tsv"

    outputC.create_detailed_healing_information_file(considered_diamond_hits_list, healing_data_path)

    #output_path = "/home/johannes/Desktop/considered_hits_len_distribution.png"
    #merged_len_dist = [[5, 71975], [25, 29666], [50, 24082], [100, 32263], [250, 35836], [500, 31516], [1000, 17666],
    #                    [5000, 9168], [float("inf"), 695]]
    #length_distribution = outputP.considered_diamond_hit_length_distribution_plot(considered_diamond_hits_list,
    #                                                                              output_path, merged_len_dist)
    #print(length_distribution)

    print(len(healing_region_list), " queries with at least one frameshift found")
    print(len(all_diamond_results), " queries had at least one Diamond hit")

    #diamond_query_list = outputP.check_overlapping_healing_positions_1(
    #    "/home/johannes/Desktop/trachinus_draco/TOGA_run_1_output/putative_false_corrected.tsv", healing_region_list)

    # Assembly file
    #assembly_file = "/share/gluster/assemblies/Tdraco/" \
    #                "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    # Create queries who might introduce putative wrong healing:
    #wrong_healing_query_dir = "/home/johannes/Desktop/trachinus_draco/putative_wrong_healed_queries/"
    #compD.custom_query_files_creation(diamond_query_list, assembly_file, 0, wrong_healing_query_dir, 1)

    # sbatch --array=1-53 CLCR_slurmarray.slurm


    old_assembly = "/share/gluster/assemblies/Tdraco/" \
                   "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    new_assembly_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/healed_assembly/"

    new_fna_file_path, simplified_distance_distribution = outputC.heal_assembly_file(healing_region_list, old_assembly,
                                                                                     new_assembly_dir)

    print(simplified_distance_distribution)
    
    #"""

    """# Evaluate the TOGA results

    HLtraDra1_file_path = "/home/johannes/Desktop/trachinus_draco/TOGA_run_1_output/loss_summ_data_HLtraDra1.tsv"
    HLtraDra3_file_path = "/home/johannes/Desktop/trachinus_draco/TOGA_run_1_output/loss_summ_data_HLtraDra3.tsv"
    toga_isoforms_tsv = "/home/johannes/Desktop/trachinus_draco/TOGA_run_1_output/toga.isoforms.tsv"
    query_annotation_gtf = "/home/johannes/Desktop/trachinus_draco/TOGA_run_1_output/query_annotation.gtf"

    results_dataframe = outputP.read_in_toga_lossgene_file(HLtraDra1_file_path, HLtraDra3_file_path, toga_isoforms_tsv,
                                                           query_annotation_gtf)

    #outputP.create_toga_result_plot(results_dataframe)

    # """

    # Evaluate the TOGA results 2

    putative_wrong_corrected_file_path = "/home/johannes/Desktop/trachinus_draco/toga_run_by_Michael_Hiller/putative_false_corrected.tsv"
    healing_data_tsv_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/storage_files/healing_data.tsv"
    output_tsv_path =  "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/storage_files/toga_result_analysis.tsv"

    outputP.check_overlapping_healing_positions_2(putative_wrong_corrected_file_path, healing_data_tsv_path,
                                                  output_tsv_path)

    diamond_output_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/output_files/"
    output_file_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/storage_files/toga_result_analysis_alignments.txt"
    outputP.search_relating_diamond_alignments(output_tsv_path, diamond_output_dir, output_file_path)


    # """


if __name__ == '__main__':
    main()

