# coding=utf-8

"""This program ... is the main file :)"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import bruteF
import compD


def main():
    print("MAIN CALLED")

    """# Short read low cov. region detection:
    short_read_coverage_file = "/share/gluster/NOTLOESUNG/freya/T_draco/t_draco.pbc.wgs_short.txt"
    short_low_cov_regions = bruteF.detect_regions(short_read_coverage_file, 15, 18)
    region_count = 0
    for scaffold in short_low_cov_regions:
        # print("##########\n" + scaffold[0] + "\n" + str(len(scaffold[1])) + "\n" + str(scaffold[1][0]))
        region_count += len(scaffold[1])
    print("Total short read low cov. regions before merging: ", region_count)
    short_low_cov_regions = bruteF.merge_close_reg(short_low_cov_regions, 250)

    region_count = 0
    for scaffold in short_low_cov_regions:
        #print("##########\n" + scaffold[0] + "\n" + str(len(scaffold[1])) + "\n" + str(scaffold[1][0]))
        region_count += len(scaffold[1])
    print("Total short read low cov. regions after merging: ", region_count)

    assembly_file = "/share/gluster/assemblies/Tdraco/" \
                    "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    # Create short read low cov. queries:
    short_read_query_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries/"
    compD.database_comparison(short_low_cov_regions, assembly_file, 500, short_read_query_dir, 5000)"""

    """# Sending the slurm jobarray to the cluster
    protein_database = "/home/johannes/Desktop/trachinus_draco/protein_db/protein_db.dmnd"
    input_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries/"
    output_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries_output/"

    compD.create_slurmarry(protein_database, input_dir, output_dir)

    # sbatch --array=1-51 CLCR_slurmarray.slurm sbatch --array=1-5 CLCR_slurmarray.slurm"""

    # Read in the diamond resutls

    output_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries_output/"

    output_region_list, healing_region_list = compD.read_in_results_3(output_dir)

    found_frameshift_count = 0
    regions_with_frameshift = 0
    queries_with_diamond_hit = 0
    for query in healing_region_list:
        queries_with_diamond_hit += 1
        if query[3]:
            found_frameshift_count += len(query[3])
            regions_with_frameshift += 1

    print(found_frameshift_count, " frameshifts found")
    print(regions_with_frameshift, " regions with at least one frameshift found")
    print(queries_with_diamond_hit, " queries had at least one Diamond hit")

    # 44639  frameshifts found
    # 13013  regions with at least one frameshift are detected
    # 44637  queries had at least one Diamond hit

    old_assembly = "/share/gluster/assemblies/Tdraco/" \
                    "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"

    new_assembly_dir = "/home/johannes/Desktop/trachinus_draco/healed_assembly/"

    compD.heal_assembly_file(healing_region_list, old_assembly, new_assembly_dir)


if __name__ == '__main__':
    main()

