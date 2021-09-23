# coding=utf-8

"""This program ... is the main file :)"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import bruteF
import compD


def main():
    print("MAIN CALLED")

    # Short read low cov. region detection:
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
    compD.database_comparison(short_low_cov_regions, assembly_file, 500, short_read_query_dir, 5000)


if __name__ == '__main__':
    main()

