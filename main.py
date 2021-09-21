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
    short_low_cov_regions = bruteF.merge_close_reg(short_low_cov_regions, 250)

    # Long read low cov. region detection:
    long_read_coverage_file = "/share/gluster/NOTLOESUNG/freya/T_draco/t_draco.pbc.wgs_pb.txt"
    long_low_cov_regions = bruteF.detect_regions(long_read_coverage_file, 25, 30)
    long_low_cov_regions = bruteF.merge_close_reg(long_low_cov_regions, 250)

    # Determine combined low cov. regions:
    combined_low_cov_regions = bruteF.combine_regions(short_low_cov_regions, long_low_cov_regions)

    assembly_file = "/share/gluster/assemblies/Tdraco/" \
                    "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"
    # Create short read low cov. queries
    short_read_query_dir = "/home/johannes/Desktop/trachinus_draco/short_read_queries/"
    compD.database_comparison(short_low_cov_regions, assembly_file, 500, short_read_query_dir, 5000)

    # Create long read low cov. queries
    long_read_query_dir = "/home/johannes/Desktop/trachinus_draco/long_read_queries/"
    compD.database_comparison(long_low_cov_regions, assembly_file, 500, long_read_query_dir, 5000)

    # Create combined low cov. queries
    combined_query_dir = "/home/johannes/Desktop/trachinus_draco/combined_queries/"
    compD.database_comparison(combined_low_cov_regions, assembly_file, 500, combined_query_dir, 5000)


if __name__ == '__main__':
    main()

