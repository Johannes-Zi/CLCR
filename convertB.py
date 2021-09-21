"""Function for extracting the coverage data out of the bam files by using samtools"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os


def convert_bam(bam_file, fna_file):
    """
    Calculating the per base coverage out of the raw .bam file
    :param bam_file: file path of the bam file
    :param fna_file: file path of the fna file
    :return: file path of the created coverage file
    """

    # Using the path tail for new filename
    file_name = os.path.split(bam_file)[1][:-4] + "_coverage.txt"

    # Index the assembly file that has been used as a reference
    command_sam = "samtools faidx " + fna_file

    # Compute per base coverage
    command_bed = "bedtools genomecov -ibam " + bam_file + " -d -g my_assembly.fna.fai> " + file_name

    # Execute the commands
    os.system(command_sam)
    os.system(command_bed)

    # Return the path of the created coverage .txt file
    return os.path.split(__file__)[0] + "/" + file_name


def main():
    print("Convert bam files main executed")


if __name__ == '__main__':
    main()
