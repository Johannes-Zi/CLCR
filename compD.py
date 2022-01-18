#!/usr/bin/python3

"""This file contains the functions for the protein database, slurmjob and queries creation"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os
import glob
import copy


def query_files_creation(list_with_scaffold_specific_low_cov_reg_lists, fna_file_path, min_length, query_dir,
                         seq_per_fasta):

    """
    This function creates the .fasta query files for a diamond slurm job, out of a input list with the low coverage
    regions of each and a path to the matching fna file to read out the sequences of the low coverage regions.
    The saved positions of the queries/low cov. regions in each scaffold are the normal positions, not python list
    positions(-1) !!!
    The query header consists of:
    >the scaffold # the start position of the merged low cov region in scaffold # end pos of the merged.. # start pos of
     the query in scaff. # end pos of the query..
    :param list_with_scaffold_specific_low_cov_reg_lists: list consisting of scaffold specific sublists with tuples
                                                          containing the start and he end of a region for database
                                                          comparison (output of detect_regions, merge_regions or
                                                          combine_regions_multiple_scaffolds).
    :param fna_file_path: path of the fna file for receiving the base sequence
    :param min_length: expanding the detected regions to 500 if they are smaller, for a better working on diamond
    :param query_dir: directory where the query files are stored in
    :param seq_per_fasta: max amount of queries stored per fasta, The fisrt and the last file could contain in special
                          cases more queries
    :return no return
    """

    # Read in the .fna file for fast access
    input_fna_file = open(fna_file_path)
    scaffold_list = []  # Filled with the sequences of the fna file
    temp_scaffold = []  # Contains the current scaffold
    current_header = ""  # Containing the header of the current scaffold

    # Filling the .fna region list/ reading in the original assembly
    for line in input_fna_file:

        # Sequence line case
        if (line[0] != ">") and (line[0] != ";"):
            for x in line.strip():  # filling up a new sequence
                temp_scaffold.append(x)

        # Sequence header case
        elif line[0] == ">":
            scaffold_list.append([current_header, temp_scaffold])  # appending the previous region
            current_header = line.strip().split()[0][1:]  # reading in the new header
            temp_scaffold = []  # resetting the list/ ready for new region

    # Append last remaining region
    scaffold_list.append([current_header, temp_scaffold])
    input_fna_file.close()
    # Delete the first empty initialising element
    scaffold_list.pop(0)

    # Initialise some parameters
    fasta_count = 1     # for individual naming at creating of the .fasta files
    count_added_queries = 0     # Counts how many queries are created overall
    already_added_queries = 0   # Counts how many queries are already added to the fasta file

    # Creating the first fasta file
    new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
    current_fasta = open(new_fasta, "w")

    # Creating the query files for all scaffolds, searching a the matching scaffold sequence for each
    # current_low_cov_list with the low cov. regions of each scaffold
    for scaffold in list_with_scaffold_specific_low_cov_reg_lists:
        scaffold_name = scaffold[0]
        current_low_cov_list = scaffold[1]
        current_tuple = 0  # saves the position of the current region tuple in current_low_cov_list

        for current_scaffold in scaffold_list:
            if scaffold_name == current_scaffold[0]:

                # Possible first cases where the region cant be expanded 250 to the left,
                # because the (startposition(of the region) - 250) < 0
                while (current_low_cov_list[current_tuple][0] - int(min_length / 2)) <= 0:

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

                    # Calculating the expand length
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]
                    if region_length < min_length:
                        expand_left = expand_right = int((min_length - region_length) / 2)
                    else:
                        expand_left = expand_right = 0

                    # Case if the starting point is moved in the negative by expand_left
                    if (current_low_cov_list[current_tuple][0] - expand_left) < 0:
                        region_start = 0
                        region_end = current_low_cov_list[current_tuple][1] + expand_right + \
                                     abs(current_low_cov_list[current_tuple][0] - expand_left)
                    # Normal case
                    else:
                        region_start = current_low_cov_list[current_tuple][0] - expand_left
                        region_end = current_low_cov_list[current_tuple][1] + expand_right

                    # Writing the query header in the fasta file
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) +
                                         "#" + str(region_start) + "#" + str(region_end) + "\n"))

                    # Writing the query sequence in the fasta file
                    current_fasta.write(("".join(current_scaffold[1][region_start:region_end]) + "\n"))
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                    # Break condition
                    if current_tuple >= len(current_low_cov_list):
                        break

                break_value = False
                # Each while run fills up a new fasta file
                while True:

                    # Break condition
                    if (current_tuple >= len(current_low_cov_list)) or \
                            ((current_low_cov_list[current_tuple][1]) >
                             (len(current_scaffold[1]) - int(min_length / 2))):
                        break_value = True
                        break

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

                    # Calculation the length a region must be expanded
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]

                    if region_length < min_length:
                        expand_len = int((min_length - region_length) / 2)
                    else:
                        expand_len = 0

                    # Writes the query header
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) +
                                         "#" + str(current_low_cov_list[current_tuple][0] - expand_len) +
                                         "#" + str(current_low_cov_list[current_tuple][1] + expand_len) + "\n"))

                    # Write the query sequence in the file
                    current_fasta.write("".join(current_scaffold[1]
                                                [(current_low_cov_list[current_tuple][0] - expand_len):
                                                 (current_low_cov_list[current_tuple][1] + expand_len)])
                                        + "\n")

                    # Update parameters
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                    # Break condition
                    if break_value:
                        break

                # Possible last cases where the region cant be expanded 250 to the right,
                # because (the endposition(of the region) + 250) > sequence length
                while current_tuple < len(current_low_cov_list):

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

                    # Calculating the expand length
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]
                    if region_length < min_length:
                        expand_left = expand_right = int((min_length - region_length) / 2)
                    else:
                        expand_left = expand_right = 0

                    # Case if the endpoint is moved behind the end of the current_scaffold by expand_right
                    if (current_low_cov_list[current_tuple][1] + expand_right) > len(current_scaffold[1]):
                        region_start = current_low_cov_list[current_tuple][0] - expand_left - \
                                       ((current_low_cov_list[current_tuple][1] + expand_right) -
                                        len(current_scaffold[1]))
                        region_end = len(current_scaffold[1])
                    # Normal case
                    else:
                        region_start = current_low_cov_list[current_tuple][0] - expand_left
                        region_end = current_low_cov_list[current_tuple][1] + expand_right

                    #  Writing the query header in the fasta file
                    current_fasta.write((">" + scaffold_name + "#" + str(current_low_cov_list[current_tuple][0]) +
                                         "#" + str(current_low_cov_list[current_tuple][1]) +
                                         "#" + str(region_start) + "#" + str(region_end) + "\n"))
                    # Writing the query sequence in the fasta file
                    current_fasta.write(("".join(current_scaffold[1][region_start:region_end]) + "\n"))
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                break       # Break because there is only one matching scaffold sequence for each scaffold

    # Close the last .fasta file
    current_fasta.close()

    # User information
    print("Amount of created fasta files: ", fasta_count)
    print("Amount of created queries: ", count_added_queries)

    return None


def custom_query_files_creation(list_with_scaffold_specific_low_cov_reg_lists, fna_file_path, min_length, query_dir,
                                seq_per_fasta):

    """
    CUSTOM version where the header for each query is saved at the third positions after the start and the end position
    of the query.

    This function creates the .fasta query files for a diamond slurm job, out of a input list with the low coverage
    regions of each and a path to the matching fna file to read out the sequences of the low coverage regions.
    The saved positions of the queries/low cov. regions in each scaffold are the normal positions, not python list
    positions(-1) !!!
    The query header consists of:
    >the scaffold # the start position of the merged low cov region in scaffold # end pos of the merged.. # start pos of
     the query in scaff. # end pos of the query..
    :param list_with_scaffold_specific_low_cov_reg_lists: list consisting of scaffold specific sublists with tuples
                                                          containing the start and he end of a region for database
                                                          comparison (output of detect_regions, merge_regions or
                                                          combine_regions_multiple_scaffolds).
    :param fna_file_path: path of the fna file for receiving the base sequence
    :param min_length: expanding the detected regions to 500 if they are smaller, for a better working on diamond
    :param query_dir: directory where the query files are stored in
    :param seq_per_fasta: max amount of queries stored per fasta, The fisrt and the last file could contain in special
                          cases more queries
    :return no return
    """

    # Read in the .fna file for fast access
    input_fna_file = open(fna_file_path)
    scaffold_list = []  # Filled with the sequences of the fna file
    temp_scaffold = []  # Contains the current scaffold
    current_header = ""  # Containing the header of the current scaffold

    # Filling the .fna region list/ reading in the original assembly
    for line in input_fna_file:

        # Sequence line case
        if (line[0] != ">") and (line[0] != ";"):
            for x in line.strip():  # filling up a new sequence
                temp_scaffold.append(x)

        # Sequence header case
        elif line[0] == ">":
            scaffold_list.append([current_header, temp_scaffold])  # appending the previous region
            current_header = line.strip().split()[0][1:]  # reading in the new header
            temp_scaffold = []  # resetting the list/ ready for new region

    # Append last remaining region
    scaffold_list.append([current_header, temp_scaffold])
    input_fna_file.close()
    # Delete the first empty initialising element
    scaffold_list.pop(0)

    # Initialise some parameters
    fasta_count = 1     # for individual naming at creating of the .fasta files
    count_added_queries = 0     # Counts how many queries are created overall
    already_added_queries = 0   # Counts how many queries are already added to the fasta file

    # Creating the first fasta file
    new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
    current_fasta = open(new_fasta, "w")

    # Creating the query files for all scaffolds, searching a the matching scaffold sequence for each
    # current_low_cov_list with the low cov. regions of each scaffold
    for scaffold in list_with_scaffold_specific_low_cov_reg_lists:
        scaffold_name = scaffold[0]
        current_low_cov_list = scaffold[1]
        current_tuple = 0  # saves the position of the current region tuple in current_low_cov_list

        for current_scaffold in scaffold_list:
            if scaffold_name == current_scaffold[0]:

                # Possible first cases where the region cant be expanded 250 to the left,
                # because the (startposition(of the region) - 250) < 0
                while (current_low_cov_list[current_tuple][0] - int(min_length / 2)) <= 0:

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

                    # Calculating the expand length
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]
                    if region_length < min_length:
                        expand_left = expand_right = int((min_length - region_length) / 2)
                    else:
                        expand_left = expand_right = 0

                    # Case if the starting point is moved in the negative by expand_left
                    if (current_low_cov_list[current_tuple][0] - expand_left) < 0:
                        region_start = 0
                        region_end = current_low_cov_list[current_tuple][1] + expand_right + \
                                     abs(current_low_cov_list[current_tuple][0] - expand_left)
                    # Normal case
                    else:
                        region_start = current_low_cov_list[current_tuple][0] - expand_left
                        region_end = current_low_cov_list[current_tuple][1] + expand_right

                    # Writing the query header in the fasta file
                    current_fasta.write((">" + scaffold[1][0][2] + "\n"))

                    # Writing the query sequence in the fasta file
                    current_fasta.write(("".join(current_scaffold[1][region_start:region_end]) + "\n"))
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                    # Break condition
                    if current_tuple >= len(current_low_cov_list):
                        break

                break_value = False
                # Each while run fills up a new fasta file
                while True:

                    # Break condition
                    if (current_tuple >= len(current_low_cov_list)) or \
                            ((current_low_cov_list[current_tuple][1]) >
                             (len(current_scaffold[1]) - int(min_length / 2))):
                        break_value = True
                        break

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

                    # Calculation the length a region must be expanded
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]

                    if region_length < min_length:
                        expand_len = int((min_length - region_length) / 2)
                    else:
                        expand_len = 0

                    # Writes the query header
                    current_fasta.write((">" + scaffold[1][0][2] + "\n"))

                    # Write the query sequence in the file
                    current_fasta.write("".join(current_scaffold[1]
                                                [(current_low_cov_list[current_tuple][0] - expand_len):
                                                 (current_low_cov_list[current_tuple][1] + expand_len)])
                                        + "\n")

                    # Update parameters
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                    # Break condition
                    if break_value:
                        break

                # Possible last cases where the region cant be expanded 250 to the right,
                # because (the endposition(of the region) + 250) > sequence length
                while current_tuple < len(current_low_cov_list):

                    # Creating .fasta file, if the current file is full (defined by seq_per_fasta)
                    if already_added_queries >= seq_per_fasta:
                        fasta_count += 1
                        current_fasta.close()
                        already_added_queries = 0
                        new_fasta = query_dir + "/" + "temp_in_" + str(fasta_count) + ".fasta"
                        current_fasta = open(new_fasta, "w")

                    # Calculating the expand length
                    region_length = current_low_cov_list[current_tuple][1] - current_low_cov_list[current_tuple][0]
                    if region_length < min_length:
                        expand_left = expand_right = int((min_length - region_length) / 2)
                    else:
                        expand_left = expand_right = 0

                    # Case if the endpoint is moved behind the end of the current_scaffold by expand_right
                    if (current_low_cov_list[current_tuple][1] + expand_right) > len(current_scaffold[1]):
                        region_start = current_low_cov_list[current_tuple][0] - expand_left - \
                                       ((current_low_cov_list[current_tuple][1] + expand_right) -
                                        len(current_scaffold[1]))
                        region_end = len(current_scaffold[1])
                    # Normal case
                    else:
                        region_start = current_low_cov_list[current_tuple][0] - expand_left
                        region_end = current_low_cov_list[current_tuple][1] + expand_right

                    #  Writing the query header in the fasta file
                    current_fasta.write((">" + scaffold[1][0][2] + "\n"))
                    # Writing the query sequence in the fasta file
                    current_fasta.write(("".join(current_scaffold[1][region_start:region_end]) + "\n"))
                    current_tuple += 1
                    count_added_queries += 1
                    already_added_queries += 1

                break       # Break because there is only one matching scaffold sequence for each scaffold

    # Close the last .fasta file
    current_fasta.close()

    # User information
    print("Amount of created fasta files: ", fasta_count)
    print("Amount of created queries: ", count_added_queries)

    return None


def create_slurmarry(protein_database, input_dir, output_dir):
    """
    Simply creates the slurmarray file for a given input and output directory (each dir with an / at the end!)
    :param protein_database: trivial
    :param input_dir: trivial
    :param output_dir: trivial
    :return: CLCR_slurmarray.slurm file at the cwd
    """

    # creating the .slurm file for the jobarray
    slurm_filename = os.getcwd() + "/" + "CLCR_slurmarray.slurm"
    slurm_file = open(slurm_filename, "w")

    # -d database, -q query, -o output path, -k max hits saved per query seq in output file
    blast_command = "diamond blastx -d " + protein_database + " "\
                    "-q " + input_dir + "temp_in_${SLURM_ARRAY_TASK_ID}.fasta " \
                    "-o " + output_dir + "temp_out_${SLURM_ARRAY_TASK_ID}.txt " \
                    "-k 25 --max-hsps 0 -c 1 -t /dev/shm -F 15 -f 0"

    slurm_file.write("#!/bin/bash\n\n"
                     "#SBATCH --partition=all\n"
                     "#SBATCH --account=praktikant\n"
                     "#SBATCH --cpus-per-task=4\n"
                     "#SBATCH --mem-per-cpu=16gb\n"
                     "#SBATCH --job-name=\"CLCR_run\"\n\n")

    slurm_file.write(blast_command)
    slurm_file.close()

    # sending the job via slurm to the cluster
    # os.system(("sbatch --array=1-" + str(fasta_count) + slurm_filename))

    return None


def main():
    print("database comparison main executed")


if __name__ == '__main__':
    main()
