#!python3

"""Second Benchmark for the CLCR program (correct frameshift healing)"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os
import glob
import random
import copy
import matplotlib.pyplot as plt
import numpy as np


def combine_faa(input_dir):
    """
    Combines the *_protein.faa files to a combined_protein_file.faa file, which is later used to create the main
    database.
    :param input_dir: Path to the input directory
    :return: None, but the combined_protein_file.faa is created at the cwd
    """

    # Returns the file path of all files in the input dir
    input_file_list = glob.glob(input_dir + "/*_protein.faa")

    combined_file_content = []  # All read lines of the input files are stored here

    print("The following files are used to create the database:")
    # Reads in the lines of each input .faa file
    for input_file in input_file_list:
        print(input_file)
        current_file = open(input_file)
        for line in current_file:
            combined_file_content.append(line)
        current_file.close()

    new_combined_file = open((os.getcwd() + "/combined_protein_file.faa"), "w")

    for line in combined_file_content:
        new_combined_file.write(line)

    new_combined_file.close()

    print("created combined_protein_file.faa at: " + os.getcwd() + "/")

    return None


def create_main_database():
    """
    Trivial function, simply creates the database. Part of the code to store the call in the code.
    :return: None
    """

    os.system("diamond makedb --in combined_protein_file.faa -d main_db ")

    return None


def create_ref_database():
    """
    Trivial function, simply creates the database. Part of the code to store the call in the code.
    :return: None
    """

    os.system("diamond makedb --in GCF_000004255.2_v.1.0_protein.faa -d ref_db ")

    return None


def create_mod_queries_old(input_cds_file_path, min_query_len, max_query_len):
    """
    Old version, where the file with the sequences, which are corrected, contains the complete CDS regions with the
    inserted frameshift. In the newer version the "short" queries are later corected and used to detect/check if the
    frameshift is correctly healed.

    Takes the CDS regions, given by the input file, and modifies them by inserting random frameshifts. Afterwards two
    output files are created at the cwd of the program. The mod_short_cds.fasta contains the modified (modified by
    inserting one frameshift per sequence) CDS subsequences with a length from 150 to 450, which are later used as
    queries for Diamond to detect the inserted frameshifts. The mod_complete_cds.fasta file contains the complete CDS
    sequences, which contain exacly the same frameshifts like the corresponding short subregion in the
    mod_short_cds.fasta file. This list is later used to check if the frameshift healing is working.
    :param input_cds_file_path: Path to the input file with all CDS regions of the main benchmark organism
    :param min_query_len: minimal query length
    :param max_query_len: maximal query length
    :return: None, But the previous defined files
    """

    input_cds_file = open(input_cds_file_path)

    cds_region_list = []       # Stores each CDS region header of the input file, together with the corresponding region
    current_cds_header = ""    # Stores the header of the current CDS region
    current_cds_seq = ""       # Stores the sequence of the current CDS region

    # Read in the input with the CDS regions
    for line in input_cds_file:
        if line[0] != "#":      # Excludes comments
            if line[0] == ">":      # Header of new region found
                # Append the previous region to the final list
                cds_region_list.append([current_cds_header, current_cds_seq])
                current_cds_seq = ""
                current_cds_header = line.strip()   # Set the header of the now following sequence

            else:
                # Creates the complete sequence by combining the subsequences in each line together
                current_cds_seq += line.strip()

    # Append the last remaining region
    cds_region_list.append([current_cds_header, current_cds_seq])

    cds_region_list.pop(0)      # Removes the initialising element

    # Initialise the list for the later query file creation
    short_mod_queries_list = []         # Stores the modified CDS subsequences with random length
    comp_mod_queries_list = []      # Stores the complete modified CDS regions, which are later used for correction

    # Creation of the modified CDS regions which are later used as queries and for the reading frame healing
    for cds_region in cds_region_list:

        # Choose a random subsequence of the current CDS sequence
        seq_len = random.randint(min_query_len, max_query_len)      # Random subsequence length from 100 to 400

        # Case when the random sequence length is bigger than the length of a very short CDS sequence
        if seq_len >= len(cds_region[1]):
            start_pos, seq_len = 0, len(cds_region[1])

        # Normal Case, where a random subsequence of the selected length is taken
        else:
            # Start pos of the subsequence in the CDS sequence
            start_pos = random.randint(0, (len(cds_region[1]) - seq_len))

        # Choose a random position in the subsequence, where the frameshift is inserted
        frameshift_pos = random.randint(start_pos, (start_pos + seq_len))

        # Create the modified sequences
        insertion = random.choice([True, False])

        # Insertion case
        if insertion:
            random_base = random.choice(["A", "T", "G", "C"])

            # Appends the complete new created modified CDS region to the output creation list
            comp_mod_queries_list.append([(cds_region[0].split()[0] + "#" + str(start_pos) + "#" +
                                           str(start_pos+seq_len) + "#" + str(frameshift_pos) + "#" + random_base),
                                          (cds_region[1][:frameshift_pos] + random_base +
                                           cds_region[1][frameshift_pos:])])

            # Appends the new created modified CDS subsequence to the corresponding output list
            short_mod_queries_list.append([(cds_region[0].split()[0] + "#" + str(start_pos) + "#" +
                                            str(start_pos+seq_len) + "#" + str(frameshift_pos) + "#" + random_base),
                                           (cds_region[1][start_pos:frameshift_pos] + random_base +
                                            cds_region[1][frameshift_pos:(start_pos + seq_len)])])

        # Deletion case
        else:
            # Appends the complete new created modified CDS region to the output creation list
            comp_mod_queries_list.append([(cds_region[0].split()[0] + "#" + str(start_pos) + "#" +
                                           str(start_pos+seq_len) + "#" + str(frameshift_pos) + "#D"),
                                          (cds_region[1][:(frameshift_pos - 1)] + cds_region[1][frameshift_pos:])])

            # Appends the new created modified CDS subsequence to the corresponding output list
            short_mod_queries_list.append([(cds_region[0].split()[0] + "#" + str(start_pos) + "#" +
                                            str(start_pos+seq_len) + "#" + str(frameshift_pos) + "#D"),
                                           (cds_region[1][start_pos:(frameshift_pos - 1)] +
                                            cds_region[1][frameshift_pos:(start_pos + seq_len)])])

    # Create the output files
    mod_comp_file_path = os.getcwd() + "/mod_complete_cds.fasta"
    mod_short_file_path = os.getcwd() + "/mod_short_cds.fasta"

    # Create the output file with the modified complete CDS sequences, for later frameshift healing
    new_mod_comp_file = open(mod_comp_file_path, "w")

    for cds_region in comp_mod_queries_list:
        new_mod_comp_file.write((cds_region[0] + "\n"))

        # Normal case where the nucleotide sequence is longer than 80
        if len(cds_region[1]) > 80:
            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(cds_region[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_mod_comp_file.write((cds_region[1][(list_pos - 80):list_pos] + "\n"))

            # Write the last remaining line (range function doesnt include the last elements)
            new_mod_comp_file.write((cds_region[1][position_list[-1]:len(cds_region[1])] + "\n"))

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_mod_comp_file.write(cds_region[1] + "\n")

    new_mod_comp_file.close()

    # Create the output file with the modified CDS subsequences, for the frameshift detection
    new_mod_short_file = open(mod_short_file_path, "w")

    for cds_region in short_mod_queries_list:
        new_mod_short_file.write((cds_region[0] + "\n"))

        # Normal case where the nucleotide sequence is longer than 80
        if len(cds_region[1]) > 80:
            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(cds_region[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_mod_short_file.write((cds_region[1][(list_pos - 80):list_pos] + "\n"))

            # Write the last remaining line (range function doesnt include the last elements)
            new_mod_short_file.write((cds_region[1][position_list[-1]:len(cds_region[1])] + "\n"))

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_mod_short_file.write(cds_region[1] + "\n")

    new_mod_short_file.close()

    return None


def create_mod_queries_new(input_cds_file_path, min_query_len, max_query_len):
    """
    NEW VERSION with only one output file.
    Takes the CDS regions, given by the input file, and modifies them by inserting random frameshifts, the length of the
    subsequences is determined by the min/max_query_len. The modified CDS subsequences are used as queries for the
    frameshift detection and are later healed with the frameshift detection/localisation information of the Diamond
    program part. When the queries and the sequences which are healed are the same, the part, with localising the
    frameshift over the position in the query and over the position of the query in the complete sequence, falls out.
    So the queries start positions in the complete sequences are always 0.
    In summary two files are created, the one with the modified querys(with frameshifts) and the file with the
    unmodified queries, which are used as reference.

    :param input_cds_file_path: Path to the input file with all CDS regions of the main benchmark organism
    :param min_query_len: minimal query length
    :param max_query_len: maximal query length
    :return: None, But the previous defined files
    """

    input_cds_file = open(input_cds_file_path)

    cds_region_list = []       # Stores each CDS region header of the input file, together with the corresponding region
    current_cds_header = ""    # Stores the header of the current CDS region
    current_cds_seq = ""       # Stores the sequence of the current CDS region

    # Read in the input with the CDS regions
    for line in input_cds_file:
        if line[0] != "#":      # Excludes comments
            if line[0] == ">":      # Header of new region found
                # Append the previous region to the final list
                cds_region_list.append([current_cds_header, current_cds_seq])
                current_cds_seq = ""
                current_cds_header = line.strip()   # Set the header of the now following sequence

            else:
                # Creates the complete sequence by combining the subsequences in each line together
                current_cds_seq += line.strip()

    # Append the last remaining region
    cds_region_list.append([current_cds_header, current_cds_seq])

    cds_region_list.pop(0)      # Removes the initialising element

    # Initialise the list for the later query file creation
    mod_queries_list = []         # Stores the modified CDS subsequences with random length
    unmod_queries_list = []       # contains the unmodified queries fot the reference run

    # Creation of the modified CDS regions which are later used as queries and for the reading frame healing
    for cds_region in cds_region_list:

        # Choose a random subsequence of the current CDS sequence
        seq_len = random.randint(min_query_len, max_query_len)      # Random subsequence length from 100 to 400

        # Case when the random sequence length is bigger than the length of a very short CDS sequence
        if seq_len >= len(cds_region[1]):
            start_pos, seq_len = 0, len(cds_region[1])

        # Normal Case, where a random subsequence of the selected length is taken
        else:
            # Start pos of the subsequence in the CDS sequence
            start_pos = random.randint(0, (len(cds_region[1]) - seq_len))

        # Choose a random position in the subsequence, where the frameshift is inserted
        frameshift_pos = random.randint(start_pos, (start_pos + seq_len))

        # Create the modified sequences
        insertion = random.choice([True, False])

        # Insertion case
        if insertion:
            random_base = random.choice(["A", "T", "G", "C"])

            # Appends the new created modified CDS subsequence to the corresponding output list
            mod_queries_list.append([(cds_region[0].split()[0] + "#" + str(0) + "#" +
                                      str(seq_len) + "#" + str(frameshift_pos - start_pos) + "#" + random_base),
                                     (cds_region[1][start_pos:frameshift_pos] + random_base +
                                      cds_region[1][frameshift_pos:(start_pos + seq_len)])])
            # Appends the unmodified CDS subsequence to the corresponding output list
            unmod_queries_list.append([(cds_region[0].split()[0] + "#" + str(0) + "#" +
                                       str(seq_len) + "#" + str(frameshift_pos - start_pos) + "#" + random_base),
                                       cds_region[1][start_pos:(start_pos + seq_len)]])

        # Deletion case
        else:

            # Appends the new created modified CDS subsequence to the corresponding output list
            mod_queries_list.append([(cds_region[0].split()[0] + "#" + str(0) + "#" +
                                      str(seq_len) + "#" + str(frameshift_pos - start_pos) + "#D"),
                                     (cds_region[1][start_pos:(frameshift_pos - 1)] +
                                      cds_region[1][frameshift_pos:(start_pos + seq_len)])])
            # Appends the unmodified CDS subsequence to the corresponding output list
            unmod_queries_list.append([(cds_region[0].split()[0] + "#" + str(0) + "#" +
                                        str(seq_len) + "#" + str(frameshift_pos - start_pos) + "#D"),
                                       cds_region[1][start_pos:(start_pos + seq_len)]])

    # Create the output file
    mod_short_file_path = os.getcwd() + "/mod_queries.fasta"
    unmod_short_file_path = os.getcwd() + "/unmod_queries.fasta"

    # Create the output file with the modified CDS subsequences, for the frameshift detection
    new_mod_file = open(mod_short_file_path, "w")

    for cds_region in mod_queries_list:
        new_mod_file.write((cds_region[0] + "\n"))

        # Normal case where the nucleotide sequence is longer than 80
        if len(cds_region[1]) > 80:
            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(cds_region[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_mod_file.write((cds_region[1][(list_pos - 80):list_pos] + "\n"))

            # Write the last remaining line (range function doesnt include the last elements)
            new_mod_file.write((cds_region[1][position_list[-1]:len(cds_region[1])] + "\n"))

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_mod_file.write(cds_region[1] + "\n")

    new_mod_file.close()

    # Create the output file with the unmodified CDS subsequences, for the frameshift detection
    new_unmod_file = open(unmod_short_file_path, "w")

    for cds_region in unmod_queries_list:
        new_unmod_file.write((cds_region[0] + "\n"))

        # Normal case where the nucleotide sequence is longer than 80
        if len(cds_region[1]) > 80:
            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(cds_region[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_unmod_file.write((cds_region[1][(list_pos - 80):list_pos] + "\n"))

            # Write the last remaining line (range function doesnt include the last elements)
            new_unmod_file.write((cds_region[1][position_list[-1]:len(cds_region[1])] + "\n"))

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_unmod_file.write(cds_region[1] + "\n")

    new_unmod_file.close()

    return None


def create_mod_queries_new_multiple_bases(input_cds_file_path, min_query_len, max_query_len):
    """
    MODIFIED version, which inserts/deletes multiple bases per frameshift position(1, 2, 4, 5, 7, 8 bases per position).
    NEW VERSION with only one output file.
    Takes the CDS regions, given by the input file, and modifies them by inserting random frameshifts, the length of the
    subsequences is determined by the min/max_query_len. The modified CDS subsequences are used as queries for the
    frameshift detection and are later healed with the frameshift detection/localisation information of the Diamond
    program part. When the queries and the sequences which are healed are the same, the part, with localising the
    frameshift over the position in the query and over the position of the query in the complete sequence, falls out.
    So the queries start positions in the complete sequences are always 0.
    In summary two files are created, the one with the modified querys(with frameshifts) and the file with the
    unmodified queries, which are used as reference

    :param input_cds_file_path: Path to the input file with all CDS regions of the main benchmark organism
    :param min_query_len: minimal query length
    :param max_query_len: maximal query length
    :return: None, But the previous defined files
    """

    input_cds_file = open(input_cds_file_path)

    cds_region_list = []       # Stores each CDS region header of the input file, together with the corresponding region
    current_cds_header = ""    # Stores the header of the current CDS region
    current_cds_seq = ""       # Stores the sequence of the current CDS region

    # Read in the input with the CDS regions
    for line in input_cds_file:
        if line[0] != "#":      # Excludes comments
            if line[0] == ">":      # Header of new region found
                # Append the previous region to the final list
                cds_region_list.append([current_cds_header, current_cds_seq])
                current_cds_seq = ""
                current_cds_header = line.strip()   # Set the header of the now following sequence

            else:
                # Creates the complete sequence by combining the subsequences in each line together
                current_cds_seq += line.strip()

    # Append the last remaining region
    cds_region_list.append([current_cds_header, current_cds_seq])

    cds_region_list.pop(0)      # Removes the initialising element

    # Initialise the list for the later query file creation
    mod_queries_list = []         # Stores the modified CDS subsequences with random length
    unmod_queries_list = []       # contains the unmodified queries fot the reference run

    # Creation of the modified CDS regions which are later used as queries and for the reading frame healing
    for cds_region in cds_region_list:

        # Choose a random subsequence of the current CDS sequence
        seq_len = random.randint(min_query_len, max_query_len)      # Random subsequence length from 100 to 400

        # Case when the random sequence length is bigger than the length of a very short CDS sequence
        if seq_len >= len(cds_region[1]):
            start_pos, seq_len = 0, len(cds_region[1])

        # Normal Case, where a random subsequence of the selected length is taken
        else:
            # Start pos of the subsequence in the CDS sequence
            start_pos = random.randint(0, (len(cds_region[1]) - seq_len))

        # Choose a random position in the subsequence, where the frameshift is inserted
        frameshift_pos = random.randint(start_pos, (start_pos + seq_len))

        # Create the modified sequences
        insertion = random.choice([True, False])

        # Insertion case
        if insertion:

            # Create a insertion with a random length(borders set by the parameters) and with random bases
            random_bases = ""
            insertion_length = random.choice([1, 2, 4, 5, 7, 8])

            for step in range(insertion_length):
                random_base = random.choice(["A", "T", "G", "C"])
                random_bases += random_base

            # Appends the new created modified CDS subsequence to the corresponding output list
            mod_queries_list.append([(cds_region[0].split()[0] + "#" + str(0) + "#" +
                                      str(seq_len) + "#" + str(frameshift_pos - start_pos) + "#" + random_bases),
                                     (cds_region[1][start_pos:frameshift_pos] + random_bases +
                                      cds_region[1][frameshift_pos:(start_pos + seq_len)])])
            # Appends the unmodified CDS subsequence to the corresponding output list
            unmod_queries_list.append([(cds_region[0].split()[0] + "#" + str(0) + "#" +
                                       str(seq_len) + "#" + str(frameshift_pos - start_pos) + "#" + random_bases),
                                       cds_region[1][start_pos:(start_pos + seq_len)]])

        # Deletion case
        else:

            # Adapt the frameshift position if necessary, if the frameshift position is the last position, no following
            # bases could be deleted...

            # Shifts the frameshift position non random on a determined amount of bases upstream
            if (frameshift_pos + 8) > (start_pos + seq_len):

                # Shifts the frameshift position 8 positions upstream
                frameshift_pos -= 8

            deletion_length = random.choice([1, 2, 4, 5, 7, 8])

            # Appends the new created modified CDS subsequence to the corresponding output list
            mod_queries_list.append([(cds_region[0].split()[0] + "#" + str(0) + "#" +
                                      str(seq_len) + "#" + str(frameshift_pos - start_pos) + "#D" +
                                      str(deletion_length)),
                                     (cds_region[1][start_pos:frameshift_pos] +
                                      cds_region[1][(frameshift_pos + deletion_length):(start_pos + seq_len)])])
            # Appends the unmodified CDS subsequence to the corresponding output list
            unmod_queries_list.append([(cds_region[0].split()[0] + "#" + str(0) + "#" +
                                        str(seq_len) + "#" + str(frameshift_pos - start_pos) + "#D" +
                                        str(deletion_length)),
                                       cds_region[1][start_pos:(start_pos + seq_len)]])

    # Create the output file
    mod_short_file_path = os.getcwd() + "/mod_queries.fasta"
    unmod_short_file_path = os.getcwd() + "/unmod_queries.fasta"

    # Create the output file with the modified CDS subsequences, for the frameshift detection
    new_mod_file = open(mod_short_file_path, "w")

    for cds_region in mod_queries_list:
        new_mod_file.write((cds_region[0] + "\n"))

        # Normal case where the nucleotide sequence is longer than 80
        if len(cds_region[1]) > 80:
            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(cds_region[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_mod_file.write((cds_region[1][(list_pos - 80):list_pos] + "\n"))

            # Write the last remaining line (range function doesnt include the last elements)
            new_mod_file.write((cds_region[1][position_list[-1]:len(cds_region[1])] + "\n"))

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_mod_file.write(cds_region[1] + "\n")

    new_mod_file.close()

    # Create the output file with the unmodified CDS subsequences, for the frameshift detection
    new_unmod_file = open(unmod_short_file_path, "w")

    for cds_region in unmod_queries_list:
        new_unmod_file.write((cds_region[0] + "\n"))

        # Normal case where the nucleotide sequence is longer than 80
        if len(cds_region[1]) > 80:
            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(cds_region[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_unmod_file.write((cds_region[1][(list_pos - 80):list_pos] + "\n"))

            # Write the last remaining line (range function doesnt include the last elements)
            new_unmod_file.write((cds_region[1][position_list[-1]:len(cds_region[1])] + "\n"))

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_unmod_file.write(cds_region[1] + "\n")

    new_unmod_file.close()

    return None


def create_mod_queries_new_multiple_positions(input_cds_file_path, min_query_len, max_query_len):
    """
    MODIFIED version which inserts 4 frameshifts per query
    NEW VERSION with only one output file.
    Takes the CDS regions, given by the input file, and modifies them by inserting random frameshifts, the length of the
    subsequences is determined by the min/max_query_len. The modified CDS subsequences are used as queries for the
    frameshift detection and are later healed with the frameshift detection/localisation information of the Diamond
    program part. When the queries and the sequences which are healed are the same, the part, with localising the
    frameshift over the position in the query and over the position of the query in the complete sequence, falls out.
    So the queries start positions in the complete sequences are always 0.
    In summary two files are created, the one with the modified querys(with frameshifts) and the file with the
    unmodified queries, which are used as reference

    :param input_cds_file_path: Path to the input file with all CDS regions of the main benchmark organism
    :param min_query_len: minimal query length
    :param max_query_len: maximal query length
    :return: None, But the previous defined files
    """

    input_cds_file = open(input_cds_file_path)

    cds_region_list = []       # Stores each CDS region header of the input file, together with the corresponding region
    current_cds_header = ""    # Stores the header of the current CDS region
    current_cds_seq = ""       # Stores the sequence of the current CDS region

    # Read in the input with the CDS regions
    for line in input_cds_file:
        if line[0] != "#":      # Excludes comments
            if line[0] == ">":      # Header of new region found
                # Append the previous region to the final list
                cds_region_list.append([current_cds_header, current_cds_seq])
                current_cds_seq = ""
                current_cds_header = line.strip()   # Set the header of the now following sequence

            else:
                # Creates the complete sequence by combining the subsequences in each line together
                current_cds_seq += line.strip()

    # Append the last remaining region
    cds_region_list.append([current_cds_header, current_cds_seq])

    cds_region_list.pop(0)      # Removes the initialising element

    # Initialise the list for the later query file creation
    mod_queries_list = []         # Stores the modified CDS subsequences with random length
    unmod_queries_list = []       # contains the unmodified queries fot the reference run

    # Creation of the modified CDS regions which are later used as queries and for the reading frame healing
    for cds_region in cds_region_list:

        # Choose a random subsequence of the current CDS sequence
        seq_len = random.randint(min_query_len, max_query_len)      # Random subsequence length from 100 to 400

        # Case when the random sequence length is bigger than the length of a very short CDS sequence
        if seq_len >= len(cds_region[1]):
            start_pos, seq_len = 0, len(cds_region[1])

        # Normal Case, where a random subsequence of the selected length is taken
        else:
            # Start pos of the subsequence in the CDS sequence
            start_pos = random.randint(0, (len(cds_region[1]) - seq_len))

        # Holds the current query sequences
        mod_query_sequence = cds_region[1][start_pos:(start_pos + seq_len)]
        unmod_query_sequence = cds_region[1][start_pos:(start_pos + seq_len)]

        # Choose the random 4 frameshift positions
        frameshift_positions = []
        for pos in range(4):
            # Choose a random position in the subsequence, where a frameshift is inserted
            # Stores later the frameshift type at the second position, X is for initialising
            # -4 that there are no problems with list out of range, because the sequence length varies with deletion
            frameshift_positions.append([random.randint(0, (seq_len - 4)), "X"])

        # Sorts the positions for a better overview, and a correct position assignment
        frameshift_positions.sort()

        for frameshift_position in frameshift_positions:

            # Random insertion or deletion
            insertion = random.choice([True, False])

            # Insertion case
            if insertion:
                random_base = random.choice(["A", "T", "G", "C"])

                # Saves what frameshift was inserted
                frameshift_position[1] = random_base

                mod_query_sequence = mod_query_sequence[0:(frameshift_position[0] - 1)] + random_base + \
                                     mod_query_sequence[(frameshift_position[0] - 1):len(mod_query_sequence)]

            # Deletion case
            else:
                # Saves what frameshift was inserted
                frameshift_position[1] = "D"

                mod_query_sequence = mod_query_sequence[0:(frameshift_position[0] - 1)] + \
                                     mod_query_sequence[frameshift_position[0]:len(mod_query_sequence)]

            # Appends the new created queries to the corresponding output lists
            query_header = (cds_region[0].split()[0] + "#" + str(0) + "#" + str(seq_len) + "#" +
                                      str(frameshift_positions[0][0]) + "#" + str(frameshift_positions[1][0]) + "#" +
                                      str(frameshift_positions[2][0]) + "#" + str(frameshift_positions[3][0]) + "#" +
                                      str(frameshift_positions[0][1]) + str(frameshift_positions[1][1]) +
                                      str(frameshift_positions[2][1]) + str(frameshift_positions[3][1]))

        mod_queries_list.append([query_header, mod_query_sequence])
        unmod_queries_list.append([query_header, unmod_query_sequence])

    # Create the output file
    mod_short_file_path = os.getcwd() + "/mod_queries.fasta"
    unmod_short_file_path = os.getcwd() + "/unmod_queries.fasta"

    # Create the output file with the modified CDS subsequences, for the frameshift detection
    new_mod_file = open(mod_short_file_path, "w")

    for cds_region in mod_queries_list:
        new_mod_file.write((cds_region[0] + "\n"))

        # Normal case where the nucleotide sequence is longer than 80
        if len(cds_region[1]) > 80:
            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(cds_region[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_mod_file.write((cds_region[1][(list_pos - 80):list_pos] + "\n"))

            # Write the last remaining line (range function doesnt include the last elements)
            new_mod_file.write((cds_region[1][position_list[-1]:len(cds_region[1])] + "\n"))

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_mod_file.write(cds_region[1] + "\n")

    new_mod_file.close()

    # Create the output file with the unmodified CDS subsequences, for the frameshift detection
    new_unmod_file = open(unmod_short_file_path, "w")

    for cds_region in unmod_queries_list:
        new_unmod_file.write((cds_region[0] + "\n"))

        # Normal case where the nucleotide sequence is longer than 80
        if len(cds_region[1]) > 80:
            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(cds_region[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_unmod_file.write((cds_region[1][(list_pos - 80):list_pos] + "\n"))

            # Write the last remaining line (range function doesnt include the last elements)
            new_unmod_file.write((cds_region[1][position_list[-1]:len(cds_region[1])] + "\n"))

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_unmod_file.write(cds_region[1] + "\n")

    new_unmod_file.close()

    return None


def exclude_putative_transition_frameshift(query_alignment, subject_alignment, frameshift_position):

    """
    Checks if the detected frameshifts correlates with a putative intron an could be thus neglected.
    :param query_alignment: aligned query subsequence
    :param subject_alignment: aligned subject subsequence
    :param frameshift_position: position where a frameshift was detected
    :return: returns True if the frameshift correlates with an putative intron, and False if not
    """

    # checks if the detected frameshift is a result of a putative Intron
    # if case for normal, non edge positions
    if not (((frameshift_position + 18) > len(query_alignment)) or ((frameshift_position - 17) < 0)):
        # later returned, True if the frameshifts correlates with a putative intron and could be neglected
        wrong_intron_left = True
        wrong_intron_right = True

        # used range starts at i-2 to be fault tolerant in intron-exon transition positions
        for y in range(frameshift_position - 17, frameshift_position - 2):  # upstream region
            if subject_alignment[y] != "-":
                wrong_intron_left = False
                break

        for y in range(frameshift_position + 3, frameshift_position + 18):  # downstream region
            if subject_alignment[y] != "-":
                wrong_intron_right = False
                break

        wrong_intron = wrong_intron_left or wrong_intron_right

    # case for frameshifts at the edge of the alignment
    else:
        # stays True if the frameshift might be the result of a putative intron
        wrong_intron = True
        # case for a downstream edge position
        if (frameshift_position + 18) > len(query_alignment):
            for y in range(frameshift_position - 17, frameshift_position - 2):
                if subject_alignment[y] != "-":
                    wrong_intron = False
                    break

        # case for a upstream edge position
        else:
            for y in range(frameshift_position + 3, frameshift_position + 18):
                if subject_alignment[y] != "-":
                    wrong_intron = False
                    break

    return wrong_intron


def read_in_results_3_detect(output_dir):
    """
    Modified version, output region list, saves additionally the original frameshift position for each region
    (additionally at the end of each query list)
    DETECTION version for benchmark, with a different header parsing, for the frameshift correction
    The function reads in the Diamond blastx results in a given directory. The results for each hit are stored in the
    output_region_list, which is also returned. The contains for every region: Scaffold, start pos., end pos. in
    scaffold, e-value, bitscore, similarity-percentage, list with all detected frameshifts with their positions in the
    query (I = insertion, D = deletion). The second output list (healing_region_list), contains the combined frameshift
    information of all Diamond hits per query, which means the scaffold, the start position in the scaffold, the end
    position in the scaffold, and a list with all positions in the query, where frameshifts are detected and need to
    be healed later . But not all Diamond hits of each Query are included, some are excluded by our overlapping
    heuristic, more detailed information in the documentation and in the code.
    PS: frameshift positions are the python list positions.
    :param output_dir: path to the output directory
    :return: output_region_list, healing_region_list
    """

    # Returns every file in the directory with .out at the end
    output_file_list = glob.glob(output_dir + "/temp_out_*.txt")

    # Sorting the file path ascending by their file number
    dir_path_len = len(output_dir) + 9  # +9 because temp_out_ consists of 9 chars
    output_file_list = sorted(output_file_list, key=lambda current_path: int(current_path[dir_path_len:-4]))

    output_region_list = []  # Containing the Diamond hits with additional information

    # Holds the length of the queries with no diamond hit, for later analysis
    no_hit_length_distribution = []

    # Saves for each Diamond hit the distance from the original frameshift position to the detected frameshift position
    original_frameshift_position = 0    # Initialising

    # Reading out all output files
    for x in output_file_list:
        current_output_file = open(x, "r")

        # Initialising the parameters
        query_hit_list = []    # Contains all hit of the current query
        current_region = (0, 0, 0)  # Scaffold, start and end position of the current considered region
        frameshift_list = []  # Contains the positions of detected frameshifts
        query_alignment = ""
        # The following two parameters are important for the localising of the frameshift in the scaffold, and the later
        # Overlapping heuristic for decide which frameshifts are considered in each query
        query_alignment_start_pos = 0   # Saves the start position of the alignment in the query
        query_alignment_end_pos = 0     # for saving the current end position of each alignment
        subject_alignment = ""
        # The following three store the respective values
        protein_hit = ""
        e_value = ""
        bit_score = ""
        similarity_percentage = ""

        # Reads in the current output file
        for line in current_output_file:
            # Case for a new Query
            if line.startswith("Query="):

                # Case if the query_hit_list, which holds the hits of each query, is not empty,
                # which means the previous region had an Diamond hit
                if query_hit_list:

                    # For every query_hit_list per query, the first initialising element must be deleted and the last
                    # hit must be afterwards be appended

                    # Appends the last remaining protein hit of the query, because each protein hit of the query is
                    # appended to the list by its following hit,which means that the last one needs to be appended later
                    query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                           subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

                    # Reset the parameter
                    query_alignment_start_pos = 0

                    # Removes the initialising element of the query, because the first element appends the last protein
                    # hit of the previous query, which is a result of the used approach
                    query_hit_list.pop(0)

                    # Appends each protein hit of the query to the output_region_list
                    for current_protein_hit in query_hit_list:

                        # The placed gaps in the query in the alignment are counted for a correct position calculation
                        query_gap_count = 0

                        # Detect frameshifts in aligned query region
                        for query_pos, char in enumerate(current_protein_hit[4]):

                            # Counts how many gaps are inserted into the query before the frameshift position
                            if char == "-":
                                query_gap_count += 1

                            if char == "\\":  # Insertion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      query_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (query_pos - query_gap_count)), "I"))

                            elif char == "/":  # Deletion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      query_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (query_pos - query_gap_count)), "D"))

                        output_region_list.append([current_region[0], current_region[1], current_region[2],
                                                   current_protein_hit[0], current_protein_hit[1],
                                                   current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                                   current_protein_hit[6], current_protein_hit[7],
                                                   original_frameshift_position])

                        # Reset values
                        frameshift_list = []

                    # Reset the the list, the if request is inactive till a new protein hit was found and thus
                    # initialising element is created
                    query_hit_list = []

                # Case if the query had no diamond hit
                else:
                    no_hit_length_distribution.append(current_region[2])

                # Setting the new region
                current_line = line.strip().split()[1].split("#")
                current_region = (current_line[0], current_line[1], current_line[2])

                # Saves the original position of the frameshift, that was inserted
                original_frameshift_position = current_line[3]

            # Line containing the protein id
            if line.startswith(">"):

                # Appends the previous protein hit to the list with all hits of each query
                # If This is the first hit of te current query, the last protein hit of the previous query is appended,
                # but this "initialising" element is later at the Query= if request removed. In a normal case each hit
                # is appended by the following hit of the same query, and the last hit is also appended at the Query= if
                query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                      subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

                # Reset/set the new parameters
                protein_hit = line.split()[0][1:]  # containing the protein hit of the current region
                subject_alignment = ""
                query_alignment = ""
                query_alignment_start_pos = 0

            # Line containing score and e-value
            if line.startswith(" Score"):
                current_line = line.strip().split()
                bit_score = current_line[2]
                e_value = current_line[7]

            # Line similarity percentage
            if line.startswith(" Identities"):
                current_line = line.strip().split()
                similarity_percentage = current_line[3][1:-2]

            # Reading in the aligned sequences
            if line.startswith("Query "):
                # Case that a new alignment starts, and the new start position is registered
                if not query_alignment_start_pos:
                    query_alignment_start_pos = int(line.split()[1])

                query_alignment_end_pos = int(line.split()[3])   # Append the currently last alignment position
                query_alignment += line.split()[2]

            # Subject alignment line case
            if line.startswith("Sbjct "):
                subject_alignment += line.split()[2]

        # Appending the last region if necessary (same as upper part)
        if query_hit_list:
            # Appends the last remaining protein hit of the query
            query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                   subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

            # Removes the initialising element of the query
            query_hit_list.pop(0)

            # Appends each protein hit of the query to the output_region_list
            for current_protein_hit in query_hit_list:

                # The placed gaps in the query in the alignment are counted for a correct position calculation
                query_gap_count = 0

                # Detect frameshifts in aligned query region
                for query_pos, char in enumerate(current_protein_hit[4]):

                    # Counts how many gaps are inserted into the query before the frameshift position
                    if char == "-":
                        query_gap_count += 1

                    if char == "\\":  # Insertion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              query_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (query_pos - query_gap_count)), "I"))

                    elif char == "/":  # Deletion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              query_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (query_pos - query_gap_count)), "D"))

                output_region_list.append([current_region[0], current_region[1], current_region[2],
                                           current_protein_hit[0], current_protein_hit[1],
                                           current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                           current_protein_hit[6], current_protein_hit[7],
                                           original_frameshift_position])

                # Reset the list for the next protein hit
                frameshift_list = []

        current_output_file.close()

    # Determine the frameshift positions in each query, by combining the information of all hits per query.
    # Contains the information of the current query  [scaffold, startpostion_in_scaffold, endposition_in_scaffold,
    # [frameshift position list]]
    current_query = ["#", "#", "#"]
    current_query_region_coverage = []     # Contains the start and end pos. of hit alignments, that
    # healing_region_list contains the final query information for all queries, and functions as output list
    healing_region_list = []

    for diamond_hit in output_region_list:
        # If case for start of new query information, old, now complete, query information could be appended to output
        if diamond_hit[0:3] != current_query[0:3]:
            healing_region_list.append(current_query)       # Appends the final information of the previous region
            # Initialise with new information
            current_query = [diamond_hit[0], diamond_hit[1], diamond_hit[2], copy.copy(diamond_hit[7])]
            current_query_region_coverage = [[int(diamond_hit[8]), int(diamond_hit[9])]]

        # Else case for adding frameshift information of an additional protein hit if possible
        else:
            # Initialise parameters
            region_overlapping = False      # Saves if the new region is overlapping with a already included region
            alignment_start_pos = diamond_hit[8]
            alignment_end_pos = diamond_hit[9]

            # Checks for all already added regions, if the new region is overlapping with one of them
            for prev_region in current_query_region_coverage:
                if ((alignment_end_pos >= prev_region[0]) and (alignment_start_pos <= prev_region[0])) or\
                    ((alignment_end_pos >= prev_region[1]) and (alignment_start_pos <= prev_region[1])) or\
                        ((alignment_start_pos >= prev_region[0]) and (alignment_end_pos <= prev_region[1])):
                    region_overlapping = True       # overlapping region found
                    break

            # Case were no previous region is overlapping with the current region
            if not region_overlapping:
                # Append the new region to the coverage list
                current_query_region_coverage.append([alignment_start_pos, alignment_end_pos])
                current_query[3] += diamond_hit[7]      # add the frameshifts of the new hit to final list

        # Remove the unnecessary start and end positions of the alignment in the output
        diamond_hit.pop(9)
        diamond_hit.pop(8)

    # Append last remaining region
    healing_region_list.append(current_query)
    # Remove the initialising object
    healing_region_list.pop(0)
    no_hit_length_distribution.pop(0)

    return output_region_list, healing_region_list, no_hit_length_distribution


def read_in_results_3_healed(output_dir):
    """
    SPECIAL version with a different header parsing, for counting the frameshifts of the healed queries.
    The function reads in the Diamond blastx results in a given directory. The results for each hit are stored in the
    output_region_list, which is also returned. The contains for every region: Scaffold, start pos., end pos. in
    scaffold, e-value, bitscore, similarity-percentage, list with all detected frameshifts with their positions in the
    query (I = insertion, D = deletion). The second output list (healing_region_list), contains the combined frameshift
    information of all Diamond hits per query, which means the scaffold, the start position in the scaffold, the end
    position in the scaffold, and a list with all positions in the query, where frameshifts are detected and need to
    be healed later . But not all Diamond hits of each Query are included, some are excluded by our overlapping
    heuristic, more detailed information in the documentation and in the code.
    PS: frameshift positions are the python list positions.
    :param output_dir: path to the output directory
    :return: output_region_list, healing_region_list
    """
    # Returns every file in the directory with .out at the end
    output_file_list = glob.glob(output_dir + "/temp_out_*.txt")

    # Sorting the file path ascending by their file number
    dir_path_len = len(output_dir) + 9  # +9 because temp_out_ consists of 9 chars
    output_file_list = sorted(output_file_list, key=lambda current_path: int(current_path[dir_path_len:-4]))

    output_region_list = []  # Containing the Diamond hits with additional information

    # Holds the length of the queries with no diamond hit, for later analysis
    no_hit_length_distribution = []

    # Reading out all output files
    for x in output_file_list:
        current_output_file = open(x, "r")

        # Initialising the parameters
        query_hit_list = []    # Contains all hit of the current query
        current_region = (0, 0, 0)  # Scaffold, start and end position of the current considered region
        frameshift_list = []  # Contains the positions of detected frameshifts
        query_alignment = ""
        # The following two parameters are important for the localising of the frameshift in the scaffold, and the later
        # Overlapping heuristic for decide which frameshifts are considered in each query
        query_alignment_start_pos = 0   # Saves the start position of the alignment in the query
        query_alignment_end_pos = 0     # for saving the current end position of each alignment
        subject_alignment = ""
        # The following three store the respective values
        protein_hit = ""
        e_value = ""
        bit_score = ""
        similarity_percentage = ""

        # Reads in the current output file
        for line in current_output_file:
            # Case for a new Query
            if line.startswith("Query="):

                # Case if the query_hit_list, which holds the hits of each query, is not empty,
                # which means the previous region had an Diamond hit
                if query_hit_list:

                    # For every query_hit_list per query, the first initialising element must be deleted and the last
                    # hit must be afterwards be appended

                    # Appends the last remaining protein hit of the query, because each protein hit of the query is
                    # appended to the list by its following hit,which means that the last one needs to be appended later
                    query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                           subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

                    # Reset the parameter
                    query_alignment_start_pos = 0

                    # Removes the initialising element of the query, because the first element appends the last protein
                    # hit of the previous query, which is a result of the used approach
                    query_hit_list.pop(0)

                    # Appends each protein hit of the query to the output_region_list
                    for current_protein_hit in query_hit_list:

                        # The placed gaps in the query in the alignment are counted for a correct position calculation
                        query_gap_count = 0

                        # Detect frameshifts in aligned query region
                        for query_pos, char in enumerate(current_protein_hit[4]):

                            # Counts how many gaps are inserted into the query before the frameshift position
                            if char == "-":
                                query_gap_count += 1

                            if char == "\\":  # Insertion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      query_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (query_pos - query_gap_count)), "I"))

                            elif char == "/":  # Deletion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      query_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (query_pos - query_gap_count)), "D"))

                        output_region_list.append([current_region[0], current_region[1], current_region[2],
                                                   current_protein_hit[0], current_protein_hit[1],
                                                   current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                                   current_protein_hit[6], current_protein_hit[7]])

                        # Reset values
                        frameshift_list = []

                    # Reset the the list, the if request is inactive till a new protein hit was found and thus
                    # initialising element is created
                    query_hit_list = []

                # Case if the query had no diamond hit
                else:
                    no_hit_length_distribution.append(current_region[2])

                # Setting the new region
                current_line = line.strip().split()[1]
                current_region = (current_line, "1", "1")

            # Line containing the protein id
            if line.startswith(">"):

                # Appends the previous protein hit to the list with all hits of each query
                # If This is the first hit of te current query, the last protein hit of the previous query is appended,
                # but this "initialising" element is later at the Query= if request removed. In a normal case each hit
                # is appended by the following hit of the same query, and the last hit is also appended at the Query= if
                query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                      subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

                # Reset/set the new parameters
                protein_hit = line.split()[0][1:]  # containing the protein hit of the current region
                subject_alignment = ""
                query_alignment = ""
                query_alignment_start_pos = 0

            # Line containing score and e-value
            if line.startswith(" Score"):
                current_line = line.strip().split()
                bit_score = current_line[2]
                e_value = current_line[7]

            # Line similarity percentage
            if line.startswith(" Identities"):
                current_line = line.strip().split()
                similarity_percentage = current_line[3][1:-2]

            # Reading in the aligned sequences
            if line.startswith("Query "):
                # Case that a new alignment starts, and the new start position is registered
                if not query_alignment_start_pos:
                    query_alignment_start_pos = int(line.split()[1])

                query_alignment_end_pos = int(line.split()[3])   # Append the currently last alignment position
                query_alignment += line.split()[2]

            # Subject alignment line case
            if line.startswith("Sbjct "):
                subject_alignment += line.split()[2]

        # Appending the last region if necessary (same as upper part)
        if query_hit_list:
            # Appends the last remaining protein hit of the query
            query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                   subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

            # Removes the initialising element of the query
            query_hit_list.pop(0)

            # Appends each protein hit of the query to the output_region_list
            for current_protein_hit in query_hit_list:

                # The placed gaps in the query in the alignment are counted for a correct position calculation
                query_gap_count = 0

                # Detect frameshifts in aligned query region
                for query_pos, char in enumerate(current_protein_hit[4]):

                    # Counts how many gaps are inserted into the query before the frameshift position
                    if char == "-":
                        query_gap_count += 1

                    if char == "\\":  # Insertion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              query_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (query_pos - query_gap_count)), "I"))

                    elif char == "/":  # Deletion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              query_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (query_pos - query_gap_count)), "D"))

                output_region_list.append([current_region[0], current_region[1], current_region[2],
                                           current_protein_hit[0], current_protein_hit[1],
                                           current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                           current_protein_hit[6], current_protein_hit[7]])

                # Reset the list for the next protein hit
                frameshift_list = []

        current_output_file.close()

    # Determine the frameshift positions in each query, by combining the information of all hits per query.
    # Contains the information of the current query  [scaffold, startpostion_in_scaffold, endposition_in_scaffold,
    # [frameshift position list]]
    current_query = ["#", "#", "#"]
    current_query_region_coverage = []     # Contains the start and end pos. of hit alignments, that
    # healing_region_list contains the final query information for all queries, and functions as output list
    healing_region_list = []

    for diamond_hit in output_region_list:
        # If case for start of new query information, old, now complete, query information could be appended to output
        if diamond_hit[0:3] != current_query[0:3]:
            healing_region_list.append(current_query)       # Appends the final information of the previous region
            # Initialise with new information
            current_query = [diamond_hit[0], diamond_hit[1], diamond_hit[2], copy.copy(diamond_hit[7])]
            current_query_region_coverage = [[int(diamond_hit[8]), int(diamond_hit[9])]]

        # Else case for adding frameshift information of an additional protein hit if possible
        else:
            # Initialise parameters
            region_overlapping = False      # Saves if the new region is overlapping with a already included region
            alignment_start_pos = diamond_hit[8]
            alignment_end_pos = diamond_hit[9]

            # Checks for all already added regions, if the new region is overlapping with one of them
            for prev_region in current_query_region_coverage:
                if ((alignment_end_pos >= prev_region[0]) and (alignment_start_pos <= prev_region[0])) or\
                    ((alignment_end_pos >= prev_region[1]) and (alignment_start_pos <= prev_region[1])) or\
                        ((alignment_start_pos >= prev_region[0]) and (alignment_end_pos <= prev_region[1])):
                    region_overlapping = True       # overlapping region found
                    break

            # Case were no preivious region is overlapping with the current region
            if not region_overlapping:
                # Append the new region to the coverage list
                current_query_region_coverage.append([alignment_start_pos, alignment_end_pos])
                current_query[3] += diamond_hit[7]      # add the frameshifts of the new hit to final list

        # Remove the unnecesseary start and end positions of the aligment in the output
        diamond_hit.pop(9)
        diamond_hit.pop(8)

    # Append last remaining region
    healing_region_list.append(current_query)
    # Remove the initialising object
    healing_region_list.pop(0)
    if no_hit_length_distribution:
        no_hit_length_distribution.pop(0)

    return output_region_list, healing_region_list, no_hit_length_distribution


def read_in_results_3_detect_multiple_positions(output_dir):
    """
    SPECIAL version which interacts together with the create_mod_queries_new_multiple_positions function, which means
    the original positoins of all 4 frameshift are saved in the output_region_list
    Modified version, output region list, saves additionally the original frameshift position for each region
    (additionally at the end of each query list)
    DETECTION version for benchmark, with a different header parsing, for the frameshift correction
    The function reads in the Diamond blastx results in a given directory. The results for each hit are stored in the
    output_region_list, which is also returned. The contains for every region: Scaffold, start pos., end pos. in
    scaffold, e-value, bitscore, similarity-percentage, list with all detected frameshifts with their positions in the
    query (I = insertion, D = deletion). The second output list (healing_region_list), contains the combined frameshift
    information of all Diamond hits per query, which means the scaffold, the start position in the scaffold, the end
    position in the scaffold, and a list with all positions in the query, where frameshifts are detected and need to
    be healed later . But not all Diamond hits of each Query are included, some are excluded by our overlapping
    heuristic, more detailed information in the documentation and in the code.
    PS: frameshift positions are the python list positions.
    :param output_dir: path to the output directory
    :return: output_region_list, healing_region_list
    """
    # Returns every file in the directory with .out at the end
    output_file_list = glob.glob(output_dir + "/temp_out_*.txt")

    # Sorting the file path ascending by their file number
    dir_path_len = len(output_dir) + 9  # +9 because temp_out_ consists of 9 chars
    output_file_list = sorted(output_file_list, key=lambda current_path: int(current_path[dir_path_len:-4]))

    output_region_list = []  # Containing the Diamond hits with additional information

    # Holds the length of the queries with no diamond hit, for later analysis
    no_hit_length_distribution = []

    # Saves for each Diamond hit the distance from the original frameshift position to the detected frameshift position
    original_frameshift_positions = 0    # Initialising

    # Reading out all output files
    for x in output_file_list:
        current_output_file = open(x, "r")

        # Initialising the parameters
        query_hit_list = []    # Contains all hit of the current query
        current_region = (0, 0, 0)  # Scaffold, start and end position of the current considered region
        frameshift_list = []  # Contains the positions of detected frameshifts
        query_alignment = ""
        # The following two parameters are important for the localising of the frameshift in the scaffold, and the later
        # Overlapping heuristic for decide which frameshifts are considered in each query
        query_alignment_start_pos = 0   # Saves the start position of the alignment in the query
        query_alignment_end_pos = 0     # for saving the current end position of each alignment
        subject_alignment = ""
        # The following three store the respective values
        protein_hit = ""
        e_value = ""
        bit_score = ""
        similarity_percentage = ""

        # Reads in the current output file
        for line in current_output_file:
            # Case for a new Query
            if line.startswith("Query="):

                # Case if the query_hit_list, which holds the hits of each query, is not empty,
                # which means the previous region had an Diamond hit
                if query_hit_list:

                    # For every query_hit_list per query, the first initialising element must be deleted and the last
                    # hit must be afterwards be appended

                    # Appends the last remaining protein hit of the query, because each protein hit of the query is
                    # appended to the list by its following hit,which means that the last one needs to be appended later
                    query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                           subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

                    # Reset the parameter
                    query_alignment_start_pos = 0

                    # Removes the initialising element of the query, because the first element appends the last protein
                    # hit of the previous query, which is a result of the used approach
                    query_hit_list.pop(0)

                    # Appends each protein hit of the query to the output_region_list
                    for current_protein_hit in query_hit_list:

                        # The placed gaps in the query in the alignment are counted for a correct position calculation
                        query_gap_count = 0

                        # Detect frameshifts in aligned query region
                        for query_pos, char in enumerate(current_protein_hit[4]):

                            # Counts how many gaps are inserted into the query before the frameshift position
                            if char == "-":
                                query_gap_count += 1

                            if char == "\\":  # Insertion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      query_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (query_pos - query_gap_count)), "I"))

                            elif char == "/":  # Deletion detected
                                # Function returns True if the frameshift correlates with a putative intron
                                wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                                      current_protein_hit[5],
                                                                                      query_pos)

                                if not wrong_intron:
                                    # Saves the positions as the python list positions of the first nucleotide in the
                                    # triplet with the frameshift
                                    frameshift_list.append(((current_protein_hit[6] - 1) +
                                                            (3 * (query_pos - query_gap_count)), "D"))

                        output_region_list.append([current_region[0], current_region[1], current_region[2],
                                                   current_protein_hit[0], current_protein_hit[1],
                                                   current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                                   current_protein_hit[6], current_protein_hit[7],
                                                   original_frameshift_positions])

                        # Reset values
                        frameshift_list = []

                    # Reset the the list, the if request is inactive till a new protein hit was found and thus
                    # initialising element is created
                    query_hit_list = []

                # Case if the query had no diamond hit
                else:
                    no_hit_length_distribution.append(current_region[2])

                # Setting the new region
                current_line = line.strip().split()[1].split("#")
                current_region = (current_line[0], current_line[1], current_line[2])

                """Hier eine Liste mit den originalen positionen erstellen"""
                # Saves the original position of the frameshift, that was inserted
                original_frameshift_positions = current_line[3:7]

            # Line containing the protein id
            if line.startswith(">"):

                # Appends the previous protein hit to the list with all hits of each query
                # If This is the first hit of te current query, the last protein hit of the previous query is appended,
                # but this "initialising" element is later at the Query= if request removed. In a normal case each hit
                # is appended by the following hit of the same query, and the last hit is also appended at the Query= if
                query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                      subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

                # Reset/set the new parameters
                protein_hit = line.split()[0][1:]  # containing the protein hit of the current region
                subject_alignment = ""
                query_alignment = ""
                query_alignment_start_pos = 0

            # Line containing score and e-value
            if line.startswith(" Score"):
                current_line = line.strip().split()
                bit_score = current_line[2]
                e_value = current_line[7]

            # Line similarity percentage
            if line.startswith(" Identities"):
                current_line = line.strip().split()
                similarity_percentage = current_line[3][1:-2]

            # Reading in the aligned sequences
            if line.startswith("Query "):
                # Case that a new alignment starts, and the new start position is registered
                if not query_alignment_start_pos:
                    query_alignment_start_pos = int(line.split()[1])

                query_alignment_end_pos = int(line.split()[3])   # Append the currently last alignment position
                query_alignment += line.split()[2]

            # Subject alignment line case
            if line.startswith("Sbjct "):
                subject_alignment += line.split()[2]

        # Appending the last region if necessary (same as upper part)
        if query_hit_list:
            # Appends the last remaining protein hit of the query
            query_hit_list.append([protein_hit, e_value, bit_score, similarity_percentage, query_alignment,
                                   subject_alignment, query_alignment_start_pos, query_alignment_end_pos])

            # Removes the initialising element of the query
            query_hit_list.pop(0)

            # Appends each protein hit of the query to the output_region_list
            for current_protein_hit in query_hit_list:

                # The placed gaps in the query in the alignment are counted for a correct position calculation
                query_gap_count = 0

                # Detect frameshifts in aligned query region
                for query_pos, char in enumerate(current_protein_hit[4]):

                    # Counts how many gaps are inserted into the query before the frameshift position
                    if char == "-":
                        query_gap_count += 1

                    if char == "\\":  # Insertion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              query_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (query_pos - query_gap_count)), "I"))

                    elif char == "/":  # Deletion detected
                        # Function returns True if the frameshift correlates with a putative intron
                        wrong_intron = exclude_putative_transition_frameshift(current_protein_hit[4],
                                                                              current_protein_hit[5],
                                                                              query_pos)

                        if not wrong_intron:
                            # Saves the positions as the python list positions of the first nucleotide in the
                            # triplet with the frameshift
                            frameshift_list.append(((current_protein_hit[6] - 1) +
                                                    (3 * (query_pos - query_gap_count)), "D"))

                output_region_list.append([current_region[0], current_region[1], current_region[2],
                                           current_protein_hit[0], current_protein_hit[1],
                                           current_protein_hit[2], current_protein_hit[3], frameshift_list,
                                           current_protein_hit[6], current_protein_hit[7],
                                           original_frameshift_positions])

                # Reset the list for the next protein hit
                frameshift_list = []

        current_output_file.close()

    # Determine the frameshift positions in each query, by combining the information of all hits per query.
    # Contains the information of the current query  [scaffold, startpostion_in_scaffold, endposition_in_scaffold,
    # [frameshift position list]]
    current_query = ["#", "#", "#"]
    current_query_region_coverage = []     # Contains the start and end pos. of hit alignments, that
    # healing_region_list contains the final query information for all queries, and functions as output list
    healing_region_list = []

    for diamond_hit in output_region_list:
        # If case for start of new query information, old, now complete, query information could be appended to output
        if diamond_hit[0:3] != current_query[0:3]:
            healing_region_list.append(current_query)       # Appends the final information of the previous region
            # Initialise with new information
            current_query = [diamond_hit[0], diamond_hit[1], diamond_hit[2], copy.copy(diamond_hit[7])]
            current_query_region_coverage = [[int(diamond_hit[8]), int(diamond_hit[9])]]

        # Else case for adding frameshift information of an additional protein hit if possible
        else:
            # Initialise parameters
            region_overlapping = False      # Saves if the new region is overlapping with a already included region
            alignment_start_pos = diamond_hit[8]
            alignment_end_pos = diamond_hit[9]

            # Checks for all already added regions, if the new region is overlapping with one of them
            for prev_region in current_query_region_coverage:
                if ((alignment_end_pos >= prev_region[0]) and (alignment_start_pos <= prev_region[0])) or\
                    ((alignment_end_pos >= prev_region[1]) and (alignment_start_pos <= prev_region[1])) or\
                        ((alignment_start_pos >= prev_region[0]) and (alignment_end_pos <= prev_region[1])):
                    region_overlapping = True       # overlapping region found
                    break

            # Case were no preivious region is overlapping with the current region
            if not region_overlapping:
                # Append the new region to the coverage list
                current_query_region_coverage.append([alignment_start_pos, alignment_end_pos])
                current_query[3] += diamond_hit[7]      # add the frameshifts of the new hit to final list

        # Remove the unnecesseary start and end positions of the aligment in the output
        diamond_hit.pop(9)
        diamond_hit.pop(8)

    # Append last remaining region
    healing_region_list.append(current_query)
    # Remove the initialising object
    healing_region_list.pop(0)
    no_hit_length_distribution.pop(0)

    return output_region_list, healing_region_list, no_hit_length_distribution


def heal_assembly_file_mod(healing_region_list, input_fna_path):
    """
    MODIFIED version, the output file contains only sequences, which are healed by the function.
    THIS WORKS ONLY IN THIS SPECIAL BENCHMARK CASE
    Gets the healing_region_list which contains the frameshift positions in each query, and inserts N's at those
    positions, to heal the reading frame. The healed assembly is afterwards saved in the same directory as the given
    original assembly .fna file.
    The new assembly file doesnt contain comments.
    :param healing_region_list: Contains for each query: scaffold, start pos., end pos. in scaff., frameshift pos. list
    :param input_fna_path:  File path to the original assembly file
    :return: File path to the new modified assembly file
    """
    # Read in the .fna file for fast access
    input_fna_file = open(input_fna_path)
    scaffold_list = []  # Filled with the sequences of the fna file

    healed_scaffolds = []   # List which contains only the healed sequences

    temp_scaffold = []  # Contains the current scaffold
    current_header = ""  # Containing the header of the current scaffold
    correction_count = 0    # Count how many frameshifts are healed

    # Filling the .fna region list/ reading in the original assembly
    for line in input_fna_file:

        # Sequence line case
        if (line[0] != ">") and (line[0] != ";"):
            for x in line.strip():  # filling up a new sequence
                temp_scaffold.append(x)

        # Sequence header case
        elif line[0] == ">":
            scaffold_list.append([current_header, temp_scaffold])  # appending the previous region
            current_header = line.strip().split("#")[0][1:]  # reading in the new header
            temp_scaffold = []  # resetting the list/ ready for new region

    # Append last remaining region
    scaffold_list.append([current_header, temp_scaffold])
    input_fna_file.close()
    # Delete the first empty initialising element
    scaffold_list.pop(0)

    # Correct the frameshifts
    # Inserting the N's from the end to the beginning, to prevent the necessity of adapting the insertion position in
    # dependency to the previous inserted frameshifts
    # For this, in this first sort, the queries in each scaffold are sorted descending by their start position
    # but the reverse option also sorts the scaffolds descending
    sorted_healing_region_list = sorted(healing_region_list, reverse=True,  key=lambda temp_query: (temp_query[0],
                                                                                                    temp_query[1]))

    for query in sorted_healing_region_list:
        # Search the corresponding scaffold
        for scaffold in scaffold_list:
            # Case where matching scaffold is found
            if scaffold[0] == query[0]:
                query_start_pos = int(query[1])
                # Second sort step, now the frameshift positions in each query are sorted descending
                temp_sorted_query = sorted(query[3], reverse=True,  key=lambda temp_query: temp_query[0])

                for frameshift_position in temp_sorted_query:      # Same usage of reversed() as before
                    # Deletion case
                    if frameshift_position[1] == "D":
                        scaffold[1] = scaffold[1][:(query_start_pos + frameshift_position[0])] + ["N"] + \
                                      scaffold[1][(query_start_pos + frameshift_position[0]):]
                        correction_count += 1
                    # Insertion case
                    else:
                        scaffold[1] = scaffold[1][:(query_start_pos + frameshift_position[0])] + ["N"] + ["N"] + \
                                      scaffold[1][(query_start_pos + frameshift_position[0]):]
                        correction_count += 1

                # Appends the sequence to the output list, a frameshift was healed
                if temp_sorted_query:
                    healed_scaffolds.append(scaffold)

                break

    # Create the new assembly .fna file path, located in the same dir as the original assembly file
    new_fna_file_path = input_fna_path[:(input_fna_path.rfind("/") + 1)] + "healed_queries_.fasta"

    # Creating the new assembly file (with the inserted N's)
    new_fna_file = open(new_fna_file_path, "w")
    for scaffold in healed_scaffolds:

        # Write the header line for the scaffold
        new_fna_file.write(">" + scaffold[0] + "\n")

        # Normal case where the nucleotide sequence is longer than 80
        if len(scaffold[1]) > 80:

            # Create lines with a length of 80, for a more readable .fasta file
            # Create list with the end positions of each 80 nucleotide long line
            position_list = list(range(80, len(scaffold[1]), 80))
            # Writes the lines in the assembly file
            for list_pos in position_list:
                # The nucleotides are stored as chars in a list and are converted to a string with join
                new_fna_file.write(("".join(scaffold[1][(list_pos - 80):list_pos])) + "\n")

            # Write the last remaining line (range function doesnt include the last elements)
            new_fna_file.write(("".join(scaffold[1][position_list[-1]:len(scaffold[1])])) + "\n")

        # Case where the nucleotide sequence is smaller than 80
        else:
            new_fna_file.write(("".join(scaffold[1]) + "\n"))

    new_fna_file.close()

    print(len(healed_scaffolds), " queries are healed")
    print(correction_count, " frameshifts are healed")

    return new_fna_file_path


def evaluate_result_data(healing_region_list):
    """
    Trivial function, which prints how many frameshifts are corrected.
    :param healing_region_list: output list of the read_n_results_3 function
    :return: None
    """
    detecteted_frameshift_count = 0
    multiple_frameshift_count = 0
    without_frameshift_count = 0
    curable_queries_count = 0

    for query in healing_region_list:
        frameshifts_per_query = len(query[3])
        detecteted_frameshift_count += frameshifts_per_query
        if frameshifts_per_query > 1:
            multiple_frameshift_count += 1

        if frameshifts_per_query == 0:
            without_frameshift_count += 1

        if frameshifts_per_query >= 1:
            curable_queries_count += 1

    print("####################")
    print("Detected frameshift: ", detecteted_frameshift_count)
    print("####################")
    print("Curable queries: ", curable_queries_count)
    print("####################")
    # print("Queries with more than one detected frameshift: ", multiple_frameshift_count)
    # print("####################")
    # print("Diamond hit without a detected frameshift: ", without_frameshift_count)
    # print("####################")
    print("Queries with a Diamond hit: ", len(healing_region_list))
    print("####################")

    return None


def create_boxplot_data(output_region_list, exclude_frameshift_hits):
    """
    Simply saves the similarity percentage and the bit score of each diamond hit in a new file, which is created at the
    cwd.
    :param output_region_list: list with the output of the read_in_results_3 function
    :param exclude_frameshift_hits: if True, diamond hits with a frameshift are excluded
    :return: saves the results in a file
    """

    new_file_path = os.getcwd() + "/boxplot_data.txt"
    boxplot_data_file = open(new_file_path, "w")

    # Excludes diamond hits with a detected frameshift
    if exclude_frameshift_hits:
        for diamond_hit in output_region_list:
            similarity_percentage = diamond_hit[6]
            bit_score = diamond_hit[5]

            # Excludes diamond hits with a detected frameshift
            if not diamond_hit[7]:
                boxplot_data_file.write(similarity_percentage + "\t" + bit_score + "\n")

    else:
        for diamond_hit in output_region_list:
            similarity_percentage = diamond_hit[6]
            bit_score = diamond_hit[5]

            boxplot_data_file.write(similarity_percentage + "\t" + bit_score + "\n")

    boxplot_data_file.close()

    return None


def create_combined_boxplot(input_file_list, element_labels, plot_title):
    """
    Receives a list with the path to all boxplot-data files of each run, which should be displayed in a boxplot.
    The boxplots are displayed in the combined boxplot in the order where the paths are sorted in the given list.
    :param input_file_list: list with the paths to each boxplot-data file
    :param element_labels: labels for each element of the boxplot
    :param plot_title: title of the created plot
    :return: None, output image is created at the cwd
    """

    # Contain the data of all input files, used for the boxplot creation
    identity_data = []

    # Reads in all input files
    for input_file in input_file_list:
        current_file = open(input_file)
        # Lists where the data of the current file is stored
        temp_identity_file_data = []     # Stores the identitys of the current file

        # Reads in the content of the current file
        for line in current_file:
            temp_identity_file_data.append(float(line.strip().split()[0][:-1]))

        current_file.close()

        # Append the data of the current file
        identity_data.append(temp_identity_file_data)

    # Create the combined boxplot at the cwd
    fig1, ax1 = plt.subplots()
    ax1.set_title(plot_title, fontsize=13, fontweight="bold")

    # To create boxes with different colours, different subgroups of boxes are created (here the following three)
    # The flanking dots outer the borders of each boxplot are turned off, the positions-list set the positions of each
    # subdataset, in the combined graph, patch_artist is necessary to change the colour of the boxplots
    boxes_1 = ax1.boxplot(identity_data[0::3], showfliers=False, positions=[1, 4, 7, 10], patch_artist=True)

    # Change the colour of all subparts of the boxplot
    for element in ["boxes", "whiskers", "fliers", "means", "medians", "caps"]:
        plt.setp(boxes_1[element], color="black")

    # Use a different coulor for the inside of the boxplot
    for patch in boxes_1["boxes"]:
        patch.set(facecolor="darkgrey")

    plt.setp(boxes_1["medians"], color="gold")

    # Second boxplot colour group
    boxes_2 = ax1.boxplot(identity_data[1::3], showfliers=False, positions=[2, 5, 8, 11], patch_artist=True)
    for element in ["boxes", "whiskers", "fliers", "means", "medians", "caps"]:
        plt.setp(boxes_2[element], color="black")

    for patch in boxes_2["boxes"]:
        patch.set(facecolor="tomato")

    plt.setp(boxes_2["medians"], color="gold")

    # Third boxplot colour group
    boxes_3 = ax1.boxplot(identity_data[2::3], showfliers=False, positions=[3, 6, 9, 12], patch_artist=True)
    for element in ["boxes", "whiskers", "fliers", "means", "medians", "caps"]:
        plt.setp(boxes_3[element], color="black")

    # Use a different coulor for the inside of the boxplot
    for patch in boxes_3["boxes"]:
        patch.set(facecolor="limegreen")

    plt.setp(boxes_3["medians"], color="gold")

    # Sets the positions where the labels ar placed on the x axis
    name_positions = [2, 5, 8, 11]
    # Places the labels on the x axis
    plt.xticks(name_positions, element_labels)
    # Adds a legend with the colour information
    ax1.legend([boxes_1["boxes"][0], boxes_2["boxes"][0], boxes_3["boxes"][0]], ["unmodified", "modified", "healed"])
    # Set the range of the y axis
    ax1.set_ylim(ymin=90.0, ymax=105.0)
    # Set the marks/ticks of the y axis
    ax1.set_yticks(ticks=[90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100])
    # Sets the labels of the two axis
    ax1.set_ylabel("Percent Identity", fontsize=11.5)
    ax1.set_xlabel("Query Length", fontsize=11.5)
    # Saves the created figure as a .png file
    plt.savefig("identity_boxplot.png", dpi=300)

    return None


def create_combined_boxplot_position_distribution(input_file_list, element_labels, plot_title):
    """
    Receives a list with the path to all boxplot-data files of each run, which should be displayed in a boxplot.
    The boxplots are displayed in the combined boxplot in the order where the paths are sorted in the given list.
    :param input_file_list: list with the paths to each boxplot-data file
    :param element_labels: labels for each element of the boxplot
    :param plot_title: title of the created plot
    :return: None, output image is created at the cwd
    """

    # Contain the data of all input files, used for the boxplot creation
    identity_data = []

    # Reads in all input files
    for input_file in input_file_list:
        current_file = open(input_file)
        # Lists where the data of the current file is stored
        temp_identity_file_data = []     # Stores the identities of the current file

        # Reads in the content of the current file
        for line in current_file:
            temp_identity_file_data.append(float(line.strip().split()[0]))

        current_file.close()

        # Append the data of the current file
        identity_data.append(temp_identity_file_data)

    # Create the combined boxplot at the cwd
    fig1, ax1 = plt.subplots()
    ax1.set_title(plot_title, fontsize=13, fontweight="bold")

    # To create boxes with different colours, different subgroups of boxes are created (here the following three)
    # The flanking dots outer the borders of each boxplot are turned off, the positions-list set the positions of each
    # subdataset, in the combined graph, patch_artist is necessary to change the colour of the boxplots
    boxes_1 = ax1.boxplot(identity_data, showfliers=False, positions=[1, 2, 3, 4], patch_artist=True)

    # Change the colour of all subparts of the boxplot
    for element in ["boxes", "whiskers", "fliers", "means", "medians", "caps"]:
        plt.setp(boxes_1[element], color="black")

    # Use a different coulor for the inside of the boxplot
    for patch in boxes_1["boxes"]:
        patch.set(facecolor="darkgrey")

    # Sets the positions where the labels ar placed on the x axis
    name_positions = [1, 2, 3, 4]
    # Places the labels on the x axis
    plt.xticks(name_positions, element_labels)
    # Set the range of the y axis
    #ax1.set_ylim(ymin=60, ymax=115)
    # Set the marks/ticks of the y axis
    #ax1.set_yticks(ticks=[60, 70, 80, 90, 100])
    # Sets the labels of the two axis
    ax1.set_ylabel("Distance to correct position", fontsize=11.5)
    ax1.set_xlabel("Query Length", fontsize=11.5)
    # Saves the created figure as a .png file
    plt.savefig("distance_boxplot.png", dpi=1000)

    return None


def create_bar_graph(data_file_path_list, element_cutoff):
    """

    :param data_file_path_list:
    :param element_cutoff:
    :return:
    """

    # Contain the data of all input files, used for the boxplot creation
    distance_data = []
    sum_for_average = 0

    # Reads in all input files
    for input_file in data_file_path_list:
        current_file = open(input_file)

        # Reads in the content of the current file
        for line in current_file:
            current_distance = abs(int(line.strip()))
            distance_data.append(current_distance)
            sum_for_average += current_distance

        current_file.close()

    print("The average distance is: ", sum_for_average/len(distance_data))
    # Count the amount of identical elements (Here queries with the same length)
    length_list = []

    # Appends the max query length
    length_list.append([element_cutoff, 0])

    #
    for elem_len in distance_data:
        new_length = True
        elem_len = int(elem_len)

        if elem_len > element_cutoff:
            elem_len = element_cutoff

        for already_added_length in length_list:
            if elem_len == already_added_length[0]:
                already_added_length[1] += 1
                new_length = False
                break

        if new_length:
            length_list.append([elem_len, 1])

    # Sort the element list ascending by their element values
    length_list = sorted(length_list, key=lambda element: int(element[0]))

    # X-coordinates of left sides of bars
    elements = [x for x in range(1, len(length_list)+1)]

    # Heights of bars
    element_counts = [y[1] for y in length_list]

    # Labels for bars
    #tick_label = [str(z[0]) for z in length_list]
    tick_label = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '>15']


    # Plotting a bar chart
    plt.bar(elements, element_counts, tick_label=tick_label, width=0.8, color=["limegreen"], edgecolor="black")

    #plt.setp(fontsize=10, rotation='vertical')

    # Label the picture
    plt.ylabel("Number of queries", fontsize=11.5)
    plt.xlabel("Distance in bp", fontsize=11.5)
    plt.title("Detected to original frameshift position distance", fontsize=13, fontweight="bold")

    # Function saves the plot
    plt.savefig("distance_bar_graph.png")

    return None


def result_analysis_1(healed_output_region_list, mod_output_region_list):
    """
    Checks how many heald sequences with a frameshift have another Diamond hit than in the modified sequence.
    :param healed_output_region_list:
    :param mod_output_region_list:
    :return:
    """

    healed_queries_with_frameshifts = []

    # Search the healed queries with a detected frameshift
    for healed_query in healed_output_region_list:
        if healed_query[7]:    # Case that frameshift-list is not empty
            healed_queries_with_frameshifts.append([healed_query[0], healed_query[3]])

    count_different_hits = 0

    for healed_query in healed_queries_with_frameshifts:
        for mod_query in mod_output_region_list:
            if healed_query[0] == mod_query[0]:
                if random.randint(0, 250) == 250:
                    print("##############")
                    print(healed_query[0])
                    print(healed_query[1])
                    print(mod_query[3])
                """if healed_query[1] != mod_query[3]:
                    print("##############")
                    print(healed_query[0])
                    print(healed_query[1])
                    print( mod_query[3])
                    count_different_hits += 1"""
                break

    print("Queries with different hits: ", count_different_hits)

    return None


def result_analysis_2(healed_output_file, mod_output_file, searched_query):
    """
    Printing the Diamond hits of the given query in the two big files
    :param healed_output_file:
    :param mod_output_file:
    :param searched_query:
    :return:
    """

    mod_file = open(mod_output_file)

    count_1 = 0
    count_increase = False

    for line in mod_file:
        if line.startswith("Query= " + searched_query):
            print("####################################################")
            count_increase = True

        if count_increase:
            print(line)
            count_1 += 1

        if count_1 == 200:
            break

    mod_file.close()

    healed_file = open(healed_output_file)

    count_2 = 0
    count_increase_2 = False

    for line in healed_file:
        if line.startswith("Query= " + searched_query):
            print("####################################################")
            count_increase_2 = True

        if count_increase_2:
            print(line)
            count_2 += 1

        if count_2 == 200:
            break

    healed_file.close()


def calculate_original_to_detected_frameshift_position_distance(output_region_list):
    """
    Calculates for each query the distance between the original position, where the new frameshift was inserted by me to
    modify the query, to the position where the frameshift was detected by Diamond using the database.
    Queries with more than one detected frameshift are excluded from the analysis.
    Afterwards all distances are stored in a boxplot-data-file, for a combined boxplot later.
    :param output_region_list:
    :return:
    """

    position_distance_list = []

    for query in output_region_list:

        # Only consider queries with exactly one frameshift
        if len(query[7]) == 1:
            position_distance_list.append(int(query[8]) - query[7][0][0])

    new_file_path = os.getcwd() + "/position_distance_boxplot_data.txt"
    boxplot_data_file = open(new_file_path, "w")

    for distance in position_distance_list:
        boxplot_data_file.write(str(distance) + "\n")

    boxplot_data_file.close()


def calculate_original_to_detected_frameshift_position_distance_multiple_positions(output_region_list):
    """
    SPECIAL version which interacts together with the create_mod_queries_new_multiple_positions function, which means
    for each detected frameshift the distance to the closest original position is calculated, so not always the correct
    distance is calculated.
    Queries with more than four detected frameshift are excluded from the analysis.
    Afterwards all distances are stored in a boxplot-data-file, for a combined boxplot later.
    :param output_region_list: output of create_mod_queries_new_multiple_positions()
    :return: None, but a file at the cwd, with the calculated data
    """

    position_distance_list = []

    for query in output_region_list:

        # Only consider queries with exactly one frameshift
        if len(query[7]) <= 4:

            # Extract the detected frameshift positions
            detected_frameshift_positions = []
            for detected_position in query[7]:
                detected_frameshift_positions.append(detected_position[0])

            # Calculate the min distance for each position

            for detected_position in detected_frameshift_positions:
                currently_min_distance = float("inf")
                for original_position in query[8]:
                    if abs(detected_position - int(original_position)) < abs(currently_min_distance):
                        currently_min_distance = detected_position - int(original_position)

                position_distance_list.append(int(currently_min_distance))

    new_file_path = os.getcwd() + "/position_distance_boxplot_data.txt"
    boxplot_data_file = open(new_file_path, "w")

    for distance in position_distance_list:
        boxplot_data_file.write(str(distance) + "\n")

    boxplot_data_file.close()


def cds_length_distribution(cds_file_path):
    """
    Creates a Boxplot with the length distribution of all unmodified complete CDS.
    :param cds_file_path: Path to the cds file
    :return: None, but creates a boxplot as .png file in the cwd
    """

    input_cds_file = open(cds_file_path)

    cds_region_list = []  # Stores each CDS region header of the input file, together with the corresponding region
    current_cds_header = ""  # Stores the header of the current CDS region
    current_cds_seq = ""  # Stores the sequence of the current CDS region

    # Read in the input with the CDS regions
    for line in input_cds_file:
        if line[0] != "#":  # Excludes comments
            if line[0] == ">":  # Header of new region found
                # Append the previous region to the final list
                cds_region_list.append([current_cds_header, current_cds_seq])
                current_cds_seq = ""
                current_cds_header = line.strip()  # Set the header of the now following sequence

            else:
                # Creates the complete sequence by combining the subsequences in each line together
                current_cds_seq += line.strip()

    # Append the last remaining region
    cds_region_list.append([current_cds_header, current_cds_seq])

    cds_region_list.pop(0)  # Removes the initialising element

    # Stores the length of ech cds
    cds_length_distribution_list = []

    # Appends the length of each cds to the list
    for cds in cds_region_list:
        cds_length_distribution_list.append(len(cds[1]))

    # Create the boxplot
    fig1, ax1 = plt.subplots()
    ax1.set_title("CDS length distribution", fontsize=13, fontweight="bold")

    boxes_1 = ax1.boxplot(cds_length_distribution_list, showfliers=True, positions=[1], patch_artist=True)

    # Change the colour of all subparts of the boxplot
    for element in ["boxes", "whiskers", "fliers", "means", "medians", "caps"]:
        plt.setp(boxes_1[element], color="black")

    # Use a different coulor for the inside of the boxplot
    for patch in boxes_1["boxes"]:
        patch.set(facecolor="darkgrey")

    ax1.set_ylabel("CDS length in bp", fontsize=11.5)

    # Saves the created figure as a .png file
    plt.savefig("cds_length_distribution.png", dpi=1000)

    return None


def main():
    print("Benchmark Two main executed")

    # ########## Query creation ##########
    """create_mod_queries_new_multiple_positions("/home/johannes/Desktop/benchmark_two/reference_complete_cds_queries"
                           "/GCF_000004255.2_v.1.0_cds_from_genomic.fna", 150, 450)
    # """

    # ########## Cluster run of mod- and unmod- queries ##########

    # ########## Read in results of unmod queries ##########
    """dir_with_output = "/home/johannes/Desktop/benchmark_two/run_04/unmod_queries_150-450_material/"

    output_region_list, healing_region_list, no_hit_length_distribution = read_in_results_3_detect(dir_with_output)
    evaluate_result_data(healing_region_list)

    create_boxplot_data(output_region_list, False)
    # create_output_plots(no_hit_length_distribution, "unmod_050-350", 100, "Query length", "counted Queries")"""

    # ########## Read in results and heal the frameshifts ##########
    """dir_with_output = "/home/johannes/Desktop/benchmark_two/run_04/mod_queries_150-450_material/"

    output_region_list, healing_region_list, no_hit_length_distribution = read_in_results_3_detect_multiple_positions(dir_with_output)
    evaluate_result_data(healing_region_list)

    create_boxplot_data(output_region_list, False)
    
    calculate_original_to_detected_frameshift_position_distance_multiple_positions(output_region_list)

    heal_path = "/home/johannes/Desktop/benchmark_two/run_04/mod_queries_150-450_material/mod_queries_150-450.fasta"

    heal_assembly_file_mod(healing_region_list, heal_path)

    # create_output_plots(no_hit_length_distribution, "mod_050-350", 70, "Query length", "counted Queries")"""

    # ########## Healed sequences cluster run ##########

    # ########## Count frameshifts in healed sequences ##########
    """dir_with_output = "/home/johannes/Desktop/benchmark_two/run_04/healed_queries_150-450_material/"

    output_region_list, healing_region_list, no_hit_length_distribution = read_in_results_3_healed(dir_with_output)
    evaluate_result_data(healing_region_list)

    create_boxplot_data(output_region_list, True)
    # create_output_plots(no_hit_length_distribution, "mod_050-350", 70, "Query length", "counted Queries")"""

    # ########## Create the combined boxplot ##########

    boxplot_file_list = ["/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_unmod_050-350.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_mod_050-350.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_healed_050-350.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_unmod_150-450.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_mod_150-450.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_healed_150-450.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_unmod_350-650.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_mod_350-650.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_healed_350-650.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_unmod_999999-9999999.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_mod_999999-9999999.txt",
                         "/home/johannes/Desktop/benchmark_two/run_02/boxplot_data_identity/boxplot_data_healed_999999-9999999.txt"]

    labels = ["050-350", "150-450", "350-650", "FULL-FULL"]

    plotname = "Query to target sequence similarity distribution"

    create_combined_boxplot(boxplot_file_list, labels, plotname)
    
    #"""

    """boxplot_file_list = ["/home/johannes/Desktop/benchmark_two/run_04/"
                         "boxplot_data_frameshift_position_distance/position_distance_boxplot_data_050-350.txt",
                         "/home/johannes/Desktop/benchmark_two/run_04/"
                         "boxplot_data_frameshift_position_distance/position_distance_boxplot_data_150-450.txt",
                         "/home/johannes/Desktop/benchmark_two/run_04/"
                         "boxplot_data_frameshift_position_distance/position_distance_boxplot_data_350-650.txt",
                         "/home/johannes/Desktop/benchmark_two/run_04/"
                         "boxplot_data_frameshift_position_distance/position_distance_boxplot_data_999999-9999999.txt"]

    labels = ["050-350", "150-450", "350-650", "FULL-FULL"]

    plotname = "Distance to correct frameshift position distribution"

    #create_combined_boxplot_position_distribution(boxplot_file_list, labels, plotname)
    create_bar_graph(boxplot_file_list, 16)

    #"""


if __name__ == '__main__':
    main()
