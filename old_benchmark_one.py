#!python3

"""First benchmark for the CLCR program (determine average cutoff distance)"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import datetime
import os
import random
import glob
import time
import matplotlib.pyplot as plt


def exclude_proteins_with_j(input_file_path):
    """
    creates a new protein aminoacid fasta file, without proteins containing J (Leucine or Isoleucine), because EXONERATE
    cant handle that at the database creation
    :param input_file_path: file path of the original unmodified protein file
    :return: None, the new file is created at the cwd of the program
    """

    # Initialisiation
    input_file = open(input_file_path)
    new_lines = []
    current_protein = []
    append_bool = False

    for line in input_file:

        if line[0] != ">":
            for y in line:
                if y == "J":
                    append_bool = False

            if append_bool:
                current_protein.append(line)

        else:
            if append_bool:
                new_lines += current_protein

            current_protein = []
            append_bool = True
            current_protein.append(line)

    # Appending the last remaining region if possible
    if append_bool:
        new_lines += current_protein

    input_file.close()

    # Creating new file without the proteins which are containing "J" (for Leucine or Isoleucine)

    new_file_path = os.getcwd() + "/" + input_file_path.split("/")[-1]
    print("new file at: ", new_file_path)

    new_file = open(new_file_path, "w")

    for new_line in new_lines:
        new_file.write(new_line)

    new_file.close()

    return None


def create_database():
    """
    Executes the command for creating the diamond database
    :return: None
    """
    diamond_command = "diamond makedb --in /share/project/johannes/bachelor/exonerate_vs_diamond_benchmark/NCBI_files" \
                      "/GCF_000001735.4_TAIR10.1_protein.faa -d /share/project/johannes/bachelor/exonerate_vs_diamond" \
                      "_benchmark/diamond_database/diamond_db"
    os.system(diamond_command)

    return None


def create_queries(cds_seq_file_path, query_quantity, max_cds_len):
    """
    reads in the cds regions out of the cds file, afterwards shuffeling the region list and returning as many regions,
    as set with the query_quantity threshold, not more than regions than given
    :param cds_seq_file_path: file path to the cds file
    :param max_cds_len: threshold, all CDS region with a lenght over the parameter are excluded, used at the
    complete_prot_* function cases to push down the runntime
    :param query_quantity: int to set how many of the given regions should be returned
    :return: list containing the randomly choosen regions as sublists consiting of the header string and a list with all
             bases(easily converted to string with "".join(list[X:Y]))
    """

    # Reading in the queries out of the cds file
    cds_regions = []
    cds_file = open(cds_seq_file_path)
    region_header = ""
    current_cds = ["#####"]

    print("--> read in cds file")
    for line in cds_file:
        if line[0] != "#":      # To exclude comments
            if line[0] != ">":
                for base in line.strip():       # Filling up a new sequence
                    current_cds.append(base)

            else:
                if len(current_cds) <= max_cds_len:
                    cds_regions.append([region_header, current_cds])
                # Reset list/ append new first argument
                current_cds = []
                region_header = line.strip()

    # Appending the last remaining region if possible
    if len(current_cds) <= max_cds_len:
        cds_regions.append([region_header, current_cds])

    cds_file.close()
    # Remove the first initialising region
    cds_regions.pop(0)

    print("--> create cds queries")
    # Choosing random regions
    random.seed()       # Initialising random number generator with system time
    random.shuffle(cds_regions)     # Shuffling the region list
    cds_regions_output = []     # Containing the exon regions which will be modified

    # Case if the input var is to big
    if query_quantity > len(cds_regions):
        query_quantity = len(cds_regions)

    # Appending the choosen region to the output list
    for i in range(query_quantity):
        cds_regions_output.append(cds_regions[i])

    return cds_regions_output


def database_comp_exonerate_del(cds_list, protein_file_path):
    """
    (create frameshift with deletion version)
    This functions determines the first nucleotide position in a query, where a inserted frameshift is detected by
    Exonerate instead of "cutting" the alignment at this position. For this, CDS regions are used, and step by step
    frameshifts are inserted into them, starting with a deletion of the start nucleotide and stopping with a deletion
    of the first nucleotide, which is marked by exonerate as a frameshift.
    :param cds_list: list containing the unmodified cds regions, which are later used as modified queries
    :param protein_file_path:   path to the protein file containing the amino acid sequences
    :return: None, but an output file containing the results named exonerate_output is created
    """

    # Reading in the protein sequences
    # Reading in the queries out of the cds file
    protein_sequences = []
    protein_file = open(protein_file_path)
    protein_header = ""
    current_protein_seq = ["#####"]

    print("--> read in protein file")
    for line in protein_file:
        if line[0] != "#":      # To exclude comments
            if line[0] != ">":
                for amin_acid in line.strip():      # Filling up a new sequence
                    current_protein_seq.append(amin_acid)

            else:
                protein_sequences.append([protein_header, current_protein_seq])
                # Reset list/ append new first argument
                current_protein_seq = []
                protein_header = line.strip()

    # Appending the last remaining region if possible
    protein_sequences.append([protein_header, current_protein_seq])

    protein_file.close()
    # Remove the first initialising region
    protein_sequences.pop(0)

    # Searching the fitting protein for each CDS
    # Consisting of lists with the sublist with the protein data, and a matching sublist with the CDS data
    combined_cds_protein_list = []      # [[cds, protein], ...]

    print("--> search matching CDS-protein pairs")
    # Searches the matching pairs and appends them into the combined list
    for cds_region in cds_list:

        cds_id = cds_region[0].split("protein_id=")[1].split("]")[0]
        for protein in protein_sequences:
            protein_id = protein[0].split()[0][1:]
            if cds_id == protein_id:
                combined_cds_protein_list.append([cds_region, protein])
                break

    # Initialising
    list_length = len(combined_cds_protein_list)
    process_count = 0
    frameshift_detection_threshold_list = []    # Saves at what position a frameshift was detected in each query

    print("--> start with exonerate runs")
    # Run to detect the detection threshold at the beginning
    for region_pair in combined_cds_protein_list:

        # bad programmed process messages
        if process_count == int(list_length * 0.1):
            print("--> 10% finished")
        elif process_count == int(list_length * 0.25):
            print("--> 25% finished")
        elif process_count == int(list_length * 0.5):
            print("--> 50% finished")
        elif process_count == int(list_length * 0.75):
            print("--> 75% finished")

        for position in range(40):      # 40 as biggest position, to prevent unnecessary runs

            # Create current CDS query
            # Insert the frameshift on the current position as a deletion
            mod_current_cds = ("".join(region_pair[0][1][:position])) + ("".join(region_pair[0][1][position+1:]))

            # Create CDS query file, old data will be overwritten by the "w" parameter
            new_cds_query = open("cds_query.fasta", "w")
            new_cds_query.write(region_pair[0][0] + "\n")       # Appending the CDS header
            new_cds_query.write(mod_current_cds + "\n")         # Appending the CDS nucleotide sequence
            new_cds_query.close()

            # Create current protein subject
            mod_current_prot = region_pair[1][1]
            mod_current_prot = "".join(mod_current_prot)

            # Create protein query file, old data will be overwritten by the "w" parameter
            new_prot_sbjct = open("prot_sbjct.fasta", "w")
            new_prot_sbjct.write(region_pair[1][0] + "\n")
            new_prot_sbjct.write(mod_current_prot + "\n")
            new_prot_sbjct.close()

            # Run exonerate with the current data, the python script waits till the shell command is finished
            os.system("exonerate -m protein2dna --showvulgar true -Q protein -T dna --showalignment false --verbose 0 "
                      "-q prot_sbjct.fasta -t cds_query.fasta >exonerate_temp_output.txt")

            # Reset the parameter
            frameshift_detected = False

            # Read in the exonerate output, file online consisting of one vulgar-format line
            new_exonerate_output = open("exonerate_temp_output.txt")
            for line in new_exonerate_output:
                for char in line:
                    if char == "F":     # vulgar format for frameshift detected
                        frameshift_detected = True

            # Breaks the for loop if a frameshift is detected
            if frameshift_detected:
                break

        # Appends the current position, where a frameshift was detected, if no frameshift was detected and the loop was
        # Completely run, also the current last loop run position will be added (like upper threshold)
        frameshift_detection_threshold_list.append(position)

        process_count += 1      # Increase the progress count

    # Saving the results in a file and calculating the mean
    result_sum = 0
    output_file = open("exonerate_output.txt", "w")

    for x in frameshift_detection_threshold_list:
        result_sum += x
        output_file.write(str(x) + "\n")

    output_file.close()

    print("-->  The mean threshold for exonerate frameshift detection is: ",
          result_sum/len(frameshift_detection_threshold_list))

    return None


def database_comp_exonerate_ins(cds_list, protein_file_path):
    """
    (create frameshift with insertion version)
    This functions determines the first nucleotide position in a query, where a inserted frameshift is detected by
    Exonerate instead of "cutting" the alignment at this position. For this, CDS regions are used, and step by step
    frameshifts are inserted into them, starting with a deletion of the start nucleotide and stopping with a insertion
    of the first nucleotide, which is marked by exonerate as a frameshift.
    :param cds_list: list containing the unmodified cds regions, which are later used as modified queries
    :param protein_file_path:   path to the protein file containing the amino acid sequences
    :return: None, but an output file containing the results named exonerate_output is created
    """

    # Reading in the protein sequences
    # Reading in the queries out of the cds file
    protein_sequences = []
    protein_file = open(protein_file_path)
    protein_header = ""
    current_protein_seq = ["#####"]

    print("--> read in protein file")
    for line in protein_file:
        if line[0] != "#":      # To exclude comments
            if line[0] != ">":
                for amin_acid in line.strip():      # Filling up a new sequence
                    current_protein_seq.append(amin_acid)

            else:
                protein_sequences.append([protein_header, current_protein_seq])
                # Reset list/ append new first argument
                current_protein_seq = []
                protein_header = line.strip()

    # Appending the last remaining region if possible
    protein_sequences.append([protein_header, current_protein_seq])

    protein_file.close()
    # Remove the first initialising region
    protein_sequences.pop(0)

    # Searching the fitting protein for each CDS
    # Consisting of lists with the sublist with the protein data, and a matching sublist with the CDS data
    combined_cds_protein_list = []      # [[cds, protein], ...]

    print("--> search matching CDS-protein pairs")
    # Searches the matching pairs and appends them into the combined list
    for cds_region in cds_list:

        cds_id = cds_region[0].split("protein_id=")[1].split("]")[0]
        for protein in protein_sequences:
            protein_id = protein[0].split()[0][1:]
            if cds_id == protein_id:
                combined_cds_protein_list.append([cds_region, protein])
                break

    # Initialising
    list_length = len(combined_cds_protein_list)
    process_count = 0
    frameshift_detection_threshold_list = []    # Saves at what position a frameshift was detected in each query
    base_list = ["A", "T", "G", "C"]

    print("--> start with exonerate runs")
    # Run to detect the detection threshold at the beginning
    for region_pair in combined_cds_protein_list:

        # bad programmed process messages
        if process_count == int(list_length * 0.1):
            print("--> 10% finished")
        elif process_count == int(list_length * 0.25):
            print("--> 25% finished")
        elif process_count == int(list_length * 0.5):
            print("--> 50% finished")
        elif process_count == int(list_length * 0.75):
            print("--> 75% finished")

        for position in range(40):      # 40 as biggest position, to prevent unnecessary runs

            insertion_base = random.choice(base_list)       # choose a random base
            # Create current CDS query
            # Insert the frameshift on the current position as a insertion
            mod_current_cds = ("".join(region_pair[0][1][:position])) + insertion_base +\
                              ("".join(region_pair[0][1][position:]))

            # Create CDS query file, old data will be overwritten by the "w" parameter
            new_cds_query = open("cds_query.fasta", "w")
            new_cds_query.write(region_pair[0][0] + "\n")       # Appending the CDS header
            new_cds_query.write(mod_current_cds + "\n")         # Appending the CDS nucleotide sequence
            new_cds_query.close()

            # Create current protein subject
            mod_current_prot = region_pair[1][1]
            mod_current_prot = "".join(mod_current_prot)

            # Create protein query file, old data will be overwritten by the "w" parameter
            new_prot_sbjct = open("prot_sbjct.fasta", "w")
            new_prot_sbjct.write(region_pair[1][0] + "\n")
            new_prot_sbjct.write(mod_current_prot + "\n")
            new_prot_sbjct.close()

            # Run exonerate with the current data, the python script waits till the shell command is finished
            os.system("exonerate -m protein2dna --showvulgar true -Q protein -T dna --showalignment false --verbose 0 "
                      "-q prot_sbjct.fasta -t cds_query.fasta >exonerate_temp_output.txt")

            # Reset the parameter
            frameshift_detected = False

            # Read in the exonerate output, file online consisting of one vulgar-format line
            new_exonerate_output = open("exonerate_temp_output.txt")
            for line in new_exonerate_output:
                for char in line:
                    if char == "F":     # vulgar format for frameshift detected
                        frameshift_detected = True

            # Breaks the for loop if a frameshift is detected
            if frameshift_detected:
                break

        # Appends the current position, where a frameshift was detected, if no frameshift was detected and the loop was
        # Completely run, also the current last loop run position will be added (like upper threshold)
        frameshift_detection_threshold_list.append(position)

        process_count += 1      # Increase the progress count

    # Saving the results in a file and calculating the mean
    result_sum = 0
    output_file = open("exonerate_output.txt", "w")

    for x in frameshift_detection_threshold_list:
        result_sum += x
        output_file.write(str(x) + "\n")

    output_file.close()

    print("-->  The mean threshold for exonerate frameshift detection is: ",
          result_sum/len(frameshift_detection_threshold_list))

    return None


def database_comp_diamond_del(query_list, protein_database):
    """
     (create frameshift with deletion version)
    Receives a query list, which is containing all CDS region, which are used for the benchmark, and the link to the
    diamond database for the diamond blastx runs. The frameshifts are inserted as deletion like in the
    database_comp_exonerate function.
    :param query_list: List containing all used CDS regions for the benchmark
    :param protein_database: path to the protein database for the diamond blastx runs
    :return: None, but an output file containing the results named diamond_output is created
    """

    # Initialising
    list_length = len(query_list)
    process_count = 0
    frameshift_detection_threshold_list = []  # Saves at what position a frameshift was detected in each query

    print("--> start with diamond runs")
    # Run to detect the detection threshold at the beginning
    for cds_region in query_list:

        # bad programmed process messages
        if process_count == int(list_length * 0.1):
            print("--> 10% finished")
        elif process_count == int(list_length * 0.25):
            print("--> 25% finished")
        elif process_count == int(list_length * 0.5):
            print("--> 50% finished")
        elif process_count == int(list_length * 0.75):
            print("--> 75% finished")

        for position in range(40):  # 40 as biggest position, to prevent unnecessary runs
            # Create current CDS query
            # Insert the frameshift on the current position as a deletion
            mod_current_cds = ("".join(cds_region[1][:position])) + ("".join(cds_region[1][position + 1:]))

            # Create CDS query file, old data will be overwritten by the "w" parameter
            new_cds_query = open("cds_query.fasta", "w")
            new_cds_query.write(cds_region[0] + "\n")  # Appending the CDS header
            new_cds_query.write(mod_current_cds + "\n")  # Appending the CDS nucleotide sequence
            new_cds_query.close()

            # Run diamond with the current data, the python script waits till the shell command is finished
            os.system("diamond blastx -d " + protein_database + " -q cds_query.fasta -o diamond_temp_output.txt -k 1 "
                                                                "--quiet -F 15 -f 0 ")

            # Reset the parameter
            frameshift_detected = False

            # Read in the exonerate output, file online consisting of one vulgar-format line
            new_diamond_output = open("diamond_temp_output.txt")
            for line in new_diamond_output:
                for char in line.strip():
                    if len(line) > 5:       # to exclude empty lines
                        if line[0:6] == "Query ":       # case for alignment line
                            if (char == "/") or (char == "\\"):  # frameshift detected
                                frameshift_detected = True
                                break
                # Inner break if a frameshift is detected
                if frameshift_detected:
                    break

            # Breaks the for loop if a frameshift is detected
            if frameshift_detected:
                break

        # Appends the current position, where a frameshift was detected, if no frameshift was detected and the loop was
        # Completely run, also the current last loop run position will be added (like upper threshold)
        frameshift_detection_threshold_list.append(position)

        process_count += 1  # Increase the progress count

    # Saving the results in a file and calculating the mean
    result_sum = 0
    output_file = open("diamond_output.txt", "w")

    for x in frameshift_detection_threshold_list:
        result_sum += x
        output_file.write(str(x) + "\n")

    output_file.close()

    print("-->  The mean threshold for diamond frameshift detection is: ",
          result_sum / len(frameshift_detection_threshold_list))

    return None


def database_comp_diamond_ins(query_list, protein_database):
    """
    (create frameshift with insertion version)
    Receives a query list, which is containing all CDS region, which are used for the benchmark, and the link to the
    diamond database for the diamond blastx runs. The frameshifts are inserted as deletion like in the
    database_comp_exonerate function.
    :param query_list: List containing all used CDS regions for the benchmark
    :param protein_database: path to the protein database for the diamond blastx runs
    :return: None, but an output file containing the results named diamond_output is created
    """

    # Initialising
    list_length = len(query_list)
    process_count = 0
    frameshift_detection_threshold_list = []  # Saves at what position a frameshift was detected in each query
    base_list = ["A", "T", "G", "C"]

    print("--> start with diamond runs")
    # Run to detect the detection threshold at the beginning
    for cds_region in query_list:

        # bad programmed process messages
        if process_count == int(list_length * 0.1):
            print("--> 10% finished")
        elif process_count == int(list_length * 0.25):
            print("--> 25% finished")
        elif process_count == int(list_length * 0.5):
            print("--> 50% finished")
        elif process_count == int(list_length * 0.75):
            print("--> 75% finished")

        for position in range(40):  # 40 as biggest position, to prevent unnecessary runs
            insertion_base = random.choice(base_list)  # choose a random base
            # Create current CDS query
            # Insert the frameshift on the current position as a insertion
            mod_current_cds = ("".join(cds_region[1][:position])) + insertion_base + \
                              ("".join(cds_region[1][position:]))

            # Create CDS query file, old data will be overwritten by the "w" parameter
            new_cds_query = open("cds_query.fasta", "w")
            new_cds_query.write(cds_region[0] + "\n")  # Appending the CDS header
            new_cds_query.write(mod_current_cds + "\n")  # Appending the CDS nucleotide sequence
            new_cds_query.close()

            # Run diamond with the current data, the python script waits till the shell command is finished
            os.system("diamond blastx -d " + protein_database + " -q cds_query.fasta -o diamond_temp_output.txt -k 1 "
                                                                "--quiet -F 15 -f 0 ")

            # Reset the parameter
            frameshift_detected = False

            # Read in the exonerate output, file online consisting of one vulgar-format line
            new_diamond_output = open("diamond_temp_output.txt")
            for line in new_diamond_output:
                for char in line.strip():
                    if len(line) > 5:       # to exclude empty lines
                        if line[0:6] == "Query ":       # case for alignment line
                            if (char == "/") or (char == "\\"):  # frameshift detected
                                frameshift_detected = True
                                break
                # Inner break if a frameshift is detected
                if frameshift_detected:
                    break

            # Breaks the for loop if a frameshift is detected
            if frameshift_detected:
                break

        # Appends the current position, where a frameshift was detected, if no frameshift was detected and the loop was
        # Completely run, also the current last loop run position will be added (like upper threshold)
        frameshift_detection_threshold_list.append(position)

        process_count += 1  # Increase the progress count

    # Saving the results in a file and calculating the mean
    result_sum = 0
    output_file = open("diamond_output.txt", "w")

    for x in frameshift_detection_threshold_list:
        result_sum += x
        output_file.write(str(x) + "\n")

    output_file.close()

    print("-->  The mean threshold for diamond frameshift detection is: ",
          result_sum / len(frameshift_detection_threshold_list))

    return None


def create_output_plots(input_file_path, picture_name):
    """
    Creates the bar graphs, to visualize the benchmark results. All position are += 1 to convert the python list
    positions to the real positions.
    :param input_file_path: path to the output file of the database_comp_* functions.
    :param picture_name: name for the output picture
    :return: a png file with the bar graph
    """

    # list containing all alignment cutoff positions of the input file
    cutoff_list = []

    # read in the input file
    input_file = open(input_file_path)

    # appends all cutoff to the list
    for line in input_file:
        if int(line.strip()) != 0:      # to exclude the very view false detections
            cutoff_list.append(int(line.strip()))

    # count the amount of cutoff per position
    count_cutoff_list = []

    for cutoff in cutoff_list:
        new_cutoff_position = True

        for cutoff_position in count_cutoff_list:
            if cutoff == cutoff_position[0]:
                cutoff_position[1] += 1
                new_cutoff_position = False
                break

        if new_cutoff_position:
            count_cutoff_list.append([cutoff, 1])

    input_file.close()

    # sort the cutoff list ascending by their cutoff positions
    count_cutoff_list = sorted(count_cutoff_list, key=lambda current_cuttoff: int(current_cuttoff[0]))
    print(count_cutoff_list)

    # x-coordinates of left sides of bars
    cuttoff_positions = [x for x in range(1, len(count_cutoff_list)+1)]

    # heights of bars
    cuttoff_counts = [y[1] for y in count_cutoff_list]

    # labels for bars
    tick_label = [str(z[0]+1) for z in count_cutoff_list]     # +1 to correct the python list positions to the real pos

    # plotting a bar chart
    plt.bar(cuttoff_positions, cuttoff_counts, tick_label=tick_label,
            width=0.8, color=['green'])

    # naming the x-axis
    plt.xlabel('cutoff positions')
    # naming the y-axis
    plt.ylabel('cutoffs per position')
    # plot title
    plt.title('Diamond alignment cutoff distribution')

    # function saves the plot
    plt.savefig(picture_name + '.png')

    return None


def complete_prot_exonerate_del(cds_list, protein_file_path):
    """
    (create frameshift with deletion version)
    This functions determines the frameshift detection rate in CDS regions by Exonerate, for this the first 50 and last
    50 nucleotides are skipped, to exclude the alignment cuttoff problem.
    :param cds_list: list containing the unmodified cds regions, which are later used as modified queries
    :param protein_file_path:   path to the protein file containing the amino acid sequences
    :return: None, but an output file containing the results named exonerate_output is created
    """

    # Reading in the protein sequences
    # Reading in the queries out of the cds file
    protein_sequences = []
    protein_file = open(protein_file_path)
    protein_header = ""
    current_protein_seq = ["#####"]

    print("--> read in protein file")
    for line in protein_file:
        if line[0] != "#":  # To exclude comments
            if line[0] != ">":
                for amin_acid in line.strip():  # Filling up a new sequence
                    current_protein_seq.append(amin_acid)

            else:
                protein_sequences.append([protein_header, current_protein_seq])
                # Reset list/ append new first argument
                current_protein_seq = []
                protein_header = line.strip()

    # Appending the last remaining region if possible
    protein_sequences.append([protein_header, current_protein_seq])

    protein_file.close()
    # Remove the first initialising region
    protein_sequences.pop(0)

    # Searching the fitting protein for each CDS
    # Consisting of lists with the sublist with the protein data, and a matching sublist with the CDS data
    combined_cds_protein_list = []  # [[cds, protein], ...]

    print("--> search matching CDS-protein pairs")
    # Searches the matching pairs and appends them into the combined list
    for cds_region in cds_list:

        cds_id = cds_region[0].split("protein_id=")[1].split("]")[0]
        for protein in protein_sequences:
            protein_id = protein[0].split()[0][1:]
            if cds_id == protein_id:
                combined_cds_protein_list.append([cds_region, protein])
                break

    # Initialising
    list_length = len(combined_cds_protein_list)
    process_count = 0
    frameshift_detectionrate_list = []  # Saves the frameshift detection rate for each CDS

    print("--> start with exonerate runs")
    # Run to detect the detection threshold at the beginning
    for region_pair in combined_cds_protein_list:

        # bad programmed process messages
        if process_count == int(list_length * 0.1):
            print("--> 10% finished")
        elif process_count == int(list_length * 0.25):
            print("--> 25% finished")
        elif process_count == int(list_length * 0.5):
            print("--> 50% finished")
        elif process_count == int(list_length * 0.75):
            print("--> 75% finished")

        # counts how many of the inserted frameshifts are detected
        frameshift_detection_count = 0

        # skipps the first and last 50 nucleotide positions to exclude the cuttoff problems
        for position in range(50, (len(region_pair[0][1]) - 50)):

            # Create current CDS query
            # Insert the frameshift on the current position as a deletion
            mod_current_cds = ("".join(region_pair[0][1][:position])) + ("".join(region_pair[0][1][position + 1:]))

            # Create CDS query file, old data will be overwritten by the "w" parameter
            new_cds_query = open("cds_query.fasta", "w")
            new_cds_query.write(region_pair[0][0] + "\n")  # Appending the CDS header
            new_cds_query.write(mod_current_cds + "\n")  # Appending the CDS nucleotide sequence
            new_cds_query.close()

            # Create current protein subject
            mod_current_prot = region_pair[1][1]
            mod_current_prot = "".join(mod_current_prot)

            # Create protein query file, old data will be overwritten by the "w" parameter
            new_prot_sbjct = open("prot_sbjct.fasta", "w")
            new_prot_sbjct.write(region_pair[1][0] + "\n")
            new_prot_sbjct.write(mod_current_prot + "\n")
            new_prot_sbjct.close()

            # Run exonerate with the current data, the python script waits till the shell command is finished
            os.system("exonerate -m protein2dna --showvulgar true -Q protein -T dna --showalignment false --verbose 0 "
                      "-q prot_sbjct.fasta -t cds_query.fasta >exonerate_temp_output.txt")

            # Read in the exonerate output, file online consisting of one vulgar-format line
            new_exonerate_output = open("exonerate_temp_output.txt")

            for line in new_exonerate_output:
                for elem in line:
                    if elem == "F":  # vulgar format for frameshift detected
                        frameshift_detection_count += 1
                        break
                break
        detection_rate = frameshift_detection_count / (len(region_pair[0][1]) - 100)     # calculate the detection rate
        # Appends the detection rate of the protein, which means how many of the inserted frameshift are detected
        frameshift_detectionrate_list.append(detection_rate)

        process_count += 1  # Increase the progress count

    # Saving the results in a file and calculating the mean
    mean_detection_percentage = 0
    output_file = open("exonerate_output.txt", "w")

    for x in frameshift_detectionrate_list:
        mean_detection_percentage += x
        output_file.write(str(x) + "\n")

    output_file.close()

    print("-->  The average mean detection rate is: ",
          mean_detection_percentage / len(frameshift_detectionrate_list))
    print("-->  The average mean detection percentage is: ",
          (mean_detection_percentage / len(frameshift_detectionrate_list)) * 100, "%")

    return None


def complete_prot_exonerate_ins(cds_list, protein_file_path):
    """
        (create frameshift with insertion version)
        This functions determines the frameshift detection rate in CDS regions by Exonerate, for this the first 50 and last
        50 nucleotides are skipped, to exclude the alignment cuttoff problem.
        :param cds_list: list containing the unmodified cds regions, which are later used as modified queries
        :param protein_file_path:   path to the protein file containing the amino acid sequences
        :return: None, but an output file containing the results named exonerate_output is created
        """

    # Reading in the protein sequences
    # Reading in the queries out of the cds file
    protein_sequences = []
    protein_file = open(protein_file_path)
    protein_header = ""
    current_protein_seq = ["#####"]

    print("--> read in protein file")
    for line in protein_file:
        if line[0] != "#":  # To exclude comments
            if line[0] != ">":
                for amin_acid in line.strip():  # Filling up a new sequence
                    current_protein_seq.append(amin_acid)

            else:
                protein_sequences.append([protein_header, current_protein_seq])
                # Reset list/ append new first argument
                current_protein_seq = []
                protein_header = line.strip()

    # Appending the last remaining region if possible
    protein_sequences.append([protein_header, current_protein_seq])

    protein_file.close()
    # Remove the first initialising region
    protein_sequences.pop(0)

    # Searching the fitting protein for each CDS
    # Consisting of lists with the sublist with the protein data, and a matching sublist with the CDS data
    combined_cds_protein_list = []  # [[cds, protein], ...]

    print("--> search matching CDS-protein pairs")
    # Searches the matching pairs and appends them into the combined list
    for cds_region in cds_list:

        cds_id = cds_region[0].split("protein_id=")[1].split("]")[0]
        for protein in protein_sequences:
            protein_id = protein[0].split()[0][1:]
            if cds_id == protein_id:
                combined_cds_protein_list.append([cds_region, protein])
                break

    # Initialising
    list_length = len(combined_cds_protein_list)
    process_count = 0
    frameshift_detectionrate_list = []  # Saves at what position a frameshift was detected in each query
    base_list = ["A", "T", "G", "C"]

    print("--> start with exonerate runs")
    # Run to detect the detection threshold at the beginning
    for region_pair in combined_cds_protein_list:

        # bad programmed process messages
        if process_count == int(list_length * 0.1):
            print("--> 10% finished")
        elif process_count == int(list_length * 0.25):
            print("--> 25% finished")
        elif process_count == int(list_length * 0.5):
            print("--> 50% finished")
        elif process_count == int(list_length * 0.75):
            print("--> 75% finished")

        # counts how many of the inserted frameshifts are detected
        frameshift_detection_count = 0

        # excludes the first and last 50 nucleotide positions to exclude the cuttoff problems
        for position in range(50, (len(region_pair[0][1]) - 50)):
            insertion_base = random.choice(base_list)  # choose a random base
            # Create current CDS query
            # Insert the frameshift on the current position as a insertion
            mod_current_cds = ("".join(region_pair[0][1][:position])) + insertion_base + \
                              ("".join(region_pair[0][1][position:]))

            # Create CDS query file, old data will be overwritten by the "w" parameter
            new_cds_query = open("cds_query.fasta", "w")
            new_cds_query.write(region_pair[0][0] + "\n")  # Appending the CDS header
            new_cds_query.write(mod_current_cds + "\n")  # Appending the CDS nucleotide sequence
            new_cds_query.close()

            # Create current protein subject
            mod_current_prot = region_pair[1][1]
            mod_current_prot = "".join(mod_current_prot)

            # Create protein query file, old data will be overwritten by the "w" parameter
            new_prot_sbjct = open("prot_sbjct.fasta", "w")
            new_prot_sbjct.write(region_pair[1][0] + "\n")
            new_prot_sbjct.write(mod_current_prot + "\n")
            new_prot_sbjct.close()

            # Run exonerate with the current data, the python script waits till the shell command is finished
            os.system("exonerate -m protein2dna --showvulgar true -Q protein -T dna --showalignment false --verbose 0 "
                      "-q prot_sbjct.fasta -t cds_query.fasta >exonerate_temp_output.txt")

            # Read in the exonerate output, file online consisting of one vulgar-format line
            new_exonerate_output = open("exonerate_temp_output.txt")

            for line in new_exonerate_output:
                for elem in line:
                    if elem == "F":  # vulgar format for frameshift detected
                        frameshift_detection_count += 1
                        break
                break
        detection_rate = frameshift_detection_count / (len(region_pair[0][1]) - 100)  # calculate the detection rate
        # Appends the detection rate of the protein, which means how many of the inserted frameshift are detected
        frameshift_detectionrate_list.append(detection_rate)

        process_count += 1  # Increase the progress count

    # Saving the results in a file and calculating the mean
    mean_detection_percentage = 0
    output_file = open("exonerate_output.txt", "w")

    for x in frameshift_detectionrate_list:
        mean_detection_percentage += x
        output_file.write(str(x) + "\n")

    output_file.close()

    print("-->  The average mean detection rate is: ",
          mean_detection_percentage / len(frameshift_detectionrate_list))
    print("-->  The average mean detection percentage is: ",
          (mean_detection_percentage / len(frameshift_detectionrate_list)) * 100, "%")

    return None


def complete_prot_diamond_del(query_list, protein_database):
    """
    (create frameshift with deletion version)
    This functions determines the frameshift detection rate in CDS regions by Diamond, for this the first 50 and last
    50 nucleotides are skipped, to exclude the alignment cuttoff problem.
    :param query_list: List containing all used CDS regions for the benchmark
    :param protein_database: path to the protein database for the diamond blastx runs
    :return: None, but an output file containing the results named diamond_output is created
    """

    # Initialising
    list_length = len(query_list)
    process_count = 0
    frameshift_detectionrate_list = []  # Saves the frameshift detection rate for each CDS

    print("--> start with diamond runs")
    # Run to detect the detection threshold at the beginning
    for cds_region in query_list:

        # bad programmed process messages
        if process_count == int(list_length * 0.1):
            print("--> 10% finished")
        elif process_count == int(list_length * 0.25):
            print("--> 25% finished")
        elif process_count == int(list_length * 0.5):
            print("--> 50% finished")
        elif process_count == int(list_length * 0.75):
            print("--> 75% finished")

        # counts how many of the inserted frameshifts are detected
        frameshift_detection_count = 0

        # skipps the first and last 50 nucleotide positions to exclude the cuttoff problems
        for position in range(50, (len(cds_region[1]) - 50)):
            # Create current CDS query
            # Insert the frameshift on the current position as a deletion
            mod_current_cds = ("".join(cds_region[1][:position])) + ("".join(cds_region[1][position + 1:]))

            # Create CDS query file, old data will be overwritten by the "w" parameter
            new_cds_query = open("cds_query.fasta", "w")
            new_cds_query.write(cds_region[0] + "\n")  # Appending the CDS header
            new_cds_query.write(mod_current_cds + "\n")  # Appending the CDS nucleotide sequence
            new_cds_query.close()

            # Run diamond with the current data, the python script waits till the shell command is finished
            os.system("diamond blastx -d " + protein_database + " -q cds_query.fasta -o diamond_temp_output.txt -k 1 "
                                                                "--quiet -F 15 -f 0 ")

            # Reset the parameter
            frameshift_detected = False

            # Read in the exonerate output, file online consisting of one vulgar-format line
            new_diamond_output = open("diamond_temp_output.txt")
            for line in new_diamond_output:
                for char in line.strip():
                    if len(line) > 5:  # to exclude empty lines
                        if line[0:6] == "Query ":  # case for alignment line
                            if (char == "/") or (char == "\\"):  # frameshift detected
                                frameshift_detected = True
                                frameshift_detection_count += 1
                                break
                # Inner break if a frameshift is detected
                if frameshift_detected:
                    break

        detection_rate = frameshift_detection_count / (len(cds_region[1]) - 100)  # calculate the detection rate
        # Appends the detection rate of the protein, which means how many of the inserted frameshift are detected
        frameshift_detectionrate_list.append(detection_rate)

        process_count += 1  # Increase the progress count

    # Saving the results in a file and calculating the mean
    mean_detection_percentage = 0
    output_file = open("diamond_output.txt", "w")

    for x in frameshift_detectionrate_list:
        mean_detection_percentage += x
        output_file.write(str(x) + "\n")

    output_file.close()

    print("-->  The average mean detection rate is: ",
          mean_detection_percentage / len(frameshift_detectionrate_list))
    print("-->  The average mean detection percentage is: ",
          (mean_detection_percentage / len(frameshift_detectionrate_list)) * 100, "%")

    return None


def complete_prot_diamond_ins(query_list, protein_database):
    """
    (create frameshift with insertion version)
    This functions determines the frameshift detection rate in CDS regions by Diamond, for this the first 50 and last
    50 nucleotides are skipped, to exclude the alignment cuttoff problem.
    :param query_list: List containing all used CDS regions for the benchmark
    :param protein_database: path to the protein database for the diamond blastx runs
    :return: None, but an output file containing the results named diamond_output is created
    """

    # Initialising
    list_length = len(query_list)
    process_count = 0
    frameshift_detectionrate_list = []  # Saves the frameshift detection rate for each CDS
    base_list = ["A", "T", "G", "C"]

    print("--> start with diamond runs")
    # Run to detect the detection threshold at the beginning
    for cds_region in query_list:

        # bad programmed process messages
        if process_count == int(list_length * 0.1):
            print("--> 10% finished")
        elif process_count == int(list_length * 0.25):
            print("--> 25% finished")
        elif process_count == int(list_length * 0.5):
            print("--> 50% finished")
        elif process_count == int(list_length * 0.75):
            print("--> 75% finished")

        # counts how many of the inserted frameshifts are detected
        frameshift_detection_count = 0

        # skipps the first and last 50 nucleotide positions to exclude the cuttoff problems
        for position in range(50, (len(cds_region[1]) - 50)):
            insertion_base = random.choice(base_list)  # choose a random base
            # Create current CDS query
            # Insert the frameshift on the current position as a insertion
            mod_current_cds = ("".join(cds_region[1][:position])) + insertion_base + \
                              ("".join(cds_region[1][position:]))

            # Create CDS query file, old data will be overwritten by the "w" parameter
            new_cds_query = open("cds_query.fasta", "w")
            new_cds_query.write(cds_region[0] + "\n")  # Appending the CDS header
            new_cds_query.write(mod_current_cds + "\n")  # Appending the CDS nucleotide sequence
            new_cds_query.close()

            # Run diamond with the current data, the python script waits till the shell command is finished
            os.system("diamond blastx -d " + protein_database + " -q cds_query.fasta -o diamond_temp_output.txt -k 1 "
                                                                "--quiet -F 15 -f 0 ")

            # Reset the parameter
            frameshift_detected = False

            # Read in the exonerate output, file online consisting of one vulgar-format line
            new_diamond_output = open("diamond_temp_output.txt")
            for line in new_diamond_output:
                for char in line.strip():
                    if len(line) > 5:  # to exclude empty lines
                        if line[0:6] == "Query ":  # case for alignment line
                            if (char == "/") or (char == "\\"):  # frameshift detected
                                frameshift_detected = True
                                frameshift_detection_count += 1
                                break
                # Inner break if a frameshift is detected
                if frameshift_detected:
                    break

        detection_rate = frameshift_detection_count / (len(cds_region[1]) - 100)  # calculate the detection rate
        # Appends the detection rate of the protein, which means how many of the inserted frameshift are detected
        frameshift_detectionrate_list.append(detection_rate)

        process_count += 1  # Increase the progress count

    # Saving the results in a file and calculating the mean
    mean_detection_percentage = 0
    output_file = open("diamond_output.txt", "w")

    for x in frameshift_detectionrate_list:
        mean_detection_percentage += x
        output_file.write(str(x) + "\n")

    output_file.close()

    print("-->  The average mean detection rate is: ",
          mean_detection_percentage / len(frameshift_detectionrate_list))
    print("-->  The average mean detection percentage is: ",
          (mean_detection_percentage / len(frameshift_detectionrate_list)) * 100, "%")

    return None


def main():

    print("benchOne MAIN FUNCTION called")

    """start_time = time.time()
    cds_file_path = "/share/project/johannes/bachelor/exonerate_vs_diamond_benchmark/NCBI_files/" \
                    "GCF_000001735.4_TAIR10.1_cds_from_genomic.fna"

    queries = create_queries(cds_file_path, 100, 800)

    protein_file_path = "/share/project/johannes/bachelor/exonerate_vs_diamond_benchmark/NCBI_files/modified_file/" \
                        "GCF_000001735.4_TAIR10.1_protein_mod.faa"

    #database_comp_exonerate_ins(queries, protein_file_path)

    complete_prot_exonerate_ins(queries, protein_file_path)

    stop_time = time.time()
    print("runtime: ", stop_time-start_time, " seconds")"""

    """start_time = time.time()
    #prot_database = "/share/project/johannes/bachelor/exonerate_vs_diamond_benchmark/diamond_database/diamond_db.dmnd"
    prot_database = "/home/johannes/Desktop/diamond_database/diamond_db.dmnd"
    cds_file_path = "/share/project/johannes/bachelor/exonerate_vs_diamond_benchmark/NCBI_files/" \
                    "GCF_000001735.4_TAIR10.1_cds_from_genomic.fna"

    queries = create_queries(cds_file_path, 50, 800)

    #database_comp_diamond_ins(queries, prot_database)
    complete_prot_diamond_ins(queries, prot_database)

    stop_time = time.time()
    print("runtime: ", stop_time - start_time, " seconds")"""

    #input_file_1 = "/share/project/johannes/bachelor/exonerate_vs_diamond_benchmark/exonerate_del_run/exonerate_output.txt"
    #create_output_plots(input_file_1, "counts_exonerate_del")

    #input_file_2 = "/share/project/johannes/bachelor/exonerate_vs_diamond_benchmark/diamond_ins_run/diamond_output.txt"
    #create_output_plots(input_file_2, "counts_diamond_ins")


if __name__ == '__main__':
    main()
