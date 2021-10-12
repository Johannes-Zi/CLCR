"""This file contains the functions for the creation of the output files"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import os
import glob
import copy


def heal_assembly_file(healing_region_list, input_fna_path, outut_dir):
    """
    Gets the healing_region_list which contains the frameshift positions in each query, and inserts N's at those
    positions, to heal the reading frame. The healed assembly is afterwards saved in the give output directory.
    The new assembly file doesnt contain comments.
    :param healing_region_list: Contains for each query: scaffold, start pos., end pos. in scaff., frameshift pos. list
    :param outut_dir: directory where the new "healed assembly" is stored (should contain a \ at the last position)
    :param input_fna_path:  File path to the original assembly file
    :return: File path to the new modified assembly file
    """
    # Read in the .fna file for fast access
    input_fna_file = open(input_fna_path)
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

    print(scaffold_list[0][0], len(scaffold_list[0][1]), scaffold_list[0][1][0:10])

    # Correct the frameshifts
    # Inserting the N's from the end to the beginning, to prevent the necessity of adapting the insertion position in
    # dependency to the previous inserted frameshifts
    # For this, in this first sort, the queries in each scaffold are sorted descending by their start position
    # but the reverse option also sorts the scaffolds descending
    sorted_healing_region_list = sorted(healing_region_list, reverse=True,  key=lambda temp_query: (temp_query[0],
                                                                                                    temp_query[1]))
    print("sorted")
    count = 0

    for query in sorted_healing_region_list:
        # Search the corresponding scaffold
        for scaffold in scaffold_list:
            # Case where matching scaffold is found
            if scaffold[0] == query[0]:
                query_start_pos = int(query[1])
                count += 1
                if count == 1200:
                    print("10% more")
                    count = 0

                # Second sort step, now the frameshift positions in each query are sorted descending
                temp_sorted_query = sorted(query[3], reverse=True,  key=lambda temp_query: temp_query[0])

                for frameshift_position in temp_sorted_query:      # Same usage of reversed() as before
                    # Deletion case
                    if frameshift_position[1] == "D":
                        scaffold[1] = scaffold[1][:(query_start_pos + frameshift_position[0])] + ["N"] + \
                                      scaffold[1][(query_start_pos + frameshift_position[0]):]
                    # Insertion case
                    else:
                        scaffold[1] = scaffold[1][:(query_start_pos + frameshift_position[0])] + ["N"] + ["N"] + \
                                      scaffold[1][(query_start_pos + frameshift_position[0]):]
                break

    # Create the new assembly .fna file path, located in the same dir as the original assembly file
    new_fna_file_path = outut_dir + "new_fna_file.txt"

    print("create new assembly file")

    # Creating the new assembly file (with the inserted N's)
    new_fna_file = open(new_fna_file_path, "w")
    for scaffold in scaffold_list:

        # Write the header line for the scaffold
        new_fna_file.write(">" + scaffold[0] + "\n")

        # Create lines with a length of 80, for a more readable .fasta file
        # Create list with the end positions of each 80 nucleotide long line
        position_list = list(range(80, len(scaffold[1]), 80))
        # Writes the lines in the assembly file
        for list_pos in position_list:
            # The nucleotides are stored as chars in a list and are converted to a string with join
            new_fna_file.write(("".join(scaffold[1][(list_pos - 80):list_pos])) + "\n")

        # Write the last remaining line (range function doesnt include the last elements)
        new_fna_file.write(("".join(scaffold[1][position_list[-1]:len(scaffold[1])])) + "\n")

    new_fna_file.close()

    return new_fna_file_path


def create_gff_adaption_list(healing_region_list):
    """
    Creates a list for the annotation file adaption to the new assembly with the inserted N's, to speed up the process.
    The healing_region_list is used to determine how many frameshifts are inserted till a certain position. Those
    information is saved in the adaption_list, which has a sublist for each scaffold. Each scaffold sublist, contains
    all positions where the amount of inserted N's has changed, which means all positions where frameshifts where
    inserted. For each of this positions is saves how many frameshifts are inserted till this position INCLUDING the
    inserted frameshifts in this position.
    :param healing_region_list: Contains for each query [scaffold, start pos., end pos., [[frameshift position, "I" or
    "D"], [...], ...]
    :return: adaption_list
    """

    # Saves for all scaffolds the positions, when the count of the inserted frameshifts rises. This brings the
    # possibility for a fast check, on how many positions a feature needs to be shifted
    # The List has a sublist for each scaffold, consisting of lists (scaffold_pos, already_inserted_frameshift_count)
    adaption_list = []
    current_scaffold_list = ["####", []]

    # Sort the List after theirs scaffolds and ascending by their start positions in the scaffold
    sorted_healing_region_list = sorted(healing_region_list, key=lambda temp_query: (temp_query[0], temp_query[1]))

    # Parameter initialising
    frameshift_count = 0  # Counts how many frameshifts are inserted in the scaffold, till the current position

    # Fill the list
    for query in sorted_healing_region_list:

        # Case where a new scaffold starts
        if query[0] != current_scaffold_list[0]:
            # Append the old scaffold to the list
            adaption_list.append(current_scaffold_list)

            # Initialise the new scaffold
            frameshift_count = 0
            current_scaffold_list = [query[0], []]

            # Sorts the frameshift positions, of the current query, ascending
            sorted_frameshifts = sorted(query[3])
            # Appends each frameshift to the list
            for single_frameshift in sorted_frameshifts:

                if single_frameshift[1] == "I":     # Insertion case, two additional N's are inserted at this position
                    frameshift_count += 2

                else:       # Deletion case, one additional N is inserted
                    frameshift_count += 1

                # Append the frameshift position and the count of how many frameshifts are inserted till this position
                # to the list, query[1] is the start position of the query in the scaffold and single_frameshift[0] is
                # the position of the frameshift in the in the query
                current_scaffold_list[1].append([(query[1] + single_frameshift[0]), frameshift_count])

        # Case, where the current query is in the same scaffold as the previous region
        else:
            # Sorts the frameshift positions, of the current query, ascending
            sorted_frameshifts = sorted(query[3])
            # Appends each frameshift to the list
            for single_frameshift in sorted_frameshifts:

                if single_frameshift[1] == "I":  # Insertion case, two additional N's are inserted at this position
                    frameshift_count += 2

                else:  # Deletion case, one additional N is inserted
                    frameshift_count += 1

                # Append the frameshift position and the count of how many frameshifts are inserted till this position
                # to the list
                current_scaffold_list[1].append([(query[1] + single_frameshift[0]), frameshift_count])

    # Append the last remaining scaffold
    adaption_list.append(current_scaffold_list)
    # Remove the initialising element
    adaption_list.pop(0)

    return adaption_list


def adapt_gff_file(gff_file_path, adaption_list):
    """
    Adapts a given annotation file to the modified assembly(inserted N's). Shifts the start and end positions of
    annotation features downstream as far a necessary (dependent on how many N's are inserted till those positions).
    The list which is used for the adaption is preprocessed by the create_gff_adaption_list function. Please look at the
    header of this function for a closer description of the used list.
    :param gff_file_path: Path to the .gff file
    :param adaption_list: List which contains sublists for each scaffold with fixed frameshifts. Each scaffold sublist,
    contains all positions where the amount of inserted N's has changed, which means all positions where frameshifts
    where inserted. For each of this positions is saves how many frameshifts are inserted till this position INCLUDING
    the inserted frameshifts in this position.
    :return: None, but creates a new_gff_file.txt with the new annotation
    """

    feature_list = []       # contains all features of the annotation file
    annotation_file = open(gff_file_path)

    for line in annotation_file:
        if line[0] != "#":      # Excludes comments
            line_split = line.split()
            # Corrects the feature start and end positions in the scaffold if necessary
            for scaffold in adaption_list:     # Scaffolds without frameshifts doesn't occur
                if scaffold[0] == line_split[0]:       # Matching scaffold found
                    # Checks how on how many positions the start position needs to be shifted
                    # Searches for the first entry, which is downstream of the start pos., and then uses the frameshift
                    # count of the previous entry

                    # Add initialising objects to fix exception cases
                    temp_frameshift_list = copy.copy(scaffold[1])
                    # Quick solution for the case where no frameshift is downstream of the start or end pos.
                    # Append a last infinite element, which is bigger than all start/end positions
                    temp_frameshift_list.append([float("inf"), 0])
                    # Quick solution for the case where no frameshift is upstream of the start or end pos.
                    # Append a first element, which is smaller than all start/end positions
                    temp_frameshift_list = [[0, 0]] + temp_frameshift_list

                    for list_pos, frameshift in enumerate(temp_frameshift_list):
                        # Case, where the first entry, which is downstream of the start pos., is found
                        if frameshift[0] > int(line_split[3]):

                            # Adds the frameshift count to the start position
                            line_split[3] = str(int(line_split[3]) + temp_frameshift_list[(list_pos - 1)][1])
                            break

                    # Checks how on how many positions the end position needs to be shifted
                    for list_pos, frameshift in enumerate(temp_frameshift_list):
                        # Case, where the first entry, which is downstream of the start pos., is found

                        if frameshift[0] > int(line_split[4]):
                            # Adds the frameshift count to the start position
                            line_split[4] = str(int(line_split[4]) + temp_frameshift_list[(list_pos - 1)][1])
                            break

                    break

            feature_list.append("\t".join(line_split))

    annotation_file.close()

    # Create the new annotation file
    new_gff_file_path = gff_file_path[:(gff_file_path.rfind("/") + 1)] + "new_gff_file.gff"

    new_gff_file = open(new_gff_file_path, "w")

    # Write the adapted features in the new file
    for feature_line in feature_list:
        new_gff_file.write(feature_line + "\n")

    new_gff_file.close()

    return None


def add_gff_information(gff_file, output_region_list, fna_file):

    """
    adding to each detected region the overlapping feature information of the annotation file
    :param gff_file: path to gff file
    :param output_region_list: list with the output regions
    :param fna_file: path to the fna file
    :return: the output_region_list with the additional feature information
    """

    # reading out the different .fasta/.fna sequences, to assign th correct positions in the complete sequence to the
    # .gff regions, for a correct assignment from low cov. region to .gff region

    # initialising the parameters
    fna_region_list = []
    current_fna_region = "#"
    current_position = 0
    region_start = 0
    input_file_2 = open(fna_file)

    print("reading in .fna file")
    # reading in the start and the end of each region
    for line in input_file_2:

        # detecting new region, appending the old region to the list
        """ACHTUNG AUSLESEN DER SCAFFOLDS SEHR SPEZIFISCH (HEADER ZEILE KANN ANDERES FORMAT HABEN!!!!!)"""
        if line[0] == ">":
            fna_region_list.append((current_fna_region, region_start, current_position))
            region_start = current_position
            current_fna_region = line[1:(len(line))].strip()

        else:
            if line[0] != ";":
                current_position += len(line.strip())

    # appending the last remaining element
    fna_region_list.append((current_fna_region, region_start, current_position))
    # deleting the initialising element
    fna_region_list.pop(0)

    # read the regions out of the gff file and saving them in region_list
    input_file = open(gff_file)
    gff_region_list = []  # containing tuple like (region_start, region_end, region_type)
    # counting the regions where no parent seq. was found at the .fna file (there should be a parent regions, and if not
    # something went wrong!!!)
    parent_not_found = 0

    print("reading in .gff file")
    # reading the exon regions out of the gff
    for line in input_file:
        if line[0] != "#":
            # for excluding the false line (lines without gff table)
            if len(line.split()) > 7:
                if (line.split()[2] == "exon") or (line.split()[2] == "intron"):

                    # calculating the positions in the whole sequence, using the scaffold region positions
                    parent_sequence = line.split()[0]

                    # searching the corresponding parent sequence for each .gff region
                    parent_found = False    # value containing the information if a parent sequence was found (there
                    # should be a parent seq. if not there ist some thing wrong with the .gff)
                    for parent in fna_region_list:
                        if parent[0] == parent_sequence:
                            parent_reg_start = parent[1]
                            gff_reg_start = parent_reg_start + int(line.split()[3])
                            gff_reg_end = parent_reg_start + int(line.split()[4])
                            parent_found = True
                            break

                    # appending the .gff region start and end in the complete sequence
                    if parent_found:
                        gff_region_list.append((gff_reg_start, gff_reg_end, line.split()[2]))

                    # counting the problem regions where no parent seq. was found
                    else:
                        parent_not_found += 1

    input_file.close()

    print("lenght of gff_region_list: ", len(gff_region_list))

    # sorting the unsorted region_list ascending by the regions end positions
    gff_region_list = sorted(gff_region_list, key=lambda current_tuple: current_tuple[1])
    # invert the sorted list to sort them descending by their region end positions
    gff_region_list.reverse()       # no runtime differences to instead using the ascending list in a for loop

    print("combining the information")
    # assign regions of output_region_list to regions in region_list if possible
    for current_out_region in output_region_list:

        current_gff_region = 0
        out_reg_start = int(current_out_region[0])      # starting position of the current "output" region
        out_reg_end = int(current_out_region[1])       # end position of the current "output" region
        region_type = []  # for saving the output region type, if there is a match with a gff region

        # runs as long as the start of the output region is smaller than the end of the current gff_region
        while out_reg_start <= gff_region_list[current_gff_region][1]:
            gff_reg_start = gff_region_list[current_gff_region][0]
            gff_reg_end = gff_region_list[current_gff_region][1]

            # checks if either the output regions start ot the end lays in a gff region, or the gff region is
            # completely covered
            if ((out_reg_end >= gff_reg_start) and (out_reg_end <= gff_reg_end)) or \
                ((out_reg_start >= gff_reg_start) and (out_reg_start <= gff_reg_end)) or \
                    ((out_reg_start <= gff_reg_start) and (out_reg_end >= gff_reg_end)):
                # case if region type is empty
                gff_type = gff_region_list[current_gff_region][2]   # current gff region type
                if not region_type:
                    # appends the gff region type
                    region_type.append(gff_type)

                # case if the output region lays in more than one gff region
                else:
                    # checks if the region is already assigned to that gff region type
                    new_type = True
                    for x in region_type:
                        if x == gff_type:
                            new_type = False

                    if new_type:    # case if a new type was found
                        region_type.append(gff_type)

            current_gff_region += 1
            # break condition
            if current_gff_region == len(gff_region_list):
                break

        # case if the output region was assigned to a gff_region
        if region_type:
            # case if both region types are covered
            if len(region_type) == 2:
                current_out_region.append("exon_intron")
            # case if only one region type is covered
            else:
                current_out_region.append(region_type[0])

        # case if there was no assignment, so the region is noncoding
        else:
            current_out_region.append("noncoding")

    print("Regions with no parentsequences (This shouldn't happen!): ", parent_not_found)

    return output_region_list


def add_gff_information_no_introns(gff_file, output_region_list, fna_file):

    """
    function when no intron regions are included in the -gff file (instead using gene regions)
    adding to each detected region the overlapping feature information of the annotation file
    :param gff_file: path to gff file
    :param output_region_list: list with the output regions
    :param fna_file: path to the fna file
    :return: the output_region_list with the additional feature information
    """

    # reading out the different .fasta/.fna sequences, to assign th correct positions in the complete sequence to the
    # .gff regions, for a correct assignment from low cov. region to .gff region

    # initialising the parameters
    fna_region_list = []
    current_fna_region = "#"
    current_position = 0
    region_start = 0
    input_file_2 = open(fna_file)

    print("reading in .fna file")
    # reading in the start and the end of each region
    for line in input_file_2:

        # detecting new region, appending the old region to the list
        """ACHTUNG AUSLESEN DER SCAFFOLDS SEHR SPEZIFISCH (HEADER ZEILE KANN ANDERES FORMAT HABEN!!!!!)"""
        if line[0] == ">":
            fna_region_list.append((current_fna_region, region_start, current_position))
            region_start = current_position
            current_fna_region = line[1:(len(line))].strip()

        else:
            if line[0] != ";":
                current_position += len(line.strip())

    input_file_2.close()

    # appending the last remaining element
    fna_region_list.append((current_fna_region, region_start, current_position))

    # deleting the initialising element
    fna_region_list.pop(0)

    # read the regions out of the gff file and saving them in gff_region_list
    input_file = open(gff_file)
    gff_region_list = []  # containing tuple like (region_start, region_end, region_type)
    # counting the regions where no parent seq. was found at the .fna file (there should be a parent regions, and if not
    # something went wrong!!!)
    parent_not_found = 0

    print("reading in .gff file")
    # reading the exon regions out of the gff
    for line in input_file:
        if line[0] != "#":
            # for excluding the false line (lines without gff table)
            if len(line.split()) > 7:
                if (line.split()[2] == "exon") or (line.split()[2] == "gene"):

                    # calculating the positions in the whole sequence, using the scaffold region positions
                    parent_sequence = line.split()[0]

                    # searching the corresponding parent sequence for each .gff region
                    parent_found = False    # value containing the information if a parent sequence was found (there
                    # should be a parent seq. if not there ist some thing wrong with the .gff)
                    for parent in fna_region_list:
                        if parent[0] == parent_sequence:
                            parent_reg_start = parent[1]
                            gff_reg_start = parent_reg_start + int(line.split()[3])
                            gff_reg_end = parent_reg_start + int(line.split()[4])
                            parent_found = True
                            break

                    # appending the .gff region start and end in the complete sequence
                    if parent_found:

                        # CDS is changed to intron, to fit for the later assignment approach (later no type
                        # change needed)
                        if line.split()[2] == "gene":
                            gff_region_list.append((gff_reg_start, gff_reg_end, "intron/UTR"))

                        # normal exon case
                        else:
                            gff_region_list.append((gff_reg_start, gff_reg_end, line.split()[2]))

                    # counting the problem regions where no parent seq. was found
                    else:
                        parent_not_found += 1

    input_file.close()

    print("lenght of gff_region_list: ", len(gff_region_list))

    for x in gff_region_list:
        if (x[0] == 5245342) or (x[0] == 5245414) or (x[0] == 5245710):
            print(x)

    # sorting the unsorted region_list ascending by the regions end positions
    gff_region_list = sorted(gff_region_list, key=lambda current_tuple: current_tuple[1])
    # invert the sorted list to sort them descending by their region end positions
    gff_region_list.reverse()       # no runtime differences to instead using the ascending list in a for loop

    print("assingning annotation features if possible")
    # assign regions of output_region_list to regions in region_list if possible

    progress_count = 0
    len_out_reg_4 = int(len(output_region_list) / 4)

    for current_out_region in output_region_list:

        # printing the progress
        if progress_count == len_out_reg_4:
            print("25% reached")
        if progress_count == (len_out_reg_4 * 2):
            print("50% reached")
        if progress_count == (len_out_reg_4 * 3):
            print("75% reached")


        current_gff_region = 0
        out_reg_start = int(current_out_region[0])      # starting position of the current "output" region
        out_reg_end = int(current_out_region[1])       # end position of the current "output" region
        region_type = []  # for saving the output region type, if there is a match with a gff region

        # runs as long as the start of the output region is smaller than the end of the current gff_region
        while out_reg_start <= gff_region_list[current_gff_region][1]:
            gff_reg_start = gff_region_list[current_gff_region][0]
            gff_reg_end = gff_region_list[current_gff_region][1]

            # checks if either the output regions start or the end lays in a gff region, or the gff region is
            # completely covered
            if ((out_reg_end >= gff_reg_start) and (out_reg_end <= gff_reg_end)) or \
                ((out_reg_start >= gff_reg_start) and (out_reg_start <= gff_reg_end)) or \
                    ((out_reg_start <= gff_reg_start) and (out_reg_end >= gff_reg_end)):
                # case if region type is empty
                gff_type = gff_region_list[current_gff_region][2]   # current gff region type
                if not region_type:
                    # appends the gff region type
                    region_type.append((gff_type, gff_reg_start, gff_reg_end))

                # case if the output region lays in more than one gff region
                else:
                    # checks if the region is already assigned to that gff region type
                    new_type = True
                    for x in region_type:
                        if x[0] == gff_type:
                            new_type = False

                    if new_type:    # case if a new type was found
                        region_type.append((gff_type, gff_reg_start, gff_reg_end))

            current_gff_region += 1
            # break condition
            if current_gff_region == len(gff_region_list):
                break

        # case if the output region was assigned to a gff_region
        if region_type:
            # case if the region covers exon and cds
            if len(region_type) == 2:

                # assigning the gff type region lenght
                for reg_type in region_type:
                    if reg_type[0] == "exon":
                        exon_start = reg_type[1]
                        exon_end = reg_type[2]
                        break

                # case where the low cow. regions lay completely in the exon
                if (out_reg_start >= exon_start) and (out_reg_end <= exon_end):
                    current_out_region.append("exon")

                # case where an exon and an intron is covered by the low cov. region
                else:
                    current_out_region.append("exon_intron")

            # case if only one region type is covered
            else:
                current_out_region.append(region_type[0][0])

        # case if there was no assignment, so the region is noncoding
        else:
            current_out_region.append("noncoding")

        progress_count += 1

    print("Regions with no parentsequences (This shouldn't happen!): ", parent_not_found)
    print("lenght output_region_list: ", len(output_region_list))

    # adding the distance to the next gene
    cds_region_list = []  # list containing all CDS Regions out of the .gff file

    # appending all CDS regions in the new list
    for region in gff_region_list:
        if region[2] == "intron":   # remember, CDS was changed to intron
            cds_region_list.append([region[0], region[1], "CDS"])

    print("adding distance to next CDS information...")
    progress_count = 0
    len_out_reg_4 = int(len(output_region_list)/4)


    # adding the information to every output region
    for current_region in output_region_list:

        if progress_count == len_out_reg_4:
            print("25% reached")
        if progress_count == (len_out_reg_4 * 2):
            print("50% reached")
        if progress_count == (len_out_reg_4 * 3):
            print("75% reached")


        shortest_distance = cds_region_list[0][1]  # initialising the variable with the largest possible distance
        current_cds_region = 0      # count for saving the currently used cds region (tuple position in list)
        out_reg_start = int(current_region[0])  # starting position of the current "output" region

        # the regions in the cds_region_list are still sorted descending by their end position
        # runs as long as the start of the output region is smaller than the end of the current cds_region
        while out_reg_start < cds_region_list[current_cds_region][1]:
            # case that the current_region is upstream from the cds_region, AND not overlapping
            if int(current_region[1]) < cds_region_list[current_cds_region][0]:
                temp_distance = cds_region_list[current_cds_region][0] - int(current_region[1])  # distance between the considere regions
                # case that the new distance from the current_region end to the cds_region start is smaller then the
                # current one
                if temp_distance < shortest_distance:
                    shortest_distance = temp_distance
            current_cds_region += 1

            if current_cds_region == len(cds_region_list):
                break

        # consideres the cds regions with a start upstream of the output region start, the first of those regions,
        # which is not overlapping with the output region, is potentially closer than the closest downstream region
        cds_region_list_len = len(cds_region_list)
        while True:
            if current_cds_region == cds_region_list_len:
                break
            # case that the first upstream cds region, which is not overlapping with the output region is found
            if cds_region_list[current_cds_region][1] < out_reg_start:
                temp_distance = int(current_region[0]) - cds_region_list[current_cds_region][1]
                if temp_distance < shortest_distance:
                    shortest_distance = temp_distance
                break
            current_cds_region += 1

        # append the calculated distance
        current_region.append(shortest_distance)

        progress_count += 1

    return output_region_list


def create_combined_output_file(output_region_list, fna_file):

    """
    Creating the output .txt file, either from regions without or with gff information
    reading out the different .fasta/.fna sequences, to transform the complete sequence positions to the position in
    the specific scaffolds
    :param output_region_list: list containing the output regions
    :param fna_file: path to the fna_file
    :return: nothing
    """

    # initialising the parameters
    fna_region_list = []
    current_fna_region = "#"
    current_position = 0
    region_start = 0
    input_file_2 = open(fna_file)

    print("reading in .fna file...")
    # reading in the start and the end of each region
    for line in input_file_2:

        # detecting new region, appending the old region to the list
        """ACHTUNG AUSLESEN DER SCAFFOLDS SEHR SPEZIFISCH (HEADER ZEILE KANN ANDERES FORMAT HABEN!!!!!)"""
        if line[0] == ">":
            fna_region_list.append((current_fna_region, region_start, current_position))
            region_start = current_position
            current_fna_region = line[1:(len(line))].strip()

        else:
            if line[0] != ";":
                current_position += len(line.strip())

    # appending the last remaining element
    fna_region_list.append((current_fna_region, region_start, current_position))
    # deleting the initialising element
    fna_region_list.pop(0)

    # transforming the complete sequence positions to the position in the specific scaffolds
    new_output_region_list = []
    for output_region in output_region_list:
        output_region_start = int(output_region[0])
        output_region_end = int(output_region[1])
        output_region_scaffold = ""     # assigned scaffold of the current output region
        # searching the fitting  scaffold:
        for scaffolds in fna_region_list:
            # case if a fitting scaffold was found
            if (output_region_start >= scaffolds[1]) and (output_region_start <= scaffolds[2]):
                # modifying the values
                # assigned scaffold of the current output region
                output_region_scaffold = scaffolds[0]
                # start position in the assigned scaffold of the current output region
                output_region_start = output_region_start - scaffolds[1]
                # end position in the assigned scaffold of the current output region
                output_region_end = output_region_end - scaffolds[1]

        new_output_region_list.append([output_region_scaffold, output_region_start, output_region_end]
                                      + output_region[2:])

    # creating the output file
    output_file_path = os.getcwd() + "/output_file.txt"

    output_file = open(output_file_path, "w")

    for region in new_output_region_list:
        new_line = ""
        for region_info in region:
            new_line += (str(region_info) + "\t")
        output_file.write((new_line + "\n"))

    output_file.close()

    return None


def filter_output_file(output_file_path, identity_threshold, frameshift_filter):

    """
    filtering the output by throwing out the regions below a given identity threshold and if configured excluding regions
    without a frameshift
    :param output_file_path: output_file_path
    :param identity_threshold: identity threshold
    :param frameshift_filter: if TRUE regions without frameshifts are excluded
    :return: nothing
    """

    # containing the regions with an value higher than the threshold
    filtered_output_regions = []

    output_file = open(output_file_path)

    if frameshift_filter:
        for line in output_file:
            if (int(line.split()[6][:-1]) >= identity_threshold) and (int(line.split()[7]) > 0):
                filtered_output_regions.append(line)
    else:
        for line in output_file:
            if int(line.split()[6][:-1]) >= identity_threshold:
                filtered_output_regions.append(line)

    output_file.close()

    print("regions after filtering with threshold ", identity_threshold, ": ", len(filtered_output_regions))

    # creating the output file
    new_file_path = os.getcwd() + "/filtered_output_file.txt"
    new_file = open(new_file_path, "w")

    for line in filtered_output_regions:
        new_file.write(line)

    new_file.close()

    return None


def transform_txt_to_gff(txt_file_path):

    """
    creating a gff file out of the original txt output file
    :param txt_file_path: file path
    :return: nothing
    """

    txt_file = open(txt_file_path)

    # reading in the txt regions
    gff_line_list = []

    """for line in txt_file:

        line_t = line.strip().split()       # temporary list
        modified_line = line_t[0] + "\t" + "." + "\t" + "LowCov" + "\t" + line_t[1] + "\t" + line_t[2] + "\t" \
                        + ".\t.\t.\t" + "protein_hit \"" + line_t[3] + "\"; evalue \"" + line_t[4] + "\"; bitscore \"" \
                        + line_t[5] + "\"; identities \"" + line_t[6] + "\"; frameshift_count \"" + line_t[7] \
                        + "\"; assigned_gff_region_type \"" + line_t[8] + "\"; distance_to_next_CDS \"" + line_t[9] \
                        + "\";\n"

        gff_line_list.append(modified_line)"""

    for line in txt_file:

        line_t = line.strip().split()       # temporary list
        modified_line = line_t[0] + "\t" + "." + "\t" + "LowCov" + "\t" + line_t[1] + "\t" + line_t[2] + "\t" \
                        + ".\t.\t.\t" + "protein_hit=\"" + line_t[3] + "\";evalue=\"" + line_t[4] + "\";bitscore=\"" \
                        + line_t[5] + "\";identities=\"" + line_t[6] + "\";frameshift_count=\"" + line_t[7] \
                        + "\";assigned_gff_region_type=\"" + line_t[8] + "\";distance_to_next_CDS=\"" + line_t[9] \
                        + "\";\n"

        gff_line_list.append(modified_line)

    txt_file.close()

    # creating the gff file
    new_gff_file_path = os.getcwd() + "/output_gff_file.gff"
    new_gff_file = open(new_gff_file_path, "w")

    for line in gff_line_list:
        new_gff_file.write(line)

    new_gff_file.close()

    return None


def main():
    print("Output Creation main executed")


if __name__ == '__main__':
    main()
