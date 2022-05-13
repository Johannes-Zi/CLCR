#!python3

"""This file contains the functions for the Toga analysis"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import pandas as pd
import glob


def read_in_toga_lossgene_file(unhealed_file_path, healed_file_path, isoforms_file_path, query_annotation_file_path,
                               putative_false_corrected_tsv_path):
    """
    Small function for evaluating the TOGA loss_summ_data.tsv files. It filters out the genes that are marked in the
    Toga run with the unmodified/unhealed as intact and in the healed assembly as lost or missing.
    For each of those genes all isoforms and all location of each isoform is detected and stored with the relating
    position data in the saved output file.
    :param unhealed_file_path: Toga loss_summ_data.tsv file path of unhealed assembly
    :param healed_file_path: Toga loss_summ_data.tsv file path of healed assembly
    :param isoforms_file_path: Toga output file
    :param query_annotation_file_path: Toga annotation file path
    :param putative_false_corrected_tsv_path: path to the output file that will be created
    :return: results_dataframe
    """

    # Read in the tsv files with different separator
    column_names_1 = ["type", "gene_id", "status"]
    dataframe_unhealed = pd.read_csv(unhealed_file_path, sep='\t', header=None, names=column_names_1,
                                     index_col="gene_id")
    dataframe_healed = pd.read_csv(healed_file_path, sep='\t', header=None, names=column_names_1, index_col="gene_id")

    gene_df_unhealed = dataframe_unhealed[dataframe_unhealed['type'] == 'GENE']

    gene_df_healed = dataframe_healed[dataframe_healed['type'] == 'GENE']

    # Create an combined dataframe with the information of the healed and unhealed file
    combined_df = pd.merge(gene_df_unhealed, gene_df_healed, how="inner",
                           suffixes=("_unhealed", "_healed"), on="gene_id").drop(columns=["type_healed", "type_unhealed"])

    # Search for the gene entries, where the gene was intact in the unhealed assembly and lost, uncertain lost or
    # missing in the healed assembly

    # Check putative wrongly healed genes
    putative_false_corrected_df = combined_df[(combined_df['status_unhealed'] == 'I') &
                                              ((combined_df['status_healed'] == 'UL') |
                                               (combined_df['status_healed'] == 'L') |
                                               (combined_df['status_healed'] == 'M'))]

    """# Check putative correctly healed genes
    putative_false_corrected_df = combined_df[((combined_df['status_unhealed'] == 'L') |
                                              (combined_df['status_unhealed'] == 'UL')) &
                                              (combined_df['status_healed'] == 'I')]"""

    # Read in the isoforms tsv file
    gene_isoforms_df = pd.read_csv(isoforms_file_path, sep='\t', index_col="gene_id").rename(
        columns={"trascript_id": "transcript_id"})

    # Print value count before first inner merge (original amount of putative wrong healed genes)
    print(putative_false_corrected_df.value_counts())

    # Not necessary, just for checks
    """# Saves the data for later checks
    putative_false_corrected_df.to_csv("/home/johannes/Desktop/trachinus_draco/TOGA_run_1_output/"
                                       "original_putative_wrong_healed.tsv", sep='\t')"""

    # Includes the informations of the isoform-identifieres for each gene
    putative_false_corrected_df = pd.merge(putative_false_corrected_df, gene_isoforms_df, how="inner", on="gene_id")

    # Read in the complete annotation of the assembly
    column_names_2 = ["scaffold", "2", "3", "start_pos", "end_pos", "6", "7", "8", "9"]
    query_annotation_df = pd.read_csv(query_annotation_file_path, sep='\t', header=None, names=column_names_2).drop(
        columns=["2", "6", "7", "8"])

    # Filter out the entries that are not "transcript" in the 3.rd column
    query_annotation_df = query_annotation_df[query_annotation_df["3"] == "transcript"].drop(columns=["3"])
    # Extract the isoform identifier out of the last column and add those as new column
    query_annotation_df["transcript_id"] = query_annotation_df["9"].str.split("\"").str[1]
    # Removes the projection-location specific identifier number at the end of each isoform identifier
    query_annotation_df["transcript_id"] = query_annotation_df["transcript_id"].map(lambda x: x.rstrip("0123456789"
                                                                                                       ).rstrip("."))
    query_annotation_df = query_annotation_df.drop(columns=["9"]).set_index("transcript_id")

    # Print before second inner merge
    print(putative_false_corrected_df)

    # Inner merge to add the positions of each projection for each gene in the assembly
    putative_false_corrected_df = pd.merge(putative_false_corrected_df.reset_index(), query_annotation_df, how="inner",
                                           on="transcript_id")

    # Print final dataframe
    print(putative_false_corrected_df)

    # Print how many of the original regions are left after merging - some could get lost if there are no entries in the
    # query_annotation file
    print("Genes with transcrips in the annotation file:",
          putative_false_corrected_df.groupby("gene_id")["transcript_id"].nunique().count())

    putative_false_corrected_df.to_csv(putative_false_corrected_tsv_path, sep='\t')

    # Result dataframe with overall results
    results_dataframe = combined_df.value_counts()

    return results_dataframe


def check_overlapping_healing_positions_1(putative_wrong_corrected_file_path,
                                          healing_region_list):
    """

    Note! The best hit of the created queries in the database is not in all cases the hit on which the healing decision
    was based, thus the new version 2 of this function was created (searching for the exact hit)

    This function checks each gene that is marked by TOGA in the unhealed assembly as intact and in healed version
    as lost/missing. This means for all locations of all isoforms of all genes is checked, if they are overlapping with
    a healing position, thus false healing might have happened in this region. For each of those detected overlapping,
    a query is created +- 400 Bp around the healing position. Healing positions that are represented multiple times by
    different queries, are filtered out. The start and end position of each query is afterwards appended to the
    filtered_diamond_query_list, which is afterwards returned for the query file creation.
    :param putative_wrong_corrected_file_path: path to the processed toga output of read_in_toga_lossgene_file()
    :param healing_region_list: output of filter_out_relevant_results()
    :return: filtered_diamond_query_list
    """
    # Read in the file
    input_file = open(putative_wrong_corrected_file_path, "r")
    # Skipps the first line
    next(input_file)

    # Stores the information for the query creation (input of the query_files_creation() function)
    diamond_query_list = []

    # Read in current potential wrong corrected region
    for line in input_file:
        splitted_line = line.split()

        current_region = splitted_line[1]
        isoform_id = splitted_line[4]
        current_scaffold = splitted_line[5]
        start_pos_in_scaff = min(int(splitted_line[6]), int(splitted_line[7]))
        end_pos_in_scaff = max(int(splitted_line[6]), int(splitted_line[7]))

        overlapping_found = False

        # Search if there is a overlapping correction position
        for query in healing_region_list:
            # Case for same scaffold
            if current_scaffold == query[0]:
                # Checks for each healing positions, if its overlapping with the current region
                for healing_pos in query[3]:
                    healing_pos_in_scaff = int(query[1]) + int(healing_pos[0])

                    # Case, where healing in the putative wrong healed region was detected
                    if (healing_pos_in_scaff >= start_pos_in_scaff) and (healing_pos_in_scaff <= end_pos_in_scaff):

                        query_start = max(healing_pos_in_scaff - 400, 0)
                        query_end = healing_pos_in_scaff + 400

                        query_header = current_region + "#" + isoform_id + "#" + current_scaffold + "#" + \
                                       str(start_pos_in_scaff) + "#" + str(end_pos_in_scaff) + "#" + \
                                       str(query_start) + "#" + str(query_end)
                        diamond_query_list.append([current_scaffold, [(query_start, query_end,
                                                                       query_header)]])

                        overlapping_found = True
                        break

            if overlapping_found:
                break

    # Contains the filtered output
    filtered_diamond_query_list = []

    # Removes duplicates, when for a single loci multiple isoforms are predicted, thus multiple queries of the same loki
    # would be created, what makes no sense
    loci_list = []
    for query in diamond_query_list:
        already_added = False
        start_pos = query[1][0][0]
        end_pos = query[1][0][1]
        # Checks if the region of the current query is already represented by another query
        for already_added_loci in loci_list:
            if (already_added_loci[0] == start_pos) and (already_added_loci[1] == end_pos):
                already_added = True
                break

        # Appends the query to the loci list
        if not already_added:
            loci_list.append([start_pos, end_pos])
            filtered_diamond_query_list.append(query)

    for x in filtered_diamond_query_list:
        print(x)

    print(len(filtered_diamond_query_list))

    return filtered_diamond_query_list


def check_overlapping_healing_positions_2(putative_wrong_corrected_file_path,
                                          healing_data_tsv_path, output_tsv_path):
    """
    Improved version, that uses the healing_data.tsv file in which for each healing position the underlying Diamond hit
    is saved.

    This function checks each gene that is marked by TOGA in the unhealed assembly as intact and in healed version
    as lost/missing. This means for all locations of all isoforms of all genes is checked, if they are overlapping with
    a healing position in the .tsv, thus false healing might have happened in this region.

    For each of those detected overlapping, a instance is saved in the output file, double representation of healing
    positions is prevented. The output .tsv could then be used to find the healing related Diamond hits in the output
    file.

    :param putative_wrong_corrected_file_path: path to the processed toga output of read_in_toga_lossgene_file()
    :param healing_data_tsv_path: output of create_detailed_healing_information_file()
    :param output_tsv_path: path to the output file, including the file name
    :return: None
    """

    # Read in the healing_data.tsv file
    healing_tsv = open(healing_data_tsv_path, "r")
    healing_postion_data = []    # Stores the data
    # Read in the lines
    for line in healing_tsv:
        healing_postion_data.append(line.strip().split())
    healing_tsv.close()

    # Read in the file with the putative wrong corrected positions
    input_file = open(putative_wrong_corrected_file_path, "r")
    # Skipps the first line
    next(input_file)

    # Stores the data for the ourput
    output_tsv_data = []

    # Read in current potential wrong corrected region
    for isoform_location_data in input_file:
        split_line = isoform_location_data.split()

        current_region = split_line[1]
        isoform_id = split_line[4]
        current_scaffold = split_line[5]
        start_pos_in_scaff = min(int(split_line[6]), int(split_line[7]))
        end_pos_in_scaff = max(int(split_line[6]), int(split_line[7]))

        # Search if there is a overlapping correction position
        for healing_postion in healing_postion_data:
            # Case for same scaffold
            if current_scaffold == healing_postion[0]:

                healing_pos_in_scaff = int(healing_postion[1])

                # Case, where healing in the putative wrong healed region was detected
                if (healing_pos_in_scaff >= start_pos_in_scaff) and (healing_pos_in_scaff <= end_pos_in_scaff):

                    # Check if the found Diamond hit is already represented
                    protein_hit_already_considered = False

                    for entry in output_tsv_data:

                        if healing_postion[7] == entry[11]:
                            protein_hit_already_considered = True
                            break

                    if not protein_hit_already_considered:
                        output_tsv_data.append(split_line[1:5] + healing_postion)

                    # Breaks, thus no more entries with different healing positions in the same query are saved
                    break

    input_file.close()

    # Create the output_file
    output_tsv_file = open(output_tsv_path, "w")

    for data_point in output_tsv_data:
        output_tsv_file.write("\t".join(data_point + ["\n"]))

    output_tsv_file.close()

    return None


def search_relating_diamond_alignments(tsv_file_path, diamond_output_dir,  output_file_path):
    """
    This function searches for given entries in the tsv file the matching Diamond alignments, and save the together in
    the several output files for easier checking . The first file contains all queries with their top hit, the other
    files contain one query each, with all corresponding protein hits.
    :param tsv_file_path: path to a output file of check_overlapping_healing_positions_2
    :param diamond_output_dir: path to the dir with the diamond output files
    :param output_file_path: path to the output directory,
    :return: None
    """

    # Read in the .tsv file
    tsv_file = open(tsv_file_path, "r")
    tsv_data = []  # Stores the data
    # Read in the lines
    for line in tsv_file:
        tsv_data.append(line.strip().split())
    tsv_file.close()

    # Returns every file in the directory with .out at the end
    output_file_list = glob.glob(diamond_output_dir + "/temp_out_*.txt")

    # Sorting the file path ascending by their file number
    dir_path_len = len(diamond_output_dir) + 9  # +9 because temp_out_ consists of 9 chars
    output_file_list = sorted(output_file_list, key=lambda current_path: int(current_path[dir_path_len:-4]))

    hit_count = 0   # Saves how many hits were found

    # Saves all results in form of the relevant lines in the file, contains one sublist for each Query, which then also
    # contains sublists with the line of each protein hit of that query
    found_results_storage = []

    # Iterate threw each datapoint
    for tsv_datapoint in tsv_data:

        # Scaffold, start pos. and end pos. of the currently searched query
        query_location = [tsv_datapoint[4]] + tsv_datapoint[7:9]
        # Name of the search protein hit
        protein_hit = tsv_datapoint[11]

        break_value = False     # Gets True when the extracted data was found and no new file needs to be read in

        # List which contains a sublist for each protein hit with the corresponding lines, first elem is query name
        query_diamond_hits = []
        protein_hit_lines = []      # Saves all lines of the current protein hit

        # Search for the corresponding Diamond hit/alignment
        for output_file in output_file_list:
            # Saves when the searched query was found, thus only the following hits are viewed
            searched_query_found = False
            current_output_file = open(output_file, "r")

            # Read the complete file, breaks after all hits of the searched query are considered
            for line in current_output_file:

                if searched_query_found:

                    # Case that the searched protein hit was found
                    if line.startswith(">"):
                        query_diamond_hits.append(protein_hit_lines)
                        protein_hit_lines = [line]      # Reset the list

                    # Case that the complete results of the current query are read in, thus the current loop is finished
                    elif line.startswith("Query= "):
                        found_results_storage.append(query_diamond_hits)
                        hit_count += 1
                        break_value = True
                        break

                    else:
                        protein_hit_lines.append(line)

                # Case that the matching query data was found
                if line.startswith("Query= "):
                    split_line = line[7:].strip().split("#")
                    current_location = [split_line[0], split_line[3], split_line[4]]
                    if current_location == query_location:
                        searched_query_found = True

                        # Append the query information line
                        protein_hit_lines = protein_hit_lines + ["\n", line]

            current_output_file.close()

            if break_value:
                break

    # Write the collected lines in the output .txt file
    output_txt_file = open((output_file_path + "queries_result_summary.txt"), "w")

    # Creates the first output file with the corresponding first protein hit alignment of each found query
    for query_data in found_results_storage:

        # Writes down first element (Query header) and second element (first protein hit)
        for protein_hit in range(min(2, len(query_data))):      # min() for queries without a protein hit

            for line in query_data[protein_hit]:
                output_txt_file.write(line)

    output_txt_file.close()

    # Create the other result files, with all protein hits of one query in each file
    name_count = 1
    for query_data in found_results_storage:
        output_file = open((output_file_path + "query_results_" + str(name_count)) + ".txt", "w")

        for diamond_hit_data in query_data:

            for line in diamond_hit_data:
                output_file.write(line)

        name_count += 1
        output_file.close()

    print("Found Diamond hits: ", hit_count)

    return None


def main():
    print("Toga analysis main executed!")

    """# Evaluate the TOGA results

    HLtraDra1_file_path = "/home/johannes/Desktop/trachinus_draco/toga_run_by_Michael_Hiller/" \
                          "loss_summ_data_HLtraDra1.tsv"
    HLtraDra3_file_path = "/home/johannes/Desktop/trachinus_draco/toga_run_by_Michael_Hiller/" \
                          "loss_summ_data_HLtraDra3.tsv"
    toga_isoforms_tsv = "/home/johannes/Desktop/trachinus_draco/toga_run_by_Michael_Hiller/toga.isoforms.tsv"
    query_annotation_gtf = "/home/johannes/Desktop/trachinus_draco/toga_run_by_Michael_Hiller/query_annotation.gtf"
    # putative_wrong_corrected_file_path = "/home/johannes/Desktop/trachinus_draco/toga_run_by_Michael_Hiller/" \
    #                                     "putative_false_corrected.tsv"
    putative_wrong_corrected_file_path = "/home/johannes/Desktop/trachinus_draco/toga_run_by_Michael_Hiller/" \
                                         "putative_rightly_corrected.tsv"

    results_dataframe = output_processing.read_in_toga_lossgene_file(HLtraDra1_file_path, HLtraDra3_file_path,
                                                                     toga_isoforms_tsv,
                                                                     query_annotation_gtf,
                                                                     putative_wrong_corrected_file_path)

    # outputP.create_toga_result_plot(results_dataframe)

    # search for overlapping
    healing_data_tsv_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                            "storage_files/healing_data.tsv"
    output_tsv_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                      "storage_files/toga_result_analysis.tsv"

    output_processing.check_overlapping_healing_positions_2(putative_wrong_corrected_file_path, healing_data_tsv_path,
                                                            output_tsv_path)

    diamond_output_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                         "output_files/"
    # output_file_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
    #                   "storage_files/putative_wrong_corrected_diamond_hits/"
    output_file_path = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_10.01.2022/" \
                       "storage_files/putative_right_corrected_diamond_hits/"
    output_processing.search_relating_diamond_alignments(output_tsv_path, diamond_output_dir, output_file_path)

    # """

if __name__ == '__main__':
    main()
