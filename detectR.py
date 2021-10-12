"""This file contains the functions for the detection of low coverage regions in the pbc file of an assembly"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"


def detect_regions(cov_file_path, cov_start, cov_end):
    """
    Detecting the low coverage regions in the per base coverage file. Saves for each low coverage region the scaffold
    and the start and end position in the scaffold. Low cov. regions that are spanning two or more scaffolds are
    excluded. Returns the low cov regions as python list positions(-1). A low cov. region is detected when the cov. is
    equal or lower than cov_start. And a low cov. region is ended when the coverage is higher than cov_end.
    :param cov_file_path: path of input coverage file
    :param cov_start: threshold for detecting a low coverage area
    :param cov_end: threshold for leaving a low coverage area
    :return: list with the low coverage regions
    """

    import_file = open(cov_file_path)  # Opens the coverage file
    current_position = 0  # Saves the current position
    low_cov = False  # Saves the current coverage situation
    low_cov_regions = []  # Low coverage regions List with the current_scaffold list of each scaffold
    region_start = int  # Saves the start point for a low coverage region

    # Saves the scaffold of the first low cov position of the current low coverage region
    start_scaffold = str    # Parameter is for excluding low cov. regions that are spanning multiple scaffolds
    current_scaffold = ["scaffoldname", [(0, 10)]]    # Saves the name of the scaffold and all low cov. regions in them.

    # Detects the low coverage regions and saves them in low_coverage_regions list
    for line in import_file:
        splitted_line = line.split()

        # New scaffold starts, old scaffold information is appended to the final output list
        if splitted_line[0] != current_scaffold[0]:
            print("New scaffold: ", splitted_line[0])

            # Appends the old scaffold to the final output list if low cov. reg. are found in the finished scaffold
            if current_scaffold[1]:
                low_cov_regions.append(current_scaffold)
            current_scaffold = [splitted_line[0], []]   # Resets the list and initialise the List for the new scaffold
            current_position = 0    # Resets the current position to 0, because a new scaffold starts

        if low_cov:  # Case for being in a low coverage region
            if int(splitted_line[2]) > cov_end:     # Checks if the low cov. region ends
                low_cov = False

                # Checks that the start and end position of the low cov. regions lays in the same scaffold
                if current_scaffold[0] == start_scaffold:
                    current_scaffold[1].append((region_start, (current_position - 1)))

        else:  # Case for being in a good coverage region
            if int(splitted_line[2]) <= cov_start:      # Checks if a new low cov. region starts
                low_cov = True
                region_start = current_position     # Saves the start position of the current low cov. region
                start_scaffold = splitted_line[0]   # Saves the scaffold of the start position

        current_position += 1

    import_file.close()  # Closes the input file

    # Appends the last remaining scaffold if low coverage regions where detected in them
    if current_scaffold[1]:
        low_cov_regions.append(current_scaffold)

    low_cov_regions.pop(0)      # Removes the initialising object

    return low_cov_regions


def merge_close_reg(input_scaffold_list, merge_distance):

    """
    Receives the output list of the detect regions function and merges all low cov. regions in each scaffold together if
    the are closer together than the given merge_distance parameter.
    :param input_scaffold_list: input list, which is containing the scaffold lists with the low coverage regions
    :param merge_distance: distance threshold, for merging regions which are closer together than the threshold
    :return: region list with the merged regions
    """

    output_list = []    # Contains the modified scaffold lists with low cov. regions of each scaffold

    # Iterates trough the scaffolds
    for scaffold_list in input_scaffold_list:

        new_scaffold_list = [scaffold_list[0], []]      # Saves the new merged low cov. regions

        # Initialising the parameters with the first low cov region
        current_start = scaffold_list[1][0][0]
        current_end = scaffold_list[1][0][1]

        # Merging the regions of the current scaffold
        for low_cov_reg in scaffold_list[1]:

            # Case in which the start position of the current region ist closer to the end position of the previous
            # region than the merge distance, so those two regions must be merged
            if (low_cov_reg[0] - merge_distance) <= current_end:
                current_end = low_cov_reg[1]

            # Case where the current region don't have to be merged with the previous region, which means no following
            # region has to be merged with the previous either. So the previous region is complete and could be appended
            # to the scaffold list
            else:
                new_scaffold_list[1].append((current_start, current_end))

                # Reset the parameters and initialise them with the new low cov. region
                current_start = low_cov_reg[0]
                current_end = low_cov_reg[1]

        # Append the last region
        new_scaffold_list[1].append((current_start, current_end))

        # Append the final scaffold list to the output list
        output_list.append(new_scaffold_list)

    return output_list


def combine_regions(input_list_short, input_list_long):

    """
    Detecting combined low coverage regions, which means regions where the short read and long read coverage is below
    the thresholds. Has to be called for each scaffold.
    :param input_list_short: list containing the short read low coverage regions
    :param input_list_long: list containing the long read low coverage regions
    :return: list which is containing the regions where both coverages are low
    """

    combined_list = []      # For saving the created combined regions

    # Saving the current tuple positions in the lists
    short_list_pos = 0
    long_list_pos = 0

    combined_region = False    # Saves the new "combined regions", which are not already represented in the input lists

    # Set starting point
    if input_list_short[0][0] <= input_list_long[0][0]:
        startpoint = input_list_short[0][0]
        active_region = "S"     # for short
    else:
        startpoint = input_list_long[0][0]
        active_region = "L"  # for long

    # Goes gradually through the regions, for each active region (short or long) the while-loop considers if the region
    # could be expanded by a region from the other input_list. If a region is expanded the "combined_region" is set TRUE
    # and if the regions could not be extended more the starting point of the region and the current end is appended in
    # combined_list.
    while True:

        # case if the current "active/considered" region is represented by an short read region
        if active_region == "S":

            # break condition
            if (long_list_pos >= len(input_list_long)) or (short_list_pos >= len(input_list_short)):
                break

            # tries every potential region which could expand the current active region
            # while-loop represents every long read region which protrudes in the current active region
            while input_list_long[long_list_pos][0] <= input_list_short[short_list_pos][1]:

                # There are two different possible cases

                # 1. the considered region expands the current active region
                if input_list_long[long_list_pos][1] > input_list_short[short_list_pos][1]:
                    active_region = "L"
                    short_list_pos += 1
                    combined_region = True      # combine argument is only set TRUE if the region is indeed expanded
                    break

                # 2. the considered region is completely in the current active region and could be skipped
                else:
                    long_list_pos += 1
                    if long_list_pos >= len(input_list_long):
                        break

            # case if no region for expanding was found and the current active region represents an end
            if active_region == "S":

                if combined_region:     # the current region is a merged region and will be appended
                    combined_list.append((startpoint, input_list_short[short_list_pos][1]))
                    combined_region = False

                short_list_pos += 1
                if short_list_pos >= len(input_list_short):
                    break

                # setting a new starting point
                if input_list_short[short_list_pos][0] <= input_list_long[long_list_pos][0]:
                    startpoint = input_list_short[short_list_pos][0]
                    active_region = "S"  # for short
                else:
                    startpoint = input_list_long[long_list_pos][0]
                    active_region = "L"  # for long

        # case if the current "active/considered" region is represented by an long read region
        else:

            # break condition
            if (long_list_pos >= len(input_list_long)) or (short_list_pos >= len(input_list_short)):
                break

            # tries every potential region which could expand the current active region
            # while-loop represents every long read region which protrudes in the current active region
            while input_list_short[short_list_pos][0] <= input_list_long[long_list_pos][1]:

                # There are two different possible cases

                # 1. the considered region expands the current active region
                if input_list_short[short_list_pos][1] > input_list_long[long_list_pos][1]:
                    active_region = "S"
                    long_list_pos += 1
                    combined_region = True  # combine argument is only set TRUE if the region is indeed expanded
                    break

                # 2. the considered region is completely in the current active region and could be skipped
                else:
                    short_list_pos += 1
                    if short_list_pos >= len(input_list_short):
                        break

            # case if no region for expanding was found and the current active region represents an end
            if active_region == "L":

                if combined_region:     # the current region is a merged region and will be appended
                    combined_list.append((startpoint, input_list_long[long_list_pos][1]))
                    combined_region = False

                long_list_pos += 1
                if long_list_pos >= len(input_list_long):
                    break

                # setting a new starting point
                if input_list_short[short_list_pos][0] <= input_list_long[long_list_pos][0]:
                    startpoint = input_list_short[short_list_pos][0]
                    active_region = "S"  # for short
                else:
                    startpoint = input_list_long[long_list_pos][0]
                    active_region = "L"  # for long

    # append the last  combined region if the last combined region is at the end of one of the input lists
    if combined_region:
        if active_region == "L":
            combined_list.append((startpoint, input_list_long[long_list_pos][1]))
        else:
            combined_list.append((startpoint, input_list_short[short_list_pos][1]))

    return combined_list


def combine_regions_multiple_scaffolds(input_scaffold_list_short, input_scaffold_list_long):
    """
    Searches the matching scaffold_list pairs with the short and long read low coverage regions of each scaffold. Calls
    the combine region for each  scaffold pair.
    :param input_scaffold_list_short: output list of detect_regions or merge_close regions with the short read low cov.
                                      regions
    :param input_scaffold_list_long: output list of detect_regions or merge_close regions with the long read low cov.
                                      regions
    :return: Combined list with sublists that contain the combined low cov. regions of each scaffold
    """

    combined_list = []      # Contains the scaffold lists with all combined low cov. regions of each scaffold

    # Search the matching long read scaffold_list for each short read scaffold_list, if possible
    for scaffold_list_short in input_scaffold_list_short:

        for scaffold_list_long in input_scaffold_list_long:

            # Case where the matching scaffold was found
            if scaffold_list_short[0] == scaffold_list_long[0]:

                combined_list.append([scaffold_list_short[0],
                                      combine_regions(scaffold_list_short[1], scaffold_list_long[1])])
                break

    return combined_list


def low_cov_length_distribution():

    low_cov_length_distribution_list = []

    return low_cov_length_distribution_list


def main():
    print("bruteforce main executed")


if __name__ == '__main__':
    main()
