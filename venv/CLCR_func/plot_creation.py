#!python3

"""This file contains the functions for the plot creation"""
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap


def low_cov_length_distribution_plot(input_scaffold_list, output_path):

    """
    This functions simply calculates the length distribution of all low coverage regions that are handed over. And
    creates an boxplot out of it, that is saved in the handed over directory.
    :param input_scaffold_list: output of one of the detect_/ or merge_close_regions function
    :param output_path: complete path to location where the plot should be saved, INCLUDING the plotname and .png
    :return: length_distribution list and  creates a boxplot at the handed over directory
    """

    # Labels for bars
    tick_labels = ["1-5", "6-25", "25-50", "51-100", "101-250", "250-500", "501-1000", "1001-5000", ">5000"]

    length_distribution = [[5, 0], [25, 0], [50, 0], [100, 0], [250, 0], [500, 0], [1000, 0], [5000, 0],
                           [float("inf"), 0]]

    # Calculating the length distribution
    for scaffold_data in input_scaffold_list:
        for low_cov_region in scaffold_data[1]:
            region_length = int(low_cov_region[1]) - int(low_cov_region[0])

            for length_range in length_distribution:
                if region_length <= length_range[0]:
                    length_range[1] += 1
                    break

    # X-coordinates of left sides of bars
    elements = [x for x in range(1, len(length_distribution) + 1)]

    # Heights of bars
    element_counts = [y[1] for y in length_distribution]

    # Plotting a bar chart
    plt.barh(elements, element_counts, tick_label=tick_labels, height=0.8, color=["limegreen"], edgecolor="black")

    # Label the picture
    plt.xlabel("Number of low coverage regions", fontsize=11.5)
    plt.ylabel("Length in bp", fontsize=11.5)
    plt.title("Low coverage regions length distribution", fontsize=13, fontweight="bold")
    plt.xticks(ticks=[20000, 40000, 60000, 80000, 100000, 120000, 140000])
    plt.xlim(xmax=145000)

    plt.tight_layout()

    # Function saves the plot
    plt.savefig(output_path, dpi=250)

    return length_distribution


def considered_diamond_hit_length_distribution_plot(considered_diamond_hits_list, output_path,
                                                    original_len_distribution):
    """
    This functions simply calculates the length distribution of all low coverage regions that are used  in the healing
    process. And creates an boxplot out of it, that is saved in the handed over directory.
    The displayed percentage in the bars, is the percentage of the original created queries in that length range, that
    were later considered in the healing process/ healed.
    :param considered_diamond_hits_list: output of one of the filter_out_relevant_results function
    :param output_path: complete path to location where the plot should be saved, INCLUDING the plotname and .png
    :param original_len_distribution: length distribution of the low cov. regions that were used as queries
    :return: length_distribution list and  creates a boxplot in the handed over directory
    """

    # Labels for bars
    tick_labels = ["1-5", "6-25", "25-50", "51-100", "101-250", "250-500", "501-1000", "1001-5000", ">5000"]

    length_distribution = [[5, 0], [25, 0], [50, 0], [100, 0], [250, 0], [500, 0], [1000, 0], [5000, 0],
                           [float("inf"), 0]]

    print(considered_diamond_hits_list[0])

    # List which saves the considered low cov regions, to prevent double counting, when more than one protein hit per
    # query was used in the healing step
    considered_low_cov_regions = []

    # Calculating the length distribution
    for protein_hit in considered_diamond_hits_list:
        low_cov_region_pos = [int(protein_hit[1]), int(protein_hit[2])]
        if low_cov_region_pos not in considered_low_cov_regions:
            considered_low_cov_regions.append(low_cov_region_pos)
            region_length = int(low_cov_region_pos[1]) - int(low_cov_region_pos[0])
            for length_range in length_distribution:
                if region_length <= length_range[0]:
                    length_range[1] += 1
                    break

    # Calculate the percentages of how many of the original low cov regions were then considered in the healing
    # regarding to their length class
    percentage_distribution = []

    for position in range(len(length_distribution)):
        considered_percentage = length_distribution[position][1] / (original_len_distribution[position][1]/100)
        percentage_distribution.append(str(round(considered_percentage, 2)) + "%")

    # X-coordinates of left sides of bars
    elements = [x for x in range(1, len(length_distribution) + 1)]

    # Heights of bars
    element_counts = [y[1] for y in length_distribution]

    # Plotting a bar chart
    plt.barh(elements, element_counts, tick_label=tick_labels, height=0.8, color=["yellow"], edgecolor="black")

    percentagelabel_position = max(element_counts) * 0.035

    # Label the bars
    for index, data in enumerate(element_counts):
        plt.text(y=index + 0.9, x=((data/2) - percentagelabel_position), s=f"{percentage_distribution[index]}",
                 fontdict=dict(fontsize=8))

    # Label the picture
    plt.xlabel("Number of low coverage regions", fontsize=11.5)
    plt.ylabel("Length in bp", fontsize=11.5)
    plt.title("Healed low coverage regions length distribution", fontsize=13, fontweight="bold")

    plt.tight_layout()

    print("figure saved")
    # Function saves the plot
    plt.savefig(output_path, dpi=250)

    return length_distribution


def create_toga_result_plot(results_dataframe):
    """
    Simply creates diagram, that displays the "gene-status flow" between the original and the healed assembly TOGA run
    :return:
    """

    print(results_dataframe)

    sns.set_style("whitegrid")

    status_type_list = ["I", "PI", "UL", "L", "M", "PM", "PG"]

    manual_noob_list_total = [60, 2, -61, -8, 9, -2, 0]

    plot_dataframe = pd.DataFrame(data=manual_noob_list_total, index=status_type_list)
    print(plot_dataframe)

    custom_colours = ["forestgreen", "seagreen", "goldenrod", "darkorange", "indianred", "firebrick", "grey"]

    ax_total = plot_dataframe.plot.bar(color=custom_colours, rot=0, figsize=(8, 5))
    ax_total.set_ylim(-130, 70)
    ax_total.legend(loc='upper right')
    ax_total.set_ylabel("amount of genes")
    ax_total.set_title("Total gene status shift")

    figure_total = ax_total.get_figure()
    figure_total.savefig("/home/johannes/Desktop/Toga_total.png", dpi=300)

    manual_noob_list_loss = [[0, 0, -56, -4, -1, 0, 0],
                              [0, 0, -2, 0, 0, 0, 0],
                              [-111, -4, 0, -25, 0, 0, 0],
                              [-10, 0, -21, 0, -9, 0, 0],
                              [0, 0, 0, -1, 0, 0, 0],
                              [0, 0, 0, -2, 0, 0, 0],
                              [0, 0, 0, 0, 0, 0, 0]]

    plot_dataframe = pd.DataFrame(data=manual_noob_list_loss, columns=status_type_list, index=status_type_list)
    print(plot_dataframe)

    custom_colours = ["forestgreen", "seagreen", "goldenrod", "darkorange", "indianred", "firebrick", "grey"]

    ax_total = plot_dataframe.plot.bar(stacked=True, color=custom_colours, rot=0, figsize=(8, 5))
    ax_total.set_ylim(-150, 0)
    ax_total.legend(loc='lower right')
    ax_total.set_ylabel("negative shift")
    ax_total.set_title("Gene status negative shift")

    figure_total = ax_total.get_figure()
    figure_total.savefig("/home/johannes/Desktop/Toga_loss.png", dpi=400)

    manual_noob_list_gain = [[0, 0, 111, 10, 0, 0, 0],
                             [0, 0, 4, 0, 0, 0, 0],
                             [56, 2, 0, 21, 0, 0, 0],
                             [4, 0, 25, 0, 1, 2, 0],
                             [1, 0, 0, 9, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0]]

    plot_dataframe = pd.DataFrame(data=manual_noob_list_gain, columns=status_type_list, index=status_type_list)
    print(plot_dataframe)

    custom_colours = ["forestgreen", "seagreen", "goldenrod", "darkorange", "indianred", "firebrick", "grey"]

    ax_total = plot_dataframe.plot.bar(stacked=True, color=custom_colours, rot=0, figsize=(8, 5), legend=False)
    ax_total.set_ylim(0, 150)
    ax_total.set_ylabel("positive shift")
    ax_total.set_title("Gene status flow")

    figure_total = ax_total.get_figure()
    figure_total.savefig("/home/johannes/Desktop/Toga_gain.png", dpi=400)

    return None


def main():
    print("Plot creation main executed")


if __name__ == '__main__':
    main()
