# coding=utf-8
# !python3

"""This file ... is the main file :) """
__author__ = "6947325: Johannes Zieres"
__credits__ = ""
__email__ = "johannes.zieres@gmail.com"


def main(input_args):
    print("CLCR MAIN CALLED")

    # Old example call

    """# Query creation
    project_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_01.02.2022/"
    cov_file_path = "/share/gluster/NOTLOESUNG/freya/T_draco/t_draco.pbc.wgs_short.txt"
    assembly_file = "/share/gluster/assemblies/Tdraco/" \
                    "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"
    low_cov_start = 15
    low_cov_end = 18
    min_query_len = 500
    queries_per_file = 5000

    create_queries(project_dir, cov_file_path, assembly_file, low_cov_start, low_cov_end, min_query_len,
                   queries_per_file)"""

    """# Create and submit slurm array job
    project_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_01.02.2022/"
    protein_database = "/home/johannes/Desktop/trachinus_draco/protein_db/protein_db.dmnd"
    auto_run = False
    prepare_slurm_run(project_dir, protein_database, auto_run)"""

    """# Create healed assembly file
    project_dir = "/home/johannes/Desktop/trachinus_draco/healing_runs/TRAdr_healing_run_01.02.2022/"
    unhealed_assembly = "/share/gluster/assemblies/Tdraco/" \
                        "t_draco_pacbio_salsa.FINAL_gap_closed.scaff_seqs_FINAL_pilon_2.fasta"
    dynamic_threshold_dist = 10

    create_healed_assembly(project_dir, unhealed_assembly, dynamic_threshold_dist)"""




