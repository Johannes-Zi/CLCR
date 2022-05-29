# CLCR

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## The CLCR approach
To improve the completeness of draft genome assemblies we created the tool CLCR, which stands for
Correction of Low Coverage Regions. The tool provides functions for short indel detection
and correction in high coverage genome assemblies, in the form of error detection in typically 
error prone assembly regions with a significantly lower read coverage. The detection
of frameshifts is implemented by the integration of the alignment tool Diamond in combination with 
the usage of a protein database with the sequences of closely related organisms.
The frameshift detection by Diamond is followed by the localization of the frameshifts in the
underlying genome in combination with extensive filter heuristics to avoid false correction.
Reading frames that are disturbed by frameshift mutations, putatively caused by a sequencing error,
are locally healed by the insertion of N’s in the assembly.

## Installation
CLCR can be retrieved as a PyPI package, it requires a Python version of 3.7 or higher and a 
[Diamond](https://github.com/bbuchfink/diamond) installation is mandatory. 

The program can be installed with the following command: <br />
```
# Install pip, if nescessary
sudo apt update
sudo apt install python3-pip

# Update pip
pip install --upgrade pip

# Install CLCR
pip install clcr
```
If there are problems with that, the code can be directly cloned from github:
```
# Install git, if nescessary
sudo apt install git-all

# Clone CLCR code from github
git clone https://github.com/Johannes-Zi/CLCR.git

# Move to the cloned directory
pip install .
# Or manually add the program to PATH when there are problems with pip
```
## Quickstart
### Query creation with *clcr.query_creation*
The first step is the creation of the query sequences for the Diamond 
blastx runs. For this, the function detects the regions with low read coverages in a given pbc (per base coverage) file. The
sequences of the detected regions are extracted from the handed over assembly and stored as .fasta query
files in the query_files directory of the handed over project. The query_files directory will be
overwritten! A log file with run information and a original_low_cov_regions.tsv with the original detected
low coverage regions before the merging step are stored at the storage_files dir. <br /> 
<br />
*ATTENTION!*: old query directory and DIAMOND output directory of current project dir is overwritten as preparation for a new cluster run!
<br /> 
<br /> 
There are three mandatory arguments required:
* genome assembly
* the matching per base coverage file
* project directory for the CLCR run

There are the following options:
```
  -h, --help            show the help message and exit

required arguments:
  -p PROJECT_DIR, --project_dir PROJECT_DIR
                        Path of the project directory
                        
  -c COV_FILE_PATH, --cov_file_path COV_FILE_PATH
                        Path of the coverage file
                        
  -a ASSEMBLY_FILE, --assembly_file ASSEMBLY_FILE
                        Path of the assembly file

optional arguments:
  --low_cov_start LOW_COV_START
                        Threshold for detecting a low cov region
                        
  --low_cov_end LOW_COV_END
                        Threshold for ending a low cov region
                        
  --min_query_len MIN_QUERY_LEN
                        Minimum query length
                        
  --queries_per_file QUERIES_PER_FILE
                        Queries sequences per query file
                        
  --verbose             Run information is print in the command line
```

### Diamond cluster run *clcr.cluster_run*
After the query files are created, the next step is to perform the
Diamond blastx searches against a protein database with the sequences 
of closely related organisms. This can be done locally on a single computer, 
or on a computer cluster. The usage of a cluster is highly recommended, the CLCR 
workflow has to be manually adapted to this. The *clcr.cluster_run* function creates a slurm-file for the Diamond blastx 
cluster run of the handed over CLCR project. 
The jobs are started automatically on the cluster, when the --auto_run parameter is activated. 
The slurm log files are stored in the slurm_files dir, a log file for the CLCR run is stored in 
the storage_files dir.<br />
<br />
*ATTENTION!*: Use this function only when your cluster supports slurm and adapt the slurm file manually 
to your local circumstances! (In this case auto submission is not recommended!)
<br />
<br />
There are the following options:
```
  -h, --help            show this help message and exit

required arguments:
  -p PROJECT_DIR, --project_dir PROJECT_DIR
                        Path of the project directory
                        
  -c PROTEIN_DATABASE, --protein_database PROTEIN_DATABASE
                        Path of the protein database

optional arguments:
  --auto_run [AUTO_RUN]
                        Activate automatic slurm job submission, when parameter is True.
  
  --verbose             Run information is print in the command line
```

### Creation of healed assembly version *clcr.assembly_healing*
The last step in the analysis is the creation of an adapted assembly
version. For this the detected frameshifts in the Diamond
blastx output are evaluated, extensively filtered and used to created a adapted assembly version with
locally healed reading frames. The healed assembly version is stored in the healed_assembly dir, and log
file for the CLCR run is stored in the storage_files dir.
<br />
<br />
There are the following options:
```
  -h, --help            show this help message and exit

required arguments:
  -p PROJECT_DIR, --project_dir PROJECT_DIR
                        Path of the project directory
                        
  -c UNHEALED_ASSEMBLY, --unhealed_assembly UNHEALED_ASSEMBLY
                        Path of the original unhealed assembly file

optional arguments:
  --dynamic_threshold_dist DYNAMIC_THRESHOLD_DIST
                        The max_detect_distance defines the distance from a detected frameshift position to
                        the original low cov. region, where a frameshift is still considered and not
                        excluded in the further analysis.
                        
  --verbose             Run information is print in the command line
```

## Contact
For error reports, questions, comments or suggestions contact us via [email](mailto:johannes.zieres@gmail.com)
