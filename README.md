——————————————————————————————————————————
VIPER - Visualization Pipeline for RNA-seq
——————————————————————————————————————————

—————
Intro
—————
ViPeR is a comprehensive RNA-seq analysis tool built using snakemake which allows for ease of use and customization. We combine the use of etc……..

In ViPeR, there will be three distinct components that are kept purposefully separate to avoid confusion, accidental deletion or editing of essential components. The three components are the DATA, VIPER, and ANALYSIS. DATA is where the user will store all of their data for analysis. VIPER is the actual source of the code, you will never have to edit anything in here. ANALYSIS will be the output of VIPER

Although included in this README are step-by-step instructions, it is assumed that the user should have a basic understanding of a command line interface on mac or linux, 


Input:
fastqs or bam files
config file - details below
meta file - details below

Output:
Mapped STAR files
…
Sample-Sample Clustering
Sample-Feature Clustering
PCA


———————————————
Getting Started
———————————————
Before ViPeR can be used, all of the necessary packages and tools must be downloaded and installed.

Installing Python2 and Python3
——————————————————————————————
https://www.python.org/downloads/


Installing wget
———————————————
To get these packages, we will use a command line tool called wget (http://www.gnu.org/software/wget/). wget is a popular tool for pulling things off of the internet, and although it should be already be installed on most linux machines, you may need to install is if you are using a Mac (http://rudix.org/packages/wget.html)


Installing Miniconda
————————————————————
We will be using Miniconda to manage most of the packages that go into ViPeR. Miniconda3 (http://conda.pydata.org/miniconda.html) is a simple tool that include the package manager Conda (http://conda.pydata.org/docs/) Conda is used to manage all of the tools used in ViPeR including version control

Navigate to a folder where you keep your coding tools or wherever you want to store the main Miniconda3 folder. Then enter the following commands:

If you are using Linux:
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh
If you are using Mac:
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	bash Miniconda3-latest-MacOSX-x86_64.sh

-Make sure to answer yes when asked if you want to prepend Miniconda3 to PATH.
-Close your terminal, open a new one and you should now have Conda working! Test by entering:
	conda update conda


Installing VIPER and setting up your environment
————————————————————————————————————————————————
1) Open a terminal and navigate to the location of the PROJECT’s data (path/to/PROJECT means the full computer path of where your project is, ex. usr/home/project)
	cd /path/to/PROJECT/
	- If you do not have a folder for your PROJECT we recommend you make one
		cd /path/to
		mkdir PROJECT
		cd PROJECT

#2) The PROJECT folder for now should only have your data in it. Next we will acquire VIPER using git
#	within your PROJECT folder you want to initialize a git repository
#		git init
#	after git initialization and set up, you will need to clone the VIPER folder in to your PROJECT folder
#		git clone git@bitbucket.org:vangalamaheshh/snakemake.git

############### THIS NEEDS TO BE IN A PUBLIC FOLDER WITH THE PROPER NAMING, need to talk about this
2) The PROJECT folder for now should only have your data in it. Next we will acquire VIPER using wget
	wget https://bitbucket.org/vangalamaheshh/snakemake/get/master.tar.gz
	tar -xf master.tar.gz
	rm master.tar.gz

3) Download all of the necessary packages. Most of the tools needed can be downloaded via Miniconda. Navigate into the snakemake folder and run the bash script to download all of the packages
	cd VIPER
	bash build_env.bash

4) Copy the config file from within your VIPER folder into your PROJECT folder. Go back to your PROJECT folder to enter:
	cp VIPER/config.yaml ./
Edit your config file, more on this in the config section

5) Create your metadata file and move it into your PROJECT folder
More on this in the meta section

6) At this point, you should have your /path/to/PROJECT folder with the following contained within it
	VIPER
	DATA
	config.yaml
	metadata.csv


Running VIPER
—————————————
If all of the above have been followed correctly, you should be ready to run VIPER! snakemake is an extensive and powerful tool with many features that the user should learn, but for now, you will only need a few simple commands

In your PROJECT folder run the following command to see if you are ready to go:
	snakemake --snakefile snakemake/RNAseq.snakefile -n
This will return a large output which basically outlines what VIPER is about to do. If no errors come back, then you are ready to go! Type in the following command and enter
	snakemake --snakefile snakemake/RNAseq.snakefile



——————— 
CONFIG:
———————
Your config file has three main sections. PATHS, PARAMS, SAMPLES:

PATHS:
——————
In this section, the user will need to specify the location of the following packages and tools. Below is an example path and a description of what each item is.

#####################NEED EXPLANATIONS HERE!!!!
bed_file: /data/static_libraries/RefGene/refseqGenesHg19.bed
	MISSING

genome_lib_dir: /home/lentaing/tmp/weinstock/newrun/ref_files/Hg19_CTAT_resource_lib
	MISSING

gtf_file: /home/lentaing/tmp/weinstock/newrun/ref_files/Hg19_CTAT_resource_lib/ref_annot.gtf
	MISSING

metasheet: metasheet.csv
	The name and location of your metasheet. This should in the same folder as your PROJECT as described above. You need to specify the name of the file and its location

picard_path: /home/lentaing/setup_files/picard-tools-1.113
	MISSING

python2: /usr/bin/python2.7
	Some tools in VIPER use python2, but snakelike is written in python3. So for tools that are written in python2, we need to specify a location for its use

ref_fasta: /data/static_libraries/assembly/humanhg19/rawgenome/hg19.fasta
	Reference genome in a .fasta or .fa format

reference: hg19
	What genome the user is using

rseqc_path: /usr/local/bin
	MISSING

star_index: /data/static_libraries/STAR/humanhg19/
	star-index for the STAR aligner
	If you don't have a star-index, make one by running the following:
	   STAR  --runMode genomeGenerate --runThreadN 24 --genomeDir /where/you/store/reference/genomes -genomeFastaFiles /dir/to/hg19/hg19.fa

star_rRNA_index: /home/lentaing/tmp/weinstock/newrun/ref_files/humanhg38_ncrna/
	MISSING

PARAMS:
———————
This section holds parameters specific to your project design

stranded: 'true'
	MISSING
library_type: 'fr-firststrand'
	MISSING

#— Parameters for plotting
#numgenes: the number of top genes for the gene heatmaps

#for kmeans later…
#globalnumcluster: the number of clusters to output in the clustering analysis for the gene heatmap global
#pairnumcluster: the number of clusters to output in the clustering analysis for the gene heatmap global

SAMPLES:
————————
The location of your files in a zipped fast format (.fastq.gz) or bam format (.bam)

samples:
	## Cell Line 1
	A1:
		- /path/to/file/A1.fastq
	B1:
		- /path/to/file/A1.bam
	C1:
		- /path/to/stranded/file/C1_R1.fastq.gz
		- /path/to/stranded/file/C1_R2.fastq.gz


—————
META:
—————
Make the metadata file in excel, and save it as a .txt or .csv, It doesn’t matter what it is named as long as it is called in the config in the spot marked “metasheet,” see the config section if confused. The format should be something like the following:

sample    cell     condition        treatment      replicates      comp_MCF7_AvB     comp_T47D_CvD
A1	  MCF7     Full_Media       NoDOX          1               1
A2	  MCF7     Full_Media       NoDOX          2               1
B1        MCF7     Full_Media       DOX            1           	   2
B2        MCF7     Full_Media       DOX            2               2
C1        T47D     Full_Media       NoDOX          1                                 1
C2        T47D     Full_Media       NoDOX          2                                 1
D1        T47D     Full_Media       DOX            1                                 2
D2        T47D     Full_Media       DOX            2                                 2

- The first column should always be sample, this should be some kind of sample ID.
- The samples that you want to perform a Differential Expression (DE) on using limma and deseq should be marked by the “comp” columns more on this below
	- This is important! The “control” should be marked with a 1, and the “treatment” should be marked with a 2.
- It is recommended that if you should have a “replicates” column to denote different samples, it is a good idea to not only have each of the sample names be unique, but also make sure that the associated metadata is unique as well to each sample.
- The rest of the  metadata columns are up to the user to write. Sample must always be first, and you are allowed to have as many “comp_XXXX” columns as you want at the end. All of the middle columns are your metadata (for this example, this is cell, condition, treatment, replicates)

Again, make this in excel so that all of the spacing is done correctly and save it out as a .txt or .csv file. This is the most common bug, so please follow this.
Notes on meta
	— use the replicates column to make sure there are no duplicates. SAMPLE NAMES MAY NOT DISTINGUISH DUPLICATES
	- To avoid bugs, the only punctuation that should be used is the underscore “_”. Dashes, periods, etc, could cause a bug because there is a lot of table formatting and manipulation
	- It is very important that you know that samples A is what you mark with 1, and samples B is what you mark with a 2. You should name your output following this format as well "comp_cond_AvB” This will let the reader know what the output DE files refer to. 
	Deseq: ”baseMeanA” refers to samples A, which follows condition 1 and “baseMeanB” refers to samples B which follows condition 2. logfc is B/A
	Limma: Logfc refers to B/A






