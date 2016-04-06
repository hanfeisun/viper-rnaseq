
# VIPER - Visualization Pipeline for RNA-seq

##Intro

__VIPER__ is a comprehensive RNA-seq analysis tool built using snakemake which allows for ease of use and customization. We combine the use of etc……..

In __VIPER__, there will be three distinct components in your __PROJECT__ folder that are kept purposefully separate to avoid confusion, accidental deletion or editing of essential components. The three components are the __DATA__, __VIPER__, and __ANALYSIS__. __DATA__ is where the user will store all of their data for analysis. Note that having __DATA__ in your __PROJECT__ folder is optional. You can start with a blank __PROJECT__ folder as well (more below.) __VIPER__ is the actual source of the code, you will never have to edit anything in here. __ANALYSIS__ will be the output of __VIPER__.

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


## Getting Started

Before __VIPER__ can be used, all of the necessary packages and tools must be downloaded and installed.


### Installing wget

To get these packages, we will use a command line tool called [wget](http://www.gnu.org/software/wget/). *wget* is a popular tool for pulling things off of the internet, and although it should be already be installed on most linux machines, you may need to install is if you are using a Mac [Link to Download wget](http://rudix.org/packages/wget.html)

If you are unsure whether or not you have *wget* enter `wget` and if you get a response that isn't `wget: command not found`, then you should have wget.


### Installing Miniconda3

We will be using Miniconda3 to manage most of the packages that go into __VIPER__. [Miniconda3](http://conda.pydata.org/miniconda.html) is a simple tool that include the package manager [Conda](http://conda.pydata.org/docs/). Conda is used to manage all of the tools used in __VIPER__ including version control

Navigate to a folder where you keep your coding tools or wherever you want to store the main Miniconda3 folder. Then enter the following commands:

- If you are using Linux or sshed into a Linux system:  

`
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
`
`
bash Miniconda3-latest-Linux-x86_64.sh
`

- If you are using Mac locally:  

`
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
`
`
bash Miniconda3-latest-MacOSX-x86_64.sh
`

- Follow the commands listed on screen, press the _enter_ key to scroll down.
- __Make sure to answer yes when asked if you want to prepend Miniconda3 to PATH.__
- Close your terminal, open a new one and you should now have Conda working! Test by entering:  
`
conda update conda
`
	- Press `y` to confirm the conda updates

## Installing VIPER and setting up your environment

1. Open a terminal and navigate to the location of the __PROJECT’s__ __DATA__ (path/to/PROJECT means the full computer path of where your project is, ex. usr/home/project)

    `
    cd /path/to/PROJECT/
    `
  - If you do not have a folder with __DATA__ for your __PROJECT__ we recommend you make a new __PROJECT__ folder. You can put __DATA__ in it now, or leave it empty.

    `	
    cd /path/to
	`  
    `
	mkdir PROJECT
	`  
    `
	cd PROJECT
    `

2. If you are familiar with git, bitbucket, have ssh keys set up etc. and would like the latest VIPER version directly from [bitbucket](https://bitbucket.org/cfce/viper) Then refer to (2a). If you are unfamiliar with git and would like to download VIPER using a compressed folder, refer to (2b)
 -  (2a) The __PROJECT__ folder should be blank or only have your data in it. Next we will acquire __VIPER__ using git.
 -  within your PROJECT folder initialize a git repository
  
    `
    git init
    `  
    
  -  after git initialization and set up, you will need to clone the __VIPER__ folder in to your __PROJECT__ folder.  

    `
    git clone git@bitbucket.org:cfce/viper.git
    `

  -  (2B) The PROJECT folder for now should only have your data in it. Next we will acquire VIPER using wget.  

    `
	wget https://bitbucket.org/cfce/viper/get/master.tar.gz
    `  
    `
	tar -xf master.tar.gz
    `  
    `
	rm master.tar.gz
    `  
    `
    mv cfce-viper-ce8b49cbce1c viper
    `  
    __NOTE:__ The last serires of digits (cfce-viper-XXXX) will be different for everyone! There should only be one thing in your folder at this though, so use tab completion to complete.

3. Download all of the necessary packages. Most of the tools needed can be downloaded via Miniconda. Run the bash script from within the __PROJECT__ folder to download all of the packages.
  
	`
	bash viper/build_env.bash
   `
    - This will take a while, While it is loading, begin reading the *_config_* and *_meta_* sections.

4. For annotation purposes, __VIPER__ uses a static library to call from. Unzip this library using the following command:  

   `
   bunzip2 viper/static/humanhg19.annot.csv.bz2
   `

5. Copy the *__config__* file from within your __VIPER__ folder into your __PROJECT__ folder. From within your __PROJECT__ folder, enter:

   `
   cp viper/config.yaml ./
   `  
  -  Edit your *__config__* file, more on this in the *__config__* section.
  -  Note for members of Dana-Farber with access to the CFCE Server: Refer to Appendix A.

6. Create your *__metasheet__* file and move it into your __PROJECT folder__.
  - More on this in the *__meta__* section  

7. At this point, you should have your */path/to/PROJECT* folder with the following contained within it:
> VIPER  
> DATA  - *optional*   
> config.yaml  
> metasheet.csv


## Running VIPER

If all of the above have been followed correctly, you should be ready to run __VIPER__! snakemake is an extensive and powerful tool with many features that the user should learn, but for now, you will only need a few simple commands

In your __PROJECT__ folder run the following command to see if you are ready to go:  

`
	snakemake --snakefile viper/viper.snakefile -n
`  

This will return a large output which basically outlines what VIPER is about to do. If no errors come back, then you are ready to go! Type in the following command and enter:

`
	snakemake --snakefile snakemake/RNAseq.snakefile
`


### CONFIG:
Your config file has three main sections. __PATHS__, __PARAMS__, __SAMPLES__:

##### PATHS:
In this section, the user will need to specify the location of the following packages and tools. Below is an example path and a description of what each item is.


>bed_file: /data/static_libraries/RefGene/refseqGenesHg19.bed  
>  -  MISSING
>
>genome_lib_dir: /home/lentaing/tmp/weinstock/newrun/ref_files/Hg19_CTAT_resource_lib  
>  -  MISSING  
>
>gtf_file: /home/lentaing/tmp/weinstock/newrun/ref_files/Hg19_CTAT_resource_lib/ref_annot.gtf  
>  -  MISSING
>
>metasheet: metasheet.csv  
>  -  The name and location of your metasheet. This should in the same folder as your PROJECT as described above. You need to specify the name of the file and its location
>
>picard_path: /home/lentaing/setup_files/picard-tools-1.113  
>  -  MISSING
>
>python2: /usr/bin/python2.7  
>  -  Some tools in VIPER use python2, but snakelike is written in python3. So for tools that are written in python2, we need to specify a location for its use
>
>ref_fasta: /data/static_libraries/assembly/humanhg19/rawgenome/hg19.fasta  
>  -  Reference genome in a .fasta or .fa format
>
>reference: hg19  
>  -  What genome the user is using
>
>rseqc_path: /usr/local/bin  
>  -  MISSING
>
>star_index: /data/static_libraries/STAR/humanhg19/  
>  -  star-index for the STAR aligner  
>  -  If you don't have a star-index, make one by running the following:
>
>`  
>	   STAR  --runMode genomeGenerate --runThreadN 24 --genomeDir /where/you/store/reference/genomes -genomeFastaFiles /dir/to/hg19/hg19.fa
>`
>
>star_rRNA_index: /home/lentaing/tmp/weinstock/newrun/ref_files/humanhg38_ncrna/  
>  -  MISSING
> 
> varscan_jar_path: /home/lentaing/setup_files/varscan-master/VarScan.v2.4.1.jar  
>  -  MISSING  
> 
> gene_annotation: viper/static/humanhg19.annot.csv  
>  -  Path to annotation files (e.g. ENSEMBL\_ID, Gene Description, Go Terms, etc.).  This file can be generated by using [SCRIPT that Len needs to include].  Pre-made annotations for hg19 and mm9 can be found in viper/static (simply bunzip2 them).


##### PARAMS:

This section holds parameters specific to your project design

>stranded: 'true'  
>  -  MISSING
>    
>library_type: 'fr-firststrand'  
>  -  MISSING
>
>RPKM\_threshold: "2.0"  
>  -  MISSING  
>
>min\_num\_samples\_expressing\_at\_threshold: "4"  
>  -  MISSING
>
>filter\_mirna: "TRUE"  
>  -  MISSING
>
>SSnumgenes: "250"  
>  -  MISSING
>  
>SFnumgenes: "1000"  
>  -  MISSING
>
>num\_kmeans\_clust: "4"  
>  -  MISSING
> 
> snp_scan_genome: "true"   
>  -  Boolean Flag "{True | False}" on whether to perform a genome-wide snp scan *IN ADDITION* to the snp scan done on chr6.


##### SAMPLES:

The location of your files in a zipped fast format (.fastq.gz) or bam format (.bam). Note that for organization purposes, these will be in your __DATA__ folder. *But they do not have to be*. You can also just specify the full path here.

>samples:  
>\## Cell Line 1  
>&nbsp;&nbsp;&nbsp;&nbsp;A1:  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- /path/to/file/A1.fastq  
>&nbsp;&nbsp;&nbsp;&nbsp;B1:  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- /path/to/file/A1.bam  
>&nbsp;&nbsp;&nbsp;&nbsp;C1:  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- /path/to/stranded/file/C1\_R1.fastq.gz  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- /path/to/stranded/file/C1\_R2.fastq.gz



### META:

Make the *__metasheet__* file in excel, and save it as a .txt or .csv, It doesn’t matter what it is named as long as it is called in the *__config__* in the spot marked “metasheet,” see the *__config__* section if confused. The format should be something like the following:

| Sample | Cell | Condition  | Treatment | Replicates | comp_MCF7_AvB | comp_T47D_CvD |
|--------|------|------------|-----------|------------|---------------|---------------|
| A1     | MCF7 | Full_Media | NoDOX     | 1          | 1             |               |
| A2     | MCF7 | Full_Media | NoDOX     | 2          | 1             |               |
| B1     | MCF7 | Full_Media | DOX       | 1          | 2             |               |
| B2     | MCF7 | Full_Media | DOX       | 2          | 2             |               |
| C1     | T47D | Full_Media | NoDOX     | 1          |               | 1             |
| C2     | T47D | Full_Media | NoDOX     | 2          |               | 1             |
| D1     | T47D | Full_Media | DOX       | 1          |               | 2             |
| D2     | T47D | Full_Media | DOX       | 2          |               | 2             |

- The first column should always be sample, this should be some kind of sample ID.
- The samples that you want to perform a Differential Expression (DE) on using limma and deseq should be marked by the “comp” columns more on this below
	- This is important! The “control” should be marked with a 1, and the “treatment” should be marked with a 2.
- It is recommended that if you should have a “replicates” column to denote different samples, it is a good idea to not only have each of the sample names be unique, but also make sure that the associated metadata is unique as well to each sample.
- The rest of the  metadata columns are up to the user to write. Sample must always be first, and you are allowed to have as many “comp_XXXX” columns as you want at the end. All of the middle columns are your metadata (for this example, this is cell, condition, treatment, replicates)

- Again, make this in excel so that all of the spacing is done correctly and save it out as a .txt or .csv file. This is the most common bug, so please follow this.
- Common Problems with *__metasheet__*
- To avoid bugs, the only punctuation that should be used is the underscore “_”. Dashes, periods, etc, could cause a bug because there is a lot of table formatting and manipulation
	- It is very important that you know that samples A is what you mark with 1, and samples B is what you mark with a 2. You should name your output following this format as well "comp\_cond\_AvB” This will let the reader know what the output DE files refer to. 
	   -  Deseq: ”baseMeanA” refers to samples A, which follows condition 1 and “baseMeanB” refers to samples B which follows condition 2. logfc is B/A
	   -  Limma: Logfc refers to B/A


### APPENDIX A: Dana-Farber CFCE Members
If you are a member of Dana-Farber and have access to the CFCE server, you will already have many of the packages you need installed globally. There will be a config with the paths filled out located within the cfce folder of __VIPER__. Run the following command to obtain your config file with paths already filled out for the cfce.  

`
cp viper/cfce/cfce_config.yaml ./config.yaml
`  



### APPENDIX B: Specific Replotting
After you have run __VIPER__ in its entirety, you may want to go back and tweak your outputs. Maybe adding or subtracting metadata columns, differential expression columns, or maybe just doing a subset of your data. Below is a list of snakemake commands to run __VIPER__ to rerun some specifics for further downstream analysis tweaking.

To learn about how snakemake works, and some of the specifics of the following commands and others, look into the [snakemake documentation](https://bitbucket.org/snakemake/snakemake/wiki/Documentation)

The following are some useful commands for rerunning and adding to the download analysis without having to rerun the whole pipeline:

`
snakemake -s viper/viper.snakefile -n
`  

`
snakemake -s viper/viper.snakefile -j 24
`

`
snakemake -s viper/viper.snakefile analysis/plots/heatmapSF_plot.pdf -f 
`
  
`
snakemake -s viper/viper.snakefile analysis/plots/heatmapSS_plot.pdf -f 
`

`
snakemake -s viper/viper.snakefile analysis/plots/pca_plot.pdf -f 
`

Adding comp columns will automatically make it generate new differential expressions analysis and adjust figures accordingly.

`
touch metasheet.csv
`



 







