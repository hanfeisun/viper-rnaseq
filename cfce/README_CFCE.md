
# VIPER - Visualization Pipeline for RNA-seq
# ---CFCE ADDITION

##Intro

__VIPER__ is a comprehensive RNA-seq analysis tool built using snakemake which allows for ease of use and customization. 

In __VIPER__, there will be three distinct components in your __PROJECT__ folder that are kept purposefully separate to avoid confusion, accidental deletion or editing of essential components. The three components are the __DATA__, __VIPER__, and __ANALYSIS__. __DATA__ is where the user will store all of their data for analysis. Note that having __DATA__ in your __PROJECT__ folder is optional. You can start with a blank __PROJECT__ folder as well (more below.) __VIPER__ is the actual source of the code, you will never have to edit anything in here. __ANALYSIS__ will be the output of __VIPER__.

__NOTE:__ the folder ref_files contain all of the reference files that viper
needs.  You probably won't have to touch or do anything here.

Although included in this README are step-by-step instructions, it is assumed that the user should have a basic understanding of a command line interface on mac or linux.


## Getting Started
1. Open a terminal and navigate to the location of the __PROJECT__ folder. (path/to/PROJECT means the full computer path of where your project is, ex. usr/home/project)

    `cd /path/to/PROJECT/`
    
  - If you do not have a __PROJECT__ folder, we recommend you start by making a new __PROJECT__ folder.

    `cd /path/to`  
    `mkdir PROJECT`  
    `cd PROJECT`
    
2. Run the following command to begin setup of your __VIPER__environment:

	`viperSetup.sh`
	__NOTE:__ if your data is mouse data, type instead: `viperSetup.sh mm9`
	
	- This will output the following:  

	> Next steps: $ means type what follows in command line and ENTER  
	> 1. $ source viper_env.bash  
	> 2. $ source activate viper  
	> 3. [EDIT config.yaml and metasheet.csv--refer to viper/README.txt]  
   NOTE: you *SHOULD NOT* have to modify the top of config.yaml--just the 'samples:' section  
	> 4. $ ln -s /path/to/your/raw/fastq/files  
	> 
	> TO RUN: $ nohup snakemake -s viper/viper.snakefile &
	>
	> Please email Len or Tosh if you have questions.  Enjoy!
	
	You have now successfully copied in all of the starting pieces of __VIPER__. We will now need to activate the environment to enable access to all of necessary packages and reference files. Follow the steps listed in the outputted blurb. They are also explicitly stated in the following instructions.
	
3. Run the two following commands to activate your environment. Note that there will be no output after running the first command.

	`source viper_env.bash`  
	`source activate viper`
	
	- After running the second command, you will see the following output:  
	
	> discarding /home/lentaing/local/miniconda3/bin from PATH  
	> prepending /home/lentaing/local/miniconda3/envs/viper/bin to PATH  
	
	The __VIPER__ environment has now been activated. 
	
4. 	You will now need to edit your *config* and *metasheet* Details on this are in the below sections marked *__CONFIG__* and *__META__*.

5. (OPTIONAL) Part of the organization of __VIPER__ is to have your __DATA__, __ANALYSIS__, and __VIPER__ all in one place. Your __PROJECT__ folder now has a fully functional `viper` folder within it. __ANALYSIS__ will be output as __VIPER__ runs. You can also create a symlink folder of your data in your __PROJECT__ folder with the following command. Although you do not need to do this, it may be useful for organizational purposes, and it may make editing your *config* easier. (Again, path/to/your/raw/fastq/files means the full computer path of where your files are, ex. mnt/cfce-stor1/port/data/01012015/YOURFASTQS)

	`ln -s /path/to/your/raw/fastq/files`

### CONFIG:
There are two main options for editing your config. You can either download the config file, open it in a text editor, and go from there. Or, the recommended away, is to use a command line text editor. These can be a little intimidating if you are not familiar with a command line interface, but they will also reduce the number of errors that can be made by using a text editor.  

One of the easiest text editors out there is [emacs](https://www.gnu.org/software/emacs/). Although there are many commands in emacs, [emacs cheatsheet](https://www.gnu.org/software/emacs/refcards/pdf/refcard.pdf), I will explicitly give a couple simple commands you will need. You can find this further below below the *samples* section.

One final note is that the config is written as a [yaml](https://en.wikipedia.org/wiki/YAML), and therefore you *__Must only use spaces, NO TABS__*

Your config file has three main sections. __PATHS__, __PARAMS__, __SAMPLES__:

##### PARAMS:

You will probably NOT need to specify path to the *metasheet* as this is already set for you. More details on this in the *__metasheet__* section. But if you name your *metasheet* with a specific name, you will need to change it in the config in this section:

> metasheet: metasheet.csv


This section holds parameters specific to your project design. *Leave quotes in there, these are purposefully there and should not be deleted*

>stranded: 'true'  
>  -  Whether or not your data is stranded. Choices are 'true' or 'false'
>    
>library_type: 'fr-firststrand'  
>  -  If stranded is true, you must specify the library type:  
>  -  Possible values are [ff-firststrand, ff-secondstrand, ff-unstranded, fr-firststrand, fr-secondstrand, fr-unstranded (default), transfrags]  
>  - The default ('fr-firststrand') is for Illumina truseq
>
>RPKM\_threshold: "1.0"  
>  -  For inital filtering, a minimal RPKM to be considered real  
>
>min\_num\_samples\_expressing\_at\_threshold: "2"  
>  -  Number of samples the gene to be expressed in to be considered significant
>
>filter\_mirna: "TRUE"  
>  -  Whether or not to remove MiRNA such as MIR and SNO rna
>
>SSnumgenes: "250"  
>  -  The number of genes used to create the Sample-Sample Correlation Heatmap
>  
>SFnumgenes: "1000"  
>  -  The number of genes shown in the Sample-Feature Heatmap
>
>num\_kmeans\_clust: "4"  
>  -  What type of Sample-Feature heatmap you want to create.   
>  -  If you choose "0" it will do hierarchical clustering.  
>  -  If you choose "4" it will do kmeans clustering with 4 clusters.  
>  -  You can also input lists of numbers "4,0,5". This will output in a single pdf: a 4 kmeans heatmap, followed by a hierarchical heatmap, followed by a 5 kmeans heatmap.
> 
> snp_scan_genome: "False"   
>  -  Boolean Flag "{True | False}" on whether to perform a genome-wide snp scan *IN ADDITION* to the snp scan done on chr6. *__This is computationally intensive, and should only be used if really necessary__*.

##### SAMPLES:

The location of your files in a zipped fast format (.fastq.gz) or bam format (.bam). Note that for organization purposes, these can be in your __DATA__ folder. *But they do not have to be*. You can also just specify the full path here. *__The sample names here must match the ones in your metasheet__*.

>samples:  
>\## Cell Line 1  
>&nbsp;&nbsp;&nbsp;&nbsp;A1:  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- /path/to/file/A1.fastq  
>&nbsp;&nbsp;&nbsp;&nbsp;B1:  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- /path/to/file/A1.bam  
>&nbsp;&nbsp;&nbsp;&nbsp;C1:  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- /path/to/stranded/file/C1\_R1.fastq.gz  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- /path/to/stranded/file/C1\_R2.fastq.gz

###### EDITING IN EMACS
Open your config by calling it with emacs:  

`emacs config.yaml`

This will open up your config, scroll down to the *samples* section using the arrow keys. Add your samples and make changes to your params in the same format as seen in the template given to you. You can also look below in the *__samples__* section for additional help.

When you are done editing, hold *control* and press X then C. This in computer syntax is `C^x C^c` This will attempt to close your emacs session. If you edited the config at all, it will ask if you want to save your changes at the bottom of the screen. Press `y` to accept changes.


### METASHEET:
Your *__metasheet__* contains all of the information specific to your project design. This includes all of the data about your samples in addition to the differential analyses that would like to perform using __VIPER__.

Make the *__metasheet__* file in excel, and save it as a .csv file. A template is provided when you activate your __VIPER__ environment. It doesn’t actually matter what your meta is named as long as it is called in the *__config__* in the spot marked “metasheet,” see the *__config__* section if confused. The format should be something like the following:

| Sample | Cell | Condition  | Treatment | Replicates | comp\_MCF7\_AvB | comp\_T47D\_CvD |
|--------|------|------------|-----------|------------|---------------|---------------|
| A1     | MCF7 | Full_Media | NoDOX     | 1          | 1             |               |
| A2     | MCF7 | Full_Media | NoDOX     | 2          | 1             |               |
| B1     | MCF7 | Full_Media | DOX       | 1          | 2             |               |
| B2     | MCF7 | Full_Media | DOX       | 2          | 2             |               |
| C1     | T47D | Full_Media | NoDOX     | 1          |               | 1             |
| C2     | T47D | Full_Media | NoDOX     | 2          |               | 1             |
| D1     | T47D | Full_Media | DOX       | 1          |               | 2             |
| D2     | T47D | Full_Media | DOX       | 2          |               | 2             |

- The first column should always be sample, this should be some kind of sample ID. *__This sample ID must match the sample names in the config file__*
- The samples that you want to perform a Differential Expression (DE) on using limma and deseq should be marked by the “comp” columns more on this below
	- This is important! The “control” should be marked with a 1, and the “treatment” should be marked with a 2.
- It is recommended that you should have a “replicates” column to denote different samples with the same metadata. 
- The rest of the  metadata columns are up to the user to write. Sample must always be first, and you are allowed to have as many “comp_XXXX” columns as you want at the end. All of the middle columns are your metadata (for this example, this is cell, condition, treatment, replicates)

- Again, make this in excel so that all of the spacing is done correctly and save it out as a .txt or .csv file. This is the most common bug, so please follow this.
- Common Problems with *__metasheet__*:
	- To avoid bugs, the only punctuation that should be used is the underscore “_”. Periods or other punctuation, could cause a bug.
	- It is very important that you know that samples A is what you mark with 1, and samples B is what you mark with a 2. You should name your output following this format as well "comp\_cond\_AvB” This will let the reader know what the output DE files refer to. 
	   -  Deseq: ”baseMeanA” refers to samples A, which follows condition 1 and “baseMeanB” refers to samples B which follows condition 2. logfc is B/A
	   -  Limma: Logfc refers to B/A


## Running VIPER

If all of the above have been followed correctly, you should be ready to run __VIPER__! [snakemake documentation](https://bitbucket.org/snakemake/snakemake/wiki/Documentation) is an extensive and powerful tool with many features that the user should learn, but for now, you will only need a few simple commands

In your __PROJECT__ folder run the following command to see if you are ready to go:  

`snakemake --snakefile viper/viper.snakefile -n`  

This will return a large output which basically outlines what __VIPER__ is about to do. If no errors (red text) come back, then you are ready to go! Type in the following command and enter:

`nohup snakemake --snakefile viper/viper.snakefile &`

Notice that we use the [nohup](https://en.wikipedia.org/wiki/Nohup) command here. This will do two things: 1) output all of the command line output into a file call *nohup.out*. 2) More importantly, this will keep __VIPER__ running even if you close your terminal. This is important for full __VIPER__ runs that may take several hours.

#### Running in Parallel:
Running specific code sections or just the downstream analysis (see App B.) is computationally easy, and will not require multiple cores. But If you are running samples from the start. You will need to use multiple cores to utilize the full speed of __VIPER__. Parallelizing __VIPER__ is as easy as using the `-j` flag. It is recommended that if you are running 1-4 samples, use 2 cores. If you are running 5-8, use 4 cores. *Please Use the mbcf to run any more than 8 samples as to avoid clogging the server*

`nohup snakemake --snakefile viper/viper.snakefile -j 4 &`

#### Resuming your VIPER Session
When you close your terminal. Your __VIPER__ session will terminate. Note that if you ran __VIPER__ using nohup as recommended, your run will continue. To resume your __VIPER__ session, change to your __PROJECT__ folder and enter the same commands from step 3 of __Getting Started__.

`source viper_env.bash`  
`source activate viper`

This will resume your __VIPER__ environment and you should be good to go.

## Running a MBCF VIPER run on CFCE:
Situation: You've asked MBCF to analyze your data using the viper pipeline.  They have put the data and the results in the port directory (/mnt/cfce-stor1/data/port/XXX/VIPER_YYY).  You WANT to re-run viper (or some parts of viper) on a CFCE server (cfce1).

STEP 1: make a project directory as described above (see Getting Started)
STEP 2: copy the data from port:
`PROJECT $ rsync -avz /mnt/port/XXX/VIPER_YYY .`
NOTE: '$' symbolized the command-line prompt
NOTE: this might take a while
STEP 3: you might have to cd to PROJECT/VIPER_YYY.  
STEP 4: SETUP PROJECT/VIPER_YYY to run on CFCE-
`$ viperSetup.sh env`
STEP 5: LOAD the viper environment
`source viper_env.bash`
`source activate viper`

To run viper, see 'Running VIPER' above.

### APPENDIX A: Specific Replotting
After you have run __VIPER__ in its entirety, you may want to go back and tweak your outputs. Maybe adding or subtracting metadata columns, differential expression columns, or maybe just doing a subset of your data. Below is a list of snakemake commands to run __VIPER__ to rerun some specifics for further downstream analysis tweaking.

To learn about how snakemake works, and some of the specifics of the following commands and others, look into the [snakemake documentation](https://bitbucket.org/snakemake/snakemake/wiki/Documentation)

The following are some useful commands for rerunning and adding to the download analysis without having to rerun the whole pipeline:

`snakemake -s viper/viper.snakefile -n`  
- Performs a "dry run" and will tell the user what __VIPER__ is about to do. Highly recommended that you run this before every real command you enter. 

`snakemake -s viper/viper.snakefile -j 4`  
- Runs __VIPER__ using 4 cores
 
`snakemake -s viper/viper.snakefile analysis/plots/heatmapSF_plot.pdf -f`  
- Reruns the Sample-Feature heatmap componenet of __VIPER__.
  
`snakemake -s viper/viper.snakefile analysis/plots/heatmapSS_plot.pdf -f`  
- Reruns the Sample-Sample heatmap componenet of __VIPER__.

`snakemake -s viper/viper.snakefile analysis/plots/pca_plot.pdf -f`  
- Reruns the PCA componenet of __VIPER__.

Adding comp columns will automatically make __VIPER__ generate new differential expressions analysis and adjust figures accordingly.

__VIPER__ as a snakemake based program will also detect whenever the metasheet is changed and run a number of rules that use this metadata. Try this by 'touching' the metasheet, then running a dry-run of __VIPER__.

`touch metasheet.csv`

### APPENDIX B: I already ran VIPER on the server, but want to rerun downstream analysis or add a few samples locally
Snakemake as a program will look at what needs to be outputted, and generate that desired output and everything upstream of it. If the upstream work is done, it will not rerun it. To work with data you have already run using __VIPER__ on the mbcf, all you have to do is copy and paste the analysis folder from the mbcf, into your __PROJECT__ folder. Note that when specifying the location of where you are copying the __ANALYSIS__ to, you must write the full path of your project, and then the word 'analysis'. This will create a new folder by that name and copy and paste everything in the 'analysis' folder from the mbcf into the new 'analysis' folder.

` cp -rp /path/to/viper/analysis /path/to/project/folder/analysis`


### APPENDIX C: I have two separate analyses, that I now want to combine and analyze
For this, we will need to use a little bit of Linux magic, specfically using the tool [rsync](http://linux.die.net/man/1/rsync). Rsync is a tool designed exactly for this purpose, and will recursively fo through your source folder, and copy folders into their proper place in the destination folder.
`rsync -avh --progress /path/from/data/you/want/analysis/ /path/to/new/project/analysis/`


### APPENDIX D: Current Errors and their fixes (to be fixed soon!)
1. PCA plot error: 
 
>Error in job pca_plot while creating output files analysis/plots/pca_plot.pdf, analysis/plots/images/.
RuleException in line 551 of /mnt/cfce-stor1/home/mgc31/code/viperproject/viper/viper.snakefile:
Output files analysis/plots/pca_plot.pdf, analysis/plots/images/ are older than input files. Did you extract an archive? Make sure that output files have a more recent modification date than the archive, e.g. by using 'touch'.
Removing output files of failed job pca_plot since they might be corrupted:
analysis/plots/pca_plot.pdf, analysis/plots/images/
Will exit after finishing currently running jobs.
Exiting because a job execution failed. Look above for error message

FIX: Delete the pca image files located in `analysis/plots/images/pca_plot_X.png` and

### CONTACT INFO:
For MBCF Issues: Zach Herbert (zherbert@mail.dfci.harvard.edu)  
For CFCE Issues: Len Taing (len.taing@gmail.com)  
For VIPER Issues: MacIntosh Cornwell (macintoshg_cornwell@dfci.harvard.edu), Len Taing (len.taing@gmail.com), Mahesh Vangala (vangalamaheshh@gmail.com)  
For questions about snakemake: Johannes Köster (koester@jimmy.harvard.edu)


 







