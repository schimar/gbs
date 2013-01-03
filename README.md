NAME

gbs (Genotyping by Sequencing)

DESCRIPTION

Python tools for handling genotyping by sequencing data. 
The UNEAK pipeline HapMap (hmp.txt) and HapMapCount (hmc.txt) files act as input. 


USAGE NOTES 

* gbs_cli.py: Open a terminal, go to respective folder (where python script and input files are)
              and type 'python gbs_cli.py <hmp_filename> <hmc_filename>
              You should then find an output file (data_sorted4.csv) 
        

* gbs.py: Follow the script from an IDE or in the console. Starting from line 144/145, read the hmp and hmc files. 
          Rename your files appropriately or change the script.  

Note: 
      Both will filter the hmp allele positions by a default threshold of 4 and sort the output 
      (rows & columns) by frequency of alleles present. Until now, you would have to transfer 
      this (csv file) to a format of your choice in order to run simple Population Genetic analyses. 
      We are currently working at the transfer into Genepop format. 

BUILD NOTES



TODO



2) GenePop
                      - establish data transfer to GenePop file format


3) plot distributions, heat maps etc.
                      - plot count_value 1 and 2 in a scatterplot (of non_N's)
                      - plot heat map
