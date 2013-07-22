NAME

gbs (Genotyping by Sequencing)

DESCRIPTION

Python and R tools for handling genotyping by sequencing data. 
The UNEAK pipeline HapMap (hmp.txt) and HapMapCount (hmc.txt) files act as input. 


USAGE NOTES 

* gbs.py: This file contains the function calls for the main pipeline (read-depth filtering, 
  sorting according to abundance of data, HWE-exact-probabilities, expected heterozygosity and some data-conversion). 
  This script calls functions, defined in gbs_mod.py. Follow the script from an IDE or in the console. 
  Import all necessary modules and read in the hmp and hmc files. 
  (Rename your files appropriately or change the script)  

* replica.py: This file contains the function calls to calculate error rates. 
  It calls the functions that are defined in gbs_mod.py.

* plot.r: contains some barplot and violinplot calls for R. 

* plot_locus.py: contains some distribution plotting calls for python/matplotlib.



BUILD NOTES



