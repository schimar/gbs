NAME

gbs (Genotyping by Sequencing)

DESCRIPTION

Python tools for handling genotyping by sequencing data. 


BUILD NOTES

USAGE NOTES


TODO

1) clean up code, refactor.
                      - hmc - resolve tab issue of first 10 entries 
                      - hmp length still hard-coded ...
                      - hmp, allele-split-loop, incorporate into big function (?)
2) establish filter 2:
                      
                      - filter_single_col working 

                      - ambiguity codes:
                          - logic! 
                          - comparison with allele_1 and allele_2 
                          - if both count_1 and count_2 are >= 4 (threshold, respectively..) then heterozygotes.
                          

3) plot distributions, heat maps etc.
                      - plot the (filtered) individual/#of plot distribution (hist) and probably overlay on the old plot (w/out filter)