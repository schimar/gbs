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
                      
                      - filter_per_row working 
                      - get_base_values working
                      - get_count_values working
                      ---> output of filter_per_row is a huge list, how should we 
                           insert this into the output DataFrame 
                           (see unfinished loop on line 142...) 
                      - ambiguity codes, either own contribution, better yet: check out BioPython !!
                      - comparison with allele_1 and allele_2  -   HOW ??

3) plot distributions, heat maps etc.
