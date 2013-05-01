library(adegenet)
library(ape)
library(pegas)
###
vignette('ReadingFiles') 

hmp <- read.loci("hmp_sorted.csv", loci.sep= ',', col.loci = 3:26 , col.pop = 2, row.names= 1)
###

summary(hmp)
print(hmp, details= T)

heterozygosity(hmp$TP156963, variance= F)
apply(hmp, 2, heterozygosity)

by(hmp, hmp$population, hw.test)
### adegenet
hmpS4 <- as.genind(hmp)



##############################
#from Bio.PopGen import GenePop
#handle = open("sub3000_sorted4_genepop.gen")
#rec = GenePop.read(handle) 
#handle.close()

####################
# see adegenet-basics.pdf

library(adegenet)
library(genetics)
library(pegas)

setwd("/home/schimar/Desktop/lab/gbs/genepop/")
#######
data <- read.genepop("sub3000_sorted4_genepop.gen")

summ <- summary(data) # too big for interpretation (at least for me...)
# N:
140
# Pop.eff
FLFL   HSPQ     KF     MI     SFQ     WWA
22      29      6      29     26       28

..........

#########################################################
# plots
par(mfrow=c(2,2))
plot(summ$pop.eff,summ$pop.nall,xlab="Colonies sample size",ylab="Number of alleles",main="Alleles numbers and sample sizes")
text(summ$pop.eff,summ$pop.nall,lab=names(summ$pop.eff))
barplot(summ$loc.nall,ylab="Number of alleles", main="Number of alleles per locus")
barplot(summ$Hexp-summ$Hobs,main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
barplot(summ$pop.eff,main="Sample sizes per population",ylab="Number of genotypes",las=3)

#########################################################
# matrix of p-values
hw1 <- HWE.test.genind(data,res="matrix")
# warnings
# Chi-squared approximation may be incorrect
dim(hw1) # 18000 tests (per locus and population)

# which tests are highly significant?
colnames(hw1)
# length of 3000, so all of them (?)

idx <- which(hw1<0.0001, TRUE)

# Here, 182 tests indicate departure from HW. Rows give populations, columns give markers. 

# Now complete tests are returned, but the significant ones are already known.
hw2 <- HWE.test.genind(data,res="full")
# again, warnings:
# Chi-squared approximation may be incorrect

#################################
# using DatABEL and GenABEL

library(GenABEL, DatABEL, genetics)

text2databel("matr.txt", outfile = "matr1", R_matrix = TRUE, type = "UNSIGNED_INT")

#Options in effect:
#   --infile    = matr.txt
#   --outfile   = matr1
#   --skiprows  = 1
#   --skipcols  = 1
#   --cnrow     = ON, using line 1 of 'matr.txt'
#   --rncol     = ON, using column 1 of 'matr.txt'
#   --transpose = OFF
#   --Rmatrix   = ON
#   --nanString = N

#zero <- text2databel("zero_all.csv", outfile= "zero1", R_matrix= T, type= "DOUBLE", transpose= T, naString= "N")

#zero <- read.csv("zero_all.csv", header = T)
#adv <- read.csv("data_adv_sorted4.csv", header = T)
as.genotype(adv$FLFL04) # works
#genadv <- t(adv)
#zero <- read.csv("zero_all_NA.csv", header=T)

tzero <- read.csv("t_zero_all_NA.csv", header= T)
t4base <- read.csv("t_4base_all_NA.csv", header = T)
tadv <- read.csv("t_adv_all_NA.csv", header = T)
tMAF <- read.csv("t_MAF_all_NA.csv", header = T)


tzero$population <- substring(tzero$X, 0, 2)
flfl <- subset(tzero, tzero$population == 'FL')
HWE.test(as.genotype(flfl$TP156963))
HWE.exact(as.genotype(flfl$TP156963))

z <- as.loci(tzero)
Fst(z)

subtzero <- tzero[c(1, 160185, 2:1000)]

tzero[c(1, 160185, 2:5000)]


