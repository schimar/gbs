
library(adegenet)
library(genetics)
library(pegas)
library(hierfstat)

setwd("/home/mschilling/Desktop/gbs/genepop/")
#				locality (for test.between)	
# populations: 	1: FLFL		1
#		2: HSPQ		1
#		3: KFO		2
#		4: MI		1
#		5: SFQ		1
#		6: WWA		2
#######################
North		Southwest
MI		WWA
SFQ		KFO
FLF		
HSPQ		

NvSW <- c(rep(1,22), rep(1, 29), rep(2, 3), rep(1, 22), rep(1,26), rep(2, 28))
###

pops <- as.factor(c(paste("P", rep(1, 16), sep= ""), paste("P", rep(2, 20), sep= ""), paste("P", rep(3, 6), sep= ""), paste("P", rep(4, 20), sep= ""), paste("P", rep(5, 20), sep= ""), paste("P", rep(6, 19), sep= "")))

hipops <- c(rep(1, 16), rep(2, 20), rep(3, 6), rep(4, 20), rep(5, 20), rep(6, 19))

###
zero <- read.genepop("zero_select.gen")
hizero <- genind2hierfstat(zero)
fszero <- pp.fst(hizero)

# new fst table

     [,1]       [,2]      [,3]       [,4]        [,5]       [,6]
[1,]  NaN 0.01688912 0.1592397 0.01141109 0.008914716 0.15873070
[2,]  NaN        NaN 0.1500468 0.01512255 0.008248107 0.16484889
[3,]  NaN        NaN       NaN 0.14914138 0.156277103 0.05510967
[4,]  NaN        NaN       NaN        NaN 0.007059692 0.16442417
[5,]  NaN        NaN       NaN        NaN         NaN 0.16546178
[6,]  NaN        NaN       NaN        NaN         NaN        NaN


zerostats <- basic.stats(hizero, diploid= T, digits= 4)
zerostats$Ho
write.csv(zerostats$Ho, "Hobs_zero_select.csv")
#######################
base4 <- read.genepop("4base_select.gen")
hi4base <- genind2hierfstat(base4)
fs4base <- pp.fst(hi4base)

# new fst table
     [,1]         [,2]         [,3]          [,4]          [,5]          [,6]
[1,]  NaN -0.001227949 0.0017601018 -0.0002017351 -0.0011838834 -0.0007055007
[2,]  NaN          NaN 0.0001734298 -0.0008255471  0.0000784478 -0.0008330801
[3,]  NaN          NaN          NaN  0.0009709388  0.0040654781 -0.0034436769
[4,]  NaN          NaN          NaN           NaN -0.0008073182  0.0000594714
[5,]  NaN          NaN          NaN           NaN           NaN -0.0001669820
[6,]  NaN          NaN          NaN           NaN           NaN           NaN

base4stats <- basic.stats(hi4base, diploid= T, digits= 4)
base4stats$Ho
write.csv(base4stats$Ho, "Hobs_4base_select.csv")


##
loci <- hi4base[c(2:27911)]
var4base <- hierfstat::varcomp.glob(levels= hi4base$pop, loci= loci, diploid=T) # took about an hour
# test for significance of effect of N-SW divide:
test4NvSW <- test.between(loci, rand.unit= hi4base$pop, test= NvSW, nperm=1000)
# test for significance of population effect:
test4pop <- test.between(loci, rand.unit = NvSW, test= hi4base$pop, nperm= 1000) 
#######################
adv <- read.genepop("adv_select.gen")
hiadv <- genind2hierfstat(adv)

fsadv <- pp.fst(hiadv) 
# new fst table
     [,1]        [,2]       [,3]        [,4]        [,5]       [,6]
[1,]  NaN 0.003061196 0.04278085 0.002818004 0.003717721 0.03935200
[2,]  NaN         NaN 0.04540788 0.003504800 0.001355175 0.04370501
[3,]  NaN         NaN        NaN 0.047050434 0.047196920 0.01464858
[4,]  NaN         NaN        NaN         NaN 0.001275478 0.04868107
[5,]  NaN         NaN        NaN         NaN         NaN 0.04552022
[6,]  NaN         NaN        NaN         NaN         NaN        NaN

advstats <- basic.stats(hiadv, diploid= T, digits= 4)
write.csv(advstats$Ho, "Hobs_adv_select.csv")

#######################
maf <- read.genepop("MAF_select.gen")
hiMAF <- genind2hierfstat(maf)

fsmaf <- pp.fst(hiMAF)
# new table 

     [,1]       [,2]      [,3]        [,4]        [,5]       [,6]
[1,]  NaN 0.01640308 0.1606804 0.009843969 0.009130159 0.15967480
[2,]  NaN        NaN 0.1495929 0.015318897 0.007312567 0.16660696
[3,]  NaN        NaN       NaN 0.147682033 0.156963412 0.05478745
[4,]  NaN        NaN       NaN         NaN 0.005506727 0.16557485
[5,]  NaN        NaN       NaN         NaN         NaN 0.16595860
[6,]  NaN        NaN       NaN         NaN         NaN        NaN


mafstats <- basic.stats(hiMAF, diploid= T, digits= 4)
write.csv(mafstats$Ho, "Hobs_maf_select.csv")


####
# Colin's dataset

col <- read.genepop('SubsetAPG_data_5-31-12_GP.gen')
hicol <- genind2hierfstat(col)

fscol <- pp.fst(hicol)
     [,1]          [,2]      [,3]        [,4]          [,5]       [,6]
[1,]  NaN -0.0009699136 0.1185477 0.001388089  0.0096114882 0.07199823
[2,]  NaN           NaN 0.1152206 0.000536980 -0.0002834597 0.07357368
[3,]  NaN           NaN       NaN 0.123674020  0.1202985364 0.04622776
[4,]  NaN           NaN       NaN         NaN  0.0108092514 0.08118123
[5,]  NaN           NaN       NaN         NaN           NaN 0.07763793
[6,]  NaN           NaN       NaN         NaN           NaN        NaN
# FLFL, HSPQ, KFO, MI, SFQ, WWA

colstats <- basic.stats(hicol, diploid= T, digits= 4)

colstats$overall
    Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
0.7149 0.7505 0.7893 0.0388 0.7970 0.0465 0.0491 0.0583 0.0474 0.1864 


