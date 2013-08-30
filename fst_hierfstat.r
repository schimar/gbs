
library(adegenet)
library(genetics)
library(pegas)
library(hierfstat)
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


zero <- read.genepop("zero.gen")
hizero <- genind2hierfstat(zero)
fszero <- pp.fst(hizero)

     [,1]       [,2]       [,3]       [,4]       [,5]        [,6]
[1,]  NaN 0.03655102 0.08712774 0.04295521 0.02946134 0.069948397
[2,]  NaN        NaN 0.08684725 0.03844623 0.02416698 0.069258059
[3,]  NaN        NaN        NaN 0.04129497 0.09363674 0.002770562
[4,]  NaN        NaN        NaN        NaN 0.02388941 0.040033628
[5,]  NaN        NaN        NaN        NaN        NaN 0.057629733
[6,]  NaN        NaN        NaN        NaN        NaN         NaN

zerostats <- basic.stats(hizero, diploid= T, digits= 4)
zerostats$Ho
write.csv(zerostats$Ho, "Hobs_zero.csv")
#######################
base4 <- read.genepop("4base.gen")
hi4base <- genind2hierfstat(base4)
pp.fst(hi4base)
     [,1]      [,2]       [,3]       [,4]       [,5]        [,6]
[1,]  NaN 0.0296248 0.09254855 0.03964764 0.03230319 0.064452822
[2,]  NaN       NaN 0.08318049 0.04018450 0.02344656 0.066976868
[3,]  NaN       NaN        NaN 0.03635674 0.08666103 0.005072949
[4,]  NaN       NaN        NaN        NaN 0.02696786 0.042396589
[5,]  NaN       NaN        NaN        NaN        NaN 0.061897294
[6,]  NaN       NaN        NaN        NaN        NaN         NaN
There were 14 warnings (use warnings())


base4stats <- basic.stats(hi4base, diploid= T, digits= 4)
base4stats$Ho
write.csv(base4stats$Ho, "Hobs_4base.csv")


##
loci <- hi4base[c(2:27911)]
var4base <- hierfstat::varcomp.glob(levels= hi4base$pop, loci= loci, diploid=T) # took about an hour
# test for significance of effect of N-SW divide:
test4NvSW <- test.between(loci, rand.unit= hi4base$pop, test= NvSW, nperm=1000)
# test for significance of population effect:
test4pop <- test.between(loci, rand.unit = NvSW, test= hi4base$pop, nperm= 1000) 
#######################
adv <- read.genepop("adv.gen")
hiadv <- genind2hierfstat(adv)

fsadv <- pp.fst(hiadv) 
     [,1]        [,2]       [,3]        [,4]        [,5]         [,6]
[1,]  NaN 0.008137494 0.01932828 0.011592544 0.009720402  0.014444368
[2,]  NaN         NaN 0.02193204 0.011517502 0.006875923  0.016881457
[3,]  NaN         NaN        NaN 0.006287892 0.021967516 -0.002697101
[4,]  NaN         NaN        NaN         NaN 0.006809729  0.010199065
[5,]  NaN         NaN        NaN         NaN         NaN  0.015871558
[6,]  NaN         NaN        NaN         NaN         NaN          NaN



advstats <- basic.stats(hiadv, diploid= T, digits= 4)
write.csv(advstats$Ho, "Hobs_adv.csv")

#######################
maf <- read.genepop("MAF.gen")
hiMAF <- genind2hierfstat(maf)

fsmaf <- pp.fst(hiMAF)
     [,1]       [,2]       [,3]       [,4]       [,5]          [,6]
[1,]  NaN 0.03663721 0.08590356 0.04520125 0.03152035  0.0713551777
[2,]  NaN        NaN 0.08230460 0.04046002 0.02533975  0.0693412906
[3,]  NaN        NaN        NaN 0.03234230 0.08812859 -0.0006338958
[4,]  NaN        NaN        NaN        NaN 0.02441456  0.0410168487
[5,]  NaN        NaN        NaN        NaN        NaN  0.0587082733
[6,]  NaN        NaN        NaN        NaN        NaN           NaN



mafstats <- basic.stats(hiMAF, diploid= T, digits= 4)
write.csv(mafstats$Ho, "Hobs_maf.csv")



