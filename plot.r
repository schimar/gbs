library(reshape)
library(gplots)

##################
# barplot_filter_types

filt_matr <- matrix(rbind(homozygotes= c(5758665, 6191765, 2203459, 915009), heterozygotes= c(562149, 129049, 163955, 163955), ambiguous= c(0, 0, 0, 1288450), deparse.level=1), nrow= 3, ncol=4, dimnames= list(c("homozygous", "heterozygous", "ambiguous"), c("initial", "MAF", "4base", "adv")))

###
barplot2(filt_matr, beside=F, col = c("white", "gray50", "gray72"), main = "Filter results", xlab = "Filter type", ylab= "Frequency", ylim= c(0, 7000000), col.sub= "gray50", plot.grid= TRUE, grid.inc= 6, legend = rownames(filt_matr), axis.lty= 1)
text(c(0.7, 1.92, 3.11, 4.32, 0.7, 1.92, 3.11, 4.32, 4.32), 
     c(500000, 500000, 500000, 500000, 6000000, 6000000, 2500000, 1000000, 1700000), labels= c("91%", "98%", "93%", "39%", 
	     "9%", "2%", "6.9%", "6.9%", "54%"))
box()
##
# used this one! 
filt_matr2 <- matrix(rbind(homozygotes= c(5758665, 2203459, 915009, 6191765), heterozygotes= c(562149, 163955, 163955, 129049), ambiguous= c(0, 0, 1288450, 0), deparse.level=1), nrow= 3, ncol=4, dimnames= list(c("homozygous", "heterozygous", "ambiguous"), c("initial", "4base", "adv", "MAF")))

##
barplot2(filt_matr2, beside=F, col = c("white", "gray50", "gray72"), xlab = "Filter type", ylab= "Frequency", ylim= c(0, 7000000), col.sub= "gray50", plot.grid= F, grid.inc= 6, axis.lty= 1)

legend(2.3, 6800000, c("homozygous", "heterozygous", "ambiguous"), fill = c("white", "gray50", "gray72"))

text(c(0.7, 1.92, 3.11, 4.32, 0.7, 1.92, 3.11, 3.11, 4.32), c(500000, 500000, 500000, 500000, 6000000, 2500000, 1000000, 1700000, 6000000), labels= c("91%", "93%", "39%", "98%", "9%", "6.9%", "6.9%", "54%", "2%"))
box()

###
dev.copy2pdf(file="barplot_filter_types_ordered.pdf")

##################
# histograms of prob_homo

hist(homo_prob_MAF_total$V2, ylim=c(0,200000))


##################
# filter expl.  plots:

x <- c(4, 6, 20)
y <- c(0, 3, 6)
par(mfrow=c(1, 3))
plot(y~x, xlim= c(0, 80), ylim= c(0, 80), col= c("green", "blue", "red"), pch= c(1, 2, 0))
#plot(y~x, col= "white", main= "Simple filter with threshold of 4")
abline(h=3.8)#, col= 'red')
abline(v=3.8)#, col= 'red')
text(40, 0, labels= "homozygotes")
text(0, 40, labels= "homozygotes", srt=90)
text(40, 40, labels = "heterozygotes", srt=45, cex= 1.5)
lines(x= c(0, 0, 3.78, 3.78), y= c(0, 3.78, 3.78, 0), col= "red")
lines(x= c(0, 3.78, 3.78, 3.78), y= c(0, 0, 3.78, 0), col= "red")

#lines(x= c(4, 0), y= c(4, 4), col= "red")
#lines(x= c(4, 4), y= c(4, 0), col= "red")
#
#plot(y~x, col= "white", main= "Advanced filter (threshold= 4 and 2*threshold for ambiguous counts")
plot(y~x, xlim= c(0, 80), ylim= c(0, 80), col= c("green", "blue", "red"), pch= c(1, 2, 0))
#plot(y~x, col= "white", main= "Simple filter with threshold of 4")
abline(h=3.8)#, col= 'red')
abline(v=3.8)#, col= 'red')
text(40, 0, labels= "homozygotes")
text(0, 40, labels= "homozygotes", srt=90)
text(40, 40, labels = "heterozygotes", srt=45, cex= 1.5)
lines(x= c(0, 0, 3.78, 3.78), y= c(0, 3.78, 3.78, 0), col= "red")
lines(x= c(0, 3.78, 3.78, 3.78), y= c(0, 0, 3.78, 0), col= "red")
lines(x= c(0, 0, 0, 3.78, 3.78, 3.78), y= c(3.78, 7.78, 7.78, 7.78, 7.78, 3.78), col= "purple")
lines(x= c(3.78, 7.78, 7.78, 7.78, 7.78, 3.78), y= c(0, 0, 0, 3.78, 3.78, 3.78), col= "purple")
#lines(x= c(0, 8, 8), y= c(0, 0, 4), col= "lightblue")
#lines(x= c(0, 0, 4), y= c(0, 8, 8), col= "lightblue")

#
plot(y~x, xlim= c(0, 80), ylim= c(0, 80), col= c("green", "blue", "red"), pch= c(1, 2 ,0))
#plot(y~x, col= "white", main = "MAF filter (0.45)")
abline(a = 0, b= 0.818)#, col= "red")
abline(a=0, b=1.25)#, col= "red")
text(50, 50, labels = "heterozygotes", srt=45, cex= 1.5)
text(50, 10, labels= "homozygotes")
text(10, 50, labels= "homozygotes", srt=90)
#



     x <- stats::runif(12); y <- stats::rnorm(12)
     i <- order(x,y); x <- x[i]; y <- y[i]
     plot(x, y, main="arrows(.) and segments(.)")
     ## draw arrows from point to point :
     s <- seq(length(x)-1)# one shorter than data
     arrows(x[s], y[s], x[s+1], y[s+1], col= 1:3)
     s <- s[-length(s)]
     segments(x[s], y[s], x[s+2], y[s+2], col= 'pink')



###############################################################

# hwe_exact test - plots

###############################################################
library(ggplot2)
library(vioplot)
     
setwd('/home/mschilling/Desktop/gbs/hwe')

hwe_zero <- read.csv("hwe_zero_select.csv", header = T)
hwe_4base <- read.csv("hwe_4base_select.csv", header = T)
hwe_adv <- read.csv("hwe_adv_select.csv", header = T)
hwe_MAF <- read.csv("hwe_MAF_select.csv", header = T)



setwd( "/home/schimar/Desktop/lab/gbs/hwe")

pool_zero <- read.csv("pooled_hwe_zero_select.csv", header = T)
pool_4base <- read.csv("pooled_hwe_4base_select.csv", header = T)
pool_adv <- read.csv("pooled_hwe_adv_select.csv", header = T)
pool_MAF <- read.csv("pooled_hwe_MAF_select.csv", header = T)


# "violin_hwe_per_pop.pdf"
dat <- pool_adv

x1 <- dat$hwe[dat$population == 'FLFL']
x2 <- dat$hwe[dat$population == 'HSPQ']
x3 <- dat$hwe[dat$population == 'KFO']
x4 <- dat$hwe[dat$population == 'MI']
x5 <- dat$hwe[dat$population == 'SFQ']
x6 <- dat$hwe[dat$population == 'WWA']

# x3 has only 1s... 
par(mfrow= c(2,2)) # then 4 times the one below, just change dat and x1 - x6
vioplot(x1, x2, x4, x5, x6, names = c('FLFL', 'HSPQ', 'MI', 'SFQ', 'WWA'), col= 'gray80')
abline(h= 0.05)

# left-upper: hmp
# right-upper: 4base
# left-lower: MAF
# right-lower: adv)

###########################
# get the "P_hwe to total" ratio (with xn from below..)
length(which(x1 > 0.05)) / length(x1)
length(which(x2 > 0.05)) / length(x2)
length(which(x3 > 0.05)) / length(x3)
length(which(x4 > 0.05)) / length(x4)


# "hwe_pooled_all_filters.pdf" (without KFO!)
x1 <- pool_zero$hwe[pool_zero$population != 'KFO']
x2 <- pool_4base$hwe[pool_4base$population != 'KFO']
x3 <- pool_adv$hwe[pool_adv$population != 'KFO']
x4 <- pool_MAF$hwe[pool_MAF$population != 'KFO']

vioplot(x1, x2, x3, x4, h=0.05, names= c('unfiltered', 'TF', 'AF', 'MAFF'), col= 'grey80')
abline(h=0.05)
text(0.8, 0.82, label= '80%')
text(1.8, 0.82, label= '79%')
text(2.8, 0.82, label= '95%')
text(3.8, 0.82, label= '57%')

###
which(complete.cases(hwe_4base)== 'TRUE')


####
boxplot(hwe_4base[, c(2:7)])

geom_violin(
