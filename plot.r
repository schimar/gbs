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

