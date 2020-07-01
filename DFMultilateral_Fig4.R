#########################################################################
###### R Code for FIGURE 4 ##############################################
###### Reference: Diewert, W.E. and K.J. Fox (2020), 
###### "Subsititution Bias in Multilatreal Methods for CPI Construction Using Scanner Data"
#########################################################################
###### Contact: Kevin Fox (K.Fox@unsw.edu.au)
###### 30 June 2020

###### Exercises on scanner data: Bottled Juice and Analgesics
###### From Dominick's supermarket Database
###### https://www.chicagobooth.edu/research/kilts/datasets/dominicks

###### Output from this script: Figure 4
###### Combined plot of results from Dominicks data on bottled juice and analgesics
###### Uses results from DFMultilateral_Analgesics.R and DFMultilateral_BottledJuice.R

library(rstudioapi) # extract path for setting up the working directory
library(openxlsx)

### set up working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

df_a7_bj<-readWorkbook("df_a7_bottled_juice.xlsx")
df_a8_an<-readWorkbook("df_a8_analgesics.xlsx")

pdf("Figure_4.pdf")
# Bottled Juice
layout(matrix(c(1,2),nrow=1), widths=c(1,1), heights=c(1.5,1.5), TRUE)
op<-par(oma=c(4,0,0,0)) # Room for the legend
#        mfrow = c(1, 2))
plot(df_a7_bj[,c("month","f_ch")], col=2, ylim=c(0.8, 1.42),
     ylab="Index", xlab="Month", main = "Bottled Juice",type="l", lwd=2, cex.main=1)
lines(df_a7_bj[,c("month","wtpd")], col=4, lwd=2)
lines(df_a7_bj[,c("month","gk")], col=6, lwd=2)
lines(df_a7_bj[,c("month","al")], col="violetred", lwd=2)
lines(df_a7_bj[,c("month","lq")], col=3, lwd=2)
lines(df_a7_bj[,c("month","ccdi")], col="gold", lwd=2)

######

# Analgesics
plot(df_a8_an[,c("month","f_ch")], col=2, ylim=c(0.8, 1.42),
     ylab="Index", xlab="Month", main="Analgesics", type="l", lwd=2, cex.main=1)
lines(df_a8_an[,c("month","wtpd")], col=4, lwd=2)
lines(df_a8_an[,c("month","gk")], col=6, lwd=2)
lines(df_a8_an[,c("month","al")], col="violetred", lwd=2)
lines(df_a8_an[,c("month","lq")], col=3, lwd=2)
lines(df_a8_an[,c("month","ccdi")], col="gold", lwd=2)

legend_texts = c(as.expression(bquote(FCH))
                 , as.expression(bquote(WTPD))
                 , as.expression(bquote(GK))
                 , as.expression(bquote(AL))
                 , as.expression(bquote(LQ))
                 , as.expression(bquote(CCDI))
)

# add a legend
par(op) # Leave the last plot
op <- par(usr=c(0,1,0,1), # Reset the coordinates
          xpd=NA)         # Allow plotting outside the plot region
legend(-1.4,-0.12 # Find suitable coordinates by trial and error
       , legend = legend_texts
       , col = c(2,4,6,"violetred", 3, "gold")
       , lty = 1, lwd=3, seg.len=1, x.intersp=0.3, cex=1, bty="n", horiz=T)

dev.off()


#########################################################################

