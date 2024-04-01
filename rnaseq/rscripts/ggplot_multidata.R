gfg_data <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
                       y1 = c(1.1, 2.4, 3.5, 4.1, 5.9, 6.7,  
                              7.1, 8.3, 9.4, 10.0), 
                       y2 = c(7, 5, 1, 7, 4, 9, 2, 3, 1, 4), 
                       y3 = c(5, 6, 4, 5, 1, 8, 7, 4, 5, 4), 
                       y4 = c(1, 4, 8, 9, 6, 1, 1, 8, 9, 1), 
                       y5 = c(1, 1, 1, 3, 3, 7, 7, 10, 10, 10)) 

View(gfg_data)
gfg_plot <- ggplot(gfg_data, aes(x)) +   
  geom_line(aes(y = y1), color = "black") + 
  geom_line(aes(y = y2), color = "red") + 
  geom_line(aes(y = y3), color = "green") + 
  geom_line(aes(y = y4), color = "blue") + 
  geom_line(aes(y = y5), color = "purple") 
gfg_plot

gfg_data <- data.frame(x = c(1,2,3,4,5,6,7,8,9,10), 
                       y1 = c(1.1,2.4,3.5,4.1,5.9,6.7, 
                              7.1,8.3,9.4,10.0), 
                       y2 = c(7,5,1,7,4,9,2,3,1,4), 
                       y3 = c(5,6,4,5,1,8,7,4,5,4), 
                       y4 = c(1,4,8,9,6,1,1,8,9,1), 
                       y5 = c(1,1,1,3,3,7,7,10,10,10)) 

data_long <- melt(gfg_data, id = "x") 
gfg_plot <- ggplot(data_long,             
                   aes(x = x, 
                       y = value, 
                       color = variable)) +  geom_line() 
gfg_plot

pmeth <- read.csv("/Users/aron/Downloads/PercentMethyl_NHIPpromoter_Matrix.csv")
rownames(pmeth) <- pmeth$X
pmeth <- pmeth[,-c(1)]
pmeth$Y <- c("TD", "TD", "TD", "TD", "ASD", "ASD", "ASD", "ASD")
View(pmeth)
library(reshape2)
data_long <- melt(pmeth, id = "X") 
write.csv(data_long, "data_long.csv")
data_long <- read.csv("data_long.csv")
library(ggplot2)
ggplot(data_long, aes(x = variable, y = value, color = genotype, group = genotype)) + geom_point() + geom_smooth(method = "loess", se = FALSE) + ggtitle("% Methylation by CpG") + xlab("% Methylation") + ylab("CpG Coordinates") + labs(color = "Samples") + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

ggplot(data_long, aes(x = variable, y = value, color = Y, group = Y)) + geom_point() + geom_smooth(method = "loess") + ggtitle("% Methylation by CpG") + xlab("% Methylation") + ylab("CpG Coordinates") + labs(color = "Samples") + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
