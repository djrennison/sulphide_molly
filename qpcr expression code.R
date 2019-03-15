#individual 362 was removed from the dataset because they only expressed RH2-1 

setwd("/Users/DJR/Desktop/PhD Data/Molly Project")
xx <- read.csv("combined_molly_qpcr_noind362.csv")
yy <- read.csv("Molly samples.csv")
library(dplyr)
head(xx)
xx$Drainage <- yy$Drainage[match(xx$Identifier, yy$Identifier)]
xx$H2S <- yy$H2S[match(xx$Identifier, yy$Identifier)]

xx %>% group_by(Gene, Identifier) %>% mutate(meanCT = mean(Ct), CT_stdev = sd(Ct))  %>% distinct(Gene, Identifier)-> averaged_data
averaged_data$Expression <- (1/((1 + averaged_data$Efficiancy)^ averaged_data$meanCT))
averaged_data %>% group_by(Identifier) %>% mutate(total_cone = sum(Expression)) %>% ungroup() %>% 
group_by(Identifier, Gene) %>% summarize(percent_cone = Expression/total_cone) -> percent_data
percent_data $Drainage <- yy$Drainage[match(percent_data $Identifier, yy$Identifier)]
percent_data $H2S <- yy$H2S[match(percent_data $Identifier, yy$Identifier)]

xy <- subset(percent_data, percent_data $Gene == "LWS3")

#significance testing #mixed effects model 
library(nlme)
z <- lme(percent_cone ~ H2S, random = ~1|Drainage, data =xy)
summary(z)
anova(z)

xy <- subset(percent_data, percent_data $Gene == "SWS2A")
xy <- subset(percent_data, percent_data $Gene == "SWS2B")
xy <- subset(percent_data, percent_data $Gene == "RH2-1")

library(ggplot2)

a <- ggplot(percent_data,aes(y= percent_cone, x=Drainage, fill= as.factor(H2S)))
a+ geom_boxplot() + facet_wrap("Gene", scales="free_y")


write.csv(percent_data, "temp_molly.csv")


 a <- ggplot(xy,aes(y= percent_cone, x=H2S, color= as.factor(H2S)))
  a+ geom_boxplot()

  a <- ggplot(xy,aes(y= percent_cone, x=as.factor(H2S), color= as.factor(Drainage)))

  a <- ggplot(xy,aes(y= percent_cone, x=Drainage, color= as.factor(H2S)))
  a <- ggplot(xy,aes(y= percent_cone, x=as.factor(H2S)))
Proportion a+ geom_boxplot()
 a+ geom_point() 
