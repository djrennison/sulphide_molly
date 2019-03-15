setwd("/Users/DJR/Desktop/")
x <- read.csv("Standards_plate1_threshold50_#265_extra_col.csv", stringsAsFactors=FALSE)
x <- read.csv("Standards_plate1_threshold50_#265_extra_col_test.csv", stringsAsFactors=FALSE)
x <- read.csv("Standards_plate1_threshold50_#265_extra_col_test2.csv", stringsAsFactors=FALSE)
y <- read.csv("Standards_plate1_v2.csv", stringsAsFactors=FALSE)
x <- read.csv("Standards_plate2v2_threshold50.csv", stringsAsFactors=FALSE) 

x <- read.csv("Standards_plate1newdilution_threshold50.csv", stringsAsFactors=FALSE) 
x <- read.csv("LWS_standards_threshold50.csv", stringsAsFactors=FALSE) 


head(x)
x$logdilution <- log(x$Dilution)
y$logdilution <- log(y$Dilution)


#SWS1

SWS1 <- subset(x, x$Gene=="SWS1")
plot(logdilution ~ Ct, data= SWS1)
z <- lm(logdilution ~ Ct, data= SWS1)
abline(z)
summary(z)
r2 = 0.9379
Ct         0.71424
(exp(0.71424))-1
E= 1.042634



#SWS2A
SWS2 <- subset(x, x$Gene=="SWS2A")
plot(logdilution ~ Ct, data= SWS2)
z <- lm(logdilution ~ Ct, data= SWS2)
abline(z)
summary(z)
r2 = 0.9518
Ct  0.62589
(exp(0.62589))-1
E= 0.8699094


#SWS2B
SWS2 <- subset(x, x$Gene=="SWS2B")

plot(logdilution ~ Ct, data= SWS2)
z <- lm(logdilution ~ Ct, data= SWS2)
abline(z)
summary(z)
r2 = 0.9933
Ct  0.63032
(exp(0.63032))-1
E= 0.8782115


RH2 <- subset(x, x$Gene=="RH2-2")

plot(logdilution ~ Ct, data= RH2)
z <- lm(logdilution ~ Ct, data= RH2)
abline(z)
summary(z)
r2 = 0.9974
Ct  0.6224
(exp(0.6224))-1
E= 0.8633948

######
	RH2 <- subset(x, x$Gene=="RH2-1")

plot(logdilution ~ Ct, data= RH2)
z <- lm(logdilution ~ Ct, data= RH2)
abline(z)
summary(z)
r2 = 0.9908
Ct  0.74249
(exp(0.74249))-1
E= 1.101161


LWS3 <- subset(x, x$Gene=="LWS3")

plot(logdilution ~ Ct, data= LWS3)
z <- lm(logdilution ~ Ct, data= LWS3)
abline(z)
summary(z)
r2 = 0.9505
Ct  0.67380
(exp(0.67380))-1
E= 0.9616775


LWS1 <- subset(x, x$Gene =="LWS1")

plot(logdilution ~ Ct, data= LWS1)
z <- lm(logdilution ~ Ct, data= LWS1)
abline(z)
summary(z)
r2 = 0.9961
Ct  0.68920
(exp(0.68920))-1
E= 0.9921212

LWS2 <- subset(x, x$Gene =="LWS2")

plot(logdilution ~ Ct, data= LWS2)
z <- lm(logdilution ~ Ct, data= LWS2)
abline(z)
summary(z)
r2 = 0.9993
Ct  0.646044
(exp(0.646044))-1
E= 0.9079779

LWSr <- subset(x, x$Gene =="LWSR")

plot(logdilution ~ Ct, data= LWSr)
z <- lm(logdilution ~ Ct, data= LWSr)
abline(z)
summary(z)
r2 = 0.9923
Ct  0.60522
(exp(0.60522))-1
E= 0.8316551





