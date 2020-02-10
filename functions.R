library(tidyverse)
library(yarrr)
# NOTES 
# code from proceedings paper, needs to be adjusted
#

# using:
#	5) f.Einstein.conv  convert units from W/m^2 to Einstein/m^2/sec







#
# code with functions used for modelling sensory ecology
#
#	1) the absorbance of light at different wavelength
#	2) log ratio transformation for matrix with opsin value 
#	3) opsin expression linear model (performs lin mod for each of the opsins at once)
#	4) f.bootstrap: bootstrap to test for differences in spectral sensitivity between two populations
#	5) f.Einstein.conv: convert units from W/m^2 to Einstein/m^2/sec
#	6) mean + 95% confidence interval the t distribution
#	7) spectral sensitivities from opsin expression data
#	8) lambda50 (or 1/2 integral across a vector)
# 9) spectral sensitivity (opsin express * absorbance curves)




## ========================================================================================================
# ======   5   convert watt/m^2 to Einstein/m^2/sec ======================================================
# ========================================================================================================

# function to convert absolute irradiance measure from W/m^2 to Einstein/m^2/s

# h = Planck’s constant (6.626 x 10^ -34 joule/sec),
h <-   6.626 * 10^-34
# c = speed of light (2.998 x 10 8 m/sec) 
c <- 2.998 * 10^8 
# Acc surface area of cosine corrector: diameter = 0.0039 m
Acc <- pi * (0.0039/2)^2
# L = Avogadro constant  6.02214 * 10 ^ 23 particles/mole
L <- 6.02214 * 10^23

# our irradiance is in W/m^2 -> J/s/m^2

# wave_length_range in nm
# integr.time in mu_sec

f.Einstein.conv <- function(irrad.data, wavelength.range.nm, integr.time.musec){
	wavelength.m <- wavelength.range.nm * 10^-9
	integr.time.sec <- integr.time.musec * 10^-6
	
	E <- h * (c / wavelength.m) / integr.time.sec
	# E is joules * m per second
	# plot(E)
	
	return( (irrad.data / E) * (1 / Acc) * (1 / L) )
}

#////////////////////////////////////////////////



# ============================================================================================
# ============================== 1 absorbance functions ======================================
# ============================================================================================

# --- spectral sensitivity ---

# R code to caculate the integrated sensitivity grpah for the four stickleback opsins
# multiply this with the expression of each opsin and get the predicted sensitivity for 
# A1 and A2 pigments.
#
#
# REFERENCES
# Calculating the 'integrated sensitivity graph' from
# In search of the visual pigment template
# - VICTOR I. GOVARDOVSKII, NANNA FYHRQUIST, TOM REUTER, DMITRY G. KUZMIN and KRISTIAN DONNER
# Visual Neuroscience / Volume 17 / Issue 04 / July 2000, pp 509 528 DOI: null, Published online: 
# September 2000
#
# The lambda max values for four opsin genes found in stickleback from 
# - Inigo Flamarique et al (J Evol Biol) Pronounced heritable variation and limited phenotypic plasticity in visual 
# pigments and opsin expression of threespine stickleback photoreceptorsin press
#
# Combine: cone expression with cone absorbance functions (mean = lambda max) 
# two calculations, the alpha band (S_alpha) eq 1 and beta band (S_beta) eq 4
#
# This has to be done for two cone type A1 and A2

# --------- SECTION A ----------------------------------------------------------------------------------
# write functions to calculate the A and B pigments expressions for A1 and A2 pigments
# for the A1 and A2 pigments

#ratio A1 and A2 opsins
# note, marines are alwawys 100% A1
ratio.A1_A2 <- c(1,0)
if(sum(ratio.A1_A2) != 1) print("WARNING: A1-A2 ratio is not 1!!!")

# function to calculate lambda max
opsin.max.range <- matrix(NA, nrow = 9, ncol = 2)
rownames(opsin.max.range) <- opsin.names
colnames(opsin.max.range) <- c("A1_0", "A1_1")
opsin.max.range[1,] <- c(571, 571 ) #For LWS-S
opsin.max.range[2,] <- c(516, 516)
opsin.max.range[3,] <- c(519, 519)
opsin.max.range[4,] <- c(NA, NA)
opsin.max.range[5,] <- c(516, 516 )
opsin.max.range[6,] <- c(476, 476)
opsin.max.range[7,] <- c(438, 438)
opsin.max.range[8,] <- c(408, 408)
opsin.max.range[9,] <- c(353, 353)

f.calc.lamba.max <- function(A1.proportion){
  l.max.temp <- rep(NA,9)
  for(l.m in 1:9) l.max.temp[l.m] <- opsin.max.range[l.m, 1] + ((opsin.max.range[l.m, 2] - opsin.max.range[l.m, 1]) * A1.proportion)
  return(l.max.temp)
}

opsin.l.max <- f.calc.lamba.max(ratio.A1_A2[1])

# --- PARAMETERS ---
# --- A1 -----------
# A1 A band
A1 <- 69.7
B1 <- 28
b1 <- 0.922
C1 <- -14.9
c1 <- 1.104
D1 <- 0.674  

# B band
A1_beta <- 0.26

# --- A2 -----------
B2 <- 20.85
b2 <- 0.9101
C2 <- -10.37
c2 <- 1.1123
D2 <- 0.5343  

A2_beta <- 0.37

# ------------------------------------------------------------------------------------------------------
# --- FUNCTIONS FOR THE ALPHA AND BETA BANDS -----------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# --- A1 pigment -----------
# the A band (eq 1, 2)
# x is the lambda_max of the opsin/wavelength of interest takes A1_alpha_input (x = lambda_max/lambda)
A1.alpha.f <- function(x, lambda_max) { 
	a1 <- 0.8795 + 0.0459 * exp(-(lambda_max - 300) ^2 / 11940)				# eq 2
	return(1/(exp(A1 * (a1 - x)) + exp(B1 * (b1 - x)) + exp(C1 * (c1 - x)) + D1))	# eq 1
	}

# the B band (eq 4, 5a, 5b)
# NOTE lambda_max is for each opsin gene different!
# y is the wavelength range
A1.beta.f <- function(y, lambda_max) {
	lambda_mbeta <- 189 + 0.315 * lambda_max		# eq 5a
	b1beta <- -40.5 + 0.195 * lambda_max 			# eq 5b
	return(A1_beta * exp(-((y - lambda_mbeta)/ b1beta)^2))	# eq 4
	}

# --------------------------
# --- A2 pigment -----------
# the A band (eq 1, 6a, 6b)
A2.alpha.f <- function(x, lambda_max) { 
	A2 <- 62.7 + 1.834 * exp((lambda_max - 625)/54.2)	# eq 6a
	a2 <- 0.875 + 0.0268 * exp((lambda_max - 665)/40.7)	# eq 6b
	return(1/(exp(A2 * (a2 - x)) + exp(B2 * (b2 - x)) + exp(C2 * (c2 - x)) + D2))	# eq 1
	}

# the B band (eq 4, 8a, 8b)
# NOTE lambda_max is for each opsin gene different!
# y is the wavelength range
A2.beta.f <- function(y, lambda_max) {
	lambda_mbeta2 <- 216.7 + 0.287 * lambda_max						# eq 8a
	b2beta <- 317 - 1.149 * lambda_max + 0.00124 * (lambda_max^2)	# eq 8b
	return(A2_beta * exp(-((y - lambda_mbeta2)/ b2beta)^2))			# eq 4
	}

# the calculations use lambda_max/wavelength value (vectore of values for wavelength of interest range)
# specific for each opsin (each on one row)
wave.length.range <- seq(350,700)
wave.length.range <- wavelength.range
A.alpha.input <- matrix(NA, 9, length(wave.length.range))
for(i in 1:9) A.alpha.input[i, ] <- opsin.l.max[i]/wave.length.range
rownames(A.alpha.input) <- opsin.names
colnames(A.alpha.input) <- wave.length.range 

# make storage for 9 absorbance measures for all 9 opsins
function.names <- c("A1.alpha", "A1.beta","A2.alpha","A2.beta")
A1.alpha <- A1.beta <- A2.alpha <- A2.beta <- matrix(NA, 9, length(wave.length.range))
rownames(A1.alpha) <- rownames(A1.beta ) <- rownames(A2.alpha) <- rownames(A2.beta) <- opsin.names
colnames(A1.alpha) <- colnames(A1.beta ) <- colnames(A2.alpha) <- colnames(A2.beta) <- wave.length.range 
A1.alpha[,1:9]

# this can be coded more elegantly.....
# c("LWS","RH2","SWS1","SWS2")
# A1 A band
for(i in 1:length(opsin.names)) A1.alpha[i,] <- A1.alpha.f(A.alpha.input[i,], opsin.l.max[i])
# A1 B band
for(i in 1:length(opsin.names)) A1.beta[i,] <- A1.beta.f(wave.length.range, opsin.l.max[i])
# A2 A band
for(i in 1:length(opsin.names)) A2.alpha[i,] <- A2.alpha.f(A.alpha.input[i,], opsin.l.max[i])
# A2 B band
for(i in 1:length(opsin.names)) A2.beta[i,] <- A2.beta.f(wave.length.range, opsin.l.max[i])

# four matrices each for one of the bands
# multiply each by their ration (weight) and sum them into a final new matrix
opsin.absorbance <- matrix(NA, 9, length(wave.length.range))
rownames(opsin.absorbance) <- opsin.names
colnames(opsin.absorbance) <- wave.length.range
for(i in 1:length(opsin.names)) opsin.absorbance[i,] <- c((ratio.A1_A2[1] * A1.alpha[i,]) + (ratio.A1_A2[1] * A1.beta[i,]) + (ratio.A1_A2[2] * A2.alpha[i,]) + (ratio.A1_A2[2] * A2.beta[i,]))



basel_palette <- unname(piratepal(palette = "basel"))
pdf(paste(path.fig.subdir[4],"absorbance.curves.opsins.prop.A1_", ratio.A1_A2[1],".pdf",sep=""), 5, 5)
opsin.absorbance.tibble <- as.tibble(opsin.absorbance)
opsin.absorbance.tibble$opsin <- row.names(opsin.absorbance)
opsin.absorbance.tibble <- opsin.absorbance.tibble %>%
  pivot_longer(-opsin, names_to = "lambda", values_to = "absorbance") 
opsin.absorbance.tibble %>%
  ggplot(.,aes(x=as.numeric(lambda),y=absorbance,group=opsin,color=opsin)) + geom_line(size=2) +
  scale_color_manual(values=basel_palette,name="Opsin") +
  ylab("Absorbance") + xlab("Lambda")

dev.off()

write.csv(opsin.absorbance, paste("output/opsin/","absorbance.prop.A1_",ratio.A1_A2[1],".csv", sep=""))


# ============================================================================================
# ============================== 2 log ratio transformation ==================================
# ============================================================================================

# ---- log ratio transform function----
# used Kucera & Malmgren 1998 for the formulas
# take  matrix with in each row and indiv and each column is an opsin
transform.sum.constr <- function(data){
	
# data <- overview_opsins_ind[,7:ncol(overview_opsins_ind)]
  transf.data <- as.data.frame(matrix(NA, ncol = ncol(data), nrow = nrow(data)))
  colnames(transf.data) <- paste(colnames(data),".t", sep = "")
  
  f.geom.mean <- function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  nrow(transf.data)
  
  for(iv in 1:nrow(data)){
    g.m <- f.geom.mean(data[iv,])
    transf.data[iv, ] <- log(data[iv,]/g.m)
  }
  
 	return(transf.data)	
}


# ============================================================================================
# ============================== 3 opsin linear model function ===============================
# ============================================================================================

# function which takes opsin expression for two factor levels
# tests if irradiance difference can explain opsin expression
lm.opsin.wavelength <- function(opsin.exp, transm.diff, opsin.absorbance){
	# NOTE this function takes for the opsin expression as first column the species (factor)
	# and in the following columns the four opsin expression values
	# make storage data frame
	d.opsin.temp <-as.data.frame(matrix(NA, 4, length(wave_length_range)))

	# calculate the integral of the absorbance x the differences in transmission
	for(i in 1:nrow(opsin.absorbance)) d.opsin.temp[i,] <- opsin.absorbance[i,] * transm.diff
	sum.d.temp <- apply(d.opsin.temp,1,sum)
		names(sum.d.temp) <- opsin.l.max_names
	# multiply opsin expression with above measure
	adj.ops <- opsin.exp
	for(i in 1:4) adj.ops[, 1+i] <- opsin.exp[,1+i]* sum.d.temp[i]
	
	# perform linear model for each opsin to test if level of 
	# opsin expression adjusted for transmission differences differs between the two groups	
	m.LWS <- as.data.frame(summary(lm(adj.ops[,3]~ adj.ops$Species))$coefficients)
	m.RH2 <- as.data.frame(summary(lm(adj.ops[,4]~ adj.ops$Species))$coefficients)
	m.SWS1 <- as.data.frame(summary(lm(adj.ops[,5]~ adj.ops$Species))$coefficients)
	m.SWS2 <- as.data.frame(summary(lm(adj.ops[,6]~ adj.ops$Species))$coefficients)
	
	m.combined <- sapply(c("m.LWS","m.RH2", "m.SWS1", "m.SWS2"), get, environment(), simplify = FALSE)
		
	return(m.combined)
}


# function which does a linear model for all four opsin with one factor
lm.opsin <- function(opsin.exp){
	# NOTE this function takes for the opsin expression as first column the species (factor), second the lake
	# and in the following columns the four opsin expression values
	
	#opsin.exp <- opsin.t
	
	# perform linear model for each opsin to test if level of 
	# opsin expression adjusted for transmission differences differs between the two groups	
#	m.LWS <- as.data.frame(summary(lm(opsin.exp[,3]~ opsin.exp$Species))$coefficients)
	summary(lm(opsin.exp[,3]~ opsin.exp$Species))
#	m.RH2 <- as.data.frame(summary(lm(opsin.exp[,4]~ opsin.exp$Species))$coefficients)
#	m.SWS1 <- as.data.frame(summary(lm(opsin.exp[,5]~ opsin.exp$Species))$coefficients)
#	m.SWS2 <- as.data.frame(summary(lm(opsin.exp[,6]~ opsin.exp$Species))$coefficients)
	
	
	# version with the summary output
	m.LWS <- as.data.frame(summary(lm(opsin.exp[,3]~ opsin.exp$Species))$coefficients)
	m.LWS2 <- summary(lm(opsin.exp[,3]~ opsin.exp$Species))
	
	m.RH2 <- as.data.frame(summary(lm(opsin.exp[,4]~ opsin.exp$Species))$coefficients)
	m.RH22 <- summary(lm(opsin.exp[,4]~ opsin.exp$Species))
	
	m.SWS1 <- as.data.frame(summary(lm(opsin.exp[,5]~ opsin.exp$Species))$coefficients)
	m.SWS12 <- summary(lm(opsin.exp[,5]~ opsin.exp$Species))

	m.SWS2 <- as.data.frame(summary(lm(opsin.exp[,6]~ opsin.exp$Species))$coefficients)
	m.SWS22 <- summary(lm(opsin.exp[,6]~ opsin.exp$Species))
	
	m.combined <- sapply(c("m.LWS","m.LWS2","m.RH2","m.RH22", "m.SWS1","m.SWS12", "m.SWS2","m.SWS22"), get, environment(), simplify = FALSE)
		
	return(m.combined)
}



# function which takes opsin expression for two factor levels
# tests if irradiance difference can explain opsin expression
# opsin.exp = fixed effect (factor), random effect (factor), four opsins
# uses lme 4
glm.opsin <- function(opsin.exp){
	#opsin.exp <- data.mar.fresh.1
	colnames(opsin.exp) <- c("fixed.effect", "random.effect","LWS", "RH2", "SWS1", "SWS2")
	opsin.exp[,1] <- as.factor(opsin.exp[,1])
	opsin.exp[,2] <- as.factor(opsin.exp[,2])
	#str(opsin.exp)
	
	# perform generalised linear model for each opsin to test if level of 
	# opsin expression adjusted for transmission differences differs between the two groups	
	lws <- lmer(LWS ~ fixed.effect + (1|random.effect), data= opsin.exp)
	lws0 <- lmer(LWS ~ 1+(1|random.effect), data= opsin.exp)
	gm.LWS <- as.data.frame(anova(lws, lws0))
	
	rh2 <- lmer(RH2 ~ fixed.effect +(1|random.effect), data= opsin.exp)
	rh20 <- lmer(RH2 ~ 1+(1|random.effect), data= opsin.exp)
	gm.RH2 <- as.data.frame(anova(rh2, rh20))
	
	sws1 <- lmer(SWS1 ~ fixed.effect +(1|random.effect), data= opsin.exp)
	sws10 <- lmer(SWS1 ~ 1+(1|random.effect), data= opsin.exp)
	gm.SWS1 <- as.data.frame(anova(sws1, sws10))

	sws2 <- lmer(SWS2 ~ fixed.effect +(1|random.effect), data= opsin.exp)
	sws20 <- lmer(SWS2 ~ 1+(1|random.effect), data= opsin.exp)
	gm.SWS2 <- as.data.frame(anova(sws2, sws20))

	
	#gm.combined <- sapply(c("gm.SWS1", "gm.SWS2","gm.RH2", "gm.LWS"), get, environment(), simplify = FALSE)
	gm.combined <- sapply(c("gm.SWS1", "gm.SWS2", "gm.RH2", "gm.LWS"), get, environment(), simplify = FALSE)
	return(gm.combined)
}

# tests if irradiance difference can explain opsin expression
# opsin.exp = fixed effect (factor), random effect (factor), four opsins
# uses nlme 
glm.nlme.opsin <- function(opsin.exp){
	# uses nlme: check if library is loaded

  #opsin.exp <- data.mar.fresh.1
  colnames(opsin.exp) <- c("fixed.effect", "random.effect","LWS", "RH2", "SWS1", "SWS2")
  opsin.exp[,1] <- as.factor(opsin.exp[,1])
  opsin.exp[,2] <- as.factor(opsin.exp[,2])
  #str(opsin.exp)

	# perform generalised linear model for each opsin to test if level of 
	# opsin expression adjusted for transmission differences differs between the two groups	
	lws <- lme(LWS ~ fixed.effect, random = ~1|random.effect, data= opsin.exp)
	gm.LWS <- as.data.frame(anova(lws))
	gm.LWS2 <- summary(lws)
	
	rh2 <- lme(RH2 ~ fixed.effect, random = ~1|random.effect, data= opsin.exp)
	gm.RH2 <- as.data.frame(anova(rh2))
	gm.RH22 <- summary(rh2)
	
	sws1 <- lme(SWS1 ~ fixed.effect, random = ~1|random.effect, data= opsin.exp)
	gm.SWS1 <- as.data.frame(anova(sws1))
	gm.SWS12 <- summary(sws1)

	sws2 <- lme(SWS2 ~ fixed.effect, random = ~1|random.effect, data= opsin.exp)
	gm.SWS2 <- as.data.frame(anova(sws2))
	gm.SWS22 <- summary(sws2)
	
	gm.combined <- sapply(c("gm.LWS","gm.LWS2","gm.RH2","gm.RH22", "gm.SWS1","gm.SWS12", "gm.SWS2","gm.SWS22"), get, environment(), simplify = FALSE)
	gm.combined <- sapply(c( "gm.SWS1","gm.SWS12", "gm.SWS2","gm.SWS22","gm.RH2","gm.RH22","gm.LWS","gm.LWS2"), get, environment(), simplify = FALSE)
	
	return(gm.combined)
}



# ========================================================================================================
# ======  4  BOOTSTRAP FUNCTIONS =========================================================================
# ========================================================================================================

# make bootstrap function which takes two data sets
# enter: two matrices/data frames(>) to sample from and number of iterations
# returns the mean t value across all wavelengths for the A B contrasts (position 1)
# 	and the sample distributions, (positions 2:n+1)
f.bootstrap <- function(A, B, n){
	t.values <- rep(NA, c(n+1))
	names(t.values) <- c("t",paste(c(1:n))) 
	
	t.test.tmp <- rep(NA, ncol(A))	# to store t values
	for(i in 1:ncol(A))  t.test.tmp[i] <- as.numeric(t.test(A[, i], B[, i])$statistic)
	t.values[1] <- mean(t.test.tmp)
	
	# derive the sample distribution
	A.B <- rbind(A,B)	# combine both matrices to get one big one to randomly resample
	for(i in 2:(n + 1)){
		# need to reshuffle each of the columns on A.B and than split it in two and redo t test
		random.dr <- sample(c(1:nrow(A.B)),nrow(A.B),replace=FALSE )	# gives me vector of random numbers of length nrow
		# random.dr[ c(1: nrow(A))] selects the first set of random numbers
		# random.dr[ c( (nrow(A)+1): nrow(A.B))] the second set
		# use these 'reshuffled' rows for the next sequence of t tests, one for each column
		for(j in 1:ncol(A))  t.test.tmp[j] <- as.numeric(t.test(A.B[random.dr[ c(1: nrow(A))], j], A.B[random.dr[ c( (nrow(A)+1): nrow(A.B))], j])$statistic)
		t.values[i] <- mean(t.test.tmp)
	}
	return(t.values)	
}


# This one is pretty good, taking advantage of the definition:
#   Var[x] = E[x]^2 - E[x^2]
# We could do the same with
#   Var[x] = E[(x - E[x])^2]
# of course.
colVars <- function(x, ...) {
  Ex <- colMeans(x)
  Ex2 <- colMeans(x * x)

  n <- nrow(x)
  # Last term here converts from population to sample variance
  (Ex2 - Ex * Ex) * n / (n - 1)
}


# way! FASTER version by Rich
## Final version, with a little cleaning, but no more speed tweeks.
# use matrices!
# NOTE: this now returns t^2
t.test.cols <- function(X, Y) {
  if (!is.matrix(X) || !is.matrix(Y) || ncol(X) != ncol(Y))
    stop("X and Y must be matrices with equal numbers of columns")
  nx <- nrow(X)
  ny <- nrow(Y)

  mx <- colMeans(X)
  my <- colMeans(Y)

  vx <- colVars(X)
  vy <- colVars(Y)
	# 20150204
	# need to take the square => added ^2 
  ((mx - my)/sqrt(vx/nx + vy/ny))^2
}

# I changed this to return the test on the data separately from the
# boostrapped values.  Access it with ans$data and ans$bootstrap, or
# unlist to get the structure you had before.
f.bootstrap.f <- function(A, B, n) {
  # derive the sample distribution
  # combine both matrices to get one big one to randomly resample
  A.B <- rbind(A,B)
  f <- function() {
    i <- sample(nrow(A.B))
    i.A <- i[ seq_len(nrow(A))]
    i.B <- i[-seq_len(nrow(A))]
    mean(t.test.cols(A.B[i.A,], A.B[i.B,]))
  }
  list(data=mean(t.test.cols(A, B)),
       bootstrap=replicate(n, f()))
}





# ========================================================================================================
# ======   6   mean + 95% confidence interval the t distribution =========================================
# ========================================================================================================

# function to get mean and 95% confidence interval using the t distribution
f.confint.int <- function(data){
	n <- length(data)
	mean.data <- mean(data, na.rm = TRUE)
	CI <- qt(0.95, (n - 1)) * sd(data, na.rm = TRUE) / sqrt(n)
	return(c(mean.data, CI))
}

# ========================================================================================================
# ======   7   spectral sensitivities from opsin expression data =========================================
# ========================================================================================================

# should be in format:  Species       LWS       RH2      SWS1       SWS2

# function to calculate the spectral sensitivity of a set of indidivudals (opsin.data)

#f.opsin.sens <- function(opsin.data){
#	spec.sens.combined <-  as.data.frame(matrix(NA, nrow=nrow(opsin.data), ncol=length(wave_length_range)))
#	for(ctr.i in 1:nrow(opsin.data)){
#		# loop through all opsins
#		for(ctr.j in 1:nrow(opsin.absorbance)){
#			 spec.sens.temp[ctr.j,] <- opsin.data[ctr.i, ctr.j + 1] * opsin.absorbance[ctr.j,]
#			 }
#			 spec.sens.combined[ctr.i,] <- apply(spec.sens.temp, 2, sum)
#			#plot(unlist(spec.sens.combined[9,]))
#	}
#	return(spec.sens.combined)
#}


# ========================================================================================================
# ======   8   lambda 50 =================================================================================
# ========================================================================================================
# calculate half of the full integral of response variable y against vector x
# rich helped me write this

f.lambda50 <- function(response.variable){

	# make sure both variables are numeric
	x <- as.numeric(wave.length.range)
	y <- as.numeric(response.variable)
				
	dx <- diff(x)
	mx <- (x[-1] + x[-length(x)])/2
	my <- (y[-1] + y[-length(y)])/2
		
	#plot(mx, cumsum(my * dx))
	#sum(my * dx)
		
	lambda_50 <- uniroot(approxfun(mx, cumsum(my * dx) - sum(my * dx)/2), range(mx))$root
	return(lambda_50)
}


# ========================================================================================================
# ====== 9 spectral sensitivity ==========================================================================
# ========================================================================================================

# function to calculate the spectral sensitivity of a set of indidivudals (opsin.data)
# should be in format:  Species     SWS1       SWS2   LWS       RH2
# absorbance has same order but in rows, opsin names on the first row.
f.opsin.sens <- function(opsin.data, absorbance.data){
	#opsin.data <- comtr.t
  #absorbance.data <- absorbance.contrast
  
  # delete the first column of the absorbanc
  absorbance.data <- absorbance.data[,-1]
  
  spec.sens.temp <- as.data.frame(matrix(NA, nrow=nrow(absorbance.data), ncol=ncol(absorbance.data)))
	spec.sens.combined <- as.data.frame(matrix(NA, nrow=nrow(opsin.data), ncol=length(wave.length.range)))
	colnames(spec.sens.combined) <- wave.length.range
	for(ctr.i in 1:nrow(opsin.data)){
	  # loop through all opsins 
		for(ctr.j in 1:nrow(absorbance.data)){
		  	spec.sens.temp[ctr.j,] <- opsin.data[ctr.i, ctr.j + 1] * absorbance.data[ctr.j,]
		}
			 spec.sens.combined[ctr.i,] <- apply(spec.sens.temp, 2, sum)
			#plot(unlist(spec.sens.combined[9,]))
	}
	# add the species codes to check if all lines up fine => works well
	#spec.sens.combined <- cbind(opsin.data[,1], spec.sens.combined)

	return(spec.sens.combined)
}


# ====== 10 get median from  ==========================================================================
f.get.median <- function(lightdata, wavelength){
  #lightdata <- overview.k
  #wavelength <- wave.length.range.temp
  
  sample.sites <- unique(lightdata$site.name)
  median.all <- data.frame(matrix(NA, ncol = length(wavelength), nrow = length(sample.sites)))
  for(tr in 1:length(sample.sites)){
    #tr <- 1
    dt <- lightdata[lightdata$site.name == sample.sites[tr], c((ncol(lightdata) - length(wavelength) + 1):ncol(lightdata))]
    median.all[tr,] <- apply(dt, 2, median)
  }
  median.all.sites <- cbind(sample.sites,median.all)
  colnames(median.all.sites) <- c("site.name", wavelength)
  return(median.all.sites)
}




