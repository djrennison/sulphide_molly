
#

#par(mfrow = c(1,1))


# rm(list = ls())

# get working directory
workdir.path <- getwd()

# soruce the following R scripts so they can be accessed throughout
# the project
source("figures.R")
source("functions.R")

# libraries
# CHECK WHAT IS NEEDED
library(tidyverse)
library(cowplot)
library(lme4)
library(lmerTest)
library(robCompositions)
library(rgr)
library(zoo)
library(car)
library(CCA)
library(MASS)
library(LMERConvenienceFunctions)
library(Hmisc)



# ===== FILE MANAGEMENT =======
# ====- folder + path structure -----
# names of folders for output data (figures + data output)
folder.names <- c("figures","data", "output")
for(i in 1:length(folder.names)) if(file.exists(folder.names[i]) == FALSE) dir.create(folder.names[i])

# within output
output.folder <- c("irradiance", "colour","transmission","opsin")
for(i in 1:length(output.folder)) dir.create(paste(workdir.path, "/output/", output.folder[i], sep = ""))

# within figures
figures.folder <- c("irradiance", "colour","transmission","opsin")
for(i in 1:length(figures.folder)) dir.create(paste(workdir.path, "/figures/", figures.folder[i], sep = ""))



# === paths to folders ====
# paths to the data folder
path.data <- paste(workdir.path, "/data/", sep = "")
path.data.irr <- paste(workdir.path, "/data/irradiance_with_meta/", sep = "")
path.data.col <- paste(workdir.path, "/data/colour_with_meta/", sep = "")
path.data.meta.d <- paste(workdir.path, "/data/meta_data/", sep = "")
path.data.clean.expr <- paste(workdir.path, "/data/cleaned_expression_data/", sep = "")



# to the figures folder
path.fig <- paste(workdir.path, "/figures/", sep = "")
path.fig.irr <- paste(workdir.path, "/figures/irradiance/", sep = "")
path.fig.trans <- paste(workdir.path, "/figures/transmission/", sep = "")
path.fig.col <- paste(workdir.path, "/figures/colour/", sep = "")
path.fig.subdir <- c(path.fig.col, path.fig.irr,path.fig.trans)

# to the output folder
path.output <- paste(workdir.path, "/output/", sep = "")
path.output.irr <- paste(workdir.path, "/output/irradiance/", sep = "")
path.output.colour <- paste(workdir.path, "/output/colour/", sep = "")
path.output.transmission <- paste(workdir.path, "/output/transmission/", sep = "")


# =============================
# ====== global variables =====
# =============================
wavelength.range <- c(350:700)
opsin.names <- c("LWS1","LWS2","LWS3","LWSr","RH2-1","RH2-2","SWS2A","SWS2B","SWS1")
opsin.col <- c("red","red","red","red","green","green","blue","blue","purple")

# for the irradiance.R files
wavelength.irr <- c(350:700)  # wavelength range of the irradiance measurements || depends on calib files
window.w <- 10                # width of the window for the rolling mean

UV.blue <- c(301:400)         # green orange from Brock et al => they start at 300, but does not match our calib file
gre.or <- c(501:600)          # green orange from Brock et al




