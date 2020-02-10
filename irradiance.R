#
# R code to read the irradiance file, smoothen it and calculate 
# irradiance in Einsteins in 1nm bins across the wavelength.range (300-700)
#
# NOTES 
# code from proceedings paper, needs to be adjusted
# need meta data for the irradiance files. quite some code can be re-used
#
# TO DO
# make graphs of each location, side and downwelling (no upwelling)
# get integration time from original file
#
#

# CHECK WITH DIANA
# • changed irrandiance files: 
# Banos BanBan3_xxx.txt to Ban-3_xxx.txt => ok?
# LaBonita lab1_xxx.txt t0 lab-1_xxx.txt


#
# ======== irradiance analyses ================================
# 
# A reading all irradiance files, using rolling mean and spline curve to get intensity across the wavelength range
# B calculate the transmission for each location within each lake and write it to file
# C finsih up the irradiance file and write it to file
#
# NOTE:
# - medians are always given for specific light type and depth (bot irradiance and transmission)
# - irradiance overview is for all the light and depth data
#


library(zoo)


# === General variables ====
names.irr.axis <- c("wavelength", "intensity")
# rolling mean
window.w <- 5		# the size to the left/right of vocal point over which mean is taken
# irradiance median and for transmission
irr.depth <- 10
light.type.tr <- "s" # what light to use: s = sidewelling, u = upwelling, d = downwelling
exclude.bottom <- FALSE   # use samples when probe was at the bottom

# ===== A Irradiance ================
irr.pop.names <- list.files(path.data.irr)

#irr.meta.data <- read.csv(paste(path.meta.data,"irradiance.sites.info.all.csv", sep = "" ), stringsAsFactors = FALSE)

# list to store all the data sets in
irr.all <- vector("list", length(irr.pop.names))
names(irr.all) <- irr.pop.names

# names for the columns of the overview irradiance file
irr.file.names <- c("population","sample.location", "file.name", "light", "depth", "integr.time")  
  
for(i in 1:length(irr.pop.names)){
  print(i)
  # i <- 1
  # file names of file in folder i 
  file.names <- list.files(paste(path.data.irr, irr.pop.names[i], sep = ""))
  
  # make an overview file to store new irradiance measures in
  irrad.temp <- as.data.frame(matrix(NA, ncol = length(irr.file.names) + length(wavelength.range), nrow = length(file.names)))
  colnames(irrad.temp) <- c(irr.file.names, wavelength.range)
  # head(irrad.temp)
  

  # go through all the files in the population folder
  for(j in 1:length(file.names)){
    print(j)
    # j <- 1
    irr.t <- read.table(paste(path.data.irr, irr.pop.names[i], "/", file.names[j], sep = ""), fill = TRUE, stringsAsFactors = FALSE)
    
    # update the meta data in the irrad.temp file
    # population
    irrad.temp$population[j] <- gsub( "-.*", "", file.names[j])
    # sample location
    irrad.temp$sample.location[j] <- gsub(".*-", "", gsub( "_.*", "", file.names[j]))
    # file.name
    irrad.temp$file.name[j] <- file.names[j] 
    
    name.f.t <- gsub(".*_s*|.txt", "", file.names[j])
    # light type
    irrad.temp$light[j] <- substr(name.f.t, nchar(name.f.t), nchar(name.f.t))
    # depth
    irrad.temp$depth[j] <- as.integer(substr(name.f.t, 1, nchar(name.f.t) - 1))
    # integration time : obtain intration time from original irradiance file header (line 9)
    irrad.temp$integr.time[j] <- int.t <- as.integer(irr.t[9,4])
    
    # clean up the irradiance file
    if(is.numeric(irr.t[1,1]) == FALSE) {
      # data file is messy with a header section of 17 rows and only column 1 + 2 	with data neede
      irr.t.1 <- irr.t[-c(1:17, length(irr.t[,1])), c(1,2)]
      irr.t.1[,1] <- as.numeric(irr.t.1[,1])
      irr.t.1[,2] <- as.numeric(irr.t.1[,2])
    } else {
      irr.t.1 <- irr.t
      print("no text header rows")
    }
    colnames(irr.t.1) <- names.irr.axis
    # exclude > min(wavelength.range)  and > max(wavelength.range)
    irr.t.1a <- irr.t.1[ irr.t.1$wavelength >= min(wavelength.range) & irr.t.1$wavelength <= max(wavelength.range),]
    plot(irr.t.1a[,2] ~ irr.t.1a[,1], xlab = names.irr.axis[1], ylab = names.irr.axis[2], las = 1, type = "l")

    # convert to Einsteins
    irr.t.2 <- irr.t.1a
    irr.t.2[,2] <- f.Einstein.conv(irr.t.1a[,2], irr.t.1a[,1], int.t) 
    plot(irr.t.2[,2] ~ irr.t.2[,1], xlab = names.irr.axis[1], ylab = names.irr.axis[2], las = 1, type = "l")
    
    # get the rolling mean of the wavelength and the intensity
    irr.t.roll.mean <- data.frame(cbind(rollmean(irr.t.2$wavelength, window.w), rollmean(irr.t.2$intensity, window.w)))
    head(irr.t.roll.mean)
    colnames(irr.t.roll.mean) <- names.irr.axis
    # fit a spline to in the next step get intensity per nm
    fitted.spline <- splinefun(irr.t.roll.mean$wavelength, irr.t.roll.mean$intensity)
    new.intensity <- fitted.spline(wavelength.range)
    
    irrad.temp[ j, c((length(irr.file.names) + 1): ncol(irrad.temp))] <- new.intensity
    # irrad.temp[,1:12]
  }
  irr.all[[i]] <- irrad.temp
}


# combine th irradiances irradiance into on big file
irr.all.smooth <- do.call("rbind", irr.all)

#Check irradiance measures:
as.tibble(irr.all.smooth) %>%
  pivot_longer(c(-population,-sample.location,-file.name,-light,-depth,-integr.time),
               "lambda",values_to="irradiance") %>%
  filter(light == "s",depth == 10, population == "Ixt" | population == "Esp") %>% View()
  ggplot(.,aes(x=as.numeric(lambda),y=irradiance,group=file.name)) + 
  geom_line() +
  facet_wrap(~population)


# write data files as temp for now
write.csv(irr.all.smooth, paste( path.output.irr, "irr.all.smooth.w=", window.w, ".csv", sep = ""), row.names = FALSE)

# normalised version
# sum each row:
sum.per.row <- apply(irr.all.smooth[, -c(1:length(irr.file.names))], 1, sum)
irr.norm <- irr.all.smooth[, -c(1:length(irr.file.names))] / sum.per.row
apply(irr.norm, 1, sum) # correct, all sum to 1 now
irr.all.smooth.norm <- cbind.data.frame(irr.all.smooth[, c(1:length(irr.file.names))], irr.norm)
# write to file
write.csv(irr.all.smooth.norm, paste( path.output.irr, "irr.all.smooth.normalised.w=", window.w, ".csv", sep = ""), row.names = FALSE)


as.tibble(irr.all.smooth.norm)


# make figures for each location
depth.tmp <- 10
light.tmp <- "d"
pops <- unique(irr.all.smooth.norm$population)


irr.all.smooth.norm.tibble <- as.tibble(irr.all.smooth.norm) %>%
  pivot_longer(c(-population,-sample.location,-file.name,-light,-depth,-integr.time),
               "lambda",values_to="irradiance")

irr.all.smooth.norm.tibble %>%
  filter(light=="s") %>%
  ggplot(.,aes(x=as.numeric(lambda),y=irradiance,color=depth,group=file.name)) +
  geom_line() + 
  facet_grid(sample.location~population) + scale_color_distiller(palette = "RdYlBu")
for(i in 1:length(pops)){
  #i <- 1
  data.t <- irr.all.smooth.norm[irr.all.smooth.norm$population == pops[i],]
  data.t2 <- data.t[data.t$light == light.tmp & data.t$depth == depth.tmp, ]
  
  title.t <- paste(path.fig.irr, "irr.norm.light=",light.tmp,".depth=",depth.tmp,".pop=",pops[i],".pdf",sep="")
  max.t <- max(data.t2[,-c(1:length(irr.file.names))])
  ylim.t <- c(0, max.t)
  xlim.t <- range(wavelength.range)
  pdf(title.t, width = 7, height = 5)
  plot(NA, xlim = xlim.t, ylim = ylim.t, las = 1, ylab = "norm irradiance", xlab = "wavelength (nm)", main = pops[i] )
  for(j in 1:nrow(data.t2)) lines(wavelength.range, data.t2[j,-c(1:length(irr.file.names))])
  dev.off()
}




# %%%%%%% till here %%%%%%












# ===== B Transmission ================
#
# XX steps
# 1) calculate the kd from the smooth irradiance data
# 2) determine odd outlier kd values indicating wrong curve fitting and
#     replace this values with the mean of the sliding window
# 3) get the median, and compare all pops and ones without the odd ones
# 4) sliding window to smooth Kd
# 5) adjust range and normalise the smooth Kd
# 6) exclude bad recordings and save final version
#     • exclude bad recordings
#     • trim to fit the wave.length.range
#     • calc the median of the trimmed and  
#

# NOTE THE NORMALISATION DEPENDS ON WAVELENGTH RANGE!

# --- 1) calculate kd -------
# looks like the extreme variance in the lower wavelengths  kills the nls() convergence
# two solutions: skip the lower wavelegths or increase the window of the sliding window
#
#

# counter which goes through all irradiance files (may be the below function multiple times )
# This function is for a single site
params <- c(0, 10, 4000)
f.calc.transmission.old.code <- function(irrad.data, params){
   #params <- params
   #irrad.data <- irr.Ban
  trans.cntr <- 1
  
  site.names.1 <- unique(irrad.data$pop)
  loc.names.1 <- unique(irrad.data$sample.location)
  #Remove surf location
  loc.names.1 <- loc.names.1[loc.names.1 != "surf"]
  
  # make data frame for storage
  trans.names <- c(c("site.name","sample.location","light","n.measures","at.bottom","nls.converged"))
  transm.overview.k <- transm.overview.SS <- as.data.frame(matrix(NA, ncol = (length(trans.names) + length(wavelength.range)), nrow = length(loc.names.1)))
  colnames(transm.overview.k) <- colnames(transm.overview.SS) <- c(trans.names, wavelength.range)
  
  for(i in 1:length(site.names.1)){
     #i <- 1
    data.t <- irrad.data[irrad.data$pop == site.names.1[i],]
    for(j in 1:length(loc.names.1)){
      print(j)
      #j <- 1
      # select data
      data.t.1 <- data.t[data.t$sample.location == loc.names.1[j],]
      data.t.2 <- data.t.1[data.t.1$light == light.type.tr, ]
      data.t.3 <- data.t.2
      if(exclude.bottom == TRUE) data.t.3 <- data.t.2[data.t.2$at.bottom == FALSE, ] 
      
      # fill in meta data
      transm.overview.k$site.name[trans.cntr] <- transm.overview.SS$site.name[trans.cntr] <- site.names.1[i]
      transm.overview.k$sample.location[trans.cntr] <- transm.overview.SS$sample.location[trans.cntr] <-loc.names.1[j]
      transm.overview.k$light[trans.cntr] <- transm.overview.SS$light[trans.cntr] <- light.type.tr
      transm.overview.k$n.measures[trans.cntr] <- transm.overview.SS$n.measures[trans.cntr] <- nrow(data.t.3)   
      transm.overview.k$at.bottom[trans.cntr] <- transm.overview.SS$at.bottom[trans.cntr] <- exclude.bottom
      #transm.overview.k$pelagic[trans.cntr] <- transm.overview.SS$pelagic[trans.cntr] <- data.t.3$pelagic[1]
      
      # calculate the transmission coef for each wavelength
      #str(data.t.3)
      for(k in 1:length(wavelength.range)){
        print(k)
        # k <- 1
        # simplify variables: need depth (x.t) and intensity (y.t)
        x.t <- data.t.3$depth
        y.t <- data.t.3[, (length(irr.file.names) + k)]
        # nls to fit function
        # beer's law
        # k= attenuation coefficient, b = natural logarithm
        funct <- nls(y.t ~ b * exp(-k * x.t), start=list(k= params[1], b=params[2]), control=list(maxiter=params[3], warnOnly=TRUE))
        
        transm.overview.k[trans.cntr, (length(trans.names) + k)] <- coef(funct)[1]
        # to get an idea of the fit, save the SS residuals
        transm.overview.SS[trans.cntr, (length(trans.names) + k)] <- sum(resid(funct)^2)
        
        # plot the data and the fitted curve
        #pred.funct <- function(x.s) {coef(funct)[2] * exp(-coef(funct)[1] * x.s)}
        #x.s <- seq(min(x.t), max(x.t), 0.1)
        #plot(x.t, y.t, xlab = "Depth",ylab = "Intensity", pch = 19)
        #lines(x.s, pred.funct(x.s), col = "red")
      }
      
      # next sample location
      trans.cntr <- trans.cntr + 1
    }
  }
  
  # link both data frames into a list
  overview.list <- vector("list", 2)
  names(overview.list) <- c("k", "SS")
  overview.list[[1]] <- transm.overview.k
  overview.list[[2]] <- transm.overview.SS
  return(overview.list)
  # end of function
}  


populations <- unique(irr.all.smooth$population)
all_transmissions <- tibble(site.name=character(), sample.location=character(),
                            light=character(), n.measures=integer(),
                            nls.converged=character(), lambda=character(), 
                            transmission=double(),
                            transmission_ss=double())
for (pop in populations){
  tmp.irr <- irr.all.smooth[irr.all.smooth$population == pop,]
  tmp.k <- f.calc.transmission.old.code(tmp.irr, params)
  tmp.trans <- tmp.k[[1]]
  tmp.trans <- tmp.trans %>%
    dplyr::select(-at.bottom) %>%
    pivot_longer(c(-site.name,-sample.location,-light,-n.measures,-nls.converged),
                 names_to = "lambda", values_to = "transmission")
  tmp.SS <- tmp.k[[2]]
  tmp.SS <- tmp.SS %>%
    dplyr::select(-at.bottom) %>%
    pivot_longer(c(-site.name,-sample.location,-light,-n.measures,-nls.converged),
                 names_to = "lambda", values_to = "transmission_ss")
  all_transmissions <- rbind(all_transmissions, inner_join(tmp.SS,tmp.trans))
}
med.LC <- fig.light.measures(tLC, wave.length.range.temp, 3, "check.kd.all.")

pdf(paste(path.fig.subdir[path.subdir], "check.kd.all","_light.type_", light.type.tr,"_w_", window.w,".irr.depth_",irr.depth,".pdf",sep=""), height = 20, width = 20)
all_transmissions %>%
  ggplot(.,aes(x=as.numeric(lambda),y=transmission,group=sample.location)) +
  geom_line() + 
  facet_grid(sample.location~site.name) +
  ylab("Transmission") + xlab("Lambda")
dev.off()

write_tsv(all_transmissions,paste0(path.output.transmission,"light.type_", light.type.tr,"_w_", window.w,".irr.depth_",irr.depth,"txt"))



kd.all <- all_transmissions %>%
  #Remove one outlier locations
  filter(!(site.name == "VC" & sample.location == "4") ) %>%
  group_by(site.name,light,lambda) %>%
  dplyr::summarize(median_transmission = median(transmission)) 


pdf(paste(path.fig.subdir[path.subdir], "check.kd.median","_light.type_", light.type.tr,"_w_", window.w,".irr.depth_",irr.depth,".pdf",sep=""), height = 20, width = 20)
kd.all %>%
  ggplot(.,aes(x=as.numeric(lambda),y=median_transmission)) +
  geom_line() + 
  facet_wrap(~site.name) +
  ylab("Absorbance") + xlab("Lambda")
dev.off()

kd.all <- kd.all %>%
  group_by(site.name) %>%
  mutate(rolling_median_transmission=rollapply(median_transmission,(window.w*2),mean,align='left',fill=NA))


# --- 4) sliding window to smooth Kd -------
# - smoothen values and add columns to start and end to compensate for rollmean -
# taking a sliding window approach similar as done for the irradiance
# get the rolling mean of the wavelength and the intensity
k.roll.m <- overview.k[, c((ncol(overview.k) - length(wave.length.range.temp) + 1):ncol(overview.k))]
overview.k.smooth <- as.data.frame(matrix(NA, ncol = ncol(k.roll.m) - window.w + 1, nrow = nrow(k.roll.m) ))
for(i in 1:nrow(k.roll.m)) overview.k.smooth[i,] <- rollmean(unlist(k.roll.m[i,]), window.w)
head(overview.k.smooth)
overview.k.smooth[,1:10]
wave.length.range.smooth <- wave.length.range.temp[-c(1:((window.w - 1)/2),((length(wave.length.range.temp) - ((window.w - 1)/2) + 1):length(wave.length.range.temp)) )]
colnames(overview.k.smooth) <- wave.length.range.smooth
# measurements with negative values => problem with normalisation: check if this is a problem
rows.neg.Kd <- which(apply(overview.k.smooth, 1, min) < 0)
rows.neg.Kd.names <- paste(overview.k$site.name[rows.neg.Kd], overview.k$sample.location[rows.neg.Kd], sep =".")
# compare the effect of smoothing
# exclude the effect of window from start and end
w.effect <- (window.w -1 )/ 2
k.orig.t <- k.roll.m[, -c(1:w.effect, ( (ncol(k.roll.m) - w.effect + 1):ncol(k.roll.m)))]
plot(NA, xlim = range(wave.length.range.temp[-c(1:2, ((length(wave.length.range.temp) - 1):length(wave.length.range.temp)))]), ylim = c(-0.05, 0.05), las = 1)
for(i in 1:nrow(overview.k.smooth)) lines( wave.length.range.temp[-c(1:2, ((length(wave.length.range.temp) - 1):length(wave.length.range.temp)))], (k.orig.t[i,] - overview.k.smooth[i,]) )
# has quite an effect => good to keep the smoothing in

# --- 5) adjust range and normalise the smooth Kd -------
# --  (NOTE this will not go well with negative values, but we will exclude these at the end)
# first adjust the wavelength rage from wave.length.range.temp to wave.length.range
# wave.length.range.temp lost (window.w - 1)/2 from each side-> exclude the difference left with wave.length.range
overview.k.smooth.t <- overview.k.smooth[, match(wave.length.range, wave.length.range.smooth)]
ncol(overview.k.smooth.t)
# normalise the k smooth overview file
overview.k.smooth.norm <- as.data.frame(matrix(NA, nrow = nrow(overview.k.smooth.t), ncol = ncol(overview.k.smooth.t)))
for(i in 1:nrow(overview.k.smooth.t)) overview.k.smooth.norm[i,] <- overview.k.smooth.t[i,]/sum(overview.k.smooth.t[i,]) 
which(apply(overview.k.smooth.norm, 1, min) < 0)# same are negative as above, not surprising
str(overview.k.smooth.norm)

# ***standard***
overview.k.smooth.1 <- as.data.frame(cbind(overview.k[,1:length(trans.names)], overview.k.smooth.t) )
colnames(overview.k.smooth.1) <- c(trans.names, wave.length.range)
  write.csv(overview.k.smooth.1, paste(path.output.data.subdir[2],"k.all.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)

# get median
median.k.all <- fig.light.and.median (overview.k.smooth.1, wave.length.range, 3, "k.smooth.all")
  write.csv(median.k.all, paste(path.output.data.subdir[2],"median.k.all.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)


# ***normalised k***
overview.k.smooth.norm.1 <- as.data.frame(cbind(overview.k[,1:length(trans.names)], overview.k.smooth.norm) )
colnames(overview.k.smooth.norm.1) <- c(trans.names, wave.length.range)
 write.csv(overview.k.smooth.norm.1, paste(path.output.data.subdir[2],"k.all.norm.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)
# get median  
median.k.all.norm <- fig.light.and.median(overview.k.smooth.norm.1, wave.length.range, 4, "k.smooth.norm.all.")
  write.csv(median.k.all.norm, paste(path.output.data.subdir[2],"median.k.all.norm.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)

# --- 6) exclude bad recordings and save final version (wave.length.range) -------
# a) exclude the rows which are odd looking or go negative: do this for both the irradiance file and the transmission
excl.rows <- sort(unique(c(rows.weird.Kd,rows.neg.Kd)))

overview.k.smooth.s <- overview.k.smooth.1[-excl.rows,]
  write.csv(overview.k.smooth.s, paste(path.output.data.subdir[2],"k.s.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)
overview.k.smooth.norm.s <- overview.k.smooth.norm.1[-excl.rows,]
  write.csv(overview.k.smooth.norm.s, paste(path.output.data.subdir[2],"k.s.norm.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)

# median
median.k.s <- fig.light.and.median(overview.k.smooth.s, wave.length.range, 3, "k.smooth.s.")
  write.csv(median.k.s, paste(path.output.data.subdir[2],"median.k.s.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)
median.k.s.norm <- fig.light.and.median(overview.k.smooth.norm.s, wave.length.range, 4, "k.smooth.norm.s.")
  write.csv(median.k.s.norm, paste(path.output.data.subdir[2],"median.k.s.norm.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)

# check how much the medians differ
delta.m.s <- delta.m <- as.data.frame(matrix(NA, nrow=nrow(k.median), ncol = length(wave.length.range)))
for(i in 1:nrow(k.median)) {
  delta.m[i,] <- (k.median[i,-1] - k.median.s[i,-1])
  delta.m.s[i,] <- (k.median.norm[i,-1] - k.median.norm.s[i,-1])
}

# the standard one => big effects when you remove the odd ones
plot(NA, ylim = range(delta.m), xlim = range(wave.length.range), main = "difference standard - selected",las = 1)
for(i in 1:nrow(k.median)) lines(wave.length.range, delta.m[i,])
# the normalised one => small effect
plot(NA, ylim = range(delta.m.s), xlim = range(wave.length.range), main = "difference norm - norm selected",las = 1)
for(i in 1:nrow(k.median)) lines(wave.length.range, delta.m.s[i,])

# write data files as temp for now
#write.csv(overview.k.smooth, paste(path.output.data.subdir[1],"trans.all.smooth.temp.light_",light.type.tr,"_excl.bottom_",exclude.bottom,"_w=",window.w,".csv",sep = ""), row.names = FALSE)

# ==== C finalise clean up and normalisation and write irradiance final files =======================================
# prepare files to be saved
# - shorten wavelength range to wave.length.range if needed
# - calculate the median irradiance and normalised 

# --- irradiance ----
# which position is wave.length.range within  wave.length.range.temp
wave.pos <- match(wave.length.range, wave.length.range.temp)
# has all, raw irradiance data in it
irr.all.smooth.all <- irr.all.smooth[ ,c(1:length(irr.file.names), c(length(irr.file.names) + wave.pos))]
write.csv(irr.all.smooth.all, paste(path.output.data.subdir[1],"irr.all.smooth.all.excl.bottom_",exclude.bottom,".w_",window.w,".csv",sep = ""), row.names = FALSE)

# make a normalised version
irr.t.n <- irr.t <- irr.all.smooth.all[,c((length(irr.file.names) + 1):ncol(irr.all.smooth.all))]
for(i in 1:nrow(irr.s)) irr.t.n[i,] <- irr.t[i,]/sum(irr.t[i,]) 
irr.all.smooth.norm.all <- as.data.frame(cbind(irr.all.smooth.all[,c(1:length(irr.file.names))], irr.t.n))
head(irr.all.smooth.norm.all)
# check normalised => all good
apply(irr.all.smooth.norm.all[, c((ncol(irr.all.smooth.norm.all) - length(wave.length.range) + 1):ncol(irr.all.smooth.norm.all))], 1, sum)
write.csv(irr.all.smooth.norm.all, paste(path.output.data.subdir[1],"irr.all.smooth.norm.all.excl.bottom_",exclude.bottom,".w_",window.w,".csv",sep = ""), row.names = FALSE)

# make the version where the odd transmission samples are left out
names.to.be.removed <- sort(unique(c(rows.weird.Kd.names, rows.neg.Kd.names)))
names.irr.file <- paste(irr.all.smooth.all$site.name, irr.all.smooth.all$sample.location, sep=".")
# determine which rows have to go
rem.loc <- c()
for(i in 1:length(names.to.be.removed)) rem.loc <- c(rem.loc, which(names.irr.file == names.to.be.removed[i]))

# irradiance selected
irr.all.smooth.s <- irr.all.smooth.all[-rem.loc, ]
write.csv(irr.all.smooth.s, paste(path.output.data.subdir[1],"irr.all.smooth.s.excl.bottom_",exclude.bottom,".w_",window.w,".csv",sep = ""), row.names = FALSE)

irr.all.smooth.norm.s <- irr.all.smooth.norm.all[-rem.loc, ]
head(irr.all.smooth.norm.s)
write.csv(irr.all.smooth.norm.s, paste(path.output.data.subdir[1],"irr.all.smooth.norm.s.excl.bottom_",exclude.bottom,".w_",window.w,".csv",sep = ""), row.names = FALSE)


# get median irradiance for the four file types
# ** all **
# raw
irr.t1 <- irr.all.smooth.all[irr.all.smooth.all$depth == irr.depth & irr.all.smooth.all$light == light.type.tr,]
median.irr.all <- fig.light.and.median(irr.t1, wave.length.range, 1, "irr.all.smooth.")
write.csv(median.irr.all, paste(path.output.data.subdir[1],"median.irr.all.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)
# norm
irr.t2 <- irr.all.smooth.norm.all[irr.all.smooth.norm.all$depth == irr.depth & irr.all.smooth.norm.all$light == light.type.tr,]
median.irr.all.norm <- fig.light.and.median(irr.t2, wave.length.range, 2, "irr.all.smooth.norm")
write.csv(median.irr.all.norm, paste(path.output.data.subdir[1],"median.irr.all.norm.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)

# ** selection **
# raw
irr.t3 <- irr.all.smooth.s[irr.all.smooth.s$depth == irr.depth & irr.all.smooth.s$light == light.type.tr,]
median.irr.s <- fig.light.and.median(irr.t3, wave.length.range, 1, "irr.s.smooth.")
write.csv(median.irr.s, paste(path.output.data.subdir[1],"median.irr.s.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)
# norm
irr.t4 <- irr.all.smooth.norm.s[irr.all.smooth.norm.s$depth == irr.depth & irr.all.smooth.norm.s$light == light.type.tr,]
median.irr.s.norm <- fig.light.and.median(irr.t4, wave.length.range, 2, "irr.s.smooth.norm")
write.csv(median.irr.s.norm, paste(path.output.data.subdir[1],"median.irr.s.norm.light_",light.type.tr,".excl.bottom_",exclude.bottom,".depth_",irr.depth,".w_",window.w,".csv",sep = ""), row.names = FALSE)



















