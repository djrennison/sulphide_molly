#
# ======== Figures =========
#
# 1 fig.light.and.median
# 2 fig.light.measures
# 3 fig.sens.mar.fr
# 4 fig.sens.lim.ben
# 5 fig.slope.mar.fresh.in.fresh
# 6 fig.slopes.dS.dI.dK

jit.2 <- 0.75


# ==== global variabled ====
colour <- c("black", "grey","blue","red")
colour.h <- c("#000000", "#999999","#0000FF","#FF0000")
alpha <- 90

# ==== 1 median and light ======
# takes the data frame with light measurement
# uses the wave.length.range to get the light data
# uses the standard names to selection location and site
# • makes figures to fig subfolders for each location
# • returns  median measure per location to output data folder
# input:
# light.data    irradiance or transmission
# wavelength    wavelength range (typically (wave.length.range))
# irr.depth     depth at which irradiance is taken
# light.type.side    which light type to select (sam as transmission! light.type.tr)
# trans.or.irr      type of light data (e.g. irradiance.norm )
# path.subdir       number of the figures subdir (1 irr 2 irr.norm 3 trans 4 tr.norm 5 absorb)
# NOTE start the trans.or.irr with tra or irr => needed for data selection

fig.light.and.median <- function( light.data, wavelength, path.subdir, fig.name){
  #  light.data <-overview.k.norm.smooth
  #  wavelength <- wave.length.range.temp
  #  path.subdir <- 4
  #  fig.name <- "test"
  
  str(light.data)
  # go through the location
  names <- unique(light.data$site.name)
  
  # make overview medians
  median.names <- c("site.name")  
  median.overview <- as.data.frame(matrix(NA, ncol = (length(median.names) + length(wavelength)), nrow = length(names) ) )
  colnames(median.overview) <- c(median.names, wavelength)
  
  for(f in 1:length(names)){
    # f <- 1
    data.f <- light.data[light.data$site.name == names[f],]
    str(data.f)
    data.f[,1:10]
    # get only the light data
    l.d <- data.f[,c((ncol(data.f) - length(wavelength) + 1):ncol(data.f))]
    
    # write meta data to median file
    median.overview[f,1] <- names[f]
    # get the median
    median.t <- apply(l.d, 2, median)
    median.overview[f,c((length(median.names) + 1):ncol(median.overview))] <- median.t
    ylim.t <- range(l.d,na.rm = TRUE)
    pdf(paste(path.figures.subdir[path.subdir], fig.name,"_", names[f], "_light.type_", light.type.side,".irr.depth_",irr.depth,".pdf",sep=""), height = 5, width = 10)
    plot(NA, xlim = range(wavelength), ylim = ylim.t, ylab = trans.or.irr, xlab = "Wavelength (nm)", main = names[f] , las = 1)
    for(g in 1:nrow(data.f)) lines(wavelength, l.d[g,])
    lines(wavelength, median.t, col = "red")
    abline(v = c(350, 700), col = "grey", lty = "dashed")
    abline(h = 0 )
    dev.off()
  }
  return(median.overview)
}

# ==== 2 plot light per location ======
# written to check the irradiance and transmission data 
fig.light.measures <- function(light.data, wavelength, path.subdir, fig.name){
  #light.data <- tBan
  #wavelength <- wavelength.range
  # name.fig
  dt <- light.data[, c((ncol(light.data) - length(wavelength) + 1):ncol(light.data))]
  ylimt <- range(dt)
  pdf(paste(path.fig.subdir[path.subdir], fig.name,"_", light.data$site.name[1], "_light.type_", light.type.side,"_w_", window.w,".irr.depth_",irr.depth,".pdf",sep=""), height = 5, width = 10)
  plot(NA, ylim = ylimt, xli = range(wavelength), las = 1)
  for(tx in 1:nrow(dt)){
    lines(wavelength, dt[tx,])
    text(wavelength[1], dt[tx, 1],light.data$sample.location[tx])
  } 
  dev.off()
  return(apply(dt, 2, median))
}


# ==== 3 sensitivity marine fresh ======
# uses the output from bootstrap funnction to make the sensitivity graphs

pos.legend <- c(0.65, 0.6)

xlim.2 <- range(wave.length.range)
ylim.2 <- c(0, 0.7)

# marine freshwater sensitivities
fig.sens.mar.fr <- function(sens.data){
  
  # sens.data <- summ.data.t
  
  plot(NA, xlim = xlim.2, ylim = ylim.2, las = 2, xaxt="n",ylab = "Sensitivity", xlab = "")
  # marine
  mar.lme <- c(sens.data[1,] )
  lines(wave.length.range, mar.lme , col = colour[1])
  polygon( c(wave.length.range, rev(wave.length.range)), c( (mar.lme + sens.data[3,]) , rev(mar.lme - sens.data[3,])), col = paste(colour.h[1], alpha, sep=""), border = FALSE)
  
  #fresh
  fresh.lme <- c(sens.data[1,] + sens.data[2,])
  lines(wave.length.range, fresh.lme, col = colour[2])
  polygon( c(wave.length.range, rev(wave.length.range)), c(fresh.lme + sens.data[4,], rev(fresh.lme - sens.data[4,])), col = paste(colour.h[2], alpha, sep=""), border = FALSE)
  
  lines(wave.length.range, mar.lme , col = colour[1])
  lines(wave.length.range, fresh.lme, col = colour[2])
  
  
  axis(1, seq(wave.length.range[1], wave.length.range[length(wave.length.range)], 50))
  mtext("Wavelength (nm)",1, 3 )
  
  names.1 <- c("marine", "freshwater")
  
  # legend
  for(i in 1:2){
    points(wave.length.range[1], pos.legend[i], pch = 15, col = colour[i])
    text(wave.length.range[1], pos.legend[i], names.1[i],col = colour[i], pos = 4 )
  }
}

# ==== 4 sensitivity benthic limnetic ======
# uses the output from bootstrap funnction to make the sensitivity graphs
fig.sens.lim.ben <- function(sens.data){
  
  # sens.data <- sp.1[[3]]
  
  plot(NA, xlim = xlim.2, ylim = ylim.2, las = 2, xaxt="n",ylab = "Sensitivity", xlab = "")
  # lim
  lim.lme <- c(sens.data[1,] )
  lines(wave.length.range, lim.lme , col = colour[3])
  polygon( c(wave.length.range, rev(wave.length.range)), c( (lim.lme + sens.data[3,]) , rev(lim.lme - sens.data[3,])), col = paste(colour.h[3], alpha, sep=""), border = FALSE)
  
  # benthic
  ben.lme <- c(sens.data[1,] + sens.data[2,])
  lines(wave.length.range, ben.lme, col = colour[4])
  polygon( c(wave.length.range, rev(wave.length.range)), c(ben.lme + sens.data[4,], rev(ben.lme - sens.data[4,])), col = paste(colour.h[4], alpha, sep=""), border = FALSE)
  
  lines(wave.length.range, lim.lme , col = colour[3])
  lines(wave.length.range, ben.lme, col = colour[4])
  
  
  axis(1, seq(wave.length.range[1], wave.length.range[length(wave.length.range)], 50))
  mtext("Wavelength (nm)",1, 3 )
  
  names.1 <- c("limnetic", "benthic")
  
  # legend
  for(i in 1:2){
    points(wave.length.range[1], pos.legend[i], pch = 15, col = colour[i+2])
    text(wave.length.range[1], pos.legend[i], names.1[i],col = colour[i+2], pos = 4 )
  }
}


# ==== 5 fig.slope.mar.fresh.in.fresh ======
# test where marine and freshwater sensitivities are correalted against freshwater median light environment

offset.t <- 0.2
jit.1 <- 0.75
mean.l <- 0.15 # 1/2 the length of the line for the mean
lwd.1 <-  3
fig.slope.mar.fresh.in.fresh.K <- function(data.slope){
  #data.slope <- mar.fresh.O.1
  
  pop.n <- unique(data.slope$pop.light)
  xlim <- c(0.5, length(pop.n) + 0.5)
  ylim.K <- range(data.slope$slope.K)
  plot(NA, xlim = xlim, ylim = ylim.K, las = 1, ylab = "slope", xlab = "", xaxt = "n", main = "K")
  for(ix in 1:length(pop.n)){
    #ix <- 1
    pt <- data.slope[data.slope$pop.light == pop.n[ix],]
    points( jitter(rep(ix + offset.t, length(pt$slope.K[pt$pop.S == pop.n[ix]])), jit.1), pt$slope.K[pt$pop.S == pop.n[ix]], col = colour[1], pch = 19)
    lines(c((ix + offset.t - mean.l), (ix + offset.t + mean.l)), rep(mean(pt$slope.K[pt$pop.S == pop.n[ix]]), 2), col = colour[1], lwd = lwd.1)
    points( jitter(rep(ix - offset.t, length(pt$slope.K[pt$pop.S != pop.n[ix]])), jit.1), pt$slope.K[pt$pop.S != pop.n[ix]], pch = 19, col = colour[2])
    lines(c((ix - offset.t- mean.l), (ix - offset.t + mean.l)), rep(mean(pt$slope.K[pt$pop.S != pop.n[ix]]), 2), col = colour[2], lwd = lwd.1)
  }
  axis(1,c(1:length(pop.n)), pop.n)
  abline(h = 0, lty = "dashed")
}

fig.slope.mar.fresh.in.fresh.I <- function(data.slope){
  
  pop.n <- unique(data.slope$pop.light)
  xlim <- c(0.5, length(pop.n) + 0.5)
  
  ylim.I <- range(data.slope$slope.I)
  plot(NA, xlim = xlim, ylim = ylim.I, las = 1, ylab = "slope", xlab = "", xaxt = "n", main = "I")
  for(ix in 1:length(pop.n)){
    #ix <- 1
    pt <- data.slope[data.slope$pop.light == pop.n[ix],]
    points( jitter(rep(ix + offset.t, length(pt$slope.I[pt$pop.S == pop.n[ix]])), jit.1), pt$slope.I[pt$pop.S == pop.n[ix]], pch = 19, col = colour[1])
    lines(c((ix + offset.t - mean.l), (ix + offset.t + mean.l)), rep(mean(pt$slope.I[pt$pop.S == pop.n[ix]]), 2), col = colour[1], lwd = lwd.1)
    points( jitter(rep(ix - offset.t, length(pt$slope.I[pt$pop.S != pop.n[ix]])), jit.1), pt$slope.I[pt$pop.S != pop.n[ix]], pch = 19, col = colour[2] )
    lines(c((ix - offset.t- mean.l), (ix - offset.t + mean.l)), rep(mean(pt$slope.I[pt$pop.S != pop.n[ix]]), 2), col = colour[2], lwd = lwd.1)
  }
  axis(1,c(1:length(pop.n)), pop.n)
  abline(h = 0, lty = "dashed")
}



# ===== 6 slopes of dS vs dI/dK =====
# figure to plut the slopes for each populations
jit.2 <- 0.75

fig.slopes.dS.dI.dK <- function(slope.data, path, name){
  #slope.data <- m.f.O.1_0 
  #path <- path.figures.subdir[6]
  #name <- "test"
  slope.data <- slope.data[, 1:3]
  
  names.p <- unique(slope.data$site.name)
  if(length(names.p) <= 2){
    width <- 1.5 * length(names.p)
  } else{
    width <- length(names.p)
  }
  height <-  4
  #  K plot
  pdf( paste(path, "slope.K.",name,".pdf", sep = ""), width = width, height = height)
  ylim.K <- range(slope.data$K)
  if(ylim.K[1] > 0) ylim.K[1] <- 0
  if(ylim.K[2] < 0) ylim.K[2] <- 0
  #diff(xlim.K)
  xlim.K <- c(0.5, length(unique(slope.data$site.name)) + 0.5)
  plot(NA, xlim = xlim.K, ylim = ylim.K, las = 1, ylab = "", xlab = "", xaxt = "n")
  text(xlim.K[1] - 0.3 * diff(xlim.K), ylim.K[1] + 0.5 * diff(ylim.K), expression(beta), xpd = TRUE, srt = 0)
  
  for(tq in 1:length(names.p)){
    #tq <- 1
    d <- slope.data[slope.data$site.name == names.p[tq],]
    points( jitter(rep(tq, nrow(d)), jit.2), d$K, pch = 19, col = colour[1])
  }
  abline(h = 0, lty = "dashed")
  # to get species pair names right
  if(length(names.p) <= 2){
    axis(1, c(1:length(names.p)), substr(names.p, 1, 6))
  }else(axis(1, c(1:length(names.p)), names.p))
  dev.off()
  
  #  I plot
  pdf(paste(path, "slope.I.", name,".pdf", sep = ""), width = width, height = height)
  ylim.I <- range(slope.data$I)
  if(ylim.I[1] > 0) ylim.I[1] <- 0
  if(ylim.I[2] < 0) ylim.I[2] <- 0
  
  xlim.I <- c(0.5, length(unique(slope.data$site.name)) + 0.5)
  plot(NA, xlim = xlim.I, ylim = ylim.I, las = 1, ylab = "", xlab = "", xaxt = "n")
  text(xlim.I[1] - 0.3 * diff(xlim.I), ylim.I[1] + 0.5 * diff(ylim.I), expression(beta), xpd = TRUE, srt = 0)
  for(tq in 1:length(names.p)){
    #tq <- 1
    d <- slope.data[slope.data$site.name == names.p[tq],]
    points( jitter(rep(tq, nrow(d)), jit.2), d$I, pch = 19, col = colour[1])
  }
  abline(h = 0, lty = "dashed")
  # to get species pair names right
  if(length(names.p) <= 2){
    axis(1, c(1:length(names.p)), substr(names.p, 1, 6))
  }else(axis(1, c(1:length(names.p)), names.p))
  dev.off()
}


fig.slopes.dS.dI.dK.sp.pair <- function(slope.data, path, name){
  # slope.data <- bl.0
  # path <- path.figures.subdir[6]
  # name <- "Species pair"
  slope.data <- slope.data[, 1:3]
  
  names.p <- unique(slope.data$site.name)
  width <- 1.5 * length(names.p)
  height <- 4
  #  K plot
  pdf( paste(path, "slope.K.",name,".pdf", sep = ""), width = width, height = height)
  ylim.K <- range(slope.data$K)
  xlim.K <- c(0.5, length(unique(slope.data$site.name)) + 0.5)
  plot(NA, xlim = xlim.K, ylim = ylim.K, las = 1, ylab = "", xlab = "", xaxt = "n")
  text(xlim.K[1] - 0.4 * diff(xlim.K), ylim.K[1] + 0.5 * diff(ylim.K), expression(beta), xpd = TRUE, srt = 0)
  
  for(tq in 1:length(names.p)){
    #tq <- 1
    d <- slope.data[slope.data$site.name == names.p[tq],]
    points( jitter(rep(tq, nrow(d)), jit.2), d$K, pch = 19, col = colour[1])
  }
  abline(h = 0, lty = "dashed")
  axis(1, c(1:length(names.p)), substr(names.p, 1, 6))
  dev.off()
  
  #  I plot
  pdf(paste(path, "slope.I.", name,".pdf", sep = ""), width = width, height = height)
  ylim.I <- range(slope.data$I)
  xlim.I <- c(0.5, length(unique(slope.data$site.name)) + 0.5)
  plot(NA, xlim = xlim.I, ylim = ylim.I, las = 1, ylab = "", xlab = "", xaxt = "n")
  text(xlim.I[1] - 0.4 * diff(xlim.I), ylim.I[1] + 0.5 * diff(ylim.I), expression(beta), xpd = TRUE, srt = 0)
  
  for(tq in 1:length(names.p)){
    #tq <- 1
    d <- slope.data[slope.data$site.name == names.p[tq],]
    points( jitter(rep(tq, nrow(d)), jit.2), d$I, pch = 19, col = colour[1])
  }
  abline(h = 0, lty = "dashed")
  axis(1, c(1:length(names.p)), substr(names.p, 1, 6))
  dev.off()
}

