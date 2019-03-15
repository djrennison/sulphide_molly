#
#
# R script to analyse the colour files
#
# each indiviudal has multipl measurements (3x) for
# belly, eye, head and tail
# each indidividual comes from a population, has a sex and a number within pop.




# === General variables ====
names.col.axis <- c("wavelength", "intensity")
#
window.w.col <- 5		# the size to the left/right of vocal point over which mean is taken

col.pop.names <- list.files(path.data.col)

col.meta.data.t <- read.csv(paste(path.data.meta.d,"colour_indiviudals_info.csv", sep = "" ), stringsAsFactors = FALSE)
col.meta.data <- col.meta.data.t[,c(1:4)]

# list to store all the data sets in
col.all <- vector("list", length(col.pop.names))
names(col.all) <- col.pop.names

# names for the columns of the overview irradiance file
col.file.names <- c("population", "file.name", "individual","sex", "body.part", "rep")  
  
# go through all the folders (i) and within the folder all files (j)
for(i in 1:length(col.pop.names)){
  # i <- 1

  # file names of file in folder i 
  file.names <- list.files(paste(path.data.col, col.pop.names[i], sep = ""))
  
  # make an overview file to store new irradiance measures in
  col.temp <- as.data.frame(matrix(NA, ncol = length(col.file.names) + length(wavelength.range), nrow = length(file.names)))
  colnames(col.temp) <- c(col.file.names, wavelength.range)
  # head(irrad.temp)
  

  # go through all the files in the population folder
  for(j in 1:length(file.names)){
    print(paste("folder = ", i, " file = ",j,sep=""))
    # j <- 1
    col.t <- read.table(paste(path.data.col, col.pop.names[i], "/", file.names[j], sep = ""), fill = TRUE, stringsAsFactors = FALSE)
    
    # update the meta data in the irrad.temp file
    # population
    col.temp$population[j] <- gsub( "-.*", "", file.names[j])
    # file.name
    col.temp$file.name[j] <- file.names[j] 
    #indi.name
    col.temp$individual[j] <- gsub( "_.*", "", file.names[j])
    
    # sex
    col.temp$sex[j] <- col.meta.data$sex[ col.meta.data$individual == gsub( "_.*", "", file.names[j])]
    
    name.f.t <- gsub(".*_s*|.txt", "", file.names[j])
    # body.part rep
    col.temp$rep[j] <- as.integer(substr(name.f.t, nchar(name.f.t), nchar(name.f.t)))
    # body.part
    col.temp$body.part[j] <- substr(name.f.t, 1, nchar(name.f.t) - 1)
   
    
    # clean up the colour file
    if(is.numeric(col.t[1,1]) == FALSE) {
      # data file is messy with a header section of 17 rows and only column 1 + 2 	with data neede
      col.t.1 <- col.t[-c(1:17, length(col.t[,1])), c(1,2)]
      col.t.1[,1] <- as.numeric(col.t.1[,1])
      col.t.1[,2] <- as.numeric(col.t.1[,2])
    } else {
      col.t.1 <- col.t
      print("no text header rows")
    }
    colnames(col.t.1) <- names.col.axis
    # exclude > min(wavelength.range)  and > max(wavelength.range)
    col.t.1a <- col.t.1[ col.t.1$wavelength >= min(wavelength.range) & col.t.1$wavelength <= max(wavelength.range),]
    plot(col.t.1a[,2] ~ col.t.1a[,1], xlab = names.col.axis[1], ylab = names.col.axis[2], las = 1, type = "l")

    # get the rolling mean of the wavelength and the intensity
    col.t.roll.mean <- data.frame(cbind(rollmean(col.t.1a$wavelength, window.w.col), rollmean(col.t.1a$intensity, window.w.col)))
    head(col.t.roll.mean)
    colnames(col.t.roll.mean) <- names.col.axis
    # fit a spline to in the next step get intensity per nm
    fitted.spline <- splinefun(col.t.roll.mean$wavelength, col.t.roll.mean$intensity)
    new.intensity.col <- fitted.spline(wavelength.range)

    col.temp[ j, c((length(col.file.names) + 1): ncol(col.temp))] <- new.intensity.col
    # irrad.temp[,1:12]
  }
  col.all[[i]] <- col.temp
}

col.all.smooth <- do.call("rbind", col.all)
str(col.all.smooth)
col.all.smooth[1:30,1:10]

col.all.smooth$population

# write data files as temp for now
write.csv(col.all.smooth, paste( path.output.colour, "colour.all.smooth.w=", window.w.col, ".csv", sep = ""), row.names = FALSE)


# normalised version
# sum each row:
sum.per.row.c <- apply(col.all.smooth[, -c(1:length(col.file.names))], 1, sum)
col.norm <- col.all.smooth[, -c(1:length(col.file.names))] / sum.per.row.c
apply(col.norm, 1, sum) # correct, all sum to 1 now
col.all.smooth.norm <- cbind.data.frame(col.all.smooth[, c(1:length(col.file.names))], col.norm)
# write to file
write.csv(col.all.smooth.norm, paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".csv", sep = ""), row.names = FALSE)


# === 2 calculate summary stats per triplets ==================

# delete rows with missing values for sex and rep

NA.rep <- which(is.na(col.all.smooth$rep)==TRUE) 
NA.sex <- which(is.na(col.all.smooth$sex)==TRUE) 

col.all.smooth.norm.1 <- col.all.smooth.norm[-c(NA.rep, NA.sex),]
unique(col.all.smooth.norm.1$rep)

# delete any row with larger rep than 3
col.all.smooth.norm.2 <- col.all.smooth.norm.1[-which(col.all.smooth.norm.1$rep > 3) ,]
unique(col.all.smooth.norm.2$rep)

# get the mean of three replicates
# number of individuals
n.indiv <- unique(col.all.smooth.norm.2$individual)

col.all.smooth.norm.2[,-c(1,2,4,6)] %>% group_by(individual, body.part) %>% summarise_each(funs(mean)) -> col.summ.1

col.summ.2 <- as.data.frame(col.summ.1)
col.all.smooth.norm.summary <- merge(col.meta.data, col.summ.2, by = "individual")

write.csv(col.all.smooth.norm.summary, paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".summary.csv", sep = ""), row.names = FALSE)


# === 3 test for differences using permutation test =========
# data  first two columns fixed.effect and random.effect
# rest wavelength.range # columns with reflectance for

test.data <- col.all.smooth.norm.summary[col.all.smooth.norm.summary$body.part == "tail" | col.all.smooth.norm.summary$body == "head",]
test.data$body.part
str(test.data)
test.d <- test.data[, -c(1,2,4)]

f.permutation.pops <- function(data, random.effect, fixed.effect, n.permut){

    #data <- test.d
    #random.effect <- test.d$population
    #fixed.effect <- test.d$body.part
    #n.permut <- 3
    
    n.random.effect.t <- unique(random.effect)
    
    # list to store all the output
    permut.all <- vector("list", 4)
    names(permut.all) <- c("real.t", "real.t.sum", "null.t", "null.t.sum")
    
    # --- A --- the analysis with the 'real' fixed effect order
    temp.t <- rep(NA, ncol(data) - 2)
    for(j in 1:(ncol(data) - 2)){
      # j <- 50
      print(j)
      m.t <- lmer(data[,j + 2] ~ fixed.effect + (1|random.effect), REML = FALSE)
      #t value of fixed effect
      temp.t[j] <- as.numeric(summary(m.t)$coefficients[2,4])
    }
    # write to overview file
    permut.all$real.t <- temp.t
    permut.all$real.t.sum <- sum(temp.t)
    
    
    # --- B --- the permutation test
    # resample fixed effect within each random effect
    # do this n.permut times so all the 'reshuffled' fixed effects are ready
    permut.fixed.effect <- matrix(NA, ncol = n.permut, nrow = nrow(data))
    
    for(i in 1:n.permut){
      fixed.permut <- fixed.effect
      for(j in 1:length(n.random.effect.t)){
        #i <- 1
        locat.t <- which(random.n == n.random.effect.t[j]) 
        fixed.permut[ locat.t ] <- sample(fixed.effect[locat.t], length(locat.t), replace = FALSE )       
      }
      permut.fixed.effect[,i] <- fixed.permut
    }
    rownames(permut.fixed.effect) <- random.effect
    
    
    # for each new permutation, do a mixed-efefcts model to test for differnces
    # of fixed effect. Store the t value for the fixed effect for each wavelength range
    # the sum of t values is our test statistic (making null distribution here)
    sum.t.null <- rep(NA, n.permut)
    t.null <- matrix(NA, ncol = n.permut, nrow = (ncol(data) - 2))
    
    for(i in 1:n.permut){
      #i <- 1
      temp.t <- rep(NA, ncol(data) - 2)
      for(j in 1:(ncol(data) - 2)){
        # j <- 50
        print(paste("permut # ", i, " wavelength # ",j, sep = ""))
        m.t <- lmer(data[,j + 2] ~ permut.fixed.effect[,i] + (1|random.effect), REML = FALSE)
        #t value of fixed effect
        temp.t[j] <- as.numeric(summary(m.t)$coefficients[2,4])
      }
      t.null[,i] <- temp.t
      sum.t.null[i] <- sum(temp.t)
    }
    
    # prepare overview to be returned
    permut.all$null.t <- t.null
    permut.all$null.t.sum <- sum.t.null
 
    return(permut.all)   
}
