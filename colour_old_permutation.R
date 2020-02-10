
# === 3 test for differences using permutation test =========
# data  first two columns fixed.effect and random.effect
# rest wavelength.range # columns with reflectance for

test.data <- col.all.smooth.norm.summary[col.all.smooth.norm.summary$body.part == "tail" | col.all.smooth.norm.summary$body == "head",]
test.data$body.part
str(test.data)
test.d <- test.data[, -c(1,2,4)]

f.permutation.pops <- function(data, random.effect, fixed.effect, n.permut){
  
  data <- test.d
  random.effect <- test.d$population
  fixed.effect <- test.d$body.part
  n.permut <- 3
  
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
  
  
  # for each new permutation, do a mixed-effects model to test for differences
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