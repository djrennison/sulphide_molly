#
#
# R script to analyse the colour files
#
# each indiviudal has multipl measurements (3x) for
# belly, eye, head and tail
# each indidividual comes from a population, has a sex and a number within pop.




# === General variables ====
library(PNWColors)
library(ggforce)
library(patchwork)
library(photobiology)
library(colorspace)
library(tidyverse)

color_palette_1 <- pnw_palette(name="Bay",n=2,type="discrete")
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

#Remove extra PS0-20 belly measurement and unsexed samples
body.parts <- c("head","eye","tail","belly")
col.all.smooth <- col.all.smooth %>%
  filter(body.part %in% body.parts) %>%
  filter(!is.na(sex)) 



# write data files as temp for now
#write.csv(col.all.smooth, paste( path.output.colour, "colour.all.smooth.w=", window.w.col, ".csv", sep = ""), row.names = FALSE)

col.all.smooth <- read_csv( paste( path.output.colour, "colour.all.smooth.w=", window.w.col, ".csv", sep = ""))
#Check to see how many measures are negative
col.all.smooth %>%
  pivot_longer(c(-individual,-population,-sex,-body.part,-file.name,-rep),
               names_to = "lambda",values_to="reflectance") %>%
  group_by(file.name) %>%
  dplyr::summarize(count= sum(reflectance < 0)) %>%
  ggplot(.,aes(count)) + geom_histogram()

# Remove measures with more than 5 negative lambda values. Then replace negative reflectance values with 0.

col.all.smooth.filt <- col.all.smooth %>%
  pivot_longer(c(-individual,-population,-sex,-body.part,-file.name,-rep),
               names_to = "lambda",values_to="reflectance") %>%
  group_by(file.name) %>%
  dplyr::mutate(count= sum(reflectance < 0)) %>%
  filter(count <= 5) %>%
  mutate(reflectance = case_when(reflectance < 0 ~ 0,
                                 TRUE ~ reflectance)) %>%
  dplyr::select(-count)

#Write non-normalized to file
#write.csv(col.all.smooth.filt, paste( path.output.colour, "colour.all.smooth.filt.w=", window.w.col, ".csv", sep = ""), row.names = FALSE)
col.all.smooth.filt <- read.tsv(paste( path.output.colour, "colour.all.smooth.filt.w=", window.w.col, ".csv", sep = ""))
  



col.all.smooth.norm <- col.all.smooth.filt %>%
  group_by(population,file.name,individual,sex,body.part,rep) %>%
  mutate(sum_reflectance=sum(reflectance),norm_reflectance=reflectance/sum_reflectance) %>%
  dplyr::select(-reflectance,-sum_reflectance)
# normalised version

# write to file
write_tsv(col.all.smooth.norm, paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".tsv", sep = ""))

col.all.smooth.norm <- read_tsv(paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".tsv", sep = ""))
# === 2 calculate summary stats per triplets ==================
col.all.smooth.norm.summary <- col.all.smooth.norm %>%
  group_by(population,individual,sex,body.part,lambda) %>%
  dplyr::summarize(mean_reflectance = mean(norm_reflectance))
  

write_tsv(col.all.smooth.norm.summary, paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".summary.tsv", sep = ""))

population_information <- read_csv("data/meta_data/site_meta_data.csv") %>%
  rename(population = ID)
col.all.smooth.norm.summary %>%
  inner_join(population_information) %>%
  filter(sex == "Female") %>% 
  filter(lambda < 700) %>%
  group_by(population,sex,body.part,Drainage,H2S,Species,lambda) %>%
  do(data.frame(t(quantile(.$mean_reflectance, probs = c(0.025,0.50, 0.975))))) %>%
  rename(bottom = X2.5., mid = X50., top =X97.5.) %>%    
  filter(body.part != 0) %>%
  ggplot(.,aes()) + 
  geom_ribbon(aes(x=as.numeric(lambda),ymin=(bottom),ymax=(top),group=population,fill=as.factor(H2S)),alpha=0.5) +
  geom_line(aes(x=as.numeric(lambda),y=(mid),group=population),alpha=1,color="black",size=0.7) +
  facet_grid(Drainage~body.part) + 
  scale_fill_manual(values=color_palette_1) +
  ggtitle("Female Color")

col.all.smooth.norm.summary %>%
  inner_join(population_information) %>%
  filter(sex == "Male") %>% 
  filter(lambda < 700) %>%
  group_by(population,sex,body.part,Drainage,H2S,Species,lambda) %>%
  do(data.frame(t(quantile(.$mean_reflectance, probs = c(0.025,0.50, 0.975))))) %>%
  rename(bottom = X2.5., mid = X50., top =X97.5.) %>%    
  filter(body.part != 0) %>%
  ggplot(.,aes()) + 
  geom_ribbon(aes(x=as.numeric(lambda),ymin=(bottom),ymax=(top),group=population,fill=as.factor(H2S)),alpha=0.5) +
  geom_line(aes(x=as.numeric(lambda),y=(mid),group=population),alpha=1,color="black",size=0.7) +
  facet_grid(Drainage~body.part) + 
  scale_fill_manual(values=color_palette_1) +
  ggtitle("Male Color") 


col.all.smooth.norm.summary %>%
  inner_join(population_information) %>%
  filter(H2S == 0) %>% 
  filter(lambda < 700) %>%
  group_by(population,sex,body.part,Drainage,Species,lambda) %>%
  do(data.frame(t(quantile(.$mean_reflectance, probs = c(0.025,0.50, 0.975))))) %>%
  rename(bottom = X2.5., mid = X50., top =X97.5.) %>%    
  filter(body.part != 0) %>% 
  ggplot(.) + 
  geom_ribbon(aes(x=as.numeric(lambda),ymin=(bottom),ymax=(top),group=sex,fill=sex),alpha=0.5) +
  geom_line(aes(x=as.numeric(lambda),y=(mid),group=sex),alpha=1,color="black",size=0.7) +
  facet_grid(Drainage~body.part) + 
  scale_fill_manual(values=color_palette_1) +
  ggtitle("Non-sulphur Color") 
  

col.all.smooth.norm.summary %>%
  inner_join(population_information) %>%
  filter(H2S == 1) %>% 
  filter(lambda < 700) %>%
  group_by(population,sex,body.part,Drainage,Species,lambda) %>%
  do(data.frame(t(quantile(.$mean_reflectance, probs = c(0.025,0.50, 0.975))))) %>%
  rename(bottom = X2.5., mid = X50., top =X97.5.) %>%    
  filter(body.part != 0) %>% 
  mutate(popsex = paste0(population,sex)) %>%
  ggplot(.) + 
  geom_ribbon(aes(x=as.numeric(lambda),ymin=(bottom),ymax=(top),group=popsex,fill=sex),alpha=0.5) +
  geom_line(aes(x=as.numeric(lambda),y=(mid),group=popsex),alpha=1,color="black",size=0.7) +
  facet_grid(Drainage~body.part) + 
  scale_fill_manual(values=color_palette_1) +
  ggtitle("Sulphur Color") 
  

body.parts <- c("belly","eye","head","tail")
all_pcas <- tibble(population=character(),individual=character(),sex=character(),body.part=character(),
                   PC1=numeric(),PC2=numeric(),PC3=numeric(),PC4=numeric(),PC5=numeric(),PC6=numeric())
all_loadings <- tibble(pc1_loading=numeric(),lambda=character(),body_part=character())
scatter_plots <- list()
density_plots <- list()
i <- 0
for (chosen.body.part in body.parts){
  i <- i + 1
  col.all.smooth.norm.summary %>%
    filter(lambda < 700) %>%
    filter(body.part ==chosen.body.part) %>%
    dplyr::mutate(nm_bin = floor(as.numeric(lambda)/5)*5) %>%
    group_by(population,individual,sex,body.part,nm_bin) %>%
    dplyr::summarise(summed_reflectance = sum(mean_reflectance)) %>%
    ungroup() %>%
    mutate(ID= paste0(individual,".",body.part)) %>%
    dplyr::select(ID,nm_bin,summed_reflectance) %>%
    pivot_wider(names_from = nm_bin, values_from = summed_reflectance) ->col.all.smooth.norm.summary.forpca
  pca.info <- col.all.smooth.norm.summary %>%
    filter(lambda < 700) %>%
    filter(body.part ==chosen.body.part) %>%
    dplyr::select(population,individual,sex,body.part) %>%
    ungroup() %>%
    unique() 
  all.pca <- prcomp(col.all.smooth.norm.summary.forpca[,c(2:ncol(col.all.smooth.norm.summary.forpca))], 
                    center = TRUE,scale. = TRUE)
  pca.loadings <- as.tibble(all.pca$rotation)
  pca.loadings$lambda <- row.names(all.pca$rotation)
  pca.loadings <- pca.loadings %>% dplyr::select(PC1, lambda) %>% 
    mutate(body_part = chosen.body.part) %>%
    rename(pc1_loading = PC1)
  all_loadings <- rbind(all_loadings,pca.loadings)
  all.pca.tibble <- as.tibble(all.pca$x[,c(1:6)])
  
  min_sample <- cbind(pca.info,all.pca.tibble) %>%
    filter(PC1 == min(PC1)) %>% pull(individual)
  max_sample <- cbind(pca.info,all.pca.tibble) %>%
    filter(PC1 == max(PC1)) %>% pull(individual)
  
  
  min_color <- s_e_irrad2rgb(
    col.all.smooth.norm.summary %>% filter(individual == min_sample, body.part == chosen.body.part) %>% mutate(lambda = as.numeric(lambda)) %>% pull(lambda),
    col.all.smooth.norm.summary %>% filter(individual == min_sample, body.part == chosen.body.part) %>% pull(mean_reflectance),
  )
  min_color <- lighten(min_color,0.2)
  max_color <- s_e_irrad2rgb(
    col.all.smooth.norm.summary %>% filter(individual == max_sample, body.part == chosen.body.part) %>% mutate(lambda = as.numeric(lambda)) %>% pull(lambda),
    col.all.smooth.norm.summary %>% filter(individual == max_sample, body.part == chosen.body.part) %>% pull(mean_reflectance),
  )
  max_color <- lighten(max_color,0.2)
  
  min_x <- min(all.pca.tibble$PC1)
  max_x <- max(all.pca.tibble$PC1)
  min_y <- min(all.pca.tibble$PC2)
  max_y <- max(all.pca.tibble$PC2)
  
  
  scatter_plot <- cbind(pca.info,all.pca.tibble) %>%
    inner_join(population_information) %>%
    ggplot(.,aes(x=PC1,y=PC2,color=as.factor(H2S),shape=sex)) +
    geom_point() +
    facet_wrap(~Drainage,nrow=1) +
    stat_ellipse(aes(group=population)) +
    scale_color_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) +
    ggtitle(paste(chosen.body.part)) +
    annotate("point", x=min_x,y=max_y,color=min_color,shape=15,size=5) +
    annotate("point", x=max_x,y=max_y,color=max_color,shape=15,size=5) 
  scatter_plots[[i]] <- scatter_plot
  density_plot  <- cbind(pca.info,all.pca.tibble) %>%
    inner_join(population_information) %>%
    ggplot(.,aes(x=PC1,,fill=as.factor(H2S))) +
    geom_density(aes(group=population),alpha=0.5) +
    facet_wrap(~Drainage,nrow=1) +
    scale_fill_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) 
  density_plots[[i]] <- density_plot
  
  all_pcas <- rbind(all_pcas, cbind(pca.info,all.pca.tibble))
}


plots_part1  <- ((scatter_plots[[1]] / density_plots[[1]]) | (scatter_plots[[2]] / density_plots[[2]]) ) 
plots_part2  <- ((scatter_plots[[3]] / density_plots[[3]]) | (scatter_plots[[4]] / density_plots[[4]]) ) 
plots_part1 / plots_part2


###PCAs but separated by sex

body.parts <- c("belly","eye","head","tail")
all_pcas <- tibble(population=character(),individual=character(),sex=character(),body.part=character(),
                   PC1=numeric(),PC2=numeric(),PC3=numeric(),PC4=numeric(),PC5=numeric(),PC6=numeric())
all_loadings <- tibble(pc1_loading=numeric(),lambda=character(),body_part=character(),sex=character())
scatter_plots <- list()
density_plots <- list()
i <- 0
for (chosen.body.part in body.parts){
  for (chosen.sex in c("Male","Female")){
    i <- i + 1
    col.all.smooth.norm.summary %>%
      filter(sex == chosen.sex) %>%
      filter(lambda < 700) %>%
      filter(body.part ==chosen.body.part) %>%
      dplyr::mutate(nm_bin = floor(as.numeric(lambda)/5)*5) %>%
      group_by(population,individual,sex,body.part,nm_bin) %>%
      dplyr::summarise(summed_reflectance = sum(mean_reflectance)) %>%
      ungroup() %>%
      mutate(ID= paste0(individual,".",body.part)) %>%
      dplyr::select(ID,nm_bin,summed_reflectance) %>%
      pivot_wider(names_from = nm_bin, values_from = summed_reflectance) ->col.all.smooth.norm.summary.forpca
    pca.info <- col.all.smooth.norm.summary %>%
      filter(lambda < 700) %>%
      filter(sex == chosen.sex) %>%
      filter(body.part ==chosen.body.part) %>%
      dplyr::select(population,individual,sex,body.part) %>%
      ungroup() %>%
      unique() 
    all.pca <- prcomp(col.all.smooth.norm.summary.forpca[,c(2:ncol(col.all.smooth.norm.summary.forpca))], 
                      center = TRUE,scale. = TRUE)
    pca.loadings <- as.tibble(all.pca$rotation)
    pca.loadings$lambda <- row.names(all.pca$rotation)
    pca.loadings <- pca.loadings %>% select(PC1, lambda) %>% 
      mutate(body_part = chosen.body.part,sex=chosen.sex) %>%
      rename(pc1_loading = PC1)
    all_loadings <- rbind(all_loadings,pca.loadings)
    all.pca.tibble <- as.tibble(all.pca$x[,c(1:6)])
    
    min_sample <- cbind(pca.info,all.pca.tibble) %>%
      filter(PC1 == min(PC1)) %>% pull(individual)
    max_sample <- cbind(pca.info,all.pca.tibble) %>%
      filter(PC1 == max(PC1)) %>% pull(individual)
    
    
    min_color <- s_e_irrad2rgb(
      col.all.smooth.norm.summary %>% filter(individual == min_sample, body.part == chosen.body.part) %>% mutate(lambda = as.numeric(lambda)) %>% pull(lambda),
      col.all.smooth.norm.summary %>% filter(individual == min_sample, body.part == chosen.body.part) %>% pull(mean_reflectance),
    )
    min_color <- lighten(min_color,0.2)
    max_color <- s_e_irrad2rgb(
      col.all.smooth.norm.summary %>% filter(individual == max_sample, body.part == chosen.body.part) %>% mutate(lambda = as.numeric(lambda)) %>% pull(lambda),
      col.all.smooth.norm.summary %>% filter(individual == max_sample, body.part == chosen.body.part) %>% pull(mean_reflectance),
    )
    max_color <- lighten(max_color,0.2)
    
    min_x <- min(all.pca.tibble$PC1)
    max_x <- max(all.pca.tibble$PC1)
    min_y <- min(all.pca.tibble$PC2)
    max_y <- max(all.pca.tibble$PC2)
    
    
    scatter_plot <- cbind(pca.info,all.pca.tibble) %>%
      inner_join(population_information) %>%
      ggplot(.,aes(x=PC1,y=PC2,color=as.factor(H2S))) +
      geom_point() +
      facet_wrap(~Drainage,nrow=1) +
      stat_ellipse(aes(group=population)) +
      scale_color_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) +
      ggtitle(paste(chosen.body.part, chosen.sex)) +
      annotate("point", x=min_x,y=max_y,color=min_color,shape=15,size=5) +
      annotate("point", x=max_x,y=max_y,color=max_color,shape=15,size=5) 
    scatter_plots[[i]] <- scatter_plot
    density_plot  <- cbind(pca.info,all.pca.tibble) %>%
      inner_join(population_information) %>%
      ggplot(.,aes(x=PC1,,fill=as.factor(H2S))) +
      geom_density(aes(group=population),alpha=0.5) +
      facet_wrap(~Drainage,nrow=1) +
      scale_fill_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) 
    density_plots[[i]] <- density_plot
    
    all_pcas <- rbind(all_pcas, cbind(pca.info,all.pca.tibble))
  }
}


plots_part1  <- ((scatter_plots[[1]] / density_plots[[1]]) | (scatter_plots[[2]] / density_plots[[2]]) ) 
plots_part2  <- ((scatter_plots[[3]] / density_plots[[3]]) | (scatter_plots[[4]] / density_plots[[4]]) ) 
plots_part3  <- ((scatter_plots[[5]] / density_plots[[5]]) | (scatter_plots[[6]] / density_plots[[6]]) ) 
plots_part4  <- ((scatter_plots[[7]] / density_plots[[7]]) | (scatter_plots[[8]] / density_plots[[8]]) ) 
plots_part1 / plots_part2 / plots_part3 / plots_part4





all_loadings %>%
  ggplot(.,aes(x=as.numeric(lambda),y=abs(pc1_loading),color=body_part)) + geom_line()


all.pca.tibble <- as.tibble(all.pca$x[,c(1:6)])
cbind(pca.info,all.pca.tibble) %>%
  inner_join(population_information) %>%
  ggplot(.,aes(x=PC1,y=PC2,color=population,shape=as.factor(H2S))) +
  geom_point() +
  facet_wrap(~body.part)

error_bars <- all_pcas %>%
  inner_join(population_information) %>%
  group_by(Drainage, population,H2S,body.part) %>%
  dplyr::summarize(mean = mean(PC1),stdev = sd(PC1))
  

all_pcas %>%
  inner_join(population_information) %>%
  ggplot(.)+
  geom_violin(aes(y=PC1,x=population,fill=as.factor(H2S))) +
  geom_linerange(data=error_bars,
                aes(x=population, ymin=mean-(1.96*stdev),ymax=mean+(1.96*stdev))) +
  geom_point(data=error_bars,
             aes(x=population, y=mean)) +
  facet_grid(body.part~Drainage,scales="free") +
  scale_fill_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) 

  
#Comparing color difference with sensitivity difference at long and short wavelengths. 

bimodal_reflectance <- col.all.smooth.norm.summary %>% 
  inner_join(population_information)  %>%
  group_by(Drainage,H2S,body.part,lambda) %>%
  summarize(median_pop_reflectance = median(mean_reflectance)) %>%
  mutate(wavelength_type = case_when(lambda <= 450 ~ "short",
                                     lambda >= 600 ~ "long",
                                     TRUE ~ "middle")) %>%
  group_by(Drainage,H2S, body.part,wavelength_type) %>%
  filter(wavelength_type != "middle") %>%
  summarize(total_reflectance = sum(median_pop_reflectance)) %>%
  select(body.part,wavelength_type,total_reflectance,Drainage) %>%
  pivot_wider(names_from = H2S, values_from = total_reflectance) %>% 
  mutate(reflectance_dif = `1` -`0` )  %>%
  dplyr::select(-`0`,-`1`)

sensitivity <- read_tsv("output/opsin/sensitivity.A1.txt")

bimodal_sensitivity <- sensitivity %>%
  group_by(lambda,Drainage,H2S) %>%
  summarize(median_pop_sensitivity = median(lambda_sensitivity)) %>%
  
  mutate(wavelength_type = case_when(lambda <= 450 ~ "short",
                                   lambda >= 600 ~ "long",
                                   TRUE ~ "middle")) %>%
  group_by(Drainage,H2S,wavelength_type) %>%
  summarize(total_sen = sum(median_pop_sensitivity)) %>%
  filter(wavelength_type != "middle") %>%  
  pivot_wider(names_from = H2S, values_from = total_sen) %>% 
  mutate(sensitivity_dif = `1` -`0` )  %>%
  dplyr::select(-`0`,-`1`)
  
  
bimodal_sensitivity %>%
  inner_join(bimodal_reflectance) %>%
  ggplot(.,aes(x=sensitivity_dif,y=reflectance_dif,color=body.part)) +
  geom_point(aes()) + facet_wrap(~wavelength_type,scales="free") + geom_smooth(method="lm")
  
  
  
