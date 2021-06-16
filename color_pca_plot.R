library(tidyverse)
library(PNWColors)
library(photobiology)
library(colorspace)
library(cowplot)
library(patchwork)
meta_data <- read_csv("data/meta_data/site_meta_data.csv")
color <-read_tsv( paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".summary.tsv", sep = ""))


body.parts <- c("Belly","Eye","Head","Tail")
all_pcas <- tibble(population=character(),individual=character(),sex=character(),body.part=character(),
                   PC1=numeric(),PC2=numeric(),PC3=numeric(),PC4=numeric(),PC5=numeric(),PC6=numeric())
anova_all <- tibble(variable=character(),  Df=numeric(),`Sum Sq`=numeric(),
                    `Mean Sq`=numeric(), `F value`=numeric(),   `Pr(>F)`=numeric(),   
                    PctExp=numeric(),    body.part=character())
all_loadings <- tibble(pc1_loading=numeric(),lambda=character(),body_part=character())
scatter_plots <- list()
density_plots <- list()
box_plots <- list()
i <- 0
for (chosen.body.part.capital in body.parts){
  i <- i + 1
  chosen.body.part <- tolower(chosen.body.part.capital)
  color %>%
    filter(lambda < 700) %>%
    filter(body.part ==chosen.body.part) %>%
    dplyr::mutate(nm_bin = floor(as.numeric(lambda)/5)*5) %>%
    group_by(population,individual,sex,body.part,nm_bin) %>%
    dplyr::summarise(summed_reflectance = sum(mean_reflectance)) %>%
    ungroup() %>%
    mutate(ID= paste0(individual,".",body.part)) %>%
    dplyr::select(ID,nm_bin,summed_reflectance) %>%
    pivot_wider(names_from = nm_bin, values_from = summed_reflectance) ->col.all.smooth.norm.summary.forpca
  pca.info <- color %>%
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
  all.pca.tibble <- as_tibble(all.pca$x[,c(1:ncol(all.pca$x))])
  pr_var = data.frame(varExp = all.pca$sdev^2)
  pr_var = pr_var %>%
    mutate(pve = varExp / sum(varExp))
  write_tsv(cbind(pca.info,all.pca.tibble) %>%
              inner_join(meta_data %>% rename(population = ID)),
            paste0("output/colour/PCA_",chosen.body.part.capital,".txt"))
  
   min_sample <- cbind(pca.info,all.pca.tibble) %>%
    inner_join(meta_data %>% rename(population = ID)) %>%
    group_by(Drainage) %>%
    filter(PC1 == min(PC1)) %>% pull(individual)
  max_sample <- cbind(pca.info,all.pca.tibble) %>%
    inner_join(meta_data %>% rename(population = ID)) %>%
    group_by(Drainage) %>%
    filter(PC1 == max(PC1)) %>% pull(individual)
  
  reference_colors <- tibble(Drainage=character(),min_color=character(),
                             max_color=character())
  drainages <- unique(meta_data$Drainage)
  for (chosen_drainage in drainages){
    min_color <- s_e_irrad2rgb(
      color %>% inner_join(meta_data %>% rename(population = ID)) %>% 
        filter(individual %in% min_sample, body.part == chosen.body.part,
               Drainage == chosen_drainage) %>% 
        mutate(lambda = as.numeric(lambda)) %>% pull(lambda),
      color %>% inner_join(meta_data %>% rename(population = ID)) %>% 
        filter(individual %in% min_sample, body.part == chosen.body.part,
               Drainage == chosen_drainage) %>% pull(mean_reflectance),
    )
    min_color <- lighten(min_color,0.3)
    max_color <- s_e_irrad2rgb(
      color %>% inner_join(meta_data %>% rename(population = ID)) %>% 
        filter(individual %in% max_sample, body.part == chosen.body.part,
               Drainage == chosen_drainage) %>% 
        mutate(lambda = as.numeric(lambda)) %>% pull(lambda),
      color %>% inner_join(meta_data %>% rename(population = ID)) %>% 
        filter(individual %in% max_sample, body.part == chosen.body.part,
               Drainage == chosen_drainage) %>% pull(mean_reflectance),
    )
    max_color <- lighten(max_color,0.3)
    reference_colors <- rbind(reference_colors, tibble(Drainage=chosen_drainage,min_color=min_color,
                                                       max_color=max_color))
  }
  min_x <- min(all.pca.tibble$PC1)*1.1
  max_x <- max(all.pca.tibble$PC1)*1.1
  min_y <- min(all.pca.tibble$PC2)
  max_y <- max(all.pca.tibble$PC2)
  reference_colors$min_x <- min_x
  reference_colors$min_y <- min_y
  reference_colors$max_x <- max_x
  reference_colors$max_y <- max_y
  
  scatter_plot <- cbind(pca.info,all.pca.tibble) %>%
    inner_join(meta_data %>% rename(population = ID)) %>%
    ggplot(.) +
    geom_point(aes(x=PC1,y=PC2,color=as.factor(H2S),shape=sex)) +
    theme_cowplot() + 
    facet_wrap(~Drainage,nrow=1) +
    stat_ellipse(aes(x=PC1,y=PC2,color=as.factor(H2S),shape=sex,group=population)) +
    scale_color_manual(values=color_palette_1,name="Environment",labels=c("Non-sulphur","Sulphur")) +
    ggtitle(paste(chosen.body.part.capital)) +
      geom_rect(data=reference_colors,
                aes(xmin=min_x-2,xmax=min_x,ymin=min_y,ymax=max_y,fill=min_color)) +
      geom_rect(data=reference_colors,
                aes(xmin=max_x,xmax=max_x+2,ymin=min_y,ymax=max_y,fill=max_color)) +
      scale_fill_manual(values=sort(c(reference_colors$min_color, reference_colors$max_color)),
                        guide = 'none') +
    scale_shape_manual(values=c(20,3),name="Sex") +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) +
    xlab(paste0("PC1 ", round(100*pr_var$pve[1],1),"% PVE")) +
    ylab(paste0("PC2 ", round(100*pr_var$pve[2],1),"% PVE"))
  
  scatter_plots[[i]] <- scatter_plot

   box_plot <- cbind(pca.info,all.pca.tibble) %>%
    inner_join(meta_data %>% rename(population = ID)) %>%
     mutate(H2S = case_when(H2S == 0 ~ "B",
                            TRUE ~ "A")) %>%
    ggplot(.,aes(y=PC1,x=as.factor(H2S),color=as.factor(H2S))) +
    geom_jitter(aes(group=population,shape=sex),alpha=1,width=0.1) +
     geom_boxplot(fill=NA,color="black",outlier.shape = NA,width=0.2) + 
     theme_cowplot() +
    facet_wrap(~Drainage,nrow=1) +
    scale_color_manual(values=rev(color_palette_1),name="Environment",labels=c("Sulphur","Non-sulphur"))  +
    theme(axis.text=element_blank(),axis.title.x=element_blank()) +
    ylab("PC1") +
     scale_shape_manual(values=c(20,3),name="Sex") 
     
     
   box_plots[[i]] <- box_plot
  all_pcas <- rbind(all_pcas, cbind(pca.info,all.pca.tibble))
  
  pca_tmp <- cbind(pca.info,all.pca.tibble) %>%
    inner_join(meta_data %>% rename(population = ID))%>%
    mutate(H2S = case_when(H2S == 0 ~ "Non-sulfur",
                           TRUE ~ "Sulfur")) 
  anova_result_h2s <- aov(PC1 ~ Drainage + H2S + sex,  data = pca_tmp)
  
    summary(anova_result_h2s)
  TukeyHSD(anova_result_h2s)
  af <- anova(anova_result_h2s)
  afss <- af$"Sum Sq"
  as_tibble(cbind(af,PctExp=afss/sum(afss)*100),rownames = "variable") %>%
    mutate(body.part = chosen.body.part) -> anova_results
  anova_all <- rbind(anova_all,anova_results)
  
}

#write_tsv(anova_all,"output/colour/body_colour_anova.v1.txt")
pdf("figures/color_pcas_v3.pdf",height=10,width=14,useDingbats = F)
plots_part1  <- ((scatter_plots[[1]] / box_plots[[1]]) | (scatter_plots[[2]] / box_plots[[2]]) ) 
plots_part2  <- ((scatter_plots[[3]] / box_plots[[3]]) | (scatter_plots[[4]] / box_plots[[4]]) ) 
plots_part1 / plots_part2 +plot_annotation(tag_levels = 'A')
dev.off()

pdf("figures/color_pcas_boxplot_v1.pdf",height=3,width=10,useDingbats = F)
all_pcas %>%
  inner_join(meta_data %>% rename(population = ID)) %>%
  mutate(body.part = str_to_title(body.part)) %>%
  ggplot(.,aes(x=Drainage,y=PC1,,color=as.factor(H2S))) + 
  geom_point(aes(color=as.factor(H2S),shape=as.factor(sex)), position = position_jitterdodge(jitter.width=0.1),size=1) +
  geom_boxplot(outlier.shape = NA,alpha=0.6) +
  facet_wrap(~body.part,scales="free_y",nrow=1) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=9),
        text = element_text(size=10),
        axis.text.y=element_text(size=9)) +  
  scale_color_manual(values=color_palette_1,name="Environment",
                     labels=c("Non-sulphur", "Sulphur")) +
  scale_fill_manual(values=color_palette_1,name="Environment",
                    labels=c("Non-sulphur", "Sulphur")) +
  scale_shape_manual(values=c(20,21),name="Sex") +
  ylab("PC1")
dev.off()

all_pcas_data <- all_pcas %>%inner_join(meta_data %>% rename(population = ID)) %>%
  mutate(H2S = case_when(H2S == 1 ~ "Non-sulfur",
                         TRUE ~ "Sulfur"))

anova_result_h2s <- aov(PC1 ~ body.part,  data = all_pcas_data)
summary(anova_result_h2s)
TukeyHSD(anova_result_h2s)
af <- anova(anova_result_h2s)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

anova_result <- aov(PC1 ~ Drainage + H2S + sex,  data = all_pcas_data)

summary(anova_result_h2s)
###Now looking at PCA to see where spectral differences occur


plot_5 <- all_loadings %>%
  filter(body_part != "head") %>%
  ggplot(.,aes(x=as.numeric(lambda),y=pc1_loading,color=body_part)) + 
  geom_line() +
  geom_hline(yintercept=0,linetype="dotted") +
  ylab("Color PC1 loadings") +
  theme_cowplot()



sensitivity <- read_tsv("output/opsin/sensitivity.A1.txt")

plot_6 <- sensitivity %>%
  ggplot(.,aes(x=lambda,y=lambda_sensitivity,group=Identifier,color=as.factor(H2S))) + geom_line(alpha=0.5) +
  facet_wrap(~Drainage,nrow=4) +
  scale_color_manual(values=color_palette_1,name="Environment",labels=c("Non-sulfur","Sulfur")) +
  ylab("Sensitivity") +
  theme_cowplot()

reference_sensitivity <- sensitivity %>%
  inner_join(meta_data %>% dplyr::select(Fieldsite.ID,H2S)) %>%
  filter(H2S == 0) %>%
  group_by(Drainage, lambda) %>%
  dplyr::summarize(reference_median_sensitivity = median(lambda_sensitivity)) %>%
  ungroup()

sensitivity.differences <- sensitivity %>%
  inner_join(meta_data %>% dplyr::select(Fieldsite.ID,H2S)) %>%
  filter(H2S == "1") %>%
  ungroup() %>%
  mutate(lambda=as.numeric(lambda)) %>%
  inner_join(reference_sensitivity) %>%
  mutate(sensitivity_dif = reference_median_sensitivity - lambda_sensitivity ) %>%
  dplyr::select(Identifier,lambda,Fieldsite.ID,Drainage,sensitivity_dif) %>%
  group_by(lambda,Drainage) %>%
  dplyr::summarise(mean_sensitivity_dif = median(sensitivity_dif))

plot_7 <- sensitivity.differences %>%
  ggplot(.,aes(x=lambda,y=mean_sensitivity_dif,group=Drainage,color=Drainage)) + geom_line() +
  geom_hline(yintercept=0,linetype="dotted")  +
  ylab("Sensitivity difference")

pdf("figures/sensitivity_differences_prelim.v1.pdf",height=15,width=10)
plot_5 / plot_6 / plot_7 +plot_layout(heights = c(1, 3,1))
dev.off()
