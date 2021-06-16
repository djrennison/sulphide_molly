library(tidyverse)
library(PNWColors)
library(patchwork)
library(cowplot)
library(infer)
library(ungeviz)
library(ggridges)
color_palette_1 <- pnw_palette(name="Bay",n=2,type="discrete")
#This is for calculating comparisons with light absorbance only in the short wavelengths.
meta_data <- read_csv("data/meta_data/site_meta_data.csv")
drainages <- sort(unique(meta_data$Drainage))
body.parts <- c("belly","eye","head", "tail")

#Transmission should actually be thought of as absorbance
absorbance <- read_tsv("output/transmission/filtered.light.type_s_w_5.irr.depth_10.txt") %>%
  mutate(site.name = case_when(site.name == "Exp" ~ "Esp",
                               site.name == "lab" ~ "Lab",
                               site.name == "Vet" ~ "VS",
                               site.name == "VC" ~ "VG",
                               TRUE ~ site.name)) %>%
  rename(absorbance = transmission) %>%
  filter(light == "s") %>%
  inner_join(meta_data %>% rename(site.name = ID)) %>%
  mutate(ID = paste0(sample.location, "-",H2S)) %>%
  mutate(absorbance = 1-absorbance)




#For color, 
window.w.col <- 5
color <-read_tsv( paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".summary.tsv", sep = "")) %>%
  inner_join(meta_data %>% dplyr::select(ID,H2S,Drainage) %>% rename(population=ID))

body_abs_correlations <- tibble(type=character(),
              drainage=character(),
              body_part=character(),
              correlation=numeric())

#Correlation between body color and light absorbance (<400 nm)
for (chosen_body_part in body.parts){
  for (chosen_drainage in drainages){
    sulphur_samples <- color %>% 
      filter(Drainage == chosen_drainage, body.part == chosen_body_part) %>%
      filter(H2S == 1) %>% pull(individual) %>% unique()
    non_sulphur_samples <- color %>% 
      filter(Drainage == chosen_drainage, body.part == chosen_body_part) %>%
      filter(H2S == 0) %>% pull(individual) %>% unique()
    
    
    for (sulphur_sample in sulphur_samples){
      for (non_sulphur_sample in non_sulphur_samples){
        color_dif <- color %>%
          filter(body.part == chosen_body_part) %>%
          filter(individual == sulphur_sample | individual == non_sulphur_sample) %>%
          dplyr::select(H2S,lambda, mean_reflectance) %>%
          pivot_wider(names_from = H2S, values_from = mean_reflectance) %>%
          mutate(dif_color = `1` - `0`) %>%
          filter(lambda > 350) %>%
          dplyr::select(lambda, dif_color)
        absorb_dif <- absorbance %>%
          filter(Drainage == chosen_drainage) %>%
          group_by(H2S,lambda) %>%
          dplyr::summarise(median_absorbance = median(absorbance)) %>%
          pivot_wider(names_from = H2S, values_from = median_absorbance) %>%
          mutate(dif_absorb = `1` - `0`) %>%
          filter(lambda > 350) %>%
          dplyr::select(lambda, dif_absorb)
        correlation <- cor(color_dif$dif_color, absorb_dif$dif_absorb)
        tmp <- tibble(type="reflectance-absorbance",
                      drainage=chosen_drainage,
                      body_part=chosen_body_part,
                      correlation=correlation)
        body_abs_correlations <- rbind(body_abs_correlations, tmp)
        
        print(paste(chosen_body_part, chosen_drainage, sulphur_sample, non_sulphur_sample))
        
      }
    }
  }
}
write_tsv(body_abs_correlations, "output/body_absorbance_correlations.txt")




#Sensitivity

sensitivity <- read_tsv("output/opsin/sensitivity.A1.txt") %>%
  inner_join(meta_data %>% dplyr::select(Fieldsite.ID,H2S)) 





#Correlation between body color and sensitivity 
for (chosen_body_part in body.parts){
  body_sens_correlations <- tibble(type=character(),
                                   drainage=character(),
                                   body_part=character(),
                                   correlation=numeric())
  for (chosen_drainage in drainages){
    sulphur_samples <- color %>% 
      filter(Drainage == chosen_drainage, body.part == chosen_body_part) %>%
      filter(H2S == 1) %>% pull(individual) %>% unique()
    non_sulphur_samples <- color %>% 
      filter(Drainage == chosen_drainage, body.part == chosen_body_part) %>%
      filter(H2S == 0) %>% pull(individual) %>% unique()
    sulphur_senses <- sensitivity %>% 
      filter(Drainage == chosen_drainage) %>%
      filter(H2S == 1) %>% pull(Identifier) %>% unique()
    non_sulphur_senses <- sensitivity %>% 
      filter(Drainage == chosen_drainage) %>%
      filter(H2S == 0) %>% pull(Identifier) %>% unique()
    
    for (sulphur_sample in sulphur_samples){
      for (non_sulphur_sample in non_sulphur_samples){
        color_dif <- color %>%
          filter(body.part == chosen_body_part) %>%
          filter(individual == sulphur_sample | individual == non_sulphur_sample) %>%
          dplyr::select(H2S,lambda, mean_reflectance) %>%
          pivot_wider(names_from = H2S, values_from = mean_reflectance) %>%
          mutate(dif_color = `1` - `0`) %>%
          dplyr::select(lambda, dif_color)
        print(paste(chosen_body_part, chosen_drainage, sulphur_sample, non_sulphur_sample))
        for (sulphur_sense in sulphur_senses){
          for (non_sulphur_sense in non_sulphur_senses){
            sense_dif <- sensitivity %>%
              filter(Identifier == sulphur_sense | Identifier == non_sulphur_sense) %>%
              dplyr::select(H2S,lambda, lambda_sensitivity) %>%
              pivot_wider(names_from = H2S, values_from = lambda_sensitivity) %>%
              mutate(dif_sense = `1` - `0`) %>%
              dplyr::select(lambda, dif_sense)
        
            correlation <- cor(color_dif$dif_color, sense_dif$dif_sense)
            tmp <- tibble(type="reflectance-sensitivity",
                      drainage=chosen_drainage,
                      body_part=chosen_body_part,
                      correlation=correlation)
            body_sens_correlations <- rbind(body_sens_correlations, tmp)
        
          }
        }
      }
    }
  }
  write_tsv(body_sens_correlations, paste0("output/body_sensitivity_correlations.",chosen_body_part,".txt"))
}



#Correlation between absorbance and sensitivity 
abs_sens_correlations <- tibble(type=character(),
                                drainage=character(),
                                body_part=character(),
                                correlation=numeric())
for (chosen_drainage in drainages){
  sulphur_senses <- sensitivity %>% 
    filter(Drainage == chosen_drainage) %>%
    filter(H2S == 1) %>% pull(Identifier) %>% unique()
  non_sulphur_senses <- sensitivity %>% 
    filter(Drainage == chosen_drainage) %>%
    filter(H2S == 0) %>% pull(Identifier) %>% unique()
  for (sulphur_sense in sulphur_senses){
    for (non_sulphur_sense in non_sulphur_senses){
      print(paste(chosen_drainage, sulphur_sense, non_sulphur_sense))
      absorb_dif <- absorbance %>%
        filter(Drainage == chosen_drainage) %>%
        group_by(H2S,lambda) %>%
        dplyr::summarise(median_absorbance = median(absorbance)) %>%
        pivot_wider(names_from = H2S, values_from = median_absorbance) %>%
        mutate(dif_absorb = `1` - `0`) %>%
        filter(lambda > 350) %>%
        dplyr::select(lambda, dif_absorb)
      sense_dif <- sensitivity %>%
        filter(Identifier == sulphur_sense | Identifier == non_sulphur_sense) %>%
        dplyr::select(H2S,lambda, lambda_sensitivity) %>%
        filter(lambda > 350) %>%
        pivot_wider(names_from = H2S, values_from = lambda_sensitivity) %>%
        mutate(dif_sense = `1` - `0`) %>%
        dplyr::select(lambda, dif_sense)
      
      correlation <- cor(absorb_dif$dif_absorb, sense_dif$dif_sense)
      tmp <- tibble(type="absorbance-sensitivity",
                    drainage=chosen_drainage,
                    body_part="NA",
                    correlation=correlation)
      abs_sens_correlations <- rbind(abs_sens_correlations, tmp)
      
    }
  }
}
write_tsv(abs_sens_correlations, paste0("output/absorbance_sensitivity_correlations.txt"))

##########################
##Using median value for non-sulphur and individuals for sulphur
##########################
abs_sens_correlations <- tibble(type=character(),
                                drainage=character(),
                                body_part=character(),
                                correlation=numeric())
for (chosen_drainage in drainages){
  sulphur_senses <- sensitivity %>% 
    filter(Drainage == chosen_drainage) %>%
    filter(H2S == 1) %>% pull(Identifier) %>% unique()
  sulphur_absorbs <- absorbance %>%
    filter(Drainage == chosen_drainage,H2S == 1) %>% 
    mutate(id = paste0(site.name,"_",sample.location)) %>%
    pull(id) %>% unique()
  nonsulphur_sensitivity <- sensitivity %>%
    filter(Drainage == chosen_drainage,
           H2S == 0) %>%
    dplyr::select(H2S,lambda, lambda_sensitivity) %>%
    group_by(lambda) %>%
    dplyr::summarize(`0` = median(lambda_sensitivity)) %>%
    ungroup()
  nonsulphur_absorbance <- absorbance %>%
    filter(Drainage == chosen_drainage,H2S == 0) %>%
    group_by(lambda) %>%
    dplyr::summarise(`0` = median(absorbance)) %>% 
    ungroup()
  for (sulphur_sense in sulphur_senses){
    for (sulphur_absorb in sulphur_absorbs){
      print(paste(chosen_drainage, sulphur_sense,sulphur_absorb ))
      absorb_dif <- absorbance %>%
        filter(Drainage == chosen_drainage) %>%
        mutate(id = paste0(site.name,"_",sample.location)) %>%
        filter(id == sulphur_absorb) %>%
        group_by(lambda) %>%
        dplyr::summarise(`1` = median(absorbance)) %>%
        inner_join(nonsulphur_absorbance) %>%
        mutate(dif_absorb = `1` - `0`) %>%
        filter(lambda > 350) %>%
        dplyr::select(lambda, dif_absorb)
      sense_dif <- sensitivity %>%
        filter(Identifier == sulphur_sense) %>%
        dplyr::select(lambda, lambda_sensitivity) %>%
        filter(lambda > 350) %>%
        rename(`1` = lambda_sensitivity) %>%
        inner_join(nonsulphur_sensitivity) %>%
        mutate(dif_sense = `1` - `0`) %>%
        dplyr::select(lambda, dif_sense)
      
      correlation <- cor(absorb_dif$dif_absorb, sense_dif$dif_sense)
      tmp <- tibble(type="absorbance-sensitivity",
                    drainage=chosen_drainage,
                    body_part="NA",
                    correlation=correlation)
      abs_sens_correlations <- rbind(abs_sens_correlations, tmp)
      
    }
  }
}
abs_sens_correlations %>%
  group_by(drainage) %>%
  dplyr::summarize(mean = mean(correlation)) %>%
  ungroup() %>%
  dplyr::summarise(mean = mean(mean))
#####################
##Permutations
#####################
abs_sens_permutations <- tibble()
for (perm in 1:100){
  print(paste0("Permutation ",perm))
  for (chosen_drainage in drainages){
    abs_sens_perm_correlations <- tibble()
    perm_sense_samples <- sensitivity %>% 
      select(Drainage, H2S, Sample.ID) %>% 
      unique() %>% 
      filter(Drainage == chosen_drainage) %>% 
      mutate(H2S = sample(H2S, replace = F))
    perm_absorb_samples <- absorbance %>%
      filter(Drainage == chosen_drainage) %>%
      mutate(id = paste0(site.name,"_",sample.location)) %>%
      select(id, H2S) %>% unique() %>%
      mutate(H2S = sample(H2S, replace = F))
    sulphur_senses <- sensitivity %>% 
      filter(Drainage == chosen_drainage) %>%
      select(-H2S) %>%
      inner_join(perm_sense_samples, by = "Sample.ID") %>%
      filter(H2S == 1) %>% pull(Identifier) %>% unique()
    sulphur_absorbs <- absorbance %>%
      filter(Drainage == chosen_drainage) %>% 
      mutate(id = paste0(site.name,"_",sample.location)) %>%
      select(-H2S) %>%
      inner_join(perm_absorb_samples, by = "id") %>%
      filter(H2S == 1) %>%
      pull(id) %>% unique()
    nonsulphur_sensitivity <- sensitivity %>%
      filter(Drainage == chosen_drainage) %>%
      select(-H2S) %>%
      inner_join(perm_sense_samples, by = "Sample.ID") %>%
      filter(H2S == 0) %>%
      dplyr::select(H2S,lambda, lambda_sensitivity) %>%
      group_by(lambda) %>%
      dplyr::summarize(`0` = median(lambda_sensitivity)) %>%
      ungroup()
    nonsulphur_absorbance <- absorbance %>%
      filter(Drainage == chosen_drainage) %>%
      mutate(id = paste0(site.name,"_",sample.location)) %>%
      select(-H2S) %>%
      inner_join(perm_absorb_samples, by = "id") %>%
      filter(H2S == 0) %>%
      group_by(lambda) %>%
      dplyr::summarise(`0` = median(absorbance)) %>% 
      ungroup()
    for (sulphur_sense in sulphur_senses){
      for (sulphur_absorb in sulphur_absorbs){
        #print(paste(chosen_drainage, sulphur_sense,sulphur_absorb ))
        absorb_dif <- absorbance %>%
          filter(Drainage == chosen_drainage) %>%
          mutate(id = paste0(site.name,"_",sample.location)) %>%
          filter(id == sulphur_absorb) %>%
          group_by(lambda) %>%
          dplyr::summarise(`1` = median(absorbance)) %>%
          inner_join(nonsulphur_absorbance, by = "lambda") %>%
          mutate(dif_absorb = `1` - `0`) %>%
          filter(lambda > 350) %>%
          dplyr::select(lambda, dif_absorb)
        sense_dif <- sensitivity %>%
          filter(Identifier == sulphur_sense) %>%
          dplyr::select(lambda, lambda_sensitivity) %>%
          filter(lambda > 350) %>%
          rename(`1` = lambda_sensitivity) %>%
          inner_join(nonsulphur_sensitivity, by = "lambda") %>%
          mutate(dif_sense = `1` - `0`) %>%
          dplyr::select(lambda, dif_sense)
        
        correlation <- cor(absorb_dif$dif_absorb, sense_dif$dif_sense)
        tmp <- tibble(type="absorbance-sensitivity",
                      drainage=chosen_drainage,
                      body_part="NA",
                      correlation=correlation)
        abs_sens_perm_correlations <- rbind(abs_sens_correlations, tmp)
        
      }
    }
    tmp <- abs_sens_perm_correlations %>%
      group_by(drainage) %>%
      dplyr::summarize(mean = mean(correlation)) 
    abs_sens_permutations <- rbind(abs_sens_permutations, tmp)
  }

}
abs_sens_permutations %>%
  ggplot(.,aes(mean)) + geom_histogram() +
  facet_wrap(~drainage)

##########################
##Getting mean correlation values
##########################
median_correlations <- tibble()
for (chosen_body_part in body.parts){
  for (chosen_drainage in drainages){
    color_dif <- color %>%
      filter(body.part == chosen_body_part) %>%
      filter(Drainage == chosen_drainage) %>%
      group_by(H2S,lambda) %>%
      dplyr::summarize(median_reflectance = median(mean_reflectance)) %>%
      pivot_wider(names_from = H2S, values_from = median_reflectance) %>%
      mutate(dif_color = `1` - `0`) %>%
      filter(lambda > 350) %>%
      dplyr::select(lambda, dif_color)
    absorb_dif <- absorbance %>%
      filter(Drainage == chosen_drainage) %>%
      group_by(H2S,lambda) %>%
      dplyr::summarise(median_absorbance = median(absorbance)) %>%
      pivot_wider(names_from = H2S, values_from = median_absorbance) %>%
      mutate(dif_absorb = `1` - `0`) %>%
      filter(lambda > 350) %>%
      dplyr::select(lambda, dif_absorb)
    correlation <- cor(color_dif$dif_color, absorb_dif$dif_absorb)
    tmp <- tibble(type="reflectance-absorbance",
                  drainage=chosen_drainage,
                  body_part=chosen_body_part,
                  correlation=correlation)
    median_correlations <- rbind(median_correlations, tmp)
    color_dif <- color %>%
      filter(body.part == chosen_body_part) %>%
      filter(Drainage == chosen_drainage) %>%
      group_by(H2S,lambda) %>%
      dplyr::summarize(median_reflectance = median(mean_reflectance)) %>%
      pivot_wider(names_from = H2S, values_from = median_reflectance) %>%
      mutate(dif_color = `1` - `0`) %>%
      dplyr::select(lambda, dif_color)
    sense_dif <- sensitivity %>%
      filter(Drainage == chosen_drainage) %>%
      dplyr::select(H2S,lambda, lambda_sensitivity) %>%
      group_by(H2S,lambda) %>%
      dplyr::summarize(median_sensitivity = median(lambda_sensitivity)) %>%
      pivot_wider(names_from = H2S, values_from = median_sensitivity) %>%
      mutate(dif_sense = `1` - `0`) %>%
      dplyr::select(lambda, dif_sense)
    correlation <- cor(color_dif$dif_color, sense_dif$dif_sense)
    tmp <- tibble(type="reflectance-sensitivity",
                  drainage=chosen_drainage,
                  body_part=chosen_body_part,
                  correlation=correlation)
    median_correlations <- rbind(median_correlations, tmp)
  }
}

for (chosen_drainage in drainages){
  sense_dif <- sensitivity %>%
    filter(Drainage == chosen_drainage) %>%
    dplyr::select(H2S,lambda, lambda_sensitivity) %>%
    group_by(H2S,lambda) %>%
    dplyr::summarize(median_sensitivity = median(lambda_sensitivity)) %>%
    pivot_wider(names_from = H2S, values_from = median_sensitivity) %>%
    mutate(dif_sense = `1` - `0`) %>%
    filter(lambda > 350) %>%
    dplyr::select(lambda, dif_sense)
  absorb_dif <- absorbance %>%
    filter(Drainage == chosen_drainage) %>%
    group_by(H2S,lambda) %>%
    dplyr::summarise(median_absorbance = median(absorbance)) %>%
    pivot_wider(names_from = H2S, values_from = median_absorbance) %>%
    mutate(dif_absorb = `1` - `0`) %>%
    filter(lambda > 350) %>%
    dplyr::select(lambda, dif_absorb)
  correlation <- cor(sense_dif$dif_sense, absorb_dif$dif_absorb)
  tmp <- tibble(type="absorbance-sensitivity",
                drainage=chosen_drainage,
                body_part="NA",
                correlation=correlation)
  median_correlations <- rbind(median_correlations, tmp)
}
##########################
##Plotting all correlations
##########################
correlation_files <- list.files("./output/", "correlations")
all_correlations <- tibble()
for (file in correlation_files){
  tmp <- read_tsv(paste0("output/",file))
  all_correlations <- rbind(all_correlations, tmp)
}  
all_correlations  %>%
  ggplot(.,aes(x=body_part,y=correlation,fill=drainage)) +
  geom_violin() +
  theme_cowplot() +
  geom_hline(yintercept=0,linetype="dotted") +
  facet_wrap(~type,scales="free_x") +
  scale_fill_brewer(palette = "Set1",name="Drainage") 

median_correlations  %>%
  arrange(type, body_part, -correlation) %>%
  mutate(id = row_number()) %>%
  ggplot(.,aes(x=id,y=correlation,color=drainage)) +
  geom_point() +
  theme_cowplot() +
  geom_hline(yintercept=0,linetype="dotted") +
  scale_color_brewer(palette = "Set1",name="Drainage") +
  annotate("segment",x=0.5,xend=20.5,y=1,yend=1,size=3,color="#56B4E9") +
  annotate("segment",x=0.5,xend=4.5,y=1.1,yend=1.1,size=3,color="#009E73") +
  annotate("segment",x=4.5,xend=8.5,y=1.2,yend=1.2,size=3,color="#E69F00") +
  annotate("segment",x=8.5,xend=12.5,y=1.3,yend=1.3,size=3,color="#F0E442") +
  annotate("segment",x=12.5,xend=16.5,y=1.4,yend=1.4,size=3,color="#CC79A7") +
  annotate("segment",x=16.5,xend=20.5,y=1.5,yend=1.5,size=3,color="#D55E00") +
  annotate("segment",x=20.5,xend=36.5,y=1.1,yend=1.1,size=3,color="#009E73") +
  annotate("segment",x=20.5,xend=24.5,y=1.2,yend=1.2,size=3,color="#E69F00") +
  annotate("segment",x=24.5,xend=28.5,y=1.3,yend=1.3,size=3,color="#F0E442") +
  annotate("segment",x=28.5,xend=32.5,y=1.4,yend=1.4,size=3,color="#CC79A7") +
  annotate("segment",x=32.5,xend=36.5,y=1.5,yend=1.5,size=3,color="#D55E00") +
  geom_vline(xintercept=c(0.5,4.5,8.5,12.5,16.5,20.5,24.5,28.5,32.5,36.5)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Correlation") +
  coord_flip()

pdf("figures/correlation_medians.v1.pdf",useDingbats = F,height=4,width=6)
median_correlations  %>%
  arrange(type, body_part, -correlation) %>%
  mutate(id = paste(type, body_part,sep=" ")) %>%
  ggplot(.,aes(x=fct_rev(id),y=correlation,color=drainage)) +
  geom_point(size=3) +
  theme_cowplot() +
  geom_hline(yintercept=0,linetype="dotted") +
  scale_color_brewer(palette = "Set1",name="Drainage") +
  ylab("Correlation") +
  coord_flip(ylim=c(-1,1)) +
  theme(axis.title.y=element_blank(),
        panel.grid.major.y=element_line(size=0.1)) 
dev.off()
 
median_correlations %>%
  filter(type == "reflectance-absorbance",) %>%
  dplyr::select(correlation) %>% dplyr::pull() %>% t.test(.) 
  

library(rstatix)
median_correlations %>% 
  mutate(id = paste(type,body_part)) %>%
  group_by(id) %>% rstatix::t_test(correlation ~ 1) %>%
  dplyr::select(id,n,statistic,df, p) %>%
  write_tsv("output/correlation_medians_ttests.txt")

median_correlations %>% 
  mutate(id = paste(type,body_part)) %>%
  group_by(id) %>% 
  dplyr::summarize(mean = mean(correlation), sd = sd(correlation)) %>%
  write_tsv("output/correlation_medians_stats.txt")

all_correlations  %>%
  filter(type != "absorbance-sensitivity") %>%
  mutate(quantized_cor = case_when(correlation > 0.1 ~ "positive",
                                   correlation < -0.1 ~ "negative",
                                   TRUE ~ "null")) %>%
  ggplot(.,aes(x=drainage,fill=fct_rev(quantized_cor))) +
  geom_bar(position = "fill", stat="count") +
  theme_cowplot() +
  geom_hline(yintercept=0.5,linetype="dotted") +
  facet_grid(body_part~type,scales="free_x") +
  scale_fill_brewer(palette = "Set1", name="Correlation",
                    labels=c("Positive (r > 0.1)",
                             "Null (0.1 >= r >= -0.1)",
                             "Negative (r < -0.1)")) +
  ylab("Permutation count") + xlab("Drainage") +
  scale_x_discrete(guide = guide_axis(n.dodge=2))

all_correlations  %>%
  filter(type == "absorbance-sensitivity") %>%
  
  mutate(quantized_cor = case_when(correlation > 0.1 ~ "positive",
                                   correlation < -0.1 ~ "negative",
                                   TRUE ~ "null")) %>%
  ggplot(.,aes(x=drainage,fill=fct_rev(quantized_cor))) +
  geom_bar(position = "fill", stat="count") +
  theme_cowplot() +
  geom_hline(yintercept=0.5,linetype="dotted") +
  facet_grid(~type,scales="free_x") +
  scale_fill_brewer(palette = "Set1", name="Correlation",
                    labels=c("Positive (r > 0.1)",
                             "Null (0.1 >= r >= -0.1)",
                             "Negative (r < -0.1)")) +
  ylab("Permutation count") + xlab("Drainage") +
  scale_x_discrete(guide = guide_axis(n.dodge=2))

median_correlations  %>%
  arrange(type, body_part, -correlation) %>%
  mutate(body_part = case_when(body_part == "NA" ~ NA_character_,
                               TRUE ~ body_part)) %>%
  mutate(id = as.factor(row_number())) %>% 
  dplyr::select(-correlation)-> correlation_ids

pdf("figures/correlation_permutations.v1.pdf",useDingbats = F,height=6,width=6)
all_correlations  %>%
  inner_join(correlation_ids) %>%
  ggplot(.) +
  aes(x = correlation, y = fct_rev(id),fill=drainage) +
  annotate("rect",xmin=-1,xmax=1,ymin=1,ymax=5,fill="grey",alpha=0.3) +
  annotate("rect",xmin=-1,xmax=1,ymin=9,ymax=13,fill="grey",alpha=0.3) +
  annotate("rect",xmin=-1,xmax=1,ymin=17,ymax=21,fill="grey",alpha=0.3) +
  annotate("rect",xmin=-1,xmax=1,ymin=25,ymax=29,fill="grey",alpha=0.3) +
  annotate("rect",xmin=-1,xmax=1,ymin=33,ymax=37,fill="grey",alpha=0.3) +
  
  geom_density_ridges(
    scale = 3, bandwidth = 0.05,
    rel_min_height = 0.01,
    color = "black",
  ) +
  scale_x_continuous(
    name = "Correlation",
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    name = NULL,
    expand = expansion(add = c(0.2, 2.6))
  ) +
  theme_minimal_grid() +
  theme(
    axis.text.y = element_blank()
  ) +
  scale_fill_brewer(
    palette = "Set1",
    name="Drainage"
  ) +
  annotate("segment",x=-1,xend=-1,y=37,yend=17,size=3,color="#56B4E9")  +
  annotate("segment",x=-1.05,xend=-1.05,y=37,yend=33,size=3,color="#009E73") +
  annotate("segment",x=-1.1,xend=-1.1,y=33,yend=29,size=3,color="#E69F00") +
  annotate("segment",x=-1.15,xend=-1.15,y=29,yend=25,size=3,color="#F0E442") +
  annotate("segment",x=-1.2,xend=-1.2,y=25,yend=21,size=3,color="#CC79A7") +
  annotate("segment",x=-1.25,xend=-1.25,y=21,yend=17,size=3,color="#D55E00") +
  annotate("segment",x=-1.05,xend=-1.05,y=17,yend=1,size=3,color="#009E73") +
  annotate("segment",x=-1.1,xend=-1.1,y=17,yend=13,size=3,color="#E69F00") +
  annotate("segment",x=-1.15,xend=-1.15,y=13,yend=9,size=3,color="#F0E442") +
  annotate("segment",x=-1.2,xend=-1.2,y=9,yend=5,size=3,color="#CC79A7") +
  annotate("segment",x=-1.25,xend=-1.25,y=5,yend=1,size=3,color="#D55E00") +
  annotate("segment",x=-1.3,xend=-1.3,y=5,yend=1,size=3,color="WHITE")
dev.off()


######
#Example plot
######
color_palette_1 <- pnw_palette(name="Bay",n=2,type="discrete")
chosen_body_part <- "belly"
chosen_drainage <-"Pichucalco"
plot1 <- color %>%
  filter(body.part == chosen_body_part) %>%
  filter(Drainage == chosen_drainage) %>%
  group_by(H2S,lambda) %>%
  dplyr::summarize(median_reflectance = median(mean_reflectance)) %>%
  pivot_wider(names_from = H2S, values_from = median_reflectance) %>%
  filter(lambda > 350) %>%
  mutate(ymax = case_when(`0` > `1` ~ `0`,
                          TRUE ~ `1`),
         ymin = case_when(`0` > `1` ~ `1`,
                          TRUE ~ `0`),
         type = case_when(`0` > `1` ~ "pos",
                          TRUE ~ "neg")) %>%
  ggplot(aes(x=lambda,fill=type)) +
  #geom_line() +
  theme_cowplot() +
  geom_ribbon(aes(ymax = ymax, ymin = ymin), alpha = 0.4,color=NA) +
  scale_fill_manual(values=rev(color_palette_1)) +
  geom_line(aes(x=lambda, y=`0`),color="#0f85a0") +
  geom_line(aes(x=lambda, y=`1`),color="#dd4124") +
  ylab("Belly reflectance") +
  xlab("Wavelength") +
  theme(legend.position = "none")



plot2 <- sensitivity %>%
  filter(Drainage == chosen_drainage) %>%
  dplyr::select(H2S,lambda, lambda_sensitivity) %>%
  group_by(H2S,lambda) %>%
  dplyr::summarize(median_sensitivity = median(lambda_sensitivity)) %>%
  pivot_wider(names_from = H2S, values_from = median_sensitivity) %>%
  filter(lambda > 350) %>%
  mutate(ymax = case_when(`0` > `1` ~ `0`,
                          TRUE ~ `1`),
         ymin = case_when(`0` > `1` ~ `1`,
                          TRUE ~ `0`),
         type = case_when(`0` > `1` ~ "pos",
                          TRUE ~ "neg")) %>%
  ggplot(aes(x=lambda,fill=type)) +
  #geom_line() +
  theme_cowplot() +
  geom_ribbon(aes(ymax = ymax, ymin = ymin), alpha = 0.4,color=NA) +
  scale_fill_manual(values=rev(color_palette_1)) +
  geom_line(aes(x=lambda, y=`0`),color="#0f85a0") +
  geom_line(aes(x=lambda, y=`1`),color="#dd4124") +
  ylab("Visual Sensitivity") +
  xlab("Wavelength") +
  theme(legend.position = "none")

color_tmp <- color %>%
  filter(body.part == chosen_body_part) %>%
  filter(Drainage == chosen_drainage) %>%
  group_by(H2S,lambda) %>%
  dplyr::summarize(median_reflectance = median(mean_reflectance)) %>%
  pivot_wider(names_from = H2S, values_from = median_reflectance) %>%
  filter(lambda > 350) %>%
  mutate(dif_color = `1` - `0`) %>%
  dplyr::select(lambda, dif_color)
  
sense_tmp <- sensitivity %>%
  filter(Drainage == chosen_drainage) %>%
  dplyr::select(H2S,lambda, lambda_sensitivity) %>%
  group_by(H2S,lambda) %>%
  dplyr::summarize(median_sensitivity = median(lambda_sensitivity)) %>%
  pivot_wider(names_from = H2S, values_from = median_sensitivity) %>%
  filter(lambda > 350) %>%
  mutate(dif_sense = `1` - `0`) %>%
  dplyr::select(lambda, dif_sense)

plot3 <- color_tmp %>%
  inner_join(sense_tmp) %>%
  ggplot(.,aes(x=dif_color,y=dif_sense,color=lambda)) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  geom_point() +
  theme_minimal() +
  scale_color_viridis_c(name="Wavelength") +
  ylab("Sensitivity difference") +
  xlab("Colour difference") +
  theme(axis.text=element_text(colour="black"))

pdf("figures/correlations_example.v1.pdf",useDingbats = F,height=4,width=10)
plot1 + plot2 + plot3
dev.off()
