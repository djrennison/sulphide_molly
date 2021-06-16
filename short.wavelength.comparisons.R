library(tidyverse)
library(PNWColors)
library(patchwork)
library(cowplot)
library(infer)
color_palette_1 <- pnw_palette(name="Bay",n=2,type="discrete")
#This is for calculating comparisons with light absorbance only in the short wavelengths.
meta_data <- read_csv("data/meta_data/site_meta_data.csv")
drainages <- sort(unique(meta_data$Drainage))


#Transmission should actually be thought of as absorbance
absorbance <- read_tsv("output/transmission/filtered.light.type_s_w_5.irr.depth_10.txt") %>%
  mutate(site.name = case_when(site.name == "Exp" ~ "Esp",
                               site.name == "lab" ~ "Lab",
                               site.name == "Vet" ~ "VS",
                               site.name == "VC" ~ "VG",
                               TRUE ~ site.name)) %>%
  rename(absorbance = transmission)

short.wave.absorbance <- absorbance %>%
  inner_join(meta_data %>% rename(site.name = ID)) %>%
  filter(lambda <= 400) %>% 
  group_by(site.name,sample.location, Drainage, H2S,Species, Informal_Name) %>%
  dplyr::summarize(short_wave_absorbance = sum(absorbance)) 

t.test.results <- tibble(statistic=numeric(),t_df=numeric(),p_value=numeric(),
                         alternative=character(),lower_ci=numeric(),
                         upper_ci=numeric(),Drainage=character(),type=character(),
                         location=character())
for (chosen_drainage in drainages){
  t.test.results <- rbind(t.test.results, short.wave.absorbance %>%
    filter(Drainage == chosen_drainage) %>%
    mutate(H2S = as.factor(H2S)) %>%
    t_test(formula = short_wave_absorbance ~ H2S, 
           order = c("0", "1"),
           alternative = "two-sided") %>%
    mutate(Drainage = chosen_drainage, type = "absorbance",location="NA"))
}


absorbance.plot <- short.wave.absorbance %>%
  ggplot(.,aes(x=Drainage,y=short_wave_absorbance,color=as.factor(H2S))) + 
  geom_jitter(position=position_jitterdodge()) +
  geom_boxplot(outlier.shape=NA,alpha=0.7) +
  scale_color_manual(values=color_palette_1, labels=c("Non-sulphur","Sulphur"),
                     name="Environment") +
  ylab("Light absorbance ( <= 400 nm )") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) 

#For color, 
window.w.col <- 5
color <-read_tsv( paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".summary.tsv", sep = ""))

color.medians <- color %>%
  inner_join(meta_data %>% dplyr::select(ID,H2S,Drainage) %>% rename(population=ID)) %>%
  group_by(population, lambda,H2S,Drainage,sex,body.part) %>%
  dplyr::summarize(median_reflectance = median(mean_reflectance)) %>%
  ungroup()


short.wave.reflectance <- color %>%
  inner_join(meta_data %>% dplyr::select(ID,H2S,Drainage) %>% rename(population=ID)) %>%
  filter(lambda <= 400) %>%
  group_by(individual, population,H2S,Drainage, body.part) %>%
  dplyr::summarise(short_reflectance = sum(mean_reflectance)/n()) %>%
  rename(site.name = population)

for (chosen_drainage in drainages){
  for (chosen_body_part in sort(unique(short.wave.reflectance$body.part))){
    t.test.results <- rbind(t.test.results, short.wave.reflectance %>%
      filter(Drainage == chosen_drainage,
             body.part==chosen_body_part) %>%
      mutate(H2S = as.factor(H2S)) %>%
      t_test(formula = short_reflectance ~ H2S, 
             order = c("0", "1"),
             alternative = "two-sided") %>%
      mutate(Drainage = chosen_drainage, type = "reflectance",location=chosen_body_part))
  }

}

reflectance.plot <- short.wave.reflectance %>%
  ungroup() %>%
  mutate(body.part = case_when(body.part == "tail" ~ "Tail",
                               body.part == "eye" ~ "Eye",
                               body.part == "head" ~ "Head",
                               body.part == "belly" ~ "Belly")) %>%
  mutate(H2S = case_when(H2S == 1 & site.name == "Glo" ~ 2,
                         TRUE ~ H2S)) %>%
  ggplot(.,aes(x=Drainage,y=short_reflectance,color=as.factor(H2S))) + 
  geom_jitter(position=position_jitterdodge()) + 
  geom_boxplot(outlier.shape = NA,alpha=0.6) +
  facet_wrap(~body.part,scales="free_y") +
  scale_color_manual(values=color_palette_1[c(1,2,2)], labels=c("Non-sulphur","Sulphur"),
                     name="Environment",
                     breaks = c(0, 1)) +
  ylab("Skin reflectance ( <= 400 nm )") +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 60,hjust = 1))





sensitivity <- read_tsv("output/opsin/sensitivity.A1.txt")
short.wave.sensitivity <-  sensitivity %>%
  inner_join(meta_data %>% dplyr::select(Fieldsite.ID,H2S)) %>%
  filter(lambda <= 400) %>%
  group_by(Identifier,Fieldsite.ID,Drainage,H2S) %>%
  dplyr::summarize(short_sensitivity = mean(lambda_sensitivity))

for (chosen_drainage in drainages){
  t.test.results <- rbind(t.test.results,short.wave.sensitivity %>%
                              filter(Drainage == chosen_drainage) %>%
                              mutate(H2S = as.factor(H2S)) %>%
                              t_test(formula = short_sensitivity ~ H2S, 
                                     order = c("0", "1"),
                                     alternative = "two-sided") %>%
                              mutate(Drainage = chosen_drainage, type = "sensitivity",location="NA"))
}

sensitivity.plot <- short.wave.sensitivity %>%
  ggplot(.,aes(x=Drainage,y=short_sensitivity,color=as.factor(H2S))) + 
  geom_jitter(position=position_jitterdodge()) + 
  geom_boxplot(outlier.shape=NA,alpha=0.6) +
  scale_color_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
                     name="Environment",
                     breaks = c(0, 1)) +
  theme_cowplot() + 
  ylab("Visual sensitivity ( <= 400 nm )") +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) 

write_tsv(t.test.results, "output/shortwavelength_ttests.txt")
pdf("figures/short_wavelength_comparison.v1.pdf",height=12,width=14,useDingbats = F)
(absorbance.plot | sensitivity.plot) / reflectance.plot + plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect')
dev.off()



