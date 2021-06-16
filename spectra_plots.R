library(tidyverse)
library(PNWColors)
library(cowplot)

#Make a plot of light irradiance and absorbance
meta_data <- read_csv("data/meta_data/site_meta_data.csv")

#Absorbance
absorbance <- read_tsv("output/transmission/filtered.light.type_s_w_5.irr.depth_10.txt") %>%
  mutate(site.name = case_when(site.name == "Exp" ~ "Esp",
                               site.name == "lab" ~ "Lab",
                               site.name == "Vet" ~ "VS",
                               site.name == "VC" ~ "VG",
                               TRUE ~ site.name)) %>%
  rename(absorbance = transmission)

absorbance %>%
  filter(light == "s") %>%
  group_by(site.name,lambda) %>%
  dplyr::summarise(median_absorbance = median(absorbance)) %>%
  inner_join(meta_data %>% rename(site.name = ID)) %>%
  ungroup() -> absorbance.medians


plot1 <- absorbance.medians %>%
  ggplot(.,aes(x=lambda,y=1-median_absorbance,group=site.name)) + 
  facet_wrap(~Drainage,nrow=1) + 
  geom_line(aes(color=as.factor(H2S))) +
  theme_cowplot() +
  scale_color_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
                     name="Environment",
                     breaks = c(0, 1)) +
  ylab("Transmission (Ks)") + xlab("Wavelength (nm)")




irradiance <- read_csv("output/irradiance/irr.all.smooth.w=5.csv") %>%
  pivot_longer(c(-population,-sample.location,-file.name,-light,-depth,-integr.time),
               names_to = "lambda", values_to = "irradiance")

plot2 <- irradiance %>% 
  inner_join(meta_data %>% rename(population = ID)) %>%
  filter(light == "s", lambda > 350) %>%
  mutate(H2S = case_when(H2S == 0 ~ "Non-sulphur",
                         TRUE ~ "Sulphur")) %>%
  ggplot(.,aes(x=as.numeric(lambda),y=irradiance,group=file.name)) + 
  geom_line(aes(color=depth)) +
  theme_cowplot() +
  facet_grid(H2S~Drainage) +
  scale_color_viridis_c(name="Depth (cm)") +
  #scale_color_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
  #                   name="Environment",
  #                   breaks = c(0, 1)) +
  ylab("Normalized irradiance") + xlab("Wavelength (nm)") 


#Now gather the color data

window.w.col <- 5
color <-read_tsv( paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".summary.tsv", sep = ""))

plot3 <- color %>%
  mutate(body.part = case_when(body.part == "tail" ~ "Tail",
                               body.part == "eye" ~ "Eye",
                               body.part == "head" ~ "Head",
                               body.part == "belly" ~ "Belly")) %>%
  inner_join(meta_data %>% dplyr::select(ID,H2S,Drainage) %>% rename(population=ID)) %>%
  group_by(population, lambda,H2S,Drainage,body.part) %>%
  filter(lambda < 700) %>%
  dplyr::summarize(sample_size = sqrt(n()),
                   mean_reflectance_2 = mean(mean_reflectance), std_reflectance=2*sd(mean_reflectance)/sample_size) %>%
  ggplot(.,aes(x=lambda,y=mean_reflectance_2,group=population,color=as.factor(H2S),fill=as.factor(H2S))) + 
  theme_cowplot() + 
  geom_line() +
  geom_ribbon(aes(ymin=mean_reflectance_2-std_reflectance, ymax=mean_reflectance_2+std_reflectance),alpha=0.3,color=NA) +
  facet_grid(body.part~Drainage) +
  scale_color_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
                     name="Environment",
                     breaks = c(0, 1)) +
  scale_fill_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
                     name="Environment",
                     breaks = c(0, 1)) +
  ylab("Reflectance") + xlab("Wavelength (nm)")

plot4 <- color %>%
  mutate(body.part = case_when(body.part == "tail" ~ "Tail",
                               body.part == "eye" ~ "Eye",
                               body.part == "head" ~ "Head",
                               body.part == "belly" ~ "Belly")) %>%
  inner_join(meta_data %>% dplyr::select(ID,H2S,Drainage) %>% rename(population=ID)) %>%
  group_by(population,sex, lambda,H2S,Drainage,body.part) %>%
  filter(lambda < 700) %>%
  dplyr::summarize(sample_size = sqrt(n()),
    mean_reflectance_2 = mean(mean_reflectance), std_reflectance=2*sd(mean_reflectance)/sample_size) %>%
  mutate(popsex = paste0(population,sex)) %>%
  ggplot(.,aes(x=lambda,y=mean_reflectance_2,group=popsex,color=as.factor(H2S),linetype=sex,fill=as.factor(H2S))) + 
  theme_cowplot() + 
  geom_line() +
  geom_ribbon(aes(ymin=mean_reflectance_2-std_reflectance, ymax=mean_reflectance_2+std_reflectance),alpha=0.3,color=NA) +
  facet_grid(body.part~Drainage) +
  scale_color_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
                     name="Environment",
                     breaks = c(0, 1)) +
  scale_fill_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
                    name="Environment",
                    breaks = c(0, 1)) +
  scale_linetype(name="Sex") + 
  ylab("Reflectance") + xlab("Wavelength (nm)")

pdf("figures/sensitivity_spectrum.v2.pdf",height=6,width=8,useDingbats = F)
read_tsv("output/opsin/sensitivity.A1.txt") %>%
  inner_join(meta_data %>% dplyr::select(Fieldsite.ID,H2S)) %>%
  group_by(lambda, Drainage, H2S) %>%
  dplyr::summarize(sample_size = sqrt(n()),
    mean_sensitivity_2 = mean(lambda_sensitivity), se_sensitivity=2*sd(lambda_sensitivity)/sample_size) %>%
  ggplot(.,aes(x=lambda,y=mean_sensitivity_2,group=H2S,color=as.factor(H2S),fill=as.factor(H2S))) + 
  theme_cowplot() + 
  geom_line() +
  geom_ribbon(aes(ymin=mean_sensitivity_2-se_sensitivity, ymax=mean_sensitivity_2+se_sensitivity),alpha=0.3,color=NA) +
  facet_wrap(~Drainage,nrow=2) +
  scale_color_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
                     name="Environment",
                     breaks = c(0, 1)) +
  scale_fill_manual(values=color_palette_1[c(1,2)], labels=c("Non-sulphur","Sulphur"),
                    name="Environment",
                    breaks = c(0, 1)) +
  scale_linetype(name="Sex") + 
  ylab("Sensitivity") + xlab("Wavelength (nm)")
dev.off()


pdf("figures/All_spectra.v1.pdf",height=15,width=10)
(plot2 / plot1 ) /  plot4 + plot_layout(heights = c(1, 1,3)) +plot_annotation(tag_levels = 'A')
dev.off()

pdf("figures/Absorbance_spectrum.v1.pdf",height=3,width=8,useDingbats = F)
plot1
dev.off()

pdf("figures/Irradiance_spectrum.v2.pdf",height=6,width=8,useDingbats = F)
plot2
dev.off()

pdf("figures/Reflectance_spectrum.v1.pdf",height=4,width=8,useDingbats = F)
plot4
dev.off()
