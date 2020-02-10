#Comparing opsin measured sensitivity to transmission, irradiance, and body part reflectance.

library(tidyverse)

#First calculate all the differences between the reference (0) and test regions (1).
meta_data <- read_csv("data/meta_data/site_meta_data.csv")

#Transmission should actually be thought of as absorbance
absorbance <- read_tsv("output/transmission/light.type_s_w_5.irr.depth_10txt") %>%
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

absorbance.medians %>%
  ggplot(.,aes(x=lambda,y=median_absorbance,group=site.name)) + 
  facet_wrap(~Drainage,scales="free_y") + 
  geom_line(aes(color=as.factor(H2S)))

absorbance.medians.sulphuraria <- absorbance.medians %>%
  filter(site.name == "VS" | site.name == "Glo") %>%
  mutate(Drainage = "Pichucalco_2")

absorbance.differences <- absorbance.medians %>%
  filter(site.name != "Glo") %>%
  rbind(absorbance.medians.sulphuraria) %>%
  ungroup() %>%
  dplyr::select(Drainage,lambda,H2S,median_absorbance) %>%
  pivot_wider(names_from = H2S, values_from = median_absorbance) %>% 
  mutate(absorbance_dif = `0` - `1`)  %>%
  dplyr::select(-`0`,-`1`)

absorbance.differences %>%
  ggplot(.,aes(x=lambda,y=absorbance_dif)) + facet_wrap(~Drainage,scales="free_y") + geom_line()

#For absorbance, positive values means more light for the test (1) location


irradiance <- read_csv("output/irradiance/irr.all.smooth.w=5.csv") %>%
  pivot_longer(c(-population,-sample.location,-file.name,-light,-depth,-integr.time),
               names_to = "lambda", values_to = "irradiance")

irradiance %>% 
  filter(light == "s",depth == 10) %>%
  ggplot(.,aes(x=as.numeric(lambda),y=irradiance,group=file.name)) + 
  geom_line() +
  facet_wrap(~population)

irradiance %>%
  filter(light == "s",depth == 10) %>%
  group_by(population,lambda) %>%
  dplyr::summarise(median_irradiance = median(irradiance)) %>%
  inner_join(meta_data %>% rename(population = ID)) %>%
  ungroup() -> irradiance.medians

irradiance.medians %>% 
  ggplot(.,aes(x=as.numeric(lambda),y=median_irradiance,group=population)) + 
  geom_line(aes(color=as.factor(H2S))) +
  facet_wrap(~Drainage,scales="free_y")

irradiance.medians.sulphuraria <- irradiance.medians %>%
  filter(population == "VS" | population == "Glo") %>%
  mutate(Drainage = "Pichucalco_2")

irradiance.differences <- irradiance.medians %>%
  filter(population != "Glo") %>%
  rbind(irradiance.medians.sulphuraria) %>%
  ungroup() %>%
  dplyr::select(Drainage,lambda,H2S,median_irradiance) %>%
  pivot_wider(names_from = H2S, values_from = median_irradiance) %>% 
  mutate(irradiance_dif = `1` -`0` )  %>%
  dplyr::select(-`0`,-`1`)
  
irradiance.differences %>%
  mutate(lambda = as.numeric(lambda)) %>%
  ggplot(.,aes(x=lambda,y=irradiance_dif)) + facet_wrap(~Drainage,scales="free_y") + geom_line()

#For color, 
window.w.col <- 5
color <-read_tsv( paste( path.output.colour, "colour.all.smooth.normalised.w=", window.w.col, ".summary.tsv", sep = ""))

color.medians <- color %>%
  inner_join(meta_data %>% dplyr::select(ID,H2S,Drainage) %>% rename(population=ID)) %>%
  group_by(population, lambda,H2S,Drainage,sex,body.part) %>%
  dplyr::summarize(median_reflectance = median(mean_reflectance)) %>%
  ungroup()



color.medians.sulphuraria <- color.medians %>%
  filter(population == "VS" | population == "Glo") %>%
  mutate(Drainage = "Pichucalco_2")

color.differences <- color.medians %>%
  filter(population != "Glo") %>%
  rbind(color.medians.sulphuraria) %>%
  ungroup() %>%
  dplyr::select(Drainage,sex,body.part, lambda,H2S,median_reflectance) %>%
  pivot_wider(names_from = H2S, values_from = median_reflectance) %>% 
  mutate(reflectance_dif = `1` -`0` )  %>%
  dplyr::select(-`0`,-`1`)

color.differences %>%
  ggplot(.,aes(x=lambda,y=reflectance_dif,color=Drainage)) + geom_line() +
  facet_grid(sex~body.part) +
  geom_hline(aes(yintercept=0),linetype="dotted")

#Visual sensitivity measures

sensitivity <- read_tsv("output/opsin/sensitivity.A1.txt")
reference_sensitivity <- sensitivity %>%
  inner_join(meta_data %>% dplyr::select(Fieldsite.ID,H2S)) %>%
  filter(H2S == 0) %>%
  group_by(Drainage, lambda) %>%
  dplyr::summarize(reference_median_sensitivity = median(lambda_sensitivity)) %>%
  ungroup()

sensitivity %>%
  filter(lambda <= 400) %>%
  group_by(Identifier,Drainage,H2S) %>%
  summarize(total_uv_sen = sum(lambda_sensitivity)) %>%
  ggplot(.,aes(x=Drainage,y=total_uv_sen,color=as.factor(H2S))) + geom_boxplot()
sensitivity.differences <- sensitivity %>%
  inner_join(meta_data %>% dplyr::select(Fieldsite.ID,H2S)) %>%
  filter(H2S == "1") %>%
  ungroup() %>%
  mutate(lambda=as.numeric(lambda)) %>%
  inner_join(reference_sensitivity) %>%
  mutate(sensitivity_dif = lambda_sensitivity - reference_median_sensitivity ) %>%
  dplyr::select(Identifier,lambda,Fieldsite.ID,Drainage,sensitivity_dif)

sensitivity.differences %>%
  filter(lambda <= 400) %>%
  ggplot(.,aes(x=lambda,y=sensitivity_dif,group=Identifier)) + geom_line() +
  facet_wrap(~Drainage) +
  geom_hline(yintercept=0,linetype="dotted")

sensitivity %>%
  mutate(wavelength_type = case_when(lambda <= 450 ~ "short",
                                     lambda >= 600 ~ "long",
                                     TRUE ~ "middle")) %>%
  group_by(Identifier,Drainage,H2S,wavelength_type) %>%
  summarize(total_sen = sum(lambda_sensitivity)) %>%
  filter(wavelength_type != "middle") %>%
  ggplot(.,aes(x=Drainage,y=total_sen,color=as.factor(H2S))) + 
  geom_boxplot() +
  facet_wrap(~wavelength_type,scales="free_y")


#Now compare sensitivity against irradiance, tranmission and color
test_samples <- unique(sensitivity.differences$Identifier)

sensitivity.light.correlations <- tibble(Identifier=double(),
                                         Fieldsite.ID=character(),Drainage=character(),
                                         sex=character(),type=character(),correlation=numeric())
for (chosen_sample in test_samples){
  
  chosen_fieldsite <- sensitivity.differences %>%
    filter(Identifier == chosen_sample) %>% pull(Fieldsite.ID) %>% unique()
  chosen_drainage <- sensitivity.differences %>%
    filter(Identifier == chosen_sample) %>% pull(Drainage) %>% unique()
  
  
  tmp.1 <- sensitivity.differences %>%
    filter(Identifier == chosen_sample) %>%
    mutate(lambda = as.character(lambda)) %>%
    inner_join(irradiance.differences)
  irradiance.cor <- stats::cor(tmp.1$irradiance_dif, tmp.1$sensitivity_dif)
  
  tmp.tibble <- tibble(Identifier=chosen_sample,Fieldsite.ID=chosen_fieldsite,
                       Drainage=chosen_drainage,type="water_irradiance",
                       correlation=irradiance.cor,sex="NA")
  sensitivity.light.correlations <- rbind(sensitivity.light.correlations,tmp.tibble)
  
  
  tmp.2 <- sensitivity.differences %>%
    filter(Identifier == chosen_sample) %>%
    mutate(lambda = lambda) %>%
    inner_join(absorbance.differences)
  absorbance.cor <- stats::cor(tmp.2$absorbance_dif, tmp.2$sensitivity_dif)
  
  tmp.tibble <- tibble(Identifier=chosen_sample,Fieldsite.ID=chosen_fieldsite,
                       Drainage=chosen_drainage,type="water_absorbance",
                       correlation=absorbance.cor,sex="NA")
  sensitivity.light.correlations <- rbind(sensitivity.light.correlations,tmp.tibble)
  
  for (chosen_sex in c("Male","Female")){
    for (chosen_body in c("head","eye","tail","belly")){
      tmp.3 <- sensitivity.differences %>%
        filter(Identifier == chosen_sample) %>%
        inner_join(color.differences %>% filter(sex == chosen_sex, body.part == chosen_body) )
      reflectance.cor <- stats::cor(tmp.3$reflectance_dif, tmp.3$sensitivity_dif)
      tmp.tibble <- tibble(Identifier=chosen_sample,Fieldsite.ID=chosen_fieldsite,
                           Drainage=chosen_drainage,type=chosen_body,
                           correlation=reflectance.cor,sex=chosen_sex)
      sensitivity.light.correlations <- rbind(sensitivity.light.correlations,tmp.tibble)
      
    }
  }

}

sensitivity.light.correlations %>%
  filter(type == "water_irradiance") %>%
  ggplot(.) + geom_jitter(aes(x=Drainage,y=correlation),width=0.1) +
  geom_hline(yintercept = 0,linetype="dotted") +
  ggtitle("Sensitivity to Irradiance")

sensitivity.light.correlations %>%
  filter(type == "water_absorbance") %>%
  ggplot(.) + geom_jitter(aes(x=Drainage,y=correlation),width=0.1) +
  geom_hline(yintercept = 0,linetype="dotted")

sensitivity.light.correlations %>%
  filter(sex != "NA") %>%
  ggplot(.) + geom_jitter(aes(x=Drainage,y=correlation),width=0.1) +
  geom_hline(yintercept = 0,linetype="dotted") +
  facet_grid(sex~type)

