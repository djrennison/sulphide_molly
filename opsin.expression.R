#
# code for analysing the opsin expression data
#
#

library(PNWColors)
library(ggforce)
library(patchwork)
library(cowplot)
color_palette_1 <- pnw_palette(name="Bay",n=2,type="discrete")

# 1) analysing opsin as proportions
#

# apparently we are going fancy and use dplyr and freak'n ggplot2
# 
#

# load opsin expression data
#opsin.exp <- read.csv( paste(path.data.clean.expr, "combined_molly_qpcr_noind362.csv", sep = ""), stringsAsFactors = FALSE) %>%
#  mutate(Efficiency = Efficiancy)
opsin.exp <- read.csv( paste(path.data.clean.expr, "combined_molly_qpcr_noind362_ben_efficiency.csv", sep = ""), stringsAsFactors = FALSE) %>%
  rename(Efficiency = Ben_abs) 
#xx

# load meta data that comes with the expression data
opsin.m.d <- read.csv( paste(path.data.meta.d, "Molly samples.csv", sep = ""), stringsAsFactors = FALSE)
#yy 
opsin.exp <- inner_join(opsin.exp, opsin.m.d)


# with file xxx do (%>%) xxxx
# mutate
opsin.exp %>%
  group_by(Gene,Identifier,Efficiency) %>%
  dplyr::summarize(meanCT = mean(Ct), CT_stdev = sd(Ct)) %>%
  mutate(expression = (1/((1 + Efficiency)^meanCT))) %>%
  group_by(Identifier) %>%
  mutate(total_cone = sum(expression)) %>%
  ungroup() %>%
  mutate(percent_cone = expression/total_cone) %>%
  inner_join(opsin.m.d) %>%
  ggplot(.,aes(x=Gene,y=percent_cone)) + geom_boxplot(aes(fill=as.factor(H2S))) +
  facet_grid(~Drainage,scales="free") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=color_palette_1,name="Sulphur",
                    labels=c("Absent", "Present"))

opsin.exp %>%
  group_by(Gene,Identifier,Efficiency) %>%
  dplyr::summarize(meanCT = mean(Ct), CT_stdev = sd(Ct)) %>%
  mutate(expression = (1/((1 + Efficiency)^meanCT))) %>%
  group_by(Identifier) %>%
  mutate(total_cone = sum(expression)) %>%
  ungroup() %>%
  mutate(percent_cone = expression/total_cone) %>%
  inner_join(opsin.m.d) %>%
  ggplot(.,aes(x=Drainage,y=percent_cone)) + geom_boxplot(aes(fill=as.factor(H2S))) +
  facet_wrap(~Gene,scales="free_y",nrow=1) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=color_palette_1,name="Sulphur",
                    labels=c("Absent", "Present")) +
  ylab("Percent expression")
  

opsin.exp %>%
  group_by(Gene,Identifier,Efficiency) %>%
  dplyr::summarize(meanCT = mean(Ct), CT_stdev = sd(Ct)) %>%
  mutate(expression = (1/((1 + Efficiency)^meanCT))) %>%
  group_by(Identifier) %>%
  mutate(total_cone = sum(expression)) %>%
  ungroup() %>%
  mutate(percent_cone = expression/total_cone) %>%
  inner_join(opsin.m.d) -> opsin.exp.summarized


# Add in sensitivity to get overall sensitivity
opsin_sensitivity <- read_csv("output/opsin/absorbance.prop.A1_1.csv") %>%
  rename(Gene=X1) %>%
  pivot_longer(-Gene, names_to = "lambda", values_to = "absorbance")

opsin.exp.summarized %>%
  inner_join(opsin_sensitivity) %>%
  group_by(Identifier, lambda,Gene) %>%
  mutate(gene_sensitivity = sum(percent_cone*absorbance)) %>%
  group_by(Identifier,lambda) %>% 
  dplyr::summarize(lambda_sensitivity = sum(gene_sensitivity,na.rm=T)) %>%
  inner_join(opsin.m.d) %>%
  ggplot(.) +
  geom_line(aes(x=as.numeric(lambda),y=lambda_sensitivity,group=Identifier,color=as.factor(H2S))) +
  facet_wrap(~Drainage) +
  scale_color_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) +
  ylab("Visual sensitivity") + xlab("Lambda")




opsin.exp.summarized %>%
  inner_join(opsin_sensitivity) %>%
  group_by(Identifier, lambda,Gene) %>%
  mutate(gene_sensitivity = sum(percent_cone*absorbance)) %>%
  group_by(Identifier,lambda) %>% 
  dplyr::summarize(lambda_sensitivity = sum(gene_sensitivity,na.rm=T)) %>%
  inner_join(opsin.m.d) %>%
  group_by(Fieldsite.ID,Drainage,lambda,H2S) %>%
  do(data.frame(t(quantile(.$lambda_sensitivity, probs = c(0.025,0.50, 0.975))))) %>%
  rename(bottom = X2.5., mid = X50., top =X97.5.) %>% 
  ggplot(.,aes()) + 
  geom_ribbon(aes(x=as.numeric(lambda),ymin=(bottom),ymax=(top),group=Fieldsite.ID,fill=as.factor(H2S)),alpha=0.5) +
  geom_line(aes(x=as.numeric(lambda),y=(mid),group=Fieldsite.ID,color=as.factor(H2S)),alpha=1,size=0.4) +  
  facet_wrap(~Drainage) +
  ylab("Sensitivity from opsins") + xlab("nm") +
  scale_fill_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present"))  +
  scale_color_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) 
  


sensitivity <- opsin.exp.summarized %>%
  inner_join(opsin_sensitivity) %>%
  group_by(Identifier, lambda,Gene) %>%
  mutate(gene_sensitivity = sum(percent_cone*absorbance)) %>%
  group_by(Identifier,lambda) %>% 
  dplyr::summarize(lambda_sensitivity = sum(gene_sensitivity,na.rm=T)) %>%
  inner_join(opsin.m.d) %>%
  group_by(Fieldsite.ID,Drainage,lambda,H2S) %>%
  dplyr::summarize(median_sensitivity = median(lambda_sensitivity)) 

write_tsv(sensitivity,"output/opsin/sensitivity.median.A1.txt")


sensitivity <- opsin.exp.summarized %>%
  inner_join(opsin_sensitivity) %>%
  group_by(Identifier, lambda,Gene) %>%
  mutate(gene_sensitivity = sum(percent_cone*absorbance)) %>%
  group_by(Identifier,lambda) %>% 
  dplyr::summarize(lambda_sensitivity = sum(gene_sensitivity,na.rm=T)) %>%
  inner_join(opsin.m.d) 

write_tsv(sensitivity,"output/opsin/sensitivity.A1.txt")


#PCA using inferred sensitivity.
sensitivity %>%
  dplyr::mutate(nm_bin = floor(as.numeric(lambda)/5)*5) %>%
  group_by(Identifier, Drainage,H2S,Fieldsite.ID,nm_bin) %>%
  dplyr::summarise(summed_sensitivity = sum(lambda_sensitivity)) %>%
  ungroup() %>%
  dplyr::select(Identifier,nm_bin,summed_sensitivity) %>%
  pivot_wider(names_from = nm_bin, values_from = summed_sensitivity) ->sensitivity.forpca
pca.info <- sensitivity %>%
  dplyr::select(Identifier, Drainage,H2S,Fieldsite.ID) %>%
  ungroup() %>%
  unique()

sensitivity.pca <- prcomp(sensitivity.forpca[,c(2:ncol(sensitivity.forpca))], 
                  center = TRUE,scale. = TRUE)
sensitivity.pca.tibble <- as.tibble(sensitivity.pca$x[,c(1:6)])
cbind(sensitivity.pca.tibble,pca.info) %>%
  ggplot(.,aes(x=PC1,y=PC2,color=as.factor(H2S))) +
  geom_point() +
  facet_wrap(~Drainage,nrow=1) +
  scale_color_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) +
  ggtitle(paste("Inferred Sensitivity"))

#PCA of opsin expression directly
opsin.exp %>%
  group_by(Gene,Identifier,Efficiency) %>%
  dplyr::summarize(meanCT = mean(Ct), CT_stdev = sd(Ct)) %>%
  mutate(expression = (1/((1 + Efficiency)^meanCT))) %>%
  group_by(Identifier) %>%
  mutate(total_cone = sum(expression)) %>%
  ungroup() %>%
  mutate(percent_cone = expression/total_cone) %>%
  select(Gene,Identifier,percent_cone) %>%
  pivot_wider(names_from = Gene, values_from = percent_cone) -> opsin.exp.forpca

pca.info <- opsin.exp %>%
  dplyr::select(Identifier, Drainage,H2S,Fieldsite.ID) %>%
  ungroup() %>%
  unique() 

opsin.exp.pca <- prcomp(opsin.exp.forpca[,c(2:ncol(opsin.exp.forpca))], 
                          center = TRUE,scale. = TRUE)
opsin.exp.pca.tibble <- as.tibble(opsin.exp.pca$x[,c(1:6)])
cbind(opsin.exp.pca.tibble,pca.info) %>%
  ggplot(.,aes(x=PC1,y=PC2,color=as.factor(H2S))) +
  geom_point() +
  facet_wrap(~Drainage,nrow=1) +
  scale_color_manual(values=color_palette_1,name="Sulphur",labels=c("Absent","Present")) +
  ggtitle(paste("Opsin Expression"))

summary(lm(percent_cone ~ Gene + H2S, data = opsin.exp.summarized))
