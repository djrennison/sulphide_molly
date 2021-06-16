library(tidyverse)
library(cowplot)
library(PNWColors)
library(ggforce)

arc <- tibble(s=1.5*pi,e=2.5*pi)
angles <- read_csv("output/angle_analysis.csv")
pdf("figures/angles_direction.v1.pdf",height=6,width=3)
angles %>%
  mutate(rad = Angle * pi / 180) %>%
  mutate(drainage_type = case_when(Watershed2 == "Tacotalpa" ~ "Tacotalpa",
                                   TRUE ~ "Other")) %>%
  ggplot() +
  annotate("segment",x=-0.5,xend=0.5,y=0,yend=0,linetype="dashed") +
  annotate("segment",x=0,xend=0,y=0,yend=0.5,linetype="dashed") +
  geom_spoke(aes(x=0,y=0,angle=rad,color=drainage_type),radius = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_wrap(~fct_relevel(Trait,"Opsins"),ncol=1) +
  geom_arc(data=arc, aes(x0 = 0, y0 = 0, r = 0.5, start = s, end = e),
           linetype="dashed") +
  scale_color_manual(values=c("black","grey"),name="Drainage")
dev.off()
