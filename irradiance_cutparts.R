
# === Ban
irr.Ban <- irr.all.smooth[irr.all.smooth$population == "Ban",]
k.Ban <- f.calc.transmission.old.code(irr.Ban, params)
tBan <- k.Ban[[1]]
irr.Ban %>%
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = irradiance, -population, -sample.location, -file.name, -light,
                -depth,-integr.time) %>%
  #filter(light == "d",sample.location == "3") %>%
  filter(sample.location != "surf") %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=irradiance,group=depth,color=depth)) + geom_line() +
  facet_grid(light~sample.location,scales="free_y")


tBan %>% 
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = transmission, -site.name, -sample.location, -light,
                -n.measures,-at.bottom,-nls.converged) %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=transmission,color=sample.location)) + geom_line()

# some weird extreme outliers

# === Vet
irr.Vet <- irr.all.smooth[irr.all.smooth$population == "Vet",]
k.Vet <- f.calc.transmission.old.code(irr.Vet, params)
tVet <- k.Vet[[1]]
irr.Vet %>%
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = irradiance, -population, -sample.location, -file.name, -light,
                -depth,-integr.time) %>%
  #filter(light == "d",sample.location == "3") %>%
  filter(sample.location != "surf") %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=irradiance,group=depth,color=depth)) + geom_line() +
  facet_grid(light~sample.location,scales="free_y") +
  theme_bw()

tVet %>% 
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = transmission, -site.name, -sample.location, -light,
                -n.measures,-at.bottom,-nls.converged) %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=transmission,color=sample.location)) + geom_line() +
  theme_bw()


med.LC <- fig.light.measures(tLC, wave.length.range.temp, 3, "check.kd.all.")


# === PS0
irr.PS0 <- irr.all.smooth[irr.all.smooth$population == "PS0",]
k.PS0 <- f.calc.transmission.old.code(irr.PS0, params)
tPS0 <- k.PS0[[1]]
irr.PS0 %>%
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = irradiance, -population, -sample.location, -file.name, -light,
                -depth,-integr.time) %>%
  #filter(light == "d",sample.location == "3") %>%
  filter(sample.location != "surf") %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=irradiance,group=depth,color=depth)) + geom_line() +
  facet_grid(light~sample.location,scales="free_y")

tPS0 %>% 
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = transmission, -site.name, -sample.location, -light,
                -n.measures,-at.bottom,-nls.converged) %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=transmission,color=sample.location)) + geom_line() +
  theme_bw()

# === Ixt
irr.Ixt <- irr.all.smooth[irr.all.smooth$population == "Ixt",]
k.Ixt <- f.calc.transmission.old.code(irr.Ixt, params)
tIxt <- k.Ixt[[1]]
irr.Ixt %>%
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = irradiance, -population, -sample.location, -file.name, -light,
                -depth,-integr.time) %>%
  #filter(light == "d",sample.location == "3") %>%
  filter(sample.location != "surf") %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=irradiance,group=depth,color=depth)) + geom_line() +
  facet_grid(light~sample.location,scales="free_y")

tIxt %>% 
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = transmission, -site.name, -sample.location, -light,
                -n.measures,-at.bottom,-nls.converged) %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=transmission,color=sample.location)) + geom_line() +
  theme_bw()

# === Exp
irr.Exp <- irr.all.smooth[irr.all.smooth$population == "Exp",]
k.Exp <- f.calc.transmission.old.code(irr.Exp, params)
tExp <- k.Exp[[1]]
irr.Exp %>%
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = irradiance, -population, -sample.location, -file.name, -light,
                -depth,-integr.time) %>%
  #filter(light == "d",sample.location == "3") %>%
  filter(sample.location != "surf") %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=irradiance,group=depth,color=depth)) + geom_line() +
  facet_grid(light~sample.location,scales="free_y")

tExp %>% 
  as.tibble() %>%
  tidyr::gather(key = wavelength, value = transmission, -site.name, -sample.location, -light,
                -n.measures,-at.bottom,-nls.converged) %>%
  ggplot(.,aes(x=as.numeric(wavelength),y=transmission,color=sample.location)) + geom_line() +
  theme_bw()




# === "Oyster"
irr.O <- irr.all.smooth[irr.all.smooth$site.name == "Oyster",]
k.O <- f.calc.transmission.old.code(irr.O, params)
tO <- k.O[[1]]
med.O <- fig.light.measures(tO, wave.length.range.temp, 3, "check.kd.all.")
# 02 and 05 are way off

# === "Paxton"
irr.Pa <- irr.all.smooth[irr.all.smooth$site.name == "Paxton",]
k.Pa <- f.calc.transmission.old.code(irr.Pa, params)
tPa <- k.Pa[[1]]
med.Pa <- fig.light.measures(tPa, wave.length.range.temp, 3, "check.kd.all.")
# some extreme outlier in extreme low nm, otehrwise good

# === "Priest"
irr.Pr <- irr.all.smooth[irr.all.smooth$site.name == "Priest",]
k.Pr <- f.calc.transmission.old.code(irr.Pr, params)
tPr <- k.Pr[[1]]
med.Pr <- fig.light.measures(tPr, wave.length.range.temp, 3, "check.kd.all.")
# some extreme variation

# === "Trout"
irr.T <- irr.all.smooth[irr.all.smooth$site.name == "Trout",]
k.T <- f.calc.transmission.old.code(irr.T, params)
tT <- k.T[[1]]
med.T <- fig.light.measures(tT, wave.length.range.temp, 3, "check.kd.all.")


kd.all <- f.get.median(overview.k, wave.length.range.temp)

# --- 2) check for outliers calculate kd -------
# use the figures made above and output of the negative values 
# Kirk
tK[1,]
tK[6,]

K.out <- c("Kirk.01", "Kirk.06") 
# Little Campbell
# looks fine

# Oyster
tO[2,]
tO[5,]
O.out <- c("Oyster.02", "Oyster.05")
# both are very large values => remove 2 and 5

# Paxton
# looks fine

# Priest
tPr[6,]
plot(unlist(tPr[6,-c(1:7)]))
tPr[8,]
plot(unlist(tPr[8,-c(1:7)]))
out.Pr <- c("Priest.06", "Priest.08")

# Trout
# looks ok
rows.weird.Kd.names <- c(K.out, O.out, out.Pr)

# --- 3) compare both approached -------
# determine which rows have weird value and need to be excluded
overv.K.merge <- paste(overview.k$site.name,overview.k$sample.location, sep = ".")
rows.weird.Kd <- match(rows.weird.Kd.names, overv.K.merge)

# exclude these values to check the effect on the median
overview.k.select <- overview.k[-rows.weird.Kd,]
kd.select <- f.get.median(overview.k.select, wave.length.range.temp)
# compare both approaches
plot(NA, xlim = range(wave.length.range.temp), ylim = c(-0.005, 0.005), las = 1)
for(i in 1:nrow(kd.select)) lines(wave.length.range.temp, (kd.all[i, -1] - kd.select[i,-1]) )
# exlcuding does not change the median: still good to exlude the clearly odd ones I think both from the irradiance
# and transmission files


# NOT ANYMORE? There is a problem with some Inf values => first replace any inf with the mean of the value before and after
# === replace Inf values with mean =====
overview.k.t <- overview.k[,c((length(trans.names) + 1):ncol(overview.k))]
ncol(overview.k.t)
cntr.inf <- 0
for(i in 1:nrow(overview.k.t)){
  # hope the first input is not inf: so start at the second till the second last
  for(j in 2:(ncol(overview.k.t) - 1)){
    if(overview.k.t[i,j] == "Inf"){
      overview.k.t[i,j] <- mean(overview.k.t[i,(j - 1)],overview.k.t[i,(j + 1)])
      cntr.inf <- cntr.inf + 1
    }   
  }
}
cntr.inf