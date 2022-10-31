# Conduct ABC model inference from SLiM data
# Connor Murray 5.3.2022
# ijob -A berglandlab_standard --mem=10G -p standard -c 1
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(viridis)
library(gganimate)
library(abc)
library(abc.data)
library(ggrepel)
library(ggforce)

# Working directory
setwd("/project/berglandlab/connor/bottleneck/")

### Simulation data ###

# Simulation FST data
data.fst <- data.table(read.csv("data.fst.10reps.csv", header = T))

# simulation PCA data
data.pca <- data.table(read.csv("data.pca.10reps.csv", header = T))

# Simulation EUC data
data.euc <- data.table(read.csv("data.euc.10reps.csv", header = T))

# simulation AF variance data
data.afvar <- data.table(read.csv("data.afvar.10reps.csv", header = T))

# Aggregate simulation summaries
mean.sim <- data.table(data.fst %>% 
                 merge(data.pca,
                       by=c('model', 'nMax', 'nMin', 'replicate')) %>% 
                 merge(data.euc,
                       by=c('model', 'nMax', 'nMin', 'replicate')) %>% 
                 merge(data.afvar,
                       by=c('model', 'nMax', 'nMin', 'replicate')))
                      
# Load empirical data FST
load("Year_to_year_object.Rdata")

dt.test <- data.table(Out_comp_vector_samepops)

empir.fst <- data.table(Out_comp_vector_samepops %>% 
                          group_by(dd.sim = bin_date,
                                   population = pop1) %>% 
                          summarize(mean.FST = FST) %>% 
                          mutate(dd.sim = case_when(dd.sim == "2.Overwinter" ~ "2.overwinter",
                                                    dd.sim == "3.Multi-Year" ~ "3.multi",
                                                    TRUE ~ as.character(dd.sim))))

# Load empirical pca results excluding 2L
load("PCA.object.include2L.filtered.Rdata")
load("PCA.object.exclude2L.filtered.Rdata")

# Excluding 2L: pca.dt2
# Including 2L: pca.dt 

# Aggregate empirical mean FST and PCA 
mean.empir <- data.table(empir.fst %>%
                               filter(!dd.sim == "3.multi", 
                                      population == "Charlottesville") %>% 
                               group_by(dd.sim) %>% 
                               summarise(mean.FST = mean(mean.FST, na.rm = TRUE)) %>% 
                               pivot_wider(names_from = dd.sim, 
                                           names_prefix = "mean.FST_",
                                           values_from=mean.FST) %>% 
                           cbind(pca.dt2 %>% 
                               group_by(dd.sim) %>% 
                               summarise(mean.pc1 = mean(Dim.1, na.rm = TRUE),
                                         mean.pc2 = mean(Dim.2, na.rm = TRUE),
                                         mean.pc3 = mean(Dim.3, na.rm = TRUE)) %>% 
                                 pivot_wider(names_from = dd.sim,
                                             values_from=c(mean.pc1, mean.pc2, mean.pc3))) %>% 
                           cbind(pca.dt %>% 
                               distinct(mean.r2pc1 = r2pc1,
                                         mean.r2pc2 = r2pc2,
                                         mean.r2pc3 = r2pc3)))

# Model information
models <- data.table(nMax=as.numeric(tstrsplit(mean.sim$model, "_")[[1]]),
                     nMin=as.numeric(tstrsplit(mean.sim$model, "_")[[2]])) %>% 
          mutate(bottle=((nMax-nMin)/nMax)*100)

colnames(mean.empir)

# ABC stats
sim.stat <- mean.sim[, -c(1:4)] %>% 
            relocate(colnames(mean.empir))
  
# Thresholds to test in ABC
thresh <- seq(from=0.01, to = 0.15, by = 0.001)

# Statistics to use
cols.stats <- which(colnames(mean.empir) %like% "mean.pc" == FALSE)

# Go through several thresholds
tmp.out <- foreach(i=1:length(thresh), .combine = "rbind") %do% {
  
  # ABC using rejection methods
  tmp <- abc(target = mean.empir %>% 
               dplyr::select(!(mean.pc1_1.within:mean.pc3_2.overwinter)), 
             param = models %>% dplyr::select(-c(bottle)), 
             sumstat = sim.stat %>% 
               dplyr::select(!c(mean.pc1_1.within:mean.pc3_2.overwinter, 
                                mean.FST_3.multi, n.FST,
                                sd.FST_1.within, sd.FST_2.overwinter, sd.FST_3.multi)),
             tol=thresh[i], 
             method="rejection")
  
  # Compile data
  tmp <- data.table(summary(tmp), 
                    mean.empir,
                    threshold=thresh[i])
  
  # Progress message
  print(paste("Threshold:", thresh[i], sep=" "))

  # Finish
  return(tmp)
}

# Fix data
colnames(tmp.out)[1:3] <- c("data", "parameter", "population")
tmp.out[data %in% c("2.5% Perc.:","97.5% Perc.:","Mean:")]

out <- data.table(pivot_wider(tmp.out, names_from = data, values_from = population))
colnames(out)[14:20] <- c("Min.", "Perc0.25", "Median", "Mean", "Mode", "Perc0.975", "Max.") 

pdf("thresholds.abc.rejection.exclude2L.nopcdims.bottleneck.xlim.more.pdf")

# ABC and tolerance thresholds
ggplot(out[!parameter=="bottle"][threshold %in% seq(0.09, 0.11, by=0.001)], 
       aes(x=threshold, 
           y=as.numeric(Mean),
           ymin=as.numeric(Perc0.25),
           ymax=as.numeric(Perc0.975),
           color=parameter)) +
  geom_pointrange(
        #alpha=0.7,
        position = "dodge"
                  ) +
  #geom_vline(xintercept=0.1, size=1.1) +
  labs(x="Thresholds", 
       y="Bottleneck size %", 
       color="Parameter") +
  theme_classic() +
  scale_fill_viridis(option = "magma", begin=0, end=1) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=14),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16))

dev.off()

# Sum of squares error generator 
sim.stat %>% 
  dplyr::select(!(mean.pc1_1.within:mean.pc3_2.overwinter)) %>% 
  as.matrix(.) -> a

mean.empir %>%
  dplyr::select(!(mean.pc1_1.within:mean.pc3_2.overwinter)) %>% 
  as.matrix(.) -> b

# Sum of squares
dt <- data.table(data.table(models,
                 sse=as.numeric(a %>% 
                                sweep(., 2, b) %>% 
                                .^2 %>%
                                rowSums(.))) %>% 
                 group_by(nMin, nMax) %>% 
                 summarise(sse.mean=mean(sse)))

out[!parameter=="bottle"][threshold %like% 0.101]

pdf("sse.raster.100reps.new2.pdf", width = 9, height = 8)

#out[!parameter=="bottle"][threshold %in% seq(0.09, 0.11, by=0.001)][parameter=="nMax"]

# Plot Sum of squares error
ggplot() +
  geom_raster(data=dt[!nMax==100000],
              aes(x=as.factor(nMin),
                  y=as.factor(nMax),
                  fill=sse.mean)) +
  geom_point(data=expand.grid(c(1000,1500,2000),
                              seq(25000, 30000, by=1000)),
             aes(x=as.factor(Var1),
                 y=as.factor(Var2)), color="white", size=2) +
  labs(x="Minimum population size", 
       y="Maximum population size", 
       fill="Mean SSE") +
  theme_classic() +
  scale_x_discrete(breaks = c(50, 2500, 7500, 25000, 
                              50000, 75000)) +
  scale_y_discrete(breaks = c(50, 2500, 7500, 25000,
                              50000, 75000)) +
  scale_fill_viridis(option = "magma", begin=0, end=1) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()

pdf("bottle.raster.100reps.new1.pdf", width = 13, height = 8)

# Bottleneck size
dt %>% 
  mutate(bottle=(nMax-nMin)/nMax) %>% 
  ggplot(.,
       aes(x=nMin,
           y=nMax,
           #fill=bottle*100,
           color=sse.mean)
       ) +
  geom_point() +
  labs(x="Minimum population size", 
       y="Maximum population size", 
       fill="Bottleneck size (%)") +
  theme_classic() +
  #scale_x_discrete(breaks = c(50, 2500, 7500, 25000, 
  #                            50000, 100000)) +
  #scale_y_discrete(breaks = c(50, 2500, 7500, 25000,
  #                            50000, 100000)) +
  #annotation_logticks(scaled = F) +
  scale_color_viridis(option = "D", begin=0, end=1) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()