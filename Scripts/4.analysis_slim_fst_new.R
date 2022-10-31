# Compiles and analyze genomic data from SLiM and PLINK
# Connor Murray 11.19.2021
# ijob -A berglandlab_standard --mem=5G -p standard -c 10
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(doParallel)
library(viridis)

# VCF output folder
setwd("/project/berglandlab/connor/BACKUP_scratch/slim_bottleneck/data_new3")

# All frequency data
filenames1 <- list.files(pattern = ".csv$")
filenames1 <- filenames1[!filenames1 %like% "pca"]

# Register cores
#doParallel::registerDoParallel(cores = 10)

# Extract metadata from filenames
meta <- data.table(data.table(seed=tstrsplit(filenames1, ".", fixed=T)[[7]],
                   nMax = tstrsplit(filenames1, ".", fixed=T)[[4]],
                   nMin = tstrsplit(filenames1, ".", fixed=T)[[5]],
                   replicate = tstrsplit(filenames1, ".", fixed=T)[[6]],
                   file=filenames1) %>% 
            group_by(nMax, nMin) %>% 
            mutate(run_number=row_number(),
                   run=paste(nMax, nMin,sep="_")))

# Distinct models
runs <- unique(meta$run)

# Confidence interval functions
lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}
upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

# ReadVCF forloop - rbind stats
out <- foreach(i = 1:length(runs), .combine="rbind", .errorhandling="remove") %dopar% {
  
  # Progress message
  print(paste("Progress:", 
              round((i/length(runs))*100, digits = 3), 
              sep=" "))
  
  # All replicate files per run
  filenames2 <- meta[run==runs[i]]$file
  
  # Load allele frequency information
  dt <- data.table(do.call(rbind, lapply(filenames2, fread)),
                   run.number=i,
                   run=runs[i])
  
  # Add within years
  dt2 <- data.table(dt %>% 
            mutate(d.sim=abs(sim.gen2-sim.gen1)) %>% 
            mutate(dd.sim = case_when(sim.gen1 %in% c(1:16) & sim.gen2 %in% c(1:16) ~ "1.within", 
                   sim.gen1 %in% c(1:33) & sim.gen2 %in% c(17:33) | sim.gen1 %in% c(17:33) & sim.gen2 %in% c(1:33) ~ "2.overwinter", 
                   sim.gen1 %in% c(1:50) & sim.gen2 %in% c(34:50) | sim.gen1 %in% c(34:50) & sim.gen2 %in% c(1:50) ~ "3.multi")))

  # Average across replicates
  dt2 <- data.table(dt2 %>% 
           group_by(nMax, nMin, dd.sim, replicate) %>% 
           summarize(mean.FST = mean(FST, na.rm = T),
                     sd.FST = sd(FST, na.rm = T),
                     n.FST = length(unique(replicate))) %>% 
           mutate(model=paste(nMax, nMin, sep="_")) %>% 
           pivot_wider(names_from = dd.sim,
           values_from=c(mean.FST, sd.FST)))

  
  # Finish
  return(dt2)
}

# Output aggregated file
write.csv(out, file="../data.fst.100reps.aggregated.added.csv", quote = F, row.names = F)

pdf("../FST.pairwise.ribbon.100reps.1k-100k.more.added.pdf", width = 13, height = 8)

# Plot within vs yearly FST
ggplot(out[!c(nMax==50)][!c(nMax==100)],
         aes(x=as.numeric(nMin),
             y=as.numeric(mean.FST),
             ymin=as.numeric(lci.FST),
             ymax=as.numeric(uci.FST),
             color=dd.sim)) +
    geom_ribbon(alpha=0.5, 
                aes(fill=dd.sim)) +
    geom_line() +
    geom_hline(yintercept = 0.0, linetype=2, size=1.1) +
    scale_x_log10() +
    facet_wrap(~nMax) +
    labs(x="Minimum", 
         y="Mean FST", 
         fill="Year") +
    theme_classic() +
    guides(color=FALSE) +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold", size=12),
          legend.title = element_text(face="bold", size=15),
          legend.background = element_blank(),
          strip.text =element_text(face="bold", size=15),
          axis.text.x = element_text(face="bold", size=18, 
                                     angle = -30, vjust = 0.6),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20))

dev.off()

pdf("../FST.pairwise.raster.100reps.1k-100k.more.pdf", width = 13, height = 8)

# Plot within vs yearly FST Raster plot
ggplot(out[!nMax %in% c(50,100)],
       aes(x=as.factor(nMin),
           y=as.factor(nMax),
           fill=mean.FST)) +
  geom_tile() +
  facet_wrap(~dd.sim, scales="free") +
  labs(x="Minimum population size", 
       y="Maximum population size", 
       fill="Mean FST") +
  theme_classic() +
  scale_fill_viridis(option = "magma", begin=0, end=1) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=18, 
                                   angle = -30, vjust = 0.6),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()

# All frequency data
filenames3 <- list.files(pattern = "pca")

# Extract metadata from filenames
meta <- data.table(data.table(seed=tstrsplit(filenames3, ".", fixed=T)[[8]],
                              nMax = tstrsplit(filenames3, ".", fixed=T)[[5]],
                              nMin = tstrsplit(filenames3, ".", fixed=T)[[6]],
                              replicate = tstrsplit(filenames3, ".", fixed=T)[[7]],
                              file=filenames3) %>% 
                     group_by(nMax, nMin) %>% 
                     mutate(run_number=row_number(),
                            run=paste(nMax, nMin,sep="_")))

# Distinct models
runs <- unique(meta$run)

# ReadVCF forloop - calculate diversity statistics.
out.pca <- foreach(i = 1:length(runs), .combine="rbind", .errorhandling="remove") %dopar% {
  
  # Progress message
  print(paste("Progress:", 
              round((i/length(runs))*100, digits = 3), 
              sep=" "))
  
  # All replicate files per run
  filenames4 <- meta[run==runs[i]]$file
  
  # Load allele frequency information
  dt <- data.table(do.call(rbind, lapply(filenames4, fread)),
                   run.number=i,
                   run=runs[i]) %>% 
    mutate(nMax=tstrsplit(Sample, "_")[[7]],
           nMin=tstrsplit(Sample, "_")[[8]],
           replicate=tstrsplit(Sample, "_")[[9]],
           seed=tstrsplit(Sample, "_")[[10]])
  
  # Aggregate simulation mean FST and PCA 
  dt.r2 <- data.table(dt %>% 
            mutate(dd.sim = case_when(Year %in% "Year1" ~ "1.within",
                                      Year %in% "Year2" ~ "2.overwinter",
                                      Year %in% "Year3" ~ "3.multi")) %>% 
            filter(!dd.sim == "3.multi") %>%
            group_by(replicate) %>% 
            mutate(r2pc1=summary(lm(Dim.1 ~ Year))$r.squared,
                   r2pc2=summary(lm(Dim.2 ~ Year))$r.squared,
                   r2pc3=summary(lm(Dim.3 ~ Year))$r.squared)) 
  
  # Aggregate further
  dt <- data.table(dt.r2 %>% 
            group_by(nMax, nMin, dd.sim, replicate) %>% 
            summarize(mean.pc1 = mean(Dim.1, na.rm = TRUE),
                      mean.pc2 = mean(Dim.2, na.rm = TRUE), 
                      mean.pc3 = mean(Dim.3, na.rm = TRUE),
                      r2pc1 = mean(r2pc1, na.rm = TRUE),
                      r2pc2 = mean(r2pc2, na.rm = TRUE),
                      r2pc3 = mean(r2pc3, na.rm = TRUE)) %>% 
            mutate(model = paste(nMax, nMin, sep="_"),
                   dd.sim = dd.sim) %>% 
            pivot_wider(names_from = dd.sim,
                        values_from=c(mean.pc1, mean.pc2, mean.pc3,
                                      r2pc1, r2pc2, r2pc3)))
  
  # Finish
  return(dt)
}

# Output pca
#write.csv(out.pca, file="../data.pca.100reps.aggregated.more.csv", quote = F, row.names = F)

pdf("../pca.raster.100reps.1k-100k.more.pdf", width = 15, height = 10)

# Plot within vs yearly FST Raster plot
out.pca %>% 
  pivot_longer(cols = colnames(out.pca)[11:12]) %>% 
  ggplot(.,
       aes(x=reorder(as.factor(nMin), as.numeric(nMin)),
           y=reorder(as.factor(nMax), as.numeric(nMax)),
           fill=as.numeric(value))) +
  geom_raster() +
  facet_wrap(~name) +
  labs(x="Minimum population size", 
       y="Maximum population size", 
       fill="R2 PC1 ~ Year") +
  theme_classic() +
  scale_fill_viridis(option = "magma", begin=0, end=1) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=18, 
                                   angle=-30, vjust=0.6),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()
