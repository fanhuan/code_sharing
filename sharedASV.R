# This R script describes how the # of ASVs shared by the two orchid species were calculated.
# For any questions, please consult Huan Fan at fh@xtbg.org.cn
library(phyloseq)
library(phyloseqCompanion)
library(dplyr)
ps <- readRDS('fungi_ps_fil5_TMM.rds')
df <- as.data.frame(otu_table(ps))
samdf <- sample.data.frame(ps)
# we work on one species and repeat it on the other.
pd <- samdf %>% filter(Orchid == 'PD')
df.pd <- df[,pd$Sample]
# brutally, let's just look for ASV that are present in all, which is 
# counting the number of 0s in each row.
df.pd$sumPD <- rowSums(df.pd>0)
# Now PH
ph <- samdf %>% filter(Orchid == 'PH')
df.ph <- df[,ph$Sample]
df.ph$sumPH <- rowSums(df.ph>0)
# Number of ASVs 
df.all <- cbind(df.ph, df.pd)
df.all.1 <- df.all %>% mutate(type = ifelse(sumPH > 0, 
                                            ifelse(sumPD > 0, 'Shared', 'PH_only'), 'PD_only'))
df.all.1$type <- as.factor(df.all.1$type)
summary(df.all.1$type)