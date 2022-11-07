library(phyloseq)
library(dplyr)
library(readr)
library(microbiome)
library(ggplot2)
library(reshape2)
# The function prediction was done using FUNGuild.
# python ~/build/FUNGuild/FUNGuild.py guild -taxa fungi_ps_fil5_TMM_taxa.txt
guilds <- as.data.frame(read_delim("fungi_ps_fil5_TMM_taxa.guilds.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE))
guilds$trophicMode <- as.factor(guilds$trophicMode)
summary(guilds$trophicMode)
guilds$Phylum <- as.factor(guilds$Phylum)
summary(guilds$Phylum)
rownames(guilds) <- guilds$OTU
guilds <- guilds %>% dplyr::select(-c(OTU,citationSource,notes))
guilds$Kingdom <- guilds$trophicMode

saveFunctionTable <- function(ps,df){
  tax_table(ps) <- as.matrix(df)
  ps1 <-microbiome::aggregate_rare(ps, level = "Kingdom", detection = 0.000000001, prevalence = 0.00001)
  p1 <- microbiome::plot_composition(ps1) 
  longdf <- p1$data
  longdf$xlabel <- NULL
  widedf <- dcast(longdf, Tax ~ Sample, value.var="Abundance")
  rownames(widedf) <- widedf$Tax
  widedf$Tax <- NULL
  return(widedf)
}
widedf <- saveFunctionTable(ps_TMM, guilds)
saveRDS(widedf, file = 'trophicMode_matrix.rds')
# Growth form
guilds$Kingdom <- guilds$growthForm
widedf <- saveFunctionTable(ps_TMM, guilds)
saveRDS(widedf, file = 'growthForm_matrix.rds')
# Guild
guilds$Kingdom <- guilds$guild
widedf <- saveFunctionTable(ps_TMM, guilds)
saveRDS(widedf, file = 'guilds_matrix.rds')
