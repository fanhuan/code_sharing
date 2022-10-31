# Read in the complete dataset
library(phyloseq)
library(dplyr)
library(edgeR)
fungi_ps <- readRDS('fungi_ps.rds')
otu.df <- as.data.frame(otu_table(fungi_ps))
# count the number of zeros in each column
num.zero <- colSums(otu.df==0)
# count the number of singletons
singleton <- otu.df[,names(num.zero[num.zero == 119])] #2622 out of 3148
# get the abundance of those singletons
singles <- colSums(singleton)
# get the singletons with abundance lower than 5 
five <- singles[singles < 5] # 979
# get the ones that not
rest <- otu.df %>% select(!names(five))
# Remove singletons with less than 5 reads.
my_subset <- subset(otu_table(fungi_ps), select = colnames(rest))
ps_fil5 <- merge_phyloseq(otu_table(rest,taxa_are_rows = FALSE), 
                          tax_table(fungi_ps), sample_data(fungi_ps), 
                          phy_tree(fungi_ps))
saveRDS(ps_fil5,file = 'ps_fil5.rds')
otu_fil5 <- as(otu_table(ps_fil5), "matrix")
write.csv(otu_fil5, 'otu_fungi_filter5.csv')
# Normalization: TMM - edgeR transformation (Edwards 2015 PNAS)
otu_fil5 <- as(otu_table(ps_fil5), "matrix")
tax_matrix <- DGEList(counts=t(otu_fil5))
tax_matrix <- edgeR::calcNormFactors(tax_matrix, method = "TMM")
tax_matrix <- cpm(tax_matrix)
abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
ps_TMM <- ps_fil5
otu_table(ps_TMM) <- abundances
# Assign the estimated diversity to sample metadata
write.csv(abundances, 'otu_fungi_fil5_TMM.csv')
saveRDS(ps_TMM, file = 'fungi_ps_fil5_TMM.rds')