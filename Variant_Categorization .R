library(dplyr)

#### Grouping data in 14 distinct groups ####
# CADD score > 30
D25007_MAF0.01_CADD30 <- subset(D25007_MAF0.01, round(as.numeric(CADD_phred), 2) >= 30)
# CADD score > 15
D25007_MAF0.01_CADD15 <- subset(D25007_MAF0.01, round(as.numeric(CADD_phred), 2) >= 15)
# missense variants with CADD > 15 
D25007_MAF0.01_CADD15_missense <- subset(D25007_MAF0.01, round(as.numeric(CADD_phred), 2) >= 15 & ExonicFunc_refgene == "nonsynonymous SNV" )
# all missense variants
D25007_MAF0.01_missense <- subset(D25007_MAF0.01, ExonicFunc_refgene == "nonsynonymous SNV")
# synonymous variants with CADD >15
D25007_MAF0.01_synonymous_CADD15  <- subset(D25007_MAF0.01, ExonicFunc_refgene == "synonymous SNV" & round(as.numeric(CADD_phred), 2) >= 15)
# all synonymous variants
D25007_MAF0.01_synonymous  <- subset(D25007_MAF0.01, ExonicFunc_refgene == "synonymous SNV")
# all exonic variants
D25007_MAF0.01_exonic  <- subset(D25007_MAF0.01, Func_refgene == "exonic")
# all variants in UTR regions
D25007_MAF0.01_UTR  <- subset(D25007_MAF0.01, Func_refgene == "UTR3"|Func_refgene == "UTR5")
# all intronic variants
D25007_MAF0.01_intronic  <- subset(D25007_MAF0.01, Func_refgene == "intronic")
# all variants in promoter regions, that is 2kb both upstream and downstream of transcriptional start site
D25007_MAF0.01_promoter  <- subset(D25007_MAF0.01, Func_refgene == "upstream"| Func_refgene == "downstream")
# variants in double-elite enhancers (need GeneHancer Data or in MViewer)

# variants in conserved regions





# Extract the "Gene_refgene" column and output as txt file
gene_column <- D25007_MAF0.01_CADD20$Gene_refgene

# Write to a text file
write.table(gene_column, file = "D25007_MAF0.01_CADD20_gene.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)




