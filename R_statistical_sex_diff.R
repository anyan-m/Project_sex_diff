getwd()
setwd('/Users/aniamarczynska/Project_sex_diff')
#install.packages("BiocManager")
#install.packages('tidyverse')
#BiocManager::install("edgeR")
#install.packages('vioplot')
#install.packages('GOplot')
#install.packages('ggVennDiagram')
install.packages('ggvenn')
install.packages('grid')
install.packages('gt')
BiocManager::install('clusterProfiler')
BiocManager::install("enrichplot")
library(grid)
library(ggvenn)
library(enrichplot)
library(clusterProfiler)
library(R.utils)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(vioplot)
library(GOplot)
library(ggVennDiagram)
library(gt)
library(DOSE)
stomach_female_15 <- read.csv('/Users/aniamarczynska/Project_sex_diff/at_least_15%_expression/stomach_female_122_85_noncoding.csv', header = FALSE, sep=',' )
stomach_male_15 <- read.csv('/Users/aniamarczynska/Project_sex_diff/at_least_15%_expression/stomach_male_122_85_noncoding.csv', header = FALSE, sep=',' )
colnames_stomach <- read.csv('/Users/aniamarczynska/Project_sex_diff/stomach_85_colnames.csv')
stomach_data_all <- read.csv('/Users/aniamarczynska/Project_sex_diff/stomach_final_genes.csv', header = FALSE, sep=',' )

thyroid_female_15 <- read.csv('/Users/aniamarczynska/Project_sex_diff/at_least_15%_expression/thyroid_female_122_85_noncoding.csv', header = FALSE, sep=',' )
thyroid_male_15 <- read.csv('/Users/aniamarczynska/Project_sex_diff/at_least_15%_expression/thyroid_male_122_85_noncoding.csv', header = FALSE, sep=',' )
colnames_thyroid <- read.csv('/Users/aniamarczynska/Project_sex_diff/thyroid_85_colnames.csv')
thyroid_data_all <- read.csv('/Users/aniamarczynska/Project_sex_diff/thyroid_final_genes.csv', header = FALSE, sep=',' )

nervous_female_15 <- read.csv('/Users/aniamarczynska/Project_sex_diff/at_least_15%_expression/nervous_female_122_85_noncoding.csv', header = FALSE, sep=',' )
nervous_male_15 <- read.csv('/Users/aniamarczynska/Project_sex_diff/at_least_15%_expression/nervous_male_122_85_noncoding.csv', header = FALSE, sep=',' )
colnames_nervous <- read.csv('/Users/aniamarczynska/Project_sex_diff/nervous_85_colnames.csv')
nervous_data_all <- read.csv('/Users/aniamarczynska/Project_sex_diff/nervous_final_genes.csv', header = FALSE, sep=',' )
save.image('~/Project_sex_diff/R_sex_diff.RData')

#transpose data frame & redefine row and column names
transpose_df_f <- function(df){
  df <- as.data.frame(df)
  df_t <- t(df)
  rownames(df_t) <- colnames(df)
  colnames(df_t) <- rownames(df)
  return(df_t)}
stomach_female_t <- transpose_df_f(stomach_female_15)
stomach_male_t <- transpose_df_f(stomach_male_15)

thyroid_female_t <- transpose_df_f(thyroid_female_15)
thyroid_male_t <- transpose_df_f(thyroid_male_15)

nervous_female_t <- transpose_df_f(nervous_female_15)
nervous_male_t <- transpose_df_f(nervous_male_15)

#Make them compatible with DESeq2 & voomlima
remove_unused_col_add_header <- function(df_t){
  colnames(df_t) <- df_t[1, ]
  df_t <- df_t[-1,]
  #df_t <- subset(df_t, select = -c(Description,Gene_type) )
  return(df_t)
}
stomach_female_t <- remove_unused_col_add_header(stomach_female_t)
stomach_male_t <- remove_unused_col_add_header(stomach_male_t)
stomach_data_all <- as.data.frame(cbind(stomach_male_t, stomach_female_t))
stomach_data_all <- remove_unused_col_add_header(stomach_data_all)

thyroid_female_t <- remove_unused_col_add_header(thyroid_female_t)
thyroid_male_t <- remove_unused_col_add_header(thyroid_male_t)
thyroid_data_all <- as.data.frame(cbind(thyroid_male_t, thyroid_female_t))
thyroid_data_all <- remove_unused_col_add_header(thyroid_data_all)

nervous_female_t <- remove_unused_col_add_header(nervous_female_t)
nervous_male_t <- remove_unused_col_add_header(nervous_male_t)
nervous_data_all <- as.data.frame(cbind(nervous_male_t, nervous_female_t))
nervous_data_all <- remove_unused_col_add_header(nervous_data_all)


#Check the 'middle' column with gene names
nervous_data_all <- nervous_data_all[-c(124)]

header_sample_name <- function(df_colnames){
  rownames(df_colnames) <- df_colnames[,1]
  df_colnames <- df_colnames[-c(1)]
  return(df_colnames)}

colnames_stomach <- header_sample_name(colnames_stomach)
stomach_data_all <- header_sample_name(stomach_data_all)
colnames_thyroid <- header_sample_name(colnames_thyroid)
thyroid_data_all <- header_sample_name(thyroid_data_all)
colnames_nervous <- header_sample_name(colnames_nervous)
nervous_data_all <- header_sample_name(nervous_data_all)

#Removing the 'Sex' row
sex_remove <- c('Sex')
stomach_data_all <- stomach_data_all[!(rownames(stomach_data_all) %in% sex_remove), ]
thyroid_data_all <- thyroid_data_all[!(rownames(thyroid_data_all) %in% sex_remove), ]
nervous_data_all <- nervous_data_all[!(rownames(nervous_data_all) %in% sex_remove), ]

sapply(thyroid_data_all, class)
stomach_data_all <- mutate_all(stomach_data_all, function(x) as.numeric(x))
thyroid_data_all <- mutate_all(thyroid_data_all, function(x) as.numeric(x))
nervous_data_all <- mutate_all(nervous_data_all, function(x) as.numeric(x))
identical(all, stomach_data_all)


#check if columns and data.type are ok
sapply(nervous_data_all, class)
all(colnames(thyroid_data_all) == rownames(colnames_thyroid))

#Making a dds object
dds_stomach <- DESeqDataSetFromMatrix(countData = round(stomach_data_all), 
                                      colData = colnames_stomach, 
                                      design = ~ Sex)
dds_stomach <- DESeq(dds_stomach)
results_dds_stomach <- results(dds_stomach)

#Try plotting
rld <- vst(dds_stomach, blind=TRUE)
plotPCA(rld, intgroup='Sex')

#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#voom-limma:
d0_stomach <- DGEList()
d0_stomach <- calcNormFactors(d0_stomach)
dim(d0_stomach) # number of genes left
group_sto <- interaction(colnames_stomach)#interaction(rownames(colnames_stomach),colnames_stomach[,1])
plotMDS(d0_stomach, col = as.numeric(group_sto))

#mm <- model.matrix(~group) #with intercept
mm <- model.matrix(~0 + group_sto)#instead of (~0 + group)
y <- voom(d0_stomach, mm, plot = T) #Warning messages:
#1: Partial NA coefficients for 21841 probe(s) 
#2: In voom(d0, mm, plot = T) :
  #The experimental design has no replication. Setting weights to 1.

#ok, but I didn't want to filter them, so I'll use d0 - THIS IS EXTRA
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
plotMDS(d, col = as.numeric(group))
y <- voom(d, mm, plot = T) #same warning as above just 18489 probe(s) 

#lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(y, mm)
head(coef(fit))
coef(fit)["AT1G60030", "pH"]
#Comparisons between groups (log fold-changes) are obtained as contrasts 
#of these fitted linear models:
#Specify which groups to compare:
contr <- makeContrasts(group_stomale - group_stofemale, levels = colnames(coef(fit)))
contr

#estimate contrast for eachs group:
tmp <- contrasts.fit(fit, contr)

#Empirical Bayes smoothing of standard errors (shrinks standard errors that 
#are much larger or smaller than those from other genes towards the 
#average standard error)
tmp <- eBayes(tmp)

#What genes are most differentially expressed?
#Has higher expression in males than females (logFC is positive)
#Has lower expression in males than females (logFC is negative)
top.table.sto <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table.sto$Gene, 53)
length(which(top.table.thy$adj.P.Val < 0.05)) #53 genes are DE in stomach
#1896 in thyroid and 135 in nervous tibial

#Save to txt file
top.table.sto$Gene <- rownames(top.table.sto)
top.table.sto <- top.table.sto[,c("Gene", names(top.table.sto)[1:6])]
write.table(top.table.sto, file = "male_v_female_stomach.txt", row.names = F, sep = "\t", quote = F)

#Do a violin plot:
#selected gene is ENSG00000130021.13	PUDP	protein_coding (no zero values)
stomach_data_all['ENSG00000068078.17',]
test_df <- as.data.frame(stomach_data_all['ENSG00000068078.17',])
sex_row <- t(colnames_stomach)
stomach_selected_gene_and_sex <- rbind(test_df,sex_row)
stomach_selected_gene_and_sex <- t(stomach_selected_gene_and_sex)
stomach_selected_gene_and_sex <- as.data.frame(stomach_selected_gene_and_sex)
stomach_selected_gene_and_sex$Sex <- as.factor(stomach_selected_gene_and_sex$Sex)
stomach_selected_gene_and_sex$ENSG00000068078.17 <- as.numeric(stomach_selected_gene_and_sex$ENSG00000068078.17)

v <- ggplot(stomach_selected_gene_and_sex, aes(x=Sex, y= ENSG00000198286.9, fill=Sex))
v + geom_violin(trim = F) + 
  labs(title="A", x="Sex", y = "CARD11 Gene Expression Level (TPM)") +
  geom_jitter(height = 0, width = 0.2) + 
  scale_fill_brewer(palette = "Pastel1") + 
  geom_boxplot(width=0.05, fill="white") + theme_minimal(base_size = 18)

thyroid_data_all['ENSG00000068078.17',]
thy_selected_gene_df <- as.data.frame(thyroid_data_all['ENSG00000068078.17',])
sex_row <- t(colnames_thyroid)
thyroid_selected_gene_and_sex <- rbind(thy_selected_gene_df,sex_row)
thyroid_selected_gene_and_sex <- t(thyroid_selected_gene_and_sex)
thyroid_selected_gene_and_sex <- as.data.frame(thyroid_selected_gene_and_sex)
thyroid_selected_gene_and_sex$Sex <- as.factor(thyroid_selected_gene_and_sex$Sex)
thyroid_selected_gene_and_sex$ENSG00000068078.17 <- as.numeric(thyroid_selected_gene_and_sex$ENSG00000068078.17)

v <- ggplot(thyroid_selected_gene_and_sex, aes(x=Sex, y= ENSG00000198286.9, fill=Sex))
v + geom_violin(trim = F) + 
  labs(title="C", x="Sex", y = "CARD11 Gene Expression Level (TPM)") +
  geom_jitter(height = 0, width = 0.2) + 
  scale_fill_brewer(palette = "Pastel1") + 
  geom_boxplot(width=0.05, fill="white") + theme_minimal(base_size = 18)

nervous_data_all['ENSG00000068078.17',]
ner_selected_gene_df <- as.data.frame(nervous_data_all['ENSG00000068078.17',])
sex_row <- t(colnames_nervous)
nervous_selected_gene_and_sex <- rbind(ner_selected_gene_df,sex_row)
nervous_selected_gene_and_sex <- t(nervous_selected_gene_and_sex)
nervous_selected_gene_and_sex <- as.data.frame(nervous_selected_gene_and_sex)
nervous_selected_gene_and_sex$Sex <- as.factor(nervous_selected_gene_and_sex$Sex)
nervous_selected_gene_and_sex$ENSG00000068078.17 <- as.numeric(nervous_selected_gene_and_sex$ENSG00000068078.17)

v <- ggplot(nervous_selected_gene_and_sex, aes(x=Sex, y= ENSG00000068078.17, fill=Sex))
v + geom_violin(trim = F) + 
  labs(title="B", x="Sex", y = "CARD11 Gene Expression Level (TPM)") +
  geom_jitter(height = 0, width = 0.2) + 
  scale_fill_brewer(palette = "Pastel1") + 
  geom_boxplot(width=0.05, fill="white") + theme_minimal(base_size = 18)

stomach_selected_gene_and_sex$Tissue <- "Stomach"
thyroid_selected_gene_and_sex$Tissue <- "Thyroid"
nervous_selected_gene_and_sex$Tissue <- "Nervous tibial"
all_selected_gene_sex_tis_FGFR <- rbind(stomach_selected_gene_and_sex, nervous_selected_gene_and_sex, thyroid_selected_gene_and_sex)
dodge <- position_dodge(width = 0.9)
v <- ggplot(all_selected_gene_sex_tis_FGFR, aes(x=Tissue, y= ENSG00000068078.17, fill=Sex))
v + geom_violin(trim = F, position = dodge) + 
  labs(title="B", x="Tissue", y = "FGFR3 Gene Expression Level (TPM)") +
  #geom_jitter(height = 1, width = 0.05) + 
  scale_fill_brewer(palette = "Pastel1") + 
  geom_boxplot(position = dodge, width=0.1) +
  theme_minimal(base_size = 18)

#Venn diagram for tissues:
#link to venn_diagram: https://rdrr.io/cran/GOplot/man/GOVenn.html
???circ<-circular_dat(EC$david, EC$genelist)

#Selecting terms of interest
#stomach
sto_diff_exp_gene <- as.character(top.table.sto[1:53,'Gene'])

#thyroid
thy_diff_exp_gene <- as.character(top.table.thy[1:1896,'Gene'])

#nervous
ner_diff_exp_gene <- as.character(top.table.ner[1:135,'Gene'])

x<-list(stomach = sto_diff_exp_gene,
        thyroid = thy_diff_exp_gene,
        nervous = ner_diff_exp_gene)

ggVennDiagram(x, label_alpha = 0, 
              set_size = 9,
              label = c("count"),
              label_size = 8,
              edge_size = 0.5,
              category.names = c('A', 'C', 'B'),
              show.legend = T) + 
  scale_color_manual(values = c('stomach' = 'black','thyroid' = 'black', 'nervous' = 'black')) +
  scale_fill_gradient(low = 'green', high = 'red')
ggvenn(
  x, 
  fill_color = c("darkorange", "lightcyan1", "magenta4"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F, text_size = 6
)
#http://sape.inf.usi.ch/quick-reference/ggplot2/colour - colors
venn <- Venn(list(stomach = sto_diff_exp_gene,
                  thyroid = thy_diff_exp_gene,
                  nervous = ner_diff_exp_gene))

overlap_center_venn <- overlap(venn, c(1,2,3)) #all overlap between sets
overlap_center_venn <- sorter_gene_names(overlap_center_venn)
#names(overlap_center_venn)[names(overlap_center_venn) == "sorter_gene_names(overlap_center_venn)"] <- "Gene"

write.table(overlap_center_venn, file = "overlap_venn_for_gesea_decimal.txt", row.names = F, sep = "\t", quote = F)
just_stomach <- discern_overlap(venn, c(1))
just_thyroid <- discern_overlap(venn, c(2))
just_nervous <- discern_overlap(venn, c(3))

sorter_gene_names <- function(chars){
  chars <- gsub("\\..*","",chars, fixed = F)
  return(chars)
}

just_stomach_NO_suffix <- sorter_gene_names(just_stomach)
#names(just_stomach_NO_suffix)[names(just_stomach_NO_suffix) == "sorter_gene_names(just_stomach)"] <- "Gene"
just_thyroid_NO_suffix <- sorter_gene_names(just_thyroid)
#names(just_thyroid_NO_suffix)[names(just_thyroid_NO_suffix) == "sorter_gene_names(just_thyroid)"] <- "Gene"
just_nervous_NO_suffix <- sorter_gene_names(just_nervous)
#names(just_nervous_NO_suffix)[names(just_nervous_NO_suffix) == "sorter_gene_names(just_nervous)"] <- "Gene"

write.table(just_stomach_NO_suffix, file = "just_stomach_for_gesea.txt", row.names = F, sep = "\t", quote = F)
write.table(just_nervous_NO_suffix, file = "just_nervous_for_gesea.txt", row.names = F, sep = "\t", quote = F)
write.table(just_thyroid_NO_suffix, file = "just_thyroid_for_gesea.txt", row.names = F, sep = "\t", quote = F)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
#https://support.bioconductor.org/p/130385/ - enrichGO instead of gseGO
# sort the list in decreasing order (required for clusterProfiler)
input_enrichGO <- function(gene_list){
  genes = sort(gene_list, decreasing = TRUE)
  go_enrich <- enrichGO(gene = genes,
                        OrgDb = organism, 
                        keyType = 'ENSEMBL',
                        readable = T,
                        ont = "ALL",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.1)
  return(go_enrich)
}

sto_genes_enrichGO <- input_enrichGO(just_stomach_NO_suffix)
ner_genes_enrichGO <- input_enrichGO(just_nervous_NO_suffix)
thy_genes_enrichGO <- input_enrichGO(just_thyroid_NO_suffix)
overlap_venn_enrichGO <- input_enrichGO(overlap_center_venn)

#new GlobalEnv
#read files with gene clustering
ner_genesGSEA <- read.table('/Users/aniamarczynska/Project_sex_diff/GSEA/just_nervous_for_gesea.txt', header = F, sep = "\t", dec = ".", skip = 1, col.names = c("Gene"))

GeneCol_tolist <- function(gene_df){
  gene_list <- list(gene_df$Gene)
  return(gene_list)
}
ner_genesGSEA <- GeneCol_tolist(ner_genesGSEA)


#plots
dotplot(overlap_venn_enrichGO)
dotplot(ner_genes_enrichGO)
dotplot(thy_genes_enrichGO)
par(mfrow=c(2,2))

##HEATPLOT
heatplot(overlap_venn_enrichGO, 
         foldChange = NULL, 
         showCategory = 10)
heatplot(ner_genes_enrichGO, 
         foldChange = NULL, 
         showCategory = 10)
heatplot(thy_genes_enrichGO, 
         foldChange = NULL, 
         showCategory = 10)
get_gene_FC <- function(df){
  df <- subset(df,grepl(paste0(overlap_center_venn, collapse = "|"), Gene))
  df <- subset(df, select=c('Gene','logFC'))
  df <- df[order(df[,"Gene"]), , drop = FALSE]
  return(df)
}
gene_logfc_sto <- get_gene_FC(top.table.sto)
colnames(gene_logfc_sto) <- c('Gene', 'logFC_sto')

gene_logfc_ner <- get_gene_FC(top.table.ner)
colnames(gene_logfc_ner) <- c('Gene', 'logFC_ner')

gene_logfc_thy <- get_gene_FC(top.table.thy)
colnames(gene_logfc_thy) <- c('Gene', 'logFC_thy')

names_genes <- read.csv('/Users/aniamarczynska/Project_sex_diff/genes_overlap_names.csv', header = TRUE, sep=',')
overlap_matrix <- cbind(gene_logfc_sto, 
             gene_logfc_ner$logFC_ner, 
             gene_logfc_thy$logFC_thy,
             names_genes$Description)
rownames(overlap_matrix) <- overlap_matrix[,5]
test_df <- subset(overlap_matrix, select = -c(`names_genes$Description`,Gene) )

tissue_names <- data.frame(tissue = c('Stomach', 'Nervous', 'Thyroid'))
colnames(tissue_names) <- "Tissue type"
rownames(tissue_names) <- colnames(test_df)

mycolors2 <- list(tissue_names = brewer.pal(3, "Set1")[1:3]); 
names(mycolors2$tissue_names) <- levels(tissue_names$`Tissue type`)
mycolors = list(`Tissue type`= c(Nervous= "magenta4", Stomach="darkorange", Thyroid = "lightcyan1"))
pheatmap(test_df, 
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(100),
         annotation_col = tissue_names,
         annotation_colors = mycolors,
         scale = 'column',
         show_colnames = FALSE) #yes, good




##NEURON MAP
cnetplot(overlap_venn_enrichGO, 
         node_label="all",
         showCategory = 10,
         color_category = 'salmon2',
         color_gene = 'cyan3') 




barPlot_enrichGO <- function(go_enrich){
  barPlot <- barplot(go_enrich, 
          drop = TRUE, 
          showCategory = 10, 
          title = "D",
          font.size = 8)
  return(barPlot)
}

barPlot_enrichGO(thy_genes_enrichGO)
barPlot_enrichGO(ner_genes_enrichGO)
barPlot_enrichGO(overlap_venn_enrichGO)

#stomach table design
stomachGSEA_data <- read.delim('/Users/aniamarczynska/Project_sex_diff/GSEA/gesea_results_stomach.txt', header = F, sep = "\t", dec = ".")
stomachGSEA_data <- header_sample_name(stomachGSEA_data)

colnames(stomachGSEA_data) <- stomachGSEA_data[1,]
stomachGSEA_data <- stomachGSEA_data[-c(1),]
stomachGSEA_data_t <- as.data.frame(t(stomachGSEA_data))

gt_tbl <- gt(stomachGSEA_data)
gt_tbl




