#Written by Cuyler Luck
#Contact: cuyler.luck@ucsf.edu / cuylerluck@gmail.com or ross.okimoto@ucsf.edu

#This analysis is on stranded RNA-seq data from C2C12 cells which have been stably transduced with one of three conditions.
#For each condition there are two clonal cell lines that have been sequenced (once each).
#The clonal lines are designated by the letter and number at the end of the sample name.
#These are the meanings of the sample names:
# C_EVxx = empty vector (pMY-IRES-EGFP)
# C_CD4xx = CIC-DUX4 (pMYs-CICDUX4)
# C_dC1xx = CIC-DUX4 with C1 domain deleted (pMYs-CICDUX4-dC1)

#First load packages and set working directory:
library(data.table) #1.14.8
library(dplyr) #1.1.1
library(edgeR) #3.40.2
library(ggplot2) #3.4.1
library(pheatmap) #1.0.12
library(tidyr) #1.3.0
library(patchwork) #1.1.2
library(ggrepel) #0.9.3
library(gprofiler2) #0.2.3

setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_C2C12_dC1/GRCm39_alignment/") #this will change by user & location of data


#For the NIH/3T3 clones I made a nice plot at this point showing the transgenes were correct, you could do the same here
#I checked my grep results and they all look good.

#Formal analysis

#First, read in all data
#Skip the first four lines, they just give info on odd-mapping statistics
EV_C7 = fread(file="ReadCounts/C_EVC7ReadsPerGene.out.tab", skip = 4)
EV_G4 = fread(file="ReadCounts/C_EVG4ReadsPerGene.out.tab", skip = 4)
CD4_C9 = fread(file="ReadCounts/C_CD4C9ReadsPerGene.out.tab", skip = 4)
CD4_D6 = fread(file="ReadCounts/C_CD4D6ReadsPerGene.out.tab", skip = 4)
dC1_C5 = fread(file="ReadCounts/C_dC1C5ReadsPerGene.out.tab", skip = 4)
dC1_D5 = fread(file="ReadCounts/C_dC1D5ReadsPerGene.out.tab", skip = 4)


#This is stranded data generated using the NEBNext Ultra II Directional Library Prep kit for Illumina
#so, we can take column 4 from STAR as the output (equivalent to HTseq -stranded "reverse")
#see https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/
#see https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

EV_C7_strand = dplyr::select(EV_C7, c(1,4))
EV_G4_strand = dplyr::select(EV_G4, c(1,4))
CD4_C9_strand = dplyr::select(CD4_C9, c(1,4))
CD4_D6_strand = dplyr::select(CD4_D6, c(1,4))
dC1_C5_strand = dplyr::select(dC1_C5, c(1,4))
dC1_D5_strand = dplyr::select(dC1_D5, c(1,4))

#Let's put names on the gene count columns that correspond to the sample so that we can merge these into one big data frame
colnames(EV_C7_strand) = c("Gene","EV_C7")
colnames(EV_G4_strand) = c("Gene","EV_G4")
colnames(CD4_C9_strand) = c("Gene","CD4_C9")
colnames(CD4_D6_strand) = c("Gene","CD4_D6")
colnames(dC1_C5_strand) = c("Gene","dC1_C5")
colnames(dC1_D5_strand) = c("Gene","dC1_D5")


master = inner_join(EV_C7_strand, EV_G4_strand, by = "Gene")
master = inner_join(master, CD4_C9_strand, by = "Gene")
master = inner_join(master, CD4_D6_strand, by = "Gene")
master = inner_join(master, dC1_C5_strand, by = "Gene")
master = inner_join(master, dC1_D5_strand, by = "Gene")


#move gene names to rownames for edgeR syntax
master = tibble::column_to_rownames(master, "Gene")


#We can now pass this master data frame into an edgeR pipeline
#I am choosing to use a GLM for this analysis. So I will use the non-classic pipeline.
groups = c("EV","EV","CICDUX4","CICDUX4","dC1","dC1")
dg = DGEList(counts=master, group = groups)
keep = filterByExpr(dg)
dg = dg[keep, , keep.lib.sizes = FALSE]
dg = calcNormFactors(dg)
#reordering the levels of the group factor so "EV" becomes the reference
dg$samples$group = relevel(dg$samples$group, ref="EV")
#checking the order with
#dg$samples$group
#reveals that the present levels are EV CICDUX4 dC1
design = model.matrix(~dg$samples$group, data = dg$samples)
#checking the design matrix shows that now the intercept (baseline condition) is EV.
#coefficient 2 = CICDUX4
#coefficient 3 = dC1
#this information will be useful in doing pairwise comparisons soon.
dg = estimateDisp(dg, design)
#plotBCV(dg) to see what the BCV looks like if desired


#now we can do pairwise comparisons between groups of interest
#with this design matrix, providing glmQLFTest a number for coef means "compare this coefficient to the baseline [EV]"
#first we get the fit
fit = glmQLFit(dg, design)

#compare CICDUX4 to EV
qlf.CICDUX4.EV = glmQLFTest(fit, coef=2)

#compare dC1 to EV
qlf.dC1.EV = glmQLFTest(fit, coef=3)

#compare dC1 to CICDUX4
qlf.dC1.CICDUX4 = glmQLFTest(fit, contrast=c(0,-1,1)) #i.e. third coefficient with second coefficient

#now we can pull out tables of results from these tests for use in volcano plots etc.
CICDUX4.EV_results = qlf.CICDUX4.EV$table
dC1.EV_results = qlf.dC1.EV$table
dC1.CICDUX4.results = qlf.dC1.CICDUX4$table

#Let's add FDR-adjusted p-values to each of these
CICDUX4.EV_results = dplyr::mutate(CICDUX4.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))
dC1.EV_results = dplyr::mutate(dC1.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))
dC1.CICDUX4.results = dplyr::mutate(dC1.CICDUX4.results, fdr_adj = p.adjust(PValue, method = "fdr"))


#let's also pull out the TMM-normalized log(cpm) values from the edgeR pipeline
#this will be useful later for looking at expression of specific genes
logcpm = as.data.frame(cpm(dg, log = T))


#and we can ask edgeR to plot an MDS plot to see how well samples cluster
pdf(file="plots/C2C12_MDSplot_TMM_samples.pdf", width = 6, height = 4)
plotMDS(dg, labels = c("EV_C7", "EV_G4", "CD4_C9", "CD4_D6", "dC1_C5", "dC1_D5"), top = 500) #top = 500 is the default setting anyways
dev.off()

#starting with some volcano plots
#I'm going to move the gene names back into the dataframes just for ease.
CICDUX4.EV_results_volcano = tibble::rownames_to_column(CICDUX4.EV_results, "Gene")
dC1.EV_results_volcano = tibble::rownames_to_column(dC1.EV_results, "Gene")
dC1.CICDUX4.results_volcano = tibble::rownames_to_column(dC1.CICDUX4.results, "Gene")

#note that I am transforming the fdr-adjusted pvalue by -log10(fdr) in the plotting commands.
#I'm going to designate a log2FC of +/-2 as a significant cutoff, and -log10(fdr) of 1 as significant (i.e. q < 0.1)

#the following blocks of code just generate simple volcano plots to see how many genes meet significance cutoffs for each comparison.

#CICDUX4 vs EV
volcano_CICDUX4vsEV = ggplot(data = CICDUX4.EV_results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "gray", size = 3) +
  geom_point(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$logFC > 2 & (-log10(CICDUX4.EV_results_volcano$fdr_adj)) > 1,], 
             color = "#fc558f", size = 3) +
  geom_point(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$logFC < -2 & (-log10(CICDUX4.EV_results_volcano$fdr_adj)) > 1,], 
             color = "#e3ac8d", size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$Gene %in% c("Vgf", "Etv4"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, nudge_x = 3) +
  geom_text_repel(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$Gene %in% c("Etv1", "Etv5"),], 
                  aes(label = Gene), nudge_y = 0.7, box.padding = 0.5, nudge_x = -3) +
  theme_bw() +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-20,20), ylim = c(0,2)) +
  ggtitle("CICDUX4 vs EV")



#dC1 vs EV
volcano_dC1vsEV = ggplot(data = dC1.EV_results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "gray", size = 3) +
  geom_point(data = dC1.EV_results_volcano[dC1.EV_results_volcano$logFC > 2 & (-log10(dC1.EV_results_volcano$fdr_adj)) > 1,], 
             color = "#006fc3", size = 3) +
  geom_point(data = dC1.EV_results_volcano[dC1.EV_results_volcano$logFC < -2 & (-log10(dC1.EV_results_volcano$fdr_adj)) > 1,], 
             color = "#e3ac8d", size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% c("Vgf", "Etv4"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, nudge_x = 3) +
  geom_text_repel(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% c("Etv1", "Etv5"),], 
                  aes(label = Gene), nudge_y = 0.7, box.padding = 0.5, nudge_x = -3) +
  theme_bw() +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-20,20), ylim = c(0,2)) +
  ggtitle("dC1 vs EV")



#Plot these two together using patchwork
patchwork_volcano = volcano_CICDUX4vsEV / volcano_dC1vsEV
pdf(file = "plots/C2C12_volcano_both_vs_EV.pdf", width = 8, height = 12)
patchwork_volcano
dev.off()





#I'm also interested in seeing where all of the genes turned on (log2FC > 2, -log10(q) > 1) in C2C12 CIC-DUX4 go in the dC1 comparison
volcano_CICDUX4_EV_all_on_in_C2C12 = ggplot(data = CICDUX4.EV_results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "gray", size = 3) +
  geom_point(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$logFC > 2 & -log10(CICDUX4.EV_results_volcano$fdr_adj) > 1,], 
             color = "#ae1604", size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$Gene %in% c("Vgf", "Etv4"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, nudge_x = 3, size = 7) +
  geom_text_repel(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$Gene %in% c("Etv1", "Etv5"),], 
                  aes(label = Gene), nudge_y = 0.7, box.padding = 0.5, nudge_x = -3, size = 7) +
  theme_bw() +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-20,20), ylim = c(0,2)) +
  ggtitle("CICDUX4 vs EV, All C2C12 CIC-DUX4 Genes Turned on in Maroon") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 17))


#define what those genes are, which are turned on
all_C2C12_CICDUX4_on = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$logFC > 2 & -log10(CICDUX4.EV_results_volcano$fdr_adj) > 1,]$Gene

#same genes from full length CIC-DUX4 vs EV, but plotted for dC1 comparison
volcano_dC1_EV_all_on_in_C2C12 = ggplot(data = dC1.EV_results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "gray", size = 3) +
  geom_point(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% all_C2C12_CICDUX4_on,], 
             color = "#ae1604", size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% c("Vgf","Etv4"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, size = 7) +
  geom_text_repel(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% c("Etv1", "Etv5"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, nudge_x = -3, size = 7) +
  theme_bw() +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-20,20), ylim = c(0,2)) +
  ggtitle("dC1 vs EV, All C2C12 Genes Turned on in Maroon") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 17))

patchwork_volcano_all_C2C12_on_genes = volcano_CICDUX4_EV_all_on_in_C2C12 / volcano_dC1_EV_all_on_in_C2C12
pdf(file = "plots/volcano_all_C2C12_on_genes.pdf", width = 8, height = 12)
patchwork_volcano_all_C2C12_on_genes
dev.off()


#Heatmaps will be helpful in looking at this data too.

#I will make annotations for the samples & colors for the samples to help with visual clarity

annotation_col = data.frame(
  Condition = c(rep("EV",2), rep("CICDUX4",2), rep("dC1",2))
)
rownames(annotation_col) = colnames(logcpm)

annotation_colors = list(
  Condition = c(EV = "gray", CICDUX4 = "#fc558f", dC1 = "#006fc3")
)

#To start, the simplest thing is to make a heatmap looking at all genes meeting the criteria for significantly ON in either
#the CICDUX4 vs EV or the dC1 vs EV comparisons.

#Recall that all_C2C12_CICDUX4_on holds all genes turned on in the CIC-DUX4 vs EV comparison

#I want to check to make sure I'm not missing any genes which are significantly turned on in dC1 but not CICDUX4
on_dC1 = dC1.EV_results_volcano[dC1.EV_results_volcano$logFC > 2 & -log10(dC1.EV_results_volcano$fdr_adj) > 1,]$Gene
#77 genes. are there any which are not in the list of genes turned on in the CICDUX4 vs EV comparison?
on_dC1_not_CICDUX4 = on_dC1[!on_dC1 %in% all_C2C12_CICDUX4_on]
#Actually yes, 31 genes. That's the same number as were in the NIH/3T3 case. Weird

#Now we can combine the two lists to get all genes turned on in at least one of the comparisons
genes_on_either = c(all_C2C12_CICDUX4_on, on_dC1_not_CICDUX4)
#509 genes

#Now let's look at a heatmap of this
logcpm_heatmap_genes_on_either = logcpm[rownames(logcpm) %in% genes_on_either,]
#Plotting a heatmap with scaling within rows turned on
pdf(file = "plots/C2C12_heatmap_genes_on_either.pdf", width = 4, height = 6)
heatmap_genes_on_either = pheatmap(logcpm_heatmap_genes_on_either, 
         scale = "row", 
         show_rownames = F,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_cols = F)
dev.off()



#Another series of simple plots that will be useful is looking at log(cpm) expression of certain targets.
#To plot this data, I need to pivot it to longer first.

logcpm = tibble::rownames_to_column(logcpm, "Gene")
logcpm_pivot = pivot_longer(logcpm, names_to = "sample", cols = c(2:7), values_to = "logcpm")

#change factor order for sample
logcpm_pivot$sample = factor(logcpm_pivot$sample, levels = c("EV_C7", "EV_G4", "CD4_C9", "CD4_D6", "dC1_C5", "dC1_D5"))

#add a variable about what kind of sample each is
logcpm_pivot = dplyr::mutate(logcpm_pivot, Condition = ifelse(sample %in% c("EV_C7","EV_G4"), "EV", ifelse(sample %in% c("CD4_C9","CD4_D6"),"CICDUX4","dC1")))

#set the order of the condition variable
logcpm_pivot$Condition = factor(logcpm_pivot$Condition, levels = c("EV","CICDUX4","dC1"))

#This function will make a nice logcpm plot for any gene you pass into it.
make_logcpm_plot = function(name){
  gene = name
  ggplot(logcpm_pivot[logcpm_pivot$Gene %in% paste(gene),], aes(x = Condition, y = logcpm, fill = Condition)) + 
    stat_summary(fun = mean, geom = "bar") +
    geom_point() +
    theme_classic() +
    ggtitle(paste(gene)) +
    ylab("log2(cpm)") +
    coord_cartesian(ylim = c(-5,17)) + #this is based on the min/max of all logcpm values in the pivot dataframe
    theme(axis.title = element_text(size = 14), 
          axis.text = element_text(size = 12, hjust = 1, angle = 30),
          legend.position = "none",
          plot.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = c("gray","#fc558f","#006fc3"))
}

#Let's make a whole bunch of these for what are known in the literature to likely be true targets
#first pull out the list of "true targets" and sort alphabetically
true_targets = c("Etv1", "Etv4", "Etv5", "Vgf", "Ccnd2", "Dusp6")
true_targets = sort(true_targets)

#then make a list containing all the plot objects
logcpm_true_targets = lapply(true_targets, make_logcpm_plot)


#and arrange the plots to be plotted nicely using patchwork
patchwork_logcpm_true_targets=wrap_plots(logcpm_true_targets[1:length(true_targets)]) + plot_layout(ncol=5)

pdf(file="plots/C2C12_patchwork_logcpm_true_targets.pdf", width = 10, height = 5)
patchwork_logcpm_true_targets
dev.off()

#write some datasets to csv files to share with people.
#logcpm
write.table(logcpm, file = "outputs/C2C12_log2_of_cpm.txt", row.names = F, sep = "\t", quote = F)

#raw counts, stranded, combined
#put gene names back on master as a variable to simplify writing
master_for_writing = tibble::rownames_to_column(master,"gene")
write.table(master_for_writing, file = "outputs/counts_C2C12_dC1.txt", row.names = F, sep = "\t", quote = F)

#CIC::DUX4 vs EV edgeR output, with gene labels
write.table(CICDUX4.EV_results_volcano, file = "outputs/C2C12_CD4_vs_EV_DE.txt", row.names = F, sep = "\t", quote = F)

#dC1 vs EV edgeR output, with gene labels
write.table(dC1.EV_results_volcano, file = "outputs/C2C12_dC1_vs_EV_DE.txt", row.names = F, sep = "\t", quote = F)

#dC1 vs CIC::DUX4 edgeR output, with gene labels
write.table(dC1.CICDUX4.results_volcano, file = "outputs/C2C12_dC1_vs_CICDUX4_DE.txt", row.names = F, sep = "\t", quote = F)



#Gene ontology / GSEA to help determine if there are any interesting pathways 
#in the different blocks of genes regulated by full length vs dC1 CIC::DUX4

#I'm manually getting the genes in the four distinct blocks of the heatmap
#doing this with pheatmap settings, not all shown here

head(logcpm_heatmap_genes_on_either) #original order going into the heatmap
head(logcpm_heatmap_genes_on_either[heatmap_genes_on_either$tree_row$order,],25) #adjusted order after clustering

ordered_logcpm_heatmap_genes_on_either = logcpm_heatmap_genes_on_either[heatmap_genes_on_either$tree_row$order,]

#Block 1 (genes up in dC1, not in CD4) is from Synm to Blvrb
which(rownames(ordered_logcpm_heatmap_genes_on_either) == "Synm")
which(rownames(ordered_logcpm_heatmap_genes_on_either) == "Blvrb")
#Rows 1 thru 23
block1 = ordered_logcpm_heatmap_genes_on_either[1:23,]

#Block 2 (genes up in CD4, but off in dC1) is from Syt7 to Fosl1
which(rownames(ordered_logcpm_heatmap_genes_on_either) == "Syt7")
which(rownames(ordered_logcpm_heatmap_genes_on_either) == "Fosl1")
#Rows 24 thru 174
block2 = ordered_logcpm_heatmap_genes_on_either[24:174,]

#Block 3 (genes up in both CD4 and dC1 vs EV, but potentially to differing degrees) is from Meox1 to Tmprss5 [the end]
which(rownames(ordered_logcpm_heatmap_genes_on_either) == "Meox1")
which(rownames(ordered_logcpm_heatmap_genes_on_either) == "Tmprss5")
#Rows 175 thru 509
block3 = ordered_logcpm_heatmap_genes_on_either[175:509,]


#Now run gene list functional enrichment analysis, I can do all blocks at once with this nice package

gost_blocks = gost(list(block1 = rownames(block1),
                        block2 = rownames(block2),
                        block3 = rownames(block3)), 
                   organism = "mmusculus", 
                   correction_method = "g_SCS", 
                   evcodes = F)

#split out the results by block
gost_results = gost_blocks$result
gost_block1 = gost_results[gost_results$query == "block1",]
gost_block2 = gost_results[gost_results$query == "block2",]
gost_block3 = gost_results[gost_results$query == "block3",]

#sort descending by adjusted p-value (it is adjusted by the function)
gost_block1 = gost_block1[order(gost_block1$p_value, decreasing = F),]
gost_block2 = gost_block2[order(gost_block2$p_value, decreasing = F),]
gost_block3 = gost_block3[order(gost_block3$p_value, decreasing = F),]

#trimming down a bit for file writing (including a column of lists causes an error in write.table)
gost_block1 = dplyr::select(gost_block1, -parents)
gost_block2 = dplyr::select(gost_block2, -parents)
gost_block3 = dplyr::select(gost_block3, -parents)


#write these out 
write.table(gost_block1, file = "outputs/C2C12_block1_geneOntology.txt", row.names = F, sep = "\t", quote = F)
write.table(gost_block2, file = "outputs/C2C12_block2_geneOntology.txt", row.names = F, sep = "\t", quote = F)
write.table(gost_block3, file = "outputs/C2C12_block3_geneOntology.txt", row.names = F, sep = "\t", quote = F)


#also writing out the genes in blocks 1, 2, and 3, and the aggregate
write.table(rownames(block1), file = "outputs/C2C12_block1.txt", row.names = F, sep = "\t", quote = F, col.names = c("Gene"))
write.table(rownames(block2), file = "outputs/C2C12_block2.txt", row.names = F, sep = "\t", quote = F, col.names = c("Gene"))
write.table(rownames(block3), file = "outputs/C2C12_block3.txt", row.names = F, sep = "\t", quote = F, col.names = c("Gene"))
write.table(rownames(ordered_logcpm_heatmap_genes_on_either), file = "outputs/NIH3T3_genes_on_either.txt", row.names = F, sep = "\t", quote = F, col.names = c("Gene"))



#Some analysis of which genes are consistently regulated across C2C12 or NIH/3T3 clones

C2C12_block1 = rownames(block1)
C2C12_block2 = rownames(block2)
C2C12_block3 = rownames(block3)
C2C12_genes_on_either = rownames(ordered_logcpm_heatmap_genes_on_either)

#Read in the NIH/3T3 genes
setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_NIH3T3_dC1/GRCm39_alignment/outputs")

NIH3T3_block1 = fread("NIH3T3_block1.txt")$Gene
NIH3T3_block2 = fread("NIH3T3_block2.txt")$Gene
NIH3T3_block3 = fread("NIH3T3_block3.txt")$Gene
NIH3T3_genes_on_either = fread("NIH3T3_genes_on_either.txt")$Gene

setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_C2C12_dC1/GRCm39_alignment")
#Which genes are on (for either CIC::DUX4 or dC1, vs EV) for both cell lines?

shared_genes_on_either = intersect(C2C12_genes_on_either, NIH3T3_genes_on_either)
shared_genes_on_either = shared_genes_on_either[order(shared_genes_on_either)]
write.table(shared_genes_on_either, file = "outputs/shared_genes_on_either.txt", row.names = F, sep = "\t", quote = F, col.names = c("Gene"))

  
shared_block1 = intersect(C2C12_block1, NIH3T3_block1)
shared_block1 = shared_block1[order(shared_block1)]
write.table(shared_block1, file = "outputs/shared_block1.txt", row.names = F, sep = "\t", quote = F, col.names = c("Gene"))


shared_block2 = intersect(C2C12_block2, NIH3T3_block2)
shared_block2 = shared_block2[order(shared_block2)]
write.table(shared_block2, file = "outputs/shared_block2.txt", row.names = F, sep = "\t", quote = F, col.names = c("Gene"))


shared_block3 = intersect(C2C12_block3, NIH3T3_block3)
shared_block3 = shared_block3[order(shared_block3)]
write.table(shared_block3, file = "outputs/shared_block3.txt", row.names = F, sep = "\t", quote = F, col.names = c("Gene"))


#I was also asked to check p21 (Cdkn1a) expression in these clones, and I'll check other differentiation-involved genes at the same time:
#first make the list of genes and sort alphabetically
diff_genes = c("Cdkn1a", "Myog", "Myod1", "Myf5")
diff_genes = sort(diff_genes)

#then make a list containing all the plot objects
logcpm_diff_genes = lapply(diff_genes, make_logcpm_plot)


#and arrange the plots to be plotted nicely using patchwork
patchwork_logcpm_diff_genes=wrap_plots(logcpm_diff_genes[1:length(diff_genes)]) + plot_layout(ncol=5)

pdf(file="plots/C2C12_patchwork_logcpm_diff_genes.pdf", width = 10, height = 3)
patchwork_logcpm_diff_genes
dev.off()

