#Written by Cuyler Luck
#Contact: cuyler.luck@ucsf.edu / cuylerluck@gmail.com or ross.okimoto@ucsf.edu



#This analysis is on stranded RNA-seq data from NIH/3T3 cells which have been stably transduced with one of three conditions.
#For each condition there are two clonal cell lines that have been sequenced (once each).
#The clonal lines are designated by the letter and number at the end of the sample name.
#These are the meanings of the sample names:
# EVxx = empty vector (pMY-IRES-EGFP)
# CDxx = CIC-DUX4 (pMYs-CICDUX4)
# delxx = CIC-DUX4 with C1 domain deleted (pMYs-CICDUX4-dC1)

#First load packages and set working directory:
library(data.table) #1.14.8
library(dplyr) #1.1.1
library(edgeR) #3.40.2
library(ggplot2) #3.4.1
library(pheatmap) #1.0.12
library(tidyr) #1.3.0
library(patchwork) #1.1.2
library(ggrepel) #0.9.3

setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_NIH3T3_dC1/GRCm39_alignment/") #this will change by user & location of data


#Before we do any formal analysis, I want to validate that the conditions have the correct & expected transgenes in them.
#I searched each fastq file with grep for sequences specific to different key parts of the transgenes...
#In particular, this is copied from my documentation file:
#Looking for three specific sequences:
#  1) 30nt trans-breakpoint sequence for the human CIC-DUX4 sequence used here
#  2) 30nt sequence spanning the 5' side of the CIC C1 domain for the human CIC-DUX4 sequence used here
#  3) 30nt sequence spanning the CIC sequence with the C1 domain deleted in the human CIC-DUX4 dC1 used here

#These are the sequences:
#1) ctcggactctgggggtggaccccaagccgg
#2) ccatactcctccctgcggcgcaccctggac
#3) ccatactcctccctgcaggctgccactccc

#The output of taking the word count from text files containing grep results for each of these 3 sequences across all samples is in
# the file called "grep_results_072423.png"

#I will take these numbers and put them into a dataframe. Combining read1 + read2. 
#Also going to grab the total number of reads in each fastq file (from the Novogene final data report). Combining read1 + read2, as it is given there.

#Order of samples (in rows) going from top to bottom is EVB6, EVE8, CDD4, CDD8, delB5, delC10

grep_reads = data.frame(total = c(69071140, 67315008, 64816948, 64795392, 74275046, 70371934),
                        CICDUX4_breakpoint = c(0+0, 0+0, 61+2249, 16+550, 16+635, 21+880),
                        intact_C1 = c(0+0, 0+0, 8+476, 9+846, 0+0, 0+0),
                        deleted_C1 = c(0+0, 0+0, 0+0, 0+0, 10+1034, 26+1461),
                        sample = c("EVB6", "EVE8", "CDD4", "CDD8", "delB5", "delC10"),
                        transgene = c("EV","EV","CICDUX4","CICDUX4","dC1","dC1"))

#change order of some factor levels as desired for ordering on plots and legends
grep_reads$sample = factor(grep_reads$sample, levels = c("EVB6", "EVE8", "CDD4", "CDD8", "delB5", "delC10"))
grep_reads$transgene = factor(grep_reads$transgene, levels = c("EV","CICDUX4","dC1"))

#plot raw reads
raw = ggplot(data = grep_reads, aes(x = sample, y = total, fill = transgene)) +
  geom_col() +
  xlab("Sample") +
  ylab("Raw Read Counts \n(read1+read2)") +
  theme_classic() +
  scale_fill_manual(values=c("#e3ac8d","#fc558f","#006fc3")) +
  ggtitle("Total Reads") +
  theme(axis.text.x=element_text(size=12, angle = 30, hjust = 1), 
        axis.title=element_text(size=14), 
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16, face = "bold")) 

#plot CIC-DUX4 breakpoint-spanning reads
bp_span = ggplot(data = grep_reads, aes(x = sample, y = CICDUX4_breakpoint, fill = transgene)) +
  geom_col() +
  xlab("Sample") +
  ylab("CIC-DUX4 Breakpoint-\nSpanning Reads (read1+read2)") +
  theme_classic() +
  scale_fill_manual(values=c("#e3ac8d","#fc558f","#006fc3")) +
  ggtitle("CIC-DUX4 Breakpoint") +
  theme(axis.text.x=element_text(size=12, angle = 30, hjust = 1), 
        axis.title=element_text(size=14), 
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16, face = "bold")) 


#plot CIC-DUX4 C1-intact reads
C1_intact = ggplot(data = grep_reads, aes(x = sample, y = intact_C1, fill = transgene)) +
  geom_col() +
  xlab("Sample") +
  ylab("CIC-DUX4 Intact \nC1 Reads (read1+read2)") +
  theme_classic() +
  scale_fill_manual(values=c("#e3ac8d","#fc558f","#006fc3")) +
  ggtitle("Intact C1 Domain") +
  theme(axis.text.x=element_text(size=12, angle = 30, hjust = 1), 
        axis.title=element_text(size=14), 
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16, face = "bold")) 


#plot CIC-DUX4 C1-deleted reads
C1_deleted = ggplot(data = grep_reads, aes(x = sample, y = deleted_C1, fill = transgene)) +
  geom_col() +
  xlab("Sample") +
  ylab("CIC-DUX4 Deleted \nC1 Reads (read1+read2)") +
  theme_classic() +
  scale_fill_manual(values=c("#e3ac8d","#fc558f","#006fc3")) +
  ggtitle("Deleted C1 Domain") +
  theme(axis.text.x=element_text(size=12, angle = 30, hjust = 1), 
        axis.title=element_text(size=14), 
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16, face = "bold")) 

#plot all together
patchwork_reads = (raw | bp_span) / (C1_intact | C1_deleted)

pdf(file = "plots/grep_condition_checking.pdf",width = 9, height = 6)
patchwork_reads + plot_annotation(tag_levels = 'A')
dev.off()


#This looks great, so all conditions are what they should be.
#Now we can continue to the true analysis portion.

#First, read in all data
#Skip the first four lines, they just give info on odd-mapping statistics
EVB6 = fread(file="Reads/EVB6ReadsPerGene.out.tab", skip = 4)
EVE8 = fread(file="Reads/EVE8ReadsPerGene.out.tab", skip = 4)
CDD4 = fread(file="Reads/CDD4ReadsPerGene.out.tab", skip = 4)
CDD8 = fread(file="Reads/CDD8ReadsPerGene.out.tab", skip = 4)
delB5 = fread(file="Reads/delB5ReadsPerGene.out.tab", skip = 4)
delC10 = fread(file="Reads/delC10ReadsPerGene.out.tab", skip = 4)

#This is stranded data generated using the NEBNext Ultra II Directional Library Prep kit for Illumina
#so, we can take column 4 from STAR as the output (equivalent to HTseq -stranded "reverse")
#see https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/
#see https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

EVB6_strand = dplyr::select(EVB6, c(1,4))
EVE8_strand = dplyr::select(EVE8, c(1,4))
CDD4_strand = dplyr::select(CDD4, c(1,4))
CDD8_strand = dplyr::select(CDD8, c(1,4))
delB5_strand = dplyr::select(delB5, c(1,4))
delC10_strand = dplyr::select(delC10, c(1,4))

#Let's put names on the gene count columns that correspond to the sample so that we can merge these into one big data frame
colnames(EVB6_strand) = c("Gene","EVB6")
colnames(EVE8_strand) = c("Gene", "EVE8")
colnames(CDD4_strand) = c("Gene", "CDD4")
colnames(CDD8_strand) = c("Gene", "CDD8")
colnames(delB5_strand) = c("Gene", "delB5")
colnames(delC10_strand) = c("Gene", "delC10")

master = inner_join(EVB6_strand, EVE8_strand, by = "Gene")
master = inner_join(master, CDD4_strand, by = "Gene")
master = inner_join(master, CDD8_strand, by = "Gene")
master = inner_join(master, delB5_strand, by = "Gene")
master = inner_join(master, delC10_strand, by = "Gene")

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
pdf(file="plots/MDSplot_TMM_samples.pdf", width = 6, height = 4)
plotMDS(dg, labels = c("EVB6", "EVE8", "CDD4", "CDD8", "delB5", "delC10"), top = 500) #top = 500 is the default setting anyways
dev.off()


#A good place to start with data visualization here is some volcano plots.
#Let's first look at CIC-DUX4 vs. EV, to see if the gene signatures we expect are being altered.
#I'm going to move the gene names back into the dataframes just for ease.
CICDUX4.EV_results_volcano = tibble::rownames_to_column(CICDUX4.EV_results, "Gene")
dC1.EV_results_volcano = tibble::rownames_to_column(dC1.EV_results, "Gene")
dC1.CICDUX4.results_volcano = tibble::rownames_to_column(dC1.CICDUX4.results, "Gene")


#note that I am transforming the fdr-adjusted pvalue by -log10(fdr) in the plotting commands.
#I'm going to designate a log2FC of +/-2 as a significant cutoff, and -log10(fdr) of >2 as significant (i.e. q < 0.01)

#the following blocks of code just generate simple volcano plots to see how many genes meet significance cutoffs for each comparison.

#CICDUX4 vs EV
volcano_CICDUX4vsEV = ggplot(data = CICDUX4.EV_results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "gray", size = 3) +
  geom_point(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$logFC > 2 & (-log10(CICDUX4.EV_results_volcano$fdr_adj)) > 2,], 
             color = "#fc558f", size = 3) +
  geom_point(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$logFC < -2 & (-log10(CICDUX4.EV_results_volcano$fdr_adj)) > 2,], 
             color = "#e3ac8d", size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$Gene %in% c("Vgf", "Etv4"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, nudge_x = 3) +
  geom_text_repel(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$Gene %in% c("Etv1", "Etv5"),], 
                  aes(label = Gene), nudge_y = 0.7, box.padding = 0.5, nudge_x = -3) +
  theme_bw() +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-20,20), ylim = c(0,4)) +
  ggtitle("CICDUX4 vs EV")

#dC1 vs EV
volcano_dC1vsEV = ggplot(data = dC1.EV_results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "gray", size = 3) +
  geom_point(data = dC1.EV_results_volcano[dC1.EV_results_volcano$logFC > 2 & (-log10(dC1.EV_results_volcano$fdr_adj)) > 2,], 
             color = "#006fc3", size = 3) +
  geom_point(data = dC1.EV_results_volcano[dC1.EV_results_volcano$logFC < -2 & (-log10(dC1.EV_results_volcano$fdr_adj)) > 2,], 
             color = "#e3ac8d", size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% c("Vgf","Etv4"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5) +
  geom_text_repel(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% c("Etv1", "Etv5"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, nudge_x = -3) +
  theme_bw() +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-20,20), ylim = c(0,4)) +
  ggtitle("dC1 vs EV")

#Plot these two together using patchwork
patchwork_volcano = volcano_CICDUX4vsEV / volcano_dC1vsEV
pdf(file = "plots/volcano_both_vs_EV.pdf", width = 8, height = 12)
patchwork_volcano
dev.off()
  

#dC1 vs CICDUX4
#I'm keeping this code around in case I want to come back to it, but won't have it put out a plot for now.
# 
# ggplot(data = dC1.CICDUX4.results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
#   geom_point(color = "gray") +
#   geom_point(data = dC1.CICDUX4.results_volcano[dC1.CICDUX4.results_volcano$logFC > 2 & (-log10(dC1.CICDUX4.results_volcano$fdr_adj)) > 2,], 
#              color = "#006fc3") +
#   geom_point(data = dC1.CICDUX4.results_volcano[dC1.CICDUX4.results_volcano$logFC < -2 & (-log10(dC1.CICDUX4.results_volcano$fdr_adj)) > 2,], 
#              color = "#fc558f") +
#   geom_vline(xintercept = 2, linetype = "dashed") +
#   geom_vline(xintercept = -2, linetype = "dashed") +
#   geom_hline(yintercept = 2, linetype = "dashed") +
#   theme_bw() +
#   xlab("log2FC") +
#   ylab("-log10(q)") +
#   xlim(c(-20,20)) +
#   ylim(c(0,4))
# 


#I'm also interested in seeing where all of the genes turned on (log2FC > 2, -log10(q) > 2) in NIH/3T3 CIC-DUX4 go in the dC1 comparison
volcano_CICDUX4_EV_all_on_in_NIH3T3 = ggplot(data = CICDUX4.EV_results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "gray", size = 3) +
  geom_point(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$logFC > 2 & -log10(CICDUX4.EV_results_volcano$fdr_adj) > 2,], 
             color = "#ae1604", size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$Gene %in% c("Vgf", "Etv4"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, nudge_x = 3, size = 7) +
  geom_text_repel(data = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$Gene %in% c("Etv1", "Etv5"),], 
                  aes(label = Gene), nudge_y = 0.7, box.padding = 0.5, nudge_x = -3, size = 7) +
  theme_bw() +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-20,20), ylim = c(0,4)) +
  ggtitle("CICDUX4 vs EV, All NIH/3T3 CIC-DUX4 Genes Turned on in Maroon") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 17))

#define what those genes are, which are turned on
all_NIH3T3_CICDUX4_on = CICDUX4.EV_results_volcano[CICDUX4.EV_results_volcano$logFC > 2 & -log10(CICDUX4.EV_results_volcano$fdr_adj) > 2,]$Gene

#same genes from full length CIC-DUX4 vs EV, but plotted for dC1 comparison
volcano_dC1_EV_all_on_in_NIH3T3 = ggplot(data = dC1.EV_results_volcano, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "gray", size = 3) +
  geom_point(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% all_NIH3T3_CICDUX4_on,], 
             color = "#ae1604", size = 3) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% c("Vgf","Etv4"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, size = 7) +
  geom_text_repel(data = dC1.EV_results_volcano[dC1.EV_results_volcano$Gene %in% c("Etv1", "Etv5"),], 
                  aes(label = Gene), nudge_y = 0.5, box.padding = 0.5, nudge_x = -3, size = 7) +
  theme_bw() +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-20,20), ylim = c(0,4)) +
  ggtitle("dC1 vs EV, All NIH/3T3 Genes Turned on in Maroon") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 17))

patchwork_volcano_all_NIH3T3_on_genes = volcano_CICDUX4_EV_all_on_in_NIH3T3 / volcano_dC1_EV_all_on_in_NIH3T3
pdf(file = "plots/volcano_all_NIH3T3_on_genes.pdf", width = 8, height = 12)
patchwork_volcano_all_NIH3T3_on_genes
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

#Recall that all_NIH3T3_CICDUX4_on holds  all genes turned on in the CIC-DUX4 vs EV comparison

#I want to check to make sure I'm not missing any genes which are significantly turned on in dC1 but not CICDUX4
on_dC1 = dC1.EV_results_volcano[dC1.EV_results_volcano$logFC > 2 & -log10(dC1.EV_results_volcano$fdr_adj) > 2,]$Gene
#Only 75 genes. are there any which are not in the list of genes turned on in the CICDUX4 vs EV comparison?
on_dC1_not_CICDUX4 = on_dC1[!on_dC1 %in% all_NIH3T3_CICDUX4_on]
#Actually yes, 31 genes.

#Now we can combine the two lists to get all genes turned on in at least one of the comparisons
genes_on_either = c(all_NIH3T3_CICDUX4_on, on_dC1_not_CICDUX4)
#359 genes

#Now let's look at a heatmap of this
logcpm_heatmap_genes_on_either = logcpm[rownames(logcpm) %in% genes_on_either,]
#Plotting a heatmap with scaling within rows turned on
pdf(file = "plots/heatmap_genes_on_either.pdf", width = 4, height = 6)
pheatmap(logcpm_heatmap_genes_on_either, 
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
logcpm_pivot$sample = factor(logcpm_pivot$sample, levels = c("EVB6", "EVE8", "CDD4", "CDD8", "delB5", "delC10"))

#add a variable about what kind of sample each is
logcpm_pivot = dplyr::mutate(logcpm_pivot, Condition = ifelse(sample %in% c("EVB6","EVE8"), "EV", ifelse(sample %in% c("CDD4","CDD8"),"CICDUX4","dC1")))

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
    coord_cartesian(ylim = c(0,15)) +
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

pdf(file="plots/patchwork_logcpm_true_targets.pdf", width = 10, height = 5)
patchwork_logcpm_true_targets
dev.off()

#finally, write some datasets to csv files to share with people.
#logcpm
write.table(logcpm, file = "outputs/log2_of_cpm.txt", row.names = F, sep = "\t", quote = F)

#raw counts, stranded, combined
#put gene names back on master as a variable to simplify writing
master_for_writing = tibble::rownames_to_column(master,"gene")
write.table(master_for_writing, file = "outputs/counts_NIH3T3_dC1.txt", row.names = F, sep = "\t", quote = F)

#CIC::DUX4 vs EV edgeR output, with gene labels
write.table(CICDUX4.EV_results_volcano, file = "outputs/CD4_vs_EV_DE.txt", row.names = F, sep = "\t", quote = F)

#dC1 vs EV edgeR output, with gene labels
write.table(dC1.EV_results_volcano, file = "outputs/dC1_vs_EV_DE.txt", row.names = F, sep = "\t", quote = F)

