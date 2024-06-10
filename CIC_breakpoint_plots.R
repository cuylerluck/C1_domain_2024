  #Cuyler Luck, Okimoto Lab
  #Cuyler.Luck@ucsf.edu or Ross.Okimoto@ucsf.edu
  
  #load necessary libraries. install with install.packages() if needed.
  library(data.table) #1.14.8
  library(dplyr) #1.1.1
  library(ggplot2) #3.4.1
  library(ggExtra) #0.10.0
  library(circlize) #0.4.15
  library(RColorBrewer) #1.1-3
  library(ComplexHeatmap) #2.14.0
  
  
  #set working directory
  setwd("~/Desktop/UCSF/General/Okimoto Lab/Breakpoints/final_data_Apr2024/")
  #set a random number seed so that the colors are the same every time
  set.seed(92363)
  
  #clear last circos plot, in case you forgot to
  circos.clear()
  
  #read in master data .csv
  master = fread("CIC_fusion_breakpoint_database_checked_042924_exact_available.csv", header=TRUE, fill = TRUE)
  
  #clean doi column to remove underscores in the beginning if they are present
  for (doi in master$Source_doi) {
    if(substr(doi,1,1) == "_"){
      master[master$Source_doi == doi,]$Source_doi = substr(doi, 2, nchar(doi))
    }
  }
  
  #give each breakpoint a unique breakpoint_ID based on its row ID in the current dataframe, and make this variable a factor
  master = tibble::rowid_to_column(master, "breakpoint_ID")
  master$breakpoint_ID = as.factor(master$breakpoint_ID)
  
  #create a subset of master that only keeps RNA samples for downstream processing after the Circos plot
  master_RNA = master[master$DNA_or_RNA == "RNA",]
  
  ### Circos plot (Figure 1B) ------------------------------------------------------------------
  
  #open command to save plot as TIFF file
  tiff("plots/CIC_fusion_Circos.tiff", units = "in", width = 8, height = 8, res = 600, bg = "transparent")
  
  #gets all R colors except for gray/greys, thanks to https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  #sample 50 colors to use for coloring fusions
  fusion_color_options = sample(color,50)
  #extract unique fusions in dataset
  fusion_options = unique(master$partner)
  #put unique fusions into a dataframe, and assign each a color
  fusion_colors_assigned = data.frame(fusion_options)
  fusion_colors_assigned$color = fusion_color_options[1:length(fusion_options)]
  #write an easy function to return the color for a given partner
  getFusionColor = function(partner){
    return(fusion_colors_assigned[fusion_colors_assigned$fusion_options == partner,]$color)
  }
  #add a variable to master defining the color for each fusion, using getFusionColor
  master$fusion_color = sapply(master$partner, getFusionColor)
  
  #sample 200 colors to use for coloring DOI
  doi_color_options = sample(color,200)
  #extract unique doi in dataset
  doi_options = unique(master$Source_doi)
  #put unique doi into a dataframe, and assign each a color
  doi_colors_assigned = data.frame(doi_options)
  doi_colors_assigned$color = doi_color_options[1:length(doi_options)]
  #write getDOIColor
  getDOIColor = function(doi){
    return(doi_colors_assigned[doi_colors_assigned$doi_options == doi,]$color)
  }
  #add a variable to master defining the color for each doi, using getDOIColor
  master$DOI_color = sapply(master$Source_doi, getDOIColor)
  
  #figure out which sample IDs appear multiple times (this should be because they are the same patient with multiple breakpoints)
  sample_id_dups = master$Sample[duplicated(master$Sample)]
  #add a variable to master defining the color for a sample whose ID appears only once vs. one whose ID appears multiple times (i.e. a sample which has one documented breakpoint vs multiple)
  master$multiple_bp_sample_color = ifelse(master$Sample %in% sample_id_dups, "pink", "gray")
  
  #add a variable to master defining the color for DNA (navy) vs. RNA (maroon) vs. unclear (gray)
  master$XNA_color = ifelse(master$DNA_or_RNA == "DNA", "navy", ifelse(master$DNA_or_RNA == "RNA", "maroon", "gray"))
  
  #add a variable to master defining the color for patient vs. cell line (sample type)
  master$type_color = ifelse(master$Sample_type == "Patient", "green", "purple")
  
  #add a variable to master defining color if both CIC and partner breakpoints are available or not
  master$bothAvail_color = ifelse((!is.na(master$CIC_last_nt) & !is.na(master$partner_first_nt)), "gold", "gray")
  
  #set graphing parameters
  circos.par(start.degree = 90, gap.degree = 0, track.height = 0.1)
  
  #create Circos plot for sample characteristics
  circos.initialize(master$breakpoint_ID, xlim = c(0,1))
  
  #create a track for the fusion partner
  circos.track(master$breakpoint_ID, ylim = c(0,1), bg.col = master$fusion_color)
  #create a track for the paper DOI for the breakpoint
  circos.track(master$breakpoint_ID, ylim = c(0,1), bg.col = master$DOI_color)
  #create a track for if sample ID has multiple breakpoints or not
  circos.track(master$breakpoint_ID, ylim = c(0,1), bg.col = master$multiple_bp_sample_color)
  #create a track for RNA or DNA breakpoint
  circos.track(master$breakpoint_ID, ylim = c(0,1), bg.col = master$XNA_color)
  #create a track for patient vs. cell line
  circos.track(master$breakpoint_ID, ylim = c(0,1), bg.col = master$type_color)
  #create a track for if both CIC and partner breakpoints are available
  circos.track(master$breakpoint_ID, ylim = c(0,1), bg.col = master$bothAvail_color)
  
  #close TIFF device
  dev.off()
  circos.clear()
  
  #dev.off again
  dev.off()
  #now make the cut-away plot of a single sector to describe what each ring means with legends
  
  #turn on the TIFF device
  tiff("plots/CIC_fusion_Circos_legends.tiff", units = "in", width = 8, height = 8, res = 600, bg = "transparent")
  
  #a single sector takes up a number of degrees equal to 360/the number of total breakpoints
  single_sector = 360/length(master$breakpoint_ID)
  #start the canvas at the top (90 degrees) and only let it continue radially to the length of one sector the same
  #as in the full plot
  circos.par("canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1), "start.degree" = 90, "gap.after" = 360-single_sector, track.height = 0.1)
  #just using one sector, in this case the one for SBRCT3
  sectors = "1" 
  circos.initialize(sectors = sectors, xlim = c(0,1))
  #create a track for the fusion partner
  circos.track(sectors, ylim = c(0,1), bg.col = master$fusion_color)
  #create a track for the paper DOI for the breakpoint
  circos.track(sectors, ylim = c(0,1), bg.col = master$DOI_color)
  #create a track for if sample ID has multiple breakpoints or not
  circos.track(sectors, ylim = c(0,1), bg.col = master$multiple_bp_sample_color)
  #create a track for RNA or DNA breakpoint
  circos.track(sectors, ylim = c(0,1), bg.col = master$XNA_color)
  #create a track for patient vs. cell line
  circos.track(sectors, ylim = c(0,1), bg.col = master$type_color)
  #create a track for if both CIC and partner breakpoints are available
  circos.track(sectors, ylim = c(0,1), bg.col = master$bothAvail_color)
  
  #make legends and draw them onto the plot as necessary
  fusion_label = Legend(labels = fusion_colors_assigned$fusion_options, legend_gp = gpar(fill = fusion_colors_assigned$color), title = "Fusion Partner", labels_gp = gpar(fontsize=16), title_gp = gpar(fontsize=18))
  draw(fusion_label, x = unit(0.3, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))
  
  DOI_label = Legend(labels = "See full data table", title = "DOI of Publication", title_gp = gpar(fontsize=18), labels_gp = gpar(fontsize=16))
  draw(DOI_label, x = unit (0.3, "npc"), y = unit(0.785, "npc"), just = c("left", "center"))
  
  multiple_bp_label = Legend(labels = c("Yes", "No"), legend_gp = gpar(fill = c("pink", "gray")), title = "Multiple Breakpoints \nin the Same Sample?", labels_gp = gpar(fontsize=16), title_gp = gpar(fontsize=18))
  draw(multiple_bp_label, x = unit(0.3, "npc"), y = unit(0.675, "npc"), just = c("left", "center"))
  
  XNA_label = Legend(labels = c("DNA", "RNA", "Unclear"), legend_gp = gpar(fill = c("navy", "maroon", "gray")), title = "DNA or RNA Data", labels_gp = gpar(fontsize=16), title_gp = gpar(fontsize=18))
  draw(XNA_label, x = unit(0.3, "npc"), y = unit(0.5385, "npc"), just = c("left", "center"))
  
  type_label = Legend(labels = c("Patient", "Cell line"), legend_gp = gpar(fill = c("green", "purple")), title = "Sample Type", labels_gp = gpar(fontsize=16), title_gp = gpar(fontsize=18))
  draw(type_label, x = unit(0.3, "npc"), y = unit(0.4255, "npc"), just = c("left", "center"))
  
  bothAvail_label = Legend(labels = c("Yes", "No"), legend_gp = gpar(fill = c("gold", "gray")), title = "Exact Breakpoint Available \nfor Both Partners?", labels_gp = gpar(fontsize=16), title_gp = gpar(fontsize=18))
  draw(bothAvail_label, x = unit(0.3, "npc"), y = unit(0.305, "npc"), just = c("left", "center"))
  
  #turn off TIFF device
  dev.off()
  circos.clear()
  
  #dev.off again
  dev.off()
  ### End Circos plot (Figure 1A) ------------------------------------------------------------------
  
  
  ### Histogram of CIC breakpoints --------------------------------------------
  
  #remove RNA samples without information on both partner breakpoints
  master_RNA_both = master_RNA[!is.na(master_RNA$CIC_last_nt) & !is.na(master_RNA$partner_first_nt),]
  
  #also determining in- vs out-of-frame before doing more with data
  #For the CIC reference sequence, the break coordinate mod (%%) 3 result has this meaning:
  # CIC %% 3 = 0, break is in first nucleotide of codon
  # CIC %% 3 = 1, break is in second nucleotide of codon
  # CIC %% 3 = 2, break is in third nucleotide of codon
  
  #For DUX4, the CDS is shifted within the reference sequence, so it follows these rules:
  # DUX4 %% 3 = 0, break is in third nucleotide of codon
  # DUX4 %% 3 = 1, break is in first nucleotide of codon
  # DUX4 %% 3 = 2, break is in second nucleotide of codon
  
  #Therefore, an in-frame fusion for CIC-DUX4 is made up of one of the following combinations (assuming no filler):
  # CIC %% 3 = 0 & DUX4 %% 3 = 2
  # CIC %% 3 = 1 & DUX4 %% 3 = 0
  # CIC %% 3 = 2 & DUX4 %% 3 = 1
  
  #And if there's any filler, simply add the filler to CIC value prior to the mod. If filler is NA, add 0.
  #e.g. for Brcic2020-4: (5001 + 0) %% 3 = 0, 536 %% 3 = 2, 0 / 2 is TRUE
  #in-frame
  
  #Mathematically, CIC-DUX4 is in-frame if (((CIC + filler) %% 3) + 2) %% 3) == DUX4 %% 3
  
  #Following similar logic for NUTM1, LEUTX, and FOXO4 for our reference sequences, these values indicate in-frame fusions:
  #NUTM1: (((CIC + filler) %% 3) + 0) %% 3) == NUTM1 %% 3 
  #LEUTX: (((CIC + filler) %% 3) + 0) %% 3) == LEUTX %% 3
  #FOXO4: (((CIC + filler) %% 3) + 2) %% 3) == FOXO4 %% 3
  
  #With these rules, we can make a new variable saying whether a given breakpoint is in- or out-of-frame.
  #WARNING: This depends on accurate reporting of the sequences and accurate interpretation by us. Always check source pubs for detailed analyses.
  
  #For the RNA-only data where both partners are available, predict which are in-frame:
  
  isInFrame = function(partner, CIC_last_nt, partner_first_nt, Filler){
    if(is.na(Filler)){
      fill = 0
    }
    else{
      fill = nchar(Filler)
    }
    
    if(partner == "DUX4" | partner == "FOXO4"){
      return((((CIC_last_nt + fill) %% 3) + 2) %% 3 == (partner_first_nt %% 3))
    }
    
    else if(partner == "NUTM1" | partner == "LEUTX"){
      return((((CIC_last_nt + fill) %% 3) + 0) %% 3 == (partner_first_nt %% 3))
    }
    
  }
  
  master_RNA_both_frame = master_RNA_both
  master_RNA_both_frame$isInFrame = FALSE
  counter = 1
  for(counter in 1:length(master_RNA_both_frame$breakpoint_ID)){
    master_RNA_both_frame$isInFrame[counter] = isInFrame(master_RNA_both_frame$partner[counter], master_RNA_both_frame$CIC_last_nt[counter], master_RNA_both_frame$partner_first_nt[counter], master_RNA_both_frame$Filler[counter])
  }
  
  #I will call anything that's a UTR fusion NA for framing -- there's not really a concept of in-frame for things that fuse beyond the CDS
  #So any breakpoints where the 3' partner breakpoint is after the last coding codon [STOP excluded], get NA for isInFrame
  
  master_RNA_both_frame$isInFrame[(master_RNA_both_frame$partner == "DUX4") & (master_RNA_both_frame$partner_first_nt > 1272)] = NA
  master_RNA_both_frame$isInFrame[(master_RNA_both_frame$partner == "NUTM1") & (master_RNA_both_frame$partner_first_nt > 3862)] = NA
  master_RNA_both_frame$isInFrame[(master_RNA_both_frame$partner == "LEUTX") & (master_RNA_both_frame$partner_first_nt > 661)] = NA
  master_RNA_both_frame$isInFrame[(master_RNA_both_frame$partner == "FOXO4") & (master_RNA_both_frame$partner_first_nt > 2193)] = NA
  
  #make it a factor for plotting. feels easier
  master_RNA_both_frame$isInFrame = as.factor(master_RNA_both_frame$isInFrame)
  
  #turn on TIFF device
  tiff("plots/CIC_fusion_histogram.tiff", units = "in", width = 4.5, height = 4.5, res = 600, bg = "transparent")
  
  #plot histogram of CIC breakpoints, faceted by partner and with dashed lines for important domains
  ggplot(master_RNA_both_frame, aes(x = CIC_last_nt)) + 
    geom_freqpoly(binwidth=1) + 
    facet_grid(vars(partner), scales = "free_y") +
    theme_classic() +
    xlab("") + 
    ylab("Number of Breakpoints") +
    geom_vline(xintercept=4626, linetype = "dashed", color = "blue", alpha = 0.5) + #beginning of C1 domain (from PMID: 28278156)
    geom_vline(xintercept=4793, linetype = "dashed", color = "blue", alpha = 0.5) + #end of C1 domain (from PMID: 28278156)
    geom_vline(xintercept=4467, linetype = "dashed", color = "red", alpha = 0.5) + #beginning of NLS (from PMID: 21087211)
    geom_vline(xintercept=4484, linetype = "dashed", color = "red", alpha = 0.5) + #end of NLS (from PMID: 21087211)
    geom_vline(xintercept=4209, linetype = "dashed", color = "purple", alpha = 0.5) + #beginning of ERK binding domain (from PMID: 26124095)
    geom_vline(xintercept=4433, linetype = "dashed", color = "purple", alpha = 0.5) + #end of ERK binding domain (from PMID: 26124095)
    geom_vline(xintercept=5061, linetype = "dashed", color = "gray", alpha = 0.5) + #start of STOP codon
    geom_vline(xintercept=5063, linetype = "dashed", color = "gray", alpha = 0.5) #end of STOP codon
  
  #turn off device
  dev.off()
  
  
  ### End Histogram of CIC breakpoints
  ### Scatter plot of DUX4 vs. CIC breakpoints  --------------------------------------------
  
  #subset further to just RNA, both partner coordinate present, DUX4 samples
  master_RNA_both_DUX4 = master_RNA_both_frame[master_RNA_both_frame$partner == "DUX4",]
  
  #open TIFF object
  tiff("plots/CIC_DUX4_scatterplot.tiff", units = "in", width = 4.5, height = 4.5, res = 600, bg = "transparent")
  
  #make an object for a scatterplot of DUX4 breakpoint locations vs. CIC breakpoint locations for these samples
  CIC_DUX4_scatter = ggplot(data = master_RNA_both_DUX4, aes(x = CIC_last_nt, y = partner_first_nt)) + 
    geom_point(size = 2.5) + 
    geom_point(data = master_RNA_both_DUX4[master_RNA_both_DUX4$isInFrame == FALSE,], color = "#ce5037", size = 2.5) + #color out-of-frame, be aware this may overplot in-frame with identical coordinates
    geom_point(data = master_RNA_both_DUX4[is.na(master_RNA_both_DUX4$isInFrame),], color = "#16a2de", size = 2.5) + #color 3' UTR fusions, be aware this may overplot in-frame with identical coordinates
    theme_bw() + 
    xlab("CIC (NM_015125.5)") + 
    ylab("DUX4 (NM_001306068.3)") +
    geom_vline(xintercept=4626, linetype = "dashed", color = "blue", alpha = 0.5) + #beginning of C1 domain (from PMID: 28278156)
    geom_vline(xintercept=4793, linetype = "dashed", color = "blue", alpha = 0.5) + #end of C1 domain (from PMID: 28278156)
    geom_vline(xintercept=4467, linetype = "dashed", color = "red", alpha = 0.5) + #beginning of NLS (from PMID: 21087211)
    geom_vline(xintercept=4484, linetype = "dashed", color = "red", alpha = 0.5) + #end of NLS (from PMID: 21087211)
    geom_vline(xintercept=4209, linetype = "dashed", color = "purple", alpha = 0.5) + #beginning of ERK binding domain (from PMID: 26124095)
    geom_vline(xintercept=4433, linetype = "dashed", color = "purple", alpha = 0.5) + #end of ERK binding domain (from PMID: 26124095)
    geom_vline(xintercept=5061, linetype = "dashed", color = "gray", alpha = 0.5) + #start of CIC STOP codon
    geom_vline(xintercept=5063, linetype = "dashed", color = "gray", alpha = 0.5) + #end of CIC STOP codon
    geom_hline(yintercept=55, linetype = "dashed", color = "brown", alpha = 0.5) + #start of HOX1 (from PMID: 29618456)
    geom_hline(yintercept=237, linetype = "dashed", color = "brown", alpha = 0.5) + #end of HOX1 (from PMID: 29618456)
    geom_hline(yintercept=280, linetype = "dashed", color = "dark green", alpha = 0.5) + #start of HOX2 (from PMID: 29618456)
    geom_hline(yintercept=462, linetype = "dashed", color = "dark green", alpha = 0.5) + #end of HOX2 (from PMID: 29618456)
    geom_hline(yintercept=979, linetype = "dashed", color = "dark orange", alpha = 0.5) + #start of the larger window important in transactivation (from PMID: 26951377)
    geom_hline(yintercept=1272, linetype = "dashed", color = "dark orange", alpha = 0.5) + #end of the larger window important in transactivation (from PMID: 26951377)
    geom_hline(yintercept=1027, linetype = "dashed", color = "dark red", alpha = 0.5) + #start of the smaller window important in transactivation (from PMID: 29618456)
    geom_hline(yintercept=1272, linetype = "dashed", color = "dark red", alpha = 0.5) + #end of the smaller window important in transactivation (from PMID: 29618456)
    geom_hline(yintercept=1273, linetype = "dashed", color = "gray", alpha = 0.5) + #start of DUX4 STOP codon
    geom_hline(yintercept=1275, linetype = "dashed", color = "gray", alpha = 0.5)  #end of DUX4 STOP codon

  #plot the above object with histograms for each axis as well
  ggMarginal(CIC_DUX4_scatter, type = "histogram", binwidth = 1) 
  
  #terminate TIFF object
  dev.off()
  ### End scatter plot of DUX4 vs. CIC breakpoints  --------------------------------------------
  
  
  
  
  
  
  
