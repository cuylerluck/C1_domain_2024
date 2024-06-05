# C1_domain_2024

The code in this repository accompanies the manuscript "The Capicua C1 Domain is Required for Full Activity of the CIC::DUX4 Fusion Oncoprotein" by Luck, Jacobs, and Okimoto (2024). 
These files are not plug-and-play but serve to provide transparency for how we analyzed the data. They may need to be modified for use in different environments or with different data.
With questions please contact the lead author at cuyler.luck@ucsf.edu and the corresponding author at ross.okimoto@ucsf.edu.
RNA-seq raw and processed (counts, log2cpm) data files are available at _______.

There are two main analyses that the code is divided into, with the purpose of each file described below:

## RNA-seq processing and analysis

***documentation_with_job_details_dC1.txt*** : a text file describing processing of FASTQ files with STAR (and other tools) as performed locally on a MacBook Pro and on Wynton, UCSF's high-performance computing cluster.

***NIH3T3_dC1_analysis.R*** : R script detailing differential expression analysis of RNA-seq data and used for plotting figures.

## Mouse tumor size analysis

***C2C12_dC1_mouse_tumor_plots.R*** : a simple R script to clean and plot mouse tumor size data.

***mouse_tumor_sizing_C2C12_clonal_dC1_Mar2024.csv*** : raw data for mouse tumor sizing.
