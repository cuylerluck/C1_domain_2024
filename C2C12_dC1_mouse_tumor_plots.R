#Cuyler Luck, Okimoto Lab, April 2024

#This is a simple script to visualize mouse tumor sizing between groups of nude mice injected with one of the following cell lines:
# C2C12 pMY-IRES-EGFP (EV) clone G4
# C2C12 pMY-CIC-DUX4 clone C9
# C2C12 pMY-CIC-DUX4 dC1 clone D5

#Each group had five mice, with each mouse receiving one injection of 400,000 cells (50 uL, half cells in DPBS, half matrigel) per flank
#Injection date, day 0, was 3/11/24.

#set wd
setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_C2C12_dC1_mouse")

#load required libraries
library(data.table) #1.14.8
library(dplyr) #1.1.1
library(tidyr) #1.3.0
library(ggplot2) #3.4.1

#read in data from the master Excel sheet
#the original sheet is an xlsx, I saved a copy as a csv for simpler import
master_og = fread(file = "raw/mouse_tumor_sizing_C2C12_clonal_dC1_Mar2024.csv")

#let's rename columns to remove spaces and simplify for coding
master = master_og
colnames(master) = c("measure_date","group","mouse","dim_left","dim_right")

#and I am actually going to drop the data from 3/15/24, where we only had two dimensions
#this is because the volume formula I will use requires three (length, width, height)
master = master[!measure_date == "3/15/24",]

#first let's use the measure_date column to make a new variable that expresses days since injection
class(master$measure_date) #measure_date is currently a character, need it to be a date
master$measure_date = as.Date(master$measure_date, format = "%m/%d/%y")
class(master$measure_date) #great

start_date = as.Date("2024-03-11")
master = dplyr::mutate(master, days_since_inj = as.numeric(measure_date - start_date))

#now I need to divide out the data to single rows being individual injection sites, not individual mice
master = pivot_longer(master, cols = c(dim_left,dim_right), names_to = "side", names_prefix = "dim_")

#now need to split the single measurement column into three: length, width, height
#we always recorded head-to-tail first, then side-to-side (front to back) second, then depth (internal to external)
#so I will call these length / width / height respectively. though it doesn't actually matter with the formula
master = separate_wider_delim(master, value, " x ", names = c("length","width","height"))
class(master$length) = "numeric"
class(master$width) = "numeric"
class(master$height) = "numeric"

#and now we can calculate volume
#I am using (pi/6) * length * width * height here
master = dplyr::mutate(master, volume = (pi * length * width * height) / 6)
#this is in cubic mm.

#also adding a column that concatenates group + mouse number + side to get unique identifiers
master = dplyr::mutate(master, id = paste(group,mouse,side, sep=""))

#and let's remove any rows with NA [signifies mouse was sac'd]
master = master[!is.na(master$volume),]

#let's also get summary statistics at each date for each group, to overplot
summary = dplyr::group_by(master, group, days_since_inj) %>% summarise(mean = mean(volume), sd = sd(volume), n = n())

#and now we can plot
pdf(file = "plots/vol_vs_time.pdf", width = 5, height = 4)
ggplot(data = master) + 
  geom_line(data = master, aes(group = id, color = group, x = days_since_inj, y = volume), alpha = 0.5) +
  geom_errorbar(data = summary, aes(x = days_since_inj, ymin = mean-sd, ymax = mean+sd, color = group), width = 0.2) +
  geom_line(data = summary, aes(group = group, x = days_since_inj, y = mean, color = group), linewidth = 1) +
  geom_point(data = summary, aes(x = days_since_inj, y = mean, color = group, shape = group), size = 2.5) +
  geom_hline(yintercept = 50, linetype = 2, alpha = 0.5) +
  theme_classic() +
  xlab("Days Since Injection") +
  ylab(bquote("Tumor Volume "(mm^3))) +
  scale_color_manual(values = c("gray","#fc558f","#006fc3"),
                     name = "Injected Cell Line",
                     breaks = c("EV", "CD4", "dC1"),
                     labels = c("EV", "CIC::DUX4","CIC::DUX4 dC1")) +
  scale_shape_discrete(name = "Injected Cell Line",
                       breaks = c("EV", "CD4", "dC1"),
                       labels = c("EV", "CIC::DUX4","CIC::DUX4 dC1")) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250))
dev.off()

  

#this can be used to validate that summary worked properly
#test = master[master$group == "EV" & master$measure_date == as.Date("2024-04-08"),]
#mean(test$volume)
#sd(test$volume)
