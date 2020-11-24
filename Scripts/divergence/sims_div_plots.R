
# LOAD LIBRARIES

# Set working directory
setwd("C:/Users/paulc/OneDrive/Documents/Introgression_Project/Results/Divergence/new_data/")

#install.packages("ggplot2")
library("ggplot2")
#install.packages("reshape2")
library("reshape2")
#install.packages("randomcoloR")
library("randomcoloR")
#install.packages("RColorBrewer")
library("RColorBrewer")
#install.packages("pals")
library("pals")
#install.packages("ggpubr")
library("ggpubr")
#install.packages("ggforce")
library("ggforce")

#---------------------------------------------

# INPUTS (EDIT THESE)

# Provide file name details
prefix <- "all_pop_new_"
chr <- "ChrX"
window_size <- 100000
window_type <- "total"

#---------------------------------------------

# DATA PREPARATION

# Combine details into file names
file_name <- paste0(prefix,chr,'_',window_type,'_',as.character(format(window_size,scientific=F)))
full_div_file <- paste0(file_name,".div_full")
win_div_file <- paste0(file_name,".div_window")
names_file <- paste0(file_name,".names")

# Read in full divergence file
full_div <- read.csv(full_div_file, sep="\t", header=FALSE)

# Read in filename file
names <- read.csv(names_file, sep="\t", header=FALSE)
names_abr <- unlist(lapply(as.character(names), function(x) {sub(".*[_]([^.]+)[_].*", "\\1", x)}))

# Read in window-based divergence file
win_div <- read.csv(win_div_file, sep="\t", header=FALSE)

# Read in recombination rate data
rec <- read.csv("mel_recomb.csv")

# Subset data by genomic region
region_df <- subset(rec, rec$ï..Arm == chr)
region_df$ï..Arm <- NULL

# Add columns names to full divergence df
div_df <- data.frame(matrix(nrow=2, ncol=length(names_abr)))
div_df[1,] <- names_abr
div_df[2,] <- as.data.frame(full_div)

# Add column names and window indices to window-based divergence df
win_div$window <- seq.int(nrow(win_div))
colnames(win_div) <- c(names_abr, "window")

#---------------------------------------------

# PLOTTING PREPARATION

# Melt window-based divergence data to prepare for plotting
win_data <- melt(win_div, id.vars = "window")

# Order the categories in the window-based data
# win_data$variable <- factor(win_data$variable, levels=as.character(names_abr))
win_data$variable <- factor(win_data$variable, levels=c("I","RAL","W","FR","N",
                                                        "EA","EF","Kenya","RG",
                                                        "CO","GA","GU","NG","WAF","KF","MW",
                                                        "SD","SP","ZI","T","EG","B"))

# Transpose full divergence data to prepare for plotting
full_data <- t(div_df)

# Put the full divergence data into a df and order to categories
full_data_df <- as.data.frame(full_data)
full_data_df$V1 <- factor(full_data_df$V1, levels=c("I","RAL","W","FR","N",
                                                    "EA","EF","Kenya","RG",
                                                    "CO","GA","GU","NG","WAF","KF","MW",
                                                    "SD","SP","ZI","T","EG","B"))

# Make a color palette for each geographic region
pal_sa <- brewer.pal(6, "Oranges")
pal_ea <- brewer.pal(5, "Greens")
pal_wa <- brewer.pal(6, "Purples")
pal_na <- "999999"
pal_aus <- "darkkhaki"
pal_c <- "grey"
pal_eu <- (brewer.pal(3, "Blues"))
pal_us <- brewer.pal(4, "Reds")

# Manually assign colors to each sampling location
group.colors <- c(B=pal_c, CO=pal_wa[2], EA=pal_ea[2], EF=pal_ea[3], EG=pal_na, FR=pal_eu[2],
                  GA=pal_wa[3], GU=pal_wa[4], I=pal_us[2], KF=pal_sa[2], Kenya=pal_ea[4], 
                  MW=pal_sa[3], NG=pal_wa[5], N=pal_eu[3], RAL=pal_us[3], RG=pal_ea[5],
                  SD=pal_sa[4], SP=pal_sa[5], T=pal_aus, WAF=pal_wa[6], W=pal_us[4], ZI=pal_sa[6])

# # Open pdf
# pdf_name = paste0("AllPop_",chr,"_Full.pdf")
# pdf(pdf_name, width=15, height=10)

# Make bar chart showing full divergence comparison between populations (change y-axis of barplot to better visualize differences)
ggplot(data=full_data_df, aes(x=full_data_df$V1, y=as.numeric(full_data_df$V2), fill=full_data_df$V1)) +
  geom_bar(stat="identity")+
  theme(legend.position="none")+
  labs(fill = "Population")+
  theme_update(text = element_text(size=30))+
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))+
  ylab("Per-site divergence")+
  xlab("Population")+
  ggtitle("Divergence between D. mel and D. sims")+
  scale_fill_manual(values=group.colors)+
  coord_cartesian(ylim=c(0.047,0.056))

#ggsave("Pop_div2R.png", units="in", width=15, height=10, dpi=300, device = 'png')

#15x10

# dev.off()
# rm(list=ls())

#---------------------------------------------

# FOR FILTERED DATA

# Open pdf
# pdf_name = paste0("AllPop_",chr,"_Win_Filt_",as.character(format(window_size,scientific=F)),".pdf")
# pdf(pdf_name, width=20, height=8)

# Make a line chart showing window-based divergence comparison between populations
ggplot(data=win_data, aes(x=window, y=value, group=variable, color=variable)) +
  geom_line()+
  labs(color = "Population")+
  ylab("Per-site divergence")+
  xlab("Genomic window")+
  ggtitle("Divergence between D. mel and D. sims")+
  theme_update(text = element_text(size=40))+
  scale_color_manual(values=group.colors)+
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(hjust = 0.5), legend.key=element_rect(fill=NA))

#ggsave("Per_site_div_Dmel_Dsim.png", units="in", width=20, height=8, dpi=300, device = 'png')

#20X8

# dev.off()
# rm(list=ls())

#---------------------------------------------

# FOR TOTAL WINDOW DATA

# Convert window number into genomic coordinate
win_data$window <- win_data$window*window_size

# Put the data in plotting format
region_points <- melt(region_df, id.vars = "sex.avg.c..cM.Mb.")
region_points$variable <- NULL
region_points$value <- as.numeric(region_points$value)
# # Open pdf
# pdf_name = paste0("AllPop_",chr,"_Win_Tot_Rec_",as.character(format(window_size,scientific=F)),".pdf")
# pdf(pdf_name, width=22, height=9)

# Plot the divergence data on top of recombination rate
p1 <- ggplot(data=win_data, aes(x=window, y=value, group=variable, color=variable)) +
  geom_line()+
  labs(color = "Population")+
  ylab("Per-site\ndivergence")+
  ggtitle("Divergence between D. mel and D. sims")+
  theme_update(text = element_text(size=40))+
  scale_color_manual(values=group.colors)+
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(hjust = 0.5), legend.key=element_rect(fill=NA),
        axis.text.x = element_blank(), axis.title.x = element_blank()) +xlim(c(0,21500000))

p2 <- ggplot(data=region_points, aes(x=value, y=sex.avg.c..cM.Mb.)) + 
  geom_line() +
  ylab("Recombination\nrate (cM/Mb)")+
  xlab("Genomic position") + xlim(c(0,21500000))

ggarrange(p1, p2, nrow = 2, ncol = 1,
          common.legend = TRUE, legend = "right", align = "v")
#ggsave(paste0("Divergence_Data_of_",chr,".png"), units="in", width=22, height=9, dpi=300, device = 'png')

#22x9

# dev.off()
# rm(list=ls())

#---------------------------------------------

# PLOT RATIO OF OTHER POPS:ZI DIVERGENCE TO SIMS

# Divide each population's divergence by ZI divergence value
ratio_tab <- cbind((win_div[,1:21] / win_div[,22]),win_div$window)

# Convert window number into genomic coordinate
ratio_tab$`win_div$window` <- (ratio_tab$`win_div$window`)*window_size

# Make sure both tables contain the same windows (rec. rate info has couple extra windows)
stop_val <- max(ratio_tab$`win_div$window`)
rec_trunc <- subset(region_df, region_df$stop <= stop_val)


# Put the data in plotting format
region_points <- melt(region_df, id.vars = "sex.avg.c..cM.Mb.")
region_points$variable <- NULL
region_points$value <- as.numeric(region_points$value)

# Melt window-based divergence data to prepare for plotting
ratio_data <- melt(ratio_tab, id.vars = "win_div$window")

# Order the categories in the window-based data
ratio_data$variable <- factor(ratio_data$variable, levels=c("I","RAL","W","FR","N",
                                                        "EA","EF","Kenya","RG",
                                                        "CO","GA","GU","NG","WAF","KF","MW",
                                                        "SD","SP","T","EG","B"))

# Make a color palette for each geographic region
pal_sa <- brewer.pal(6, "Oranges")
pal_ea <- brewer.pal(5, "Greens")
pal_wa <- brewer.pal(6, "Purples")
pal_na <- "999999"
pal_aus <- "darkkhaki"
pal_c <- "grey"
pal_eu <- (brewer.pal(3, "Blues"))
pal_us <- brewer.pal(4, "Reds")

# Manually assign colors to each sampling location
group.colors <- c(B=pal_c, CO=pal_wa[2], EA=pal_ea[2], EF=pal_ea[3], EG=pal_na, FR=pal_eu[2],
                  GA=pal_wa[3], GU=pal_wa[4], I=pal_us[2], KF=pal_sa[2], Kenya=pal_ea[4], 
                  MW=pal_sa[3], NG=pal_wa[5], N=pal_eu[3], RAL=pal_us[3], RG=pal_ea[5],
                  SD=pal_sa[4], SP=pal_sa[5], T=pal_aus, WAF=pal_wa[6], W=pal_us[4])

# # Open pdf
# dev.off()
# pdf_name = paste0("AllPop_",chr,"_Win_Tot_Rec_Ratio_",as.character(format(window_size,scientific=F)),".pdf")
# pdf(pdf_name, width=22, height=9)


# Plot the divergence data on top of recombination rate

p1 <- ggplot(data=ratio_data, aes(x=`win_div$window`, y=value, group=variable, color=variable)) +
  geom_line()+
  labs(color = "Population")+
  ylab("Divergence\nratio")+
  ggtitle("Ratio of each D. mel population to ZI divergence to D. sims")+
  theme_update(text = element_text(size=40))+
  scale_color_manual(values=group.colors)+
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(hjust = 0.5), legend.key=element_rect(fill=NA),
        axis.text.x = element_blank(), axis.title.x = element_blank()) + xlim(c(0,21500000))

p2 <- ggplot(data=region_points, aes(x=value, y=sex.avg.c..cM.Mb.)) + 
  geom_line() +
  ylab("Recombination\nrate (cM/Mb)")+
  xlab("Genomic position") + xlim(c(0,21500000))

#data = as.numeric

ggarrange(p1, p2, nrow = 2, ncol = 1,
          common.legend = TRUE, legend = "right", align = "v")
ggsave(paste0("Dmel_pop_to_ZI_div_to_Dsims_of_",chr,".png"), units="in", width=22, height=9, dpi=300, device = 'png')

#22x9

# dev.off()
# rm(list=ls())

#---------------------------------------------

# PLOT RATIO OF CO/ZI DIVERGENCE TO SIMS

# Divide CO divergence by ZI divergence value
ratio_tab <- as.data.frame(cbind((win_div$CO / win_div$ZI),win_div$window))
colnames(ratio_tab) <- c("CO/ZI", "window")

# Convert window number into genomic coordinate
ratio_tab$window <- (ratio_tab$window)*window_size

# Make sure both tables contain the same windows (rec. rate info has couple extra windows)
stop_val <- max(ratio_tab$window)
rec_trunc <- subset(region_df, region_df$stop <= stop_val)



# Put the data in plotting format
region_points <- melt(region_df, id.vars = "sex.avg.c..cM.Mb.")
region_points$variable <- NULL
region_points$value <- as.numeric(region_points$value)

# Melt window-based divergence data to prepare for plotting
ratio_data <- melt(ratio_tab, id.vars = "window")

# Open pdf
# dev.off()
# pdf_name = paste0("CO_ZI_",chr,"_Win_Tot_Rec_Ratio_",as.character(format(window_size,scientific=F)),".pdf")
# pdf(pdf_name, width=22, height=9)

# Plot the divergence data on top of recombination rate
p1 <- ggplot(data=ratio_data, aes(x=window, y=value))+
  geom_line()+
  ylab("Divergence\nratio")+
  ggtitle("Ratio of CO to ZI divergence to D. sims")+
  theme_update(text = element_text(size=40))+
  scale_color_manual(values=group.colors)+
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(hjust = 0.5), legend.key=element_rect(fill=NA),
        axis.text.x = element_blank(), axis.title.x = element_blank())


p2 <- ggplot(data=region_points, aes(x=value, y=sex.avg.c..cM.Mb.)) + 
  geom_line() +
  ylab("Recombination\nrate (cM/Mb)")+
  xlab("Genomic position")

ggarrange(p1, p2, nrow = 2, ncol = 1,
          common.legend = TRUE, legend = "right", align = "v")
ggsave(paste0("Ratio_of_CO_to_ZI_",chr,".png"), width=22, height=9, dpi=300, device = 'png')

#22x9

# dev.off()
# rm(list=ls())

#---------------------------------------------

# ADDITIONAL PLOTTING OPTIONS

# Make bar chart showing full divergence comparison between populations with bars zoomed in
ggplot(data=full_data_df, aes(x=full_data_df$V1, y=as.numeric(full_data_df$V2), fill=full_data_df$V1)) +
  geom_bar(stat="identity")+
  theme(legend.position="none")+
  labs(fill = "Population")+
  theme_update(text = element_text(size=30))+
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))+
  ylab("Per-site divergence")+
  xlab("Population")+
  ggtitle("Divergence between D. mel and D. sims")+
  #scale_fill_manual(values=group.colors)+
  facet_zoom(ylim = c(0.045,0.05))
#ggsave("Dmel_Dsim_Div.png", units="in", width=15, height=10, dpi=300, device = 'png')

# Make bar chart showing full divergence comparison between populations
ggplot(data=full_data_df, aes(x=full_data_df$V1, y=as.numeric(full_data_df$V2), fill=full_data_df$V1)) +
  geom_bar(stat="identity")+
  theme(legend.position="none")+
  labs(fill = "Population")+
  theme_update(text = element_text(size=30))+
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))+
  ylab("Per-site divergence")+
  xlab("Population")+
  ggtitle("Divergence between D. mel and D. sims")+
  #scale_fill_manual(values=group.colors)
#15x10

