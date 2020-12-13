
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
prefix <- "outlier_"
chr <- "ChrX"
window_size <- 1000
window_type <- "filtered"

#---------------------------------------------

# DATA PREPARATION

# Combine details into file names
file_name <- paste0(prefix,chr,'_',window_type,'_',as.character(format(window_size,scientific=F)))
win_div_file <- paste0(file_name,".fine_div_window")

# Read in window-based divergence file
win_div <- read.csv(win_div_file, sep="\t", header=TRUE)

#---------------------------------------------
# Make a color palette for each geographic region
pal_sa <- brewer.pal(6, "Oranges")
pal_wa <- brewer.pal(6, "Purples")
# Manually assign colors to each sampling location
#group.colors <- c(CO=pal_wa[2], ZI=pal_sa[6], Ratio="black")
group.colors <- c(CO="purple", ZI=pal_sa[6], Ratio="black")

# PLOTTING FUNCTIONS
region_data <- subset(win_div, win_div$region_index == 0)
plot_outlier_div <- function(region_data){
  vars_of_intrest <- data.frame(region_data$window_start,region_data$CO_divergence, region_data$ZI_divergence)
  colnames(vars_of_intrest) <- c("Start","CO", "ZI")
  melted_region_data <- melt(vars_of_intrest, id = "Start")
  p1 <- ggplot(data=melted_region_data, aes(x=Start, y=value, group=variable, color=variable)) +
    geom_line()+
    labs(color = "Population")+
    ylab("Divergence")+
    ggtitle("Ratio of each D. mel population to ZI divergence to D. sims")+
    theme_update(text = element_text(size=40))+
    scale_color_manual(values=group.colors)+
    guides(colour = guide_legend(override.aes = list(size = 5)))+
    theme(plot.title = element_text(hjust = 0.5), legend.key=element_rect(fill=NA),
          axis.text.x = element_blank(), axis.title.x = element_blank())
  p2 <- ggplot(data=region_data, aes(x=window_start, y=CO_to_ZI_ratio)) + 
    geom_line() +
    ylab("Divergence\nratio")+
    xlab("Genomic position")

  
  ggarrange(p1, p2, nrow = 2, ncol = 1,
            common.legend = TRUE, legend = "right", align = "v")
  ggsave(paste0("Region_divergences_of_",chr,".png"), units="in", width=22, height=9, dpi=300, device = 'png')
}
plot_outlier_div(region_data)

