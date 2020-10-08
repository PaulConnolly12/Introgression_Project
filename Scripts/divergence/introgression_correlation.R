# LOAD LIBRARIES

# Set working directory
setwd("C:/Users/ekhow/Downloads/cluster_transfers/")

#install.packages("MASS")
library("MASS")
#install.packages("lmtest")
library("lmtest")

#---------------------------------------------

# INPUTS (EDIT THESE)

# Provide file name details
prefix <- "all_pop_"
chr <- "Chr3R"
window_size <- 100000
window_type <- "total"

#---------------------------------------------

# DATA PREPARATION

# Combine details into file names
file_name <- paste0(prefix,chr,'_',window_type,'_',as.character(format(window_size,scientific=F)))
win_div_file <- paste0(file_name,".div_window")
names_file <- paste0(file_name,".names")

# Read in filename file
names <- read.csv(names_file, sep="\t", header=FALSE)
names_abr <- unlist(lapply(as.character(names), function(x) {sub(".*[_]([^.]+)[_].*", "\\1", x)}))

# Read in window-based divergence file
win_div <- read.csv(win_div_file, sep="\t", header=FALSE)

# Read in recombination rate data
rec <- read.csv("mel_recomb.csv")

# Subset data by genomic region
region_df <- subset(rec, rec$ï..arm == chr)
region_df$ï..arm <- NULL

# Add column names and window indices to window-based divergence df
win_div$window <- seq.int(nrow(win_div))
colnames(win_div) <- c(names_abr, "window")

#---------------------------------------------

# MAKE CORRELATION DATA TABLE

# Divide CO divergence by ZI divergence value
ratio_tab <- as.data.frame(cbind((win_div$CO / win_div$ZI),win_div$window))
colnames(ratio_tab) <- c("CO/ZI", "window")

# Convert window number into genomic coordinate
ratio_tab$window <- (ratio_tab$window)*window_size

# Make sure both tables contain the same windows (rec. rate info has couple extra windows)
stop_val <- max(ratio_tab$window)
rec_trunc <- subset(region_df, region_df$stop <= stop_val)

# Combine these columns to create table with ratio of CO to ZI divergence and recombination rate
corr_tab <- data.frame("ratio_CO_ZI" = ratio_tab$`CO/ZI`, "recombination" = rec_trunc$sex_avg_c_cM_per_Mb)

#---------------------------------------------

# LINEAR REGRESSION ANALYSIS

# Plot the data points
plot(y=corr_tab$ratio_CO_ZI, x=corr_tab$recombination, xlab = 'Recombination rate (cM/Mb)', ylab = 'Ratio of CO to ZI divergence')

# Create regression model
model = lm(ratio_CO_ZI~recombination, data=corr_tab)
summary(model)

# Check assumptions of linear regression

# Plot histogram of the residuals
hist(model$resid, main="Histogram of Residuals",
     ylab="Residuals")

# Make a QQ plot of residuals
qqnorm(model$resid)
qqline(model$resid)

# Perform a Shapiro-Wilk test for normality of residuals
sresid <- studres(model) 
shapiro.test(sresid) #p-value < 0.05, not normal

# Perform Breusch-Pagan test For homoscedasticity
bptest(model) #p-value > 0.05, not homoscedastic

#---------------------------------------------

# CORRELATION ANALYSIS

# Use spearman's rank correlation
cor.test( ~ corr_tab$ratio_CO_ZI + corr_tab$recombination,
          data=corr_tab,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95) # cannot compute exact p-value with ties

# Perform kendall correlation
cor.test( ~ corr_tab$ratio_CO_ZI + corr_tab$recombination,
          data=corr_tab,
          method = "kendall",
          continuity = FALSE,
          conf.level = 0.95) # p-value = 0.8473

