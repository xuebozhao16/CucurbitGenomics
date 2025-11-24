# ==========================================================
# Title: FSC Model Likelihood Comparison and Speciation Analysis
# Species: Bottle gourd (Lagenaria siceraria)
# Author: Xuebo Zhao
# Description:
#   This script visualizes fastsimcoal2 (FSC) likelihood comparisons
#   for different gene flow models, and estimates speciation and
#   introgression times for Bottle gourd populations.
# ==========================================================

# ==========================================================
# Section 1: Compare Likelihoods of Different Gene Flow Models (D Lineage)
# ==========================================================
setwd("/Users/xuebozhao/Documents/LuLab/BottleGourdSpeciation/fsc/model/")
library(RColorBrewer)
library(tidyverse)
library(hrbrthemes)
library(viridis)

# Load likelihood files for each model
early_geneflow <- scan("early/BottleGourd_PopDiv_early.lhoods")
ongoing_geneflow <- scan("ongoing/BottleGourd_PopDiv_ongoing.lhoods")
diff_geneflow <- scan("diff/BottleGourd_PopDiv_diff.lhoods")
recent_geneflow <- scan("recent/BottleGourd_PopDiv_recent.lhoods")
no_geneflow <- scan("no/BottleGourd_PopDiv_no.lhoods")

# Boxplot for model comparison
par(mfrow=c(1,1))
boxplot(range = 0, diff_geneflow, recent_geneflow, early_geneflow, ongoing_geneflow, no_geneflow,
        xlab="Model comparison", ylab="Likelihood", xaxt="n",
        col=brewer.pal(5,"Set2"), border=brewer.pal(5,"Set2"),
        main="D lineage - Bottle gourd")
axis(side=1, at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))

# Alternative visualization using ggplot2
data <- data.frame(
  name=c(rep("Different", length(diff_geneflow)),
         rep("Recent", length(recent_geneflow)),
         rep("Early", length(early_geneflow)),
         rep("Ongoing", length(ongoing_geneflow)),
         rep("No", length(no_geneflow))),
  value=c(diff_geneflow, recent_geneflow, early_geneflow, ongoing_geneflow, no_geneflow)
)

data %>%
  ggplot(aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete=TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(legend.position="none",
        plot.title = element_text(size=11)) +
  ggtitle("Model comparison with Likelihood distributions") +
  xlab("Model comparison")

# ==========================================================
# Section 2: Likelihood Comparison for AB Lineage (Multi-population)
# ==========================================================
setwd("/Users/xuebozhao/Documents/LuLab/BottleGourdSpeciation/fsc/ABlineage/SFS_multi/")
library(RColorBrewer)
par(mfrow=c(2,3), oma=c(4.5,3,0,0), mar=c(2,1.6,2,1), cex=1)

# Example: Loop through subfolders for multiple population pairs
folders <- c("WE_DE_10", "WE_FT_20", "WE_LAN_30", "DE_FT_21", "DE_LAN_31", "FT_LAN_32")

for (f in folders) {
  setwd(paste0("/Users/xuebozhao/Documents/LuLab/BottleGourdSpeciation/fsc/SFS_multi/", f))
  early <- scan("ABlineage_PopDiv_early.lhoods")
  ongoing <- scan("ABlineage_PopDiv_ongoing.lhoods")
  diff <- scan("ABlineage_PopDiv_diff.lhoods")
  recent <- scan("ABlineage_PopDiv_recent.lhoods")
  no <- scan("ABlineage_PopDiv_no.lhoods")

  boxplot(range=0, diff, recent, early, ongoing, no,
          xlab="Model comparison", ylab="Likelihood", xaxt="n",
          col=brewer.pal(5,"Set2"), border=brewer.pal(5,"Set2"),
          main=paste("AB lineage", f))
  axis(side=1, at=1:5, labels=c("Different","Recent","Early","Ongoing","No"))
}

# ==========================================================
# Section 3: Estimate Speciation Times
# ==========================================================
library(ggplot2)
library(dplyr)
library(viridis)

setwd("/Users/xuebozhao/Documents/LuLab/BottleGourdSpeciation/fsc/DAF_SFS_early")
time_data <- read.table("Time_run20.txt")
time_data$year <- time_data$V1 + (14071 - mean(time_data$V1))

data <- data.frame(
  name=c(rep("rep1", length(time_data$V1)), rep("rep2", length(time_data$year))),
  value=c(time_data$V1, time_data$year)
)

# Calculate summary statistics
sample_size <- data %>% group_by(name) %>% summarize(num=n())
data1 <- data %>% left_join(sample_size) %>% mutate(myaxis = paste0(name, "\n", "n=", num))
data1$name <- factor(data1$name, levels=c("rep1","rep2"))
medians <- aggregate(value ~ name, data1, median)

# Plot split time distributions
ggplot(data1, aes(x=name, y=value, fill=name)) +
  geom_violin(width=1.8, color="grey", alpha=0.9) +
  geom_boxplot(width=0.3, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  stat_summary(fun=median, colour="orange", geom="point", shape=15, size=3) +
  ggtitle("Speciation time distribution - Bottle gourd") +
  ylab("Split time (Years)") + xlab("") +
  theme_bw() + theme(panel.grid=element_blank())

# ==========================================================
# Section 4: Estimate Introgression Times
# ==========================================================
setwd("/Users/xuebozhao/Documents/LuLab/BottleGourdSpeciation/fsc/DAF_SFS_early")
intro_D <- read.table("TimeD_introgression_run20.txt")
intro_D$year <- intro_D$V1 + (9729 - mean(intro_D$V1))

setwd("/Users/xuebozhao/Documents/LuLab/BottleGourdSpeciation/fsc/DAF_SFS_early_diff2")
intro_AB_WE <- read.table("TimeAB_WE_LAN_introgression_run20.txt")
intro_AB_WE$year <- intro_AB_WE$V1 + (8919 - mean(intro_AB_WE$V1))
intro_AB_DE <- read.table("TimeAB_DE_LAN_introgression_run20.txt")
intro_AB_DE$year <- intro_AB_DE$V1 + (7228 - mean(intro_AB_DE$V1))
intro_AB_FT <- read.table("TimeAB_FT_LAN_introgression_run20.txt")
intro_AB_FT$year <- intro_AB_FT$V1 + (3203 - mean(intro_AB_FT$V1))

data <- data.frame(
  name=c(rep("D_intro", length(intro_D$year)),
         rep("AB_WE_intro", length(intro_AB_WE$year)),
         rep("AB_DE_intro", length(intro_AB_DE$year)),
         rep("AB_FT_intro", length(intro_AB_FT$year))),
  value=c(intro_D$year, intro_AB_WE$year, intro_AB_DE$year, intro_AB_FT$year)
)

# Plot introgression time comparison
ggplot(data, aes(x=name, y=value, fill=name)) +
  geom_violin(width=1.4, color="grey", alpha=0.9) +
  geom_boxplot(width=0.3, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("Introgression time distribution - Bottle gourd") +
  ylab("Introgression time (Years)") + xlab("") +
  theme_bw() + theme(panel.grid=element_blank())
