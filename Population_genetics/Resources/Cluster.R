# ==========================================================
# Title: STRUCTURE Visualization and Population Mapping in Cucumber
# Author: Xuebo Zhao
# Description:
#   This script performs STRUCTURE/ADMIXTURE Q-matrix visualization,
#   geographic clustering, and population mapping using pie charts.
# ==========================================================

# Load required packages
library(pheatmap)
require(reshape)
require(rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)
require(gridExtra)
library(cluster)
library(factoextra)
require("RColorBrewer")

# ==========================================================
# 1️ Set Working Directory and Load STRUCTURE Q-matrix files
# ==========================================================
setwd("/Users/xuebozhao/Documents/CucumberProject/Structure_map/Select/")

# Read Q-matrix files from STRUCTURE/ADMIXTURE outputs
sfiles <- list.files(path="./Select_Qmatrix", full.names=TRUE)
slist <- readQ(files=sfiles)

# Summarize Q-matrix
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))

# ==========================================================
# 2️ Plot STRUCTURE results using pophelper
# ==========================================================
# Example 1: simple plot
p1 <- plotQ(slist[3:7], returnplot=TRUE, exportplot=FALSE, quiet=TRUE, basesize=11,
            showindlab=TRUE, showyaxis=TRUE, showticks=TRUE)

# Example 2: sorted by all clusters
p1 <- plotQ(slist[3:7], returnplot=TRUE, exportplot=FALSE, quiet=TRUE, basesize=11,
            sortind="all", showindlab=FALSE, showyaxis=TRUE, showticks=TRUE,
            sharedindlab=TRUE, clustercol=brewer.pal(6, "Set2"))

# Example 3: sorted by cluster 1
p1 <- plotQ(slist[3:10], returnplot=TRUE, exportplot=FALSE, quiet=TRUE, basesize=11,
            sortind="Cluster1", showindlab=TRUE, showyaxis=TRUE,
            showticks=TRUE, sharedindlab=TRUE)

# Combine multiple STRUCTURE bar plots
grid.arrange(p1$plot[[1]], p1$plot[[2]], p1$plot[[3]], p1$plot[[4]], p1$plot[[5]], nrow=5)

# ==========================================================
# 3️ Clustering samples by latitude/longitude
# ==========================================================
data <- read.table("land.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dataA <- data[, c(3,4)]  # Latitude, Longitude
data2 <- dataA[!is.na(dataA$Latitude),]
df <- scale(data2)

# Calculate pairwise distance and hierarchical clustering
result <- dist(df, method="euclidean")
result_hc <- hclust(d=result, method="ward.D2")

# Cut tree into 40 clusters
data2$type <- cutree(result_hc, k=40)

# Compute mean latitude and longitude for each cluster
lat_mean <- tapply(data2[,1], data2$type, mean, na.rm=TRUE)
lon_mean <- tapply(data2[,2], data2$type, mean, na.rm=TRUE)

# Assign cluster mean coordinates
data2$cluster1 <- NA
data2$cluster2 <- NA
for (i in 1:nrow(data2)) {
  for (j in 1:40) {
    if (data2[i,3] == j) {
      data2[i,4] <- as.numeric(lat_mean[j])
      data2[i,5] <- as.numeric(lon_mean[j])
    }
  }
}

write.table(data2, "cucumber_clusters.txt", sep="\t", row.names=FALSE, quote=FALSE)

# ==========================================================
# 4️⃣ Map visualization with STRUCTURE pie charts
# ==========================================================
data <- read.table("select_land.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dataA <- data[, c(1,7,8,24,25)]
colnames(dataA)[1:5] <- c("type","Latitude", "Longitude", "Sub.population.6", "value")

# Remove NA rows
cucumber <- dataA[!is.na(dataA$Latitude) & !is.na(dataA$Sub.population.6),]
cucumber_reshape <- cast(cucumber, Latitude + Longitude ~ Sub.population.6)
cucumber_reshape2 <- as.data.frame(cucumber_reshape)

# ==========================================================
# 5️ Draw STRUCTURE pie maps for different K values
# ==========================================================
# Example for K=8
mapPies(cucumber_reshape2, xlim=c(-120,140), ylim=c(35,40),
        nameX="Longitude", nameY="Latitude",
        nameZs=c('1','2','3','4','5','6','7','8'),
        symbolSize=1, zColours=brewer.pal(8, "Set2"),
        barOrient='vert', oceanCol="#D1EEEE",
        landCol="#FFDAB9", main="Cucumber Population (K=8)")

# Example for K=6
mapPies(cucumber_reshape2, xlim=c(-120,140), ylim=c(35,40),
        nameX="Longitude", nameY="Latitude",
        nameZs=c('1','2','3','4','5','6'),
        symbolSize=1, zColours=brewer.pal(6, "Set2"),
        barOrient='vert', oceanCol="#D1EEEE",
        landCol="#FFDAB9", main="Cucumber Population (K=6)")

# Example for K=4
mapPies(cucumber_reshape2, xlim=c(-120,140), ylim=c(35,36),
        nameX="Longitude", nameY="Latitude",
        nameZs=c('1','2','3','4'),
        symbolSize=1, zColours=brewer.pal(4, "Set2"),
        barOrient='vert', oceanCol="#D1EEEE",
        landCol="#FFDAB9", main="Cucumber Population (K=4)")

# ==========================================================
# 6️ Optional: Combine multiple maps
# ==========================================================
par(mfcol=c(5,1))
