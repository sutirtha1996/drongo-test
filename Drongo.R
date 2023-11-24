set.seed(1234)

library(ggplot2)
library(ggridges)
library(ggfortify)
library(Matrix)
library(grDevices)
library(gridExtra)
library(tidyverse)
library(mosaic)

############################
###########pca##############
############################

drongo<-read.csv("allnotes.csv")
subdrongo<-as.data.frame(drongo[,c(9,10,11,12,13,14,15,16,17)]) #for LDA
subdrongo<-cbind(subdrongo,drongo$Species)
colnames(subdrongo)<-c("BW","D.Time","Peak.freq","PFC.Max","PFC.Min","Peak.Time","Entropy","Start","End","Species")

# Load necessary libraries
library(dplyr)

# Function to calculate mean and standard error for numeric vectors
calculate_summary <- function(x) {
  if (is.numeric(x)) {
    mean_val <- mean(x, na.rm = TRUE)
    se_val <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
    return(paste0(round(mean_val, 2), " ± ", round(se_val, 2)))
  } else {
    return(NA)
  }
}

# Calculate mean and standard error for each trait
summary_table <- subdrongo %>%
  group_by(Species) %>%
  summarise(across(everything(), calculate_summary, .names = "{col}"))

# Print or save the summary table
print(summary_table)
write.csv(summary_table,file = "species_mean_traits.csv",row.names = F)

#summary stats#
peakdrongo<-favstats(~Peak.freq|Species,data = subdrongo)
peakdrongo
write.csv(peakdrongo, file = "speciestraits.csv", row.names = FALSE)

# install.packages("dplyr")
library(dplyr)
library(psych)
library(tidyr)

# summarise
summary_table <- aggregate(. ~ Species, data = subdrongo, function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))

# Print the summary table
print(summary_table)

write.csv(summary_table, file = "speciestraits.csv", row.names = FALSE)

?aov
boxplot(log(Peak.freq)~Species,data = subdrongo)
boxplot(Peak.freq~Species,data = subdrongo)
drongoav<-aov(log(Peak.freq)~Species,data = subdrongo)
drongoav
summary(drongoav)

# Create a boxplot using ggplot2
p1<-ggplot(subdrongo, aes(x = Species, y = Peak.freq,fill=Species))+
  labs(x = "Species", y = "Peak Frequency")+
  scale_y_log10()+
  geom_violin()+
  geom_jitter(size = 0.25, alpha = 0.3)+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_boxplot(width=0.1)
p1

#time
p2<-ggplot(subdrongo, aes(x = Species, y = D.Time,fill=Species))+
  labs(x = "Species", y = "Delta Time")+
  geom_violin()+
  scale_y_log10()+
  geom_jitter(size = 0.25, alpha = 0.3)+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_boxplot(width=0.1)

#entropy
p3<-ggplot(subdrongo, aes(x = Species, y = Entropy,fill=Species))+
  labs(x = "Species", y = "Entropy")+
  geom_violin()+
  scale_y_log10()+
  geom_jitter(size = 0.25, alpha = 0.3)+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_boxplot(width=0.1)
#bandwidth
p4<-ggplot(subdrongo, aes(x = Species, y = BW,fill=Species))+
  labs(x = "Species", y = "BW")+
  geom_violin()+
  scale_y_log10()+
  geom_jitter(size = 0.25, alpha = 0.3)+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_boxplot(width=0.1)

grid.arrange(p1,p2,p3,p4,nrow=2)

#do a pca on correlation matrix#
pc_drongo<-princomp(subdrongo[,1:9],cor = TRUE)
pc_drongo
summary(pc_drongo)
pc_loadings <- pc_drongo$loadings
print(pc_loadings)
write.csv(pc_loadings, "loadings.csv")

pc_scores <- as.data.frame(drongo[, c(20,21,22)])
pc_scores <- cbind(pc_scores, drongo$Species)
colnames(pc_scores) <- c("PC1","PC2","PC3", "Species")
#set the colors for each species#
cols <- setNames(c("#2c7bb6", "#fdae61","#4dac26","#d7191c"), 
                 c("BD", "LTD","HCD","GRT"))

ggplot(pc_scores) +
  aes(PC1, PC2, color = Species) + # define plot area
  geom_point(size = 2) + # adding data points
  scale_color_manual(values = cols)

#convex hull#
data_hulls_drongo <- pc_scores %>%
  group_by(Species) %>%
  do({
    hull_points <- chull(.$PC1, .$PC2)
    data.frame(PC1 = .$PC1[hull_points], PC2 = .$PC2[hull_points])
  })

data_hull_2<-pc_scores %>%
  group_by(Species) %>%
  do({
    hull_points <- chull(.$PC2, .$PC3)
    data.frame(PC2 = .$PC2[hull_points], PC3 = .$PC3[hull_points])
  })
#Create the scatter plot with convex hulls and set colors #
ggplot(pc_scores)+
  aes(x = PC2, y = PC3, color = Species) +
  geom_point(size=1) +
  geom_polygon(data = data_hull_2, aes(x = PC2, y = PC3, fill = Species),
               color = "black", size = 0.2,alpha=0.3) +
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)

#without convex hull#  
ggplot(pc_scores) +
  aes(PC2, PC3, color = Species) + # define plot area
  geom_point(size = 2) +  # adding data points
  scale_color_manual(values = cols)

############################
#######perform an LDA#######
############################
library(MASS)
library(caret)
library(pheatmap)

#lda on the notes, not the PCs#
lda_model2<-lda(Species ~ BW + D.Time + Peak.freq + PFC.Max + PFC.Min + Peak.Time + Entropy + Start + End, 
                  data = subdrongo)
str(lda_model2)

lda_model2
predicted_species <- predict(lda_model2, subdrongo)
str(predicted_species)

conf <- table(list(predicted=predicted_species$class, observed=subdrongo$Species))
conf
conf_matrix<-confusionMatrix(conf)
conf_matrix
# Calculate the confusion matrix as percentages
conf_matrix_percent <- (conf_matrix$table / sum(conf_matrix$table)) * 100
conf_matrix_percent
# Extract the confusion matrix as a matrix
conf_matrix_data <- as.matrix(conf_matrix)

conf_percent<-as.matrix(conf_matrix_percent)
pheatmap(
  conf_matrix_percent,
  main = "Confusion Matrix Heatmap",
  color = colorRampPalette(c("white", "darkgray"))(20),  # Adjust the color palette as needed
  fontsize = 12,
  cellwidth = 50,
  cellheight = 50,
  cluster_rows = FALSE,  # Set to TRUE if you want to cluster rows
  cluster_cols = FALSE,
  display_numbers = T,
)

# Create a control object for 10-fold cross-validation
ctrl <- trainControl(method = "cv", number = 10,classProbs = TRUE,
                     savePredictions = "final")
# Train the LDA model with cross-validation
lda_model <- train(Species ~ BW + D.Time + Peak.freq + PFC.Max + PFC.Min + Peak.Time + Entropy + Start + End, 
                   data = subdrongo, 
                   method = "lda",
                   trControl = ctrl)
confusionMatrix(lda_model) 
print(lda_model)

#randomiztion test#
library(dplyr)  # For data manipulation
library(proxy)  # For calculating Euclidean distances
library(boot)   # For bootstrapping

##########################
####Randomization Test####
##########################

# Expected Distributions

# Extract the "Species" column from the data matrix
Species <- pc_scores[, ncol(pc_scores)]
n <- 1000
avgdist_exp <- rep(NA, n)
# Get unique species values
sp <- unique(Species)
sp_avgdist <- matrix(NA, n, length(sp))
for (i in 1:n) {
  rows <- length(Species)
  row_new <- sample(1:rows, rows)
  sp_rand <- Species[row_new]
  intsp_dist <- matrix(NA, length(sp), length(sp))
  
  for (j in 1:(length(sp) - 1)) {
    ind <- sp_rand == sp[j]
    score_sub <- pc_scores[ind, 1:3]
    
    for (m in (j + 1):length(sp)) {
      ind2 <- sp_rand == sp[m]
      score_sub2 <- pc_scores[ind2, 1:3]
      intsp_dist[j, m] <- mean(dist(rbind(score_sub, score_sub2)))
    }
    
    sp_avgdist[i, j] <- mean(intsp_dist[j, !is.na(intsp_dist[j, ])])
  }
  
  avgdist_exp[i] <- mean(sp_avgdist[i, !is.na(sp_avgdist[i, ])])
}


# Observed
obs_avgdist <- numeric(length(sp))
intob_dist <- matrix(NA, length(sp), length(sp))

for (j in 1:(length(sp) - 1)) {
  ind <- Species == sp[j]
  score_sub <- pc_scores[ind, 1:3]
  
  for (m in (j + 1):length(sp)) {
    ind2 <- Species == sp[m]
    score_sub2 <- pc_scores[ind2, 1:3]
    intob_dist[j, m] <- mean(dist(rbind(score_sub, score_sub2)))
  }
  
  obs_avgdist[j] <- mean(intob_dist[j, !is.na(intob_dist[j, ])])
}
obs_avgdist
#z score#
mean_dist <- mean(avgdist_exp)
std_dist <- sd(avgdist_exp)

z_score <- (obs_avgdist - mean_dist) / std_dist
z_score
install.packages("BSDA")
library(BSDA)

# Define your observed values and the expected mean and standard deviation
observed_values <- obs_avgdist
expected_mean <- mean_dist
expected_sd <- std_dist

# Perform the Z-test
result <- z.test(observed_values, mu = expected_mean, sigma.x = expected_sd)
result
# Extract the results
z_score <- result$statistic
p_value <- result$p.value
conf_interval <- result$conf.int

z_score
p_value
conf_interval


# Create a histogram of expected distances
par(mar=c(5.1,5.1,2.1,2.1))
hist(avgdist_exp, main="Histogram of Expected Distances", xlab="Expected Distance")
# Add a point for the observed interspecific value
points(obs_avgdist,0, col="red", pch=16)

# Determine the range of the data
min_value <- min(avgdist_exp)
max_value <- max(avgdist_exp)

# Determine the amount to extend the x-axis on both sides
extension <- 1  # You can adjust this value as needed

# Create a histogram of the expected distances with an extended x-axis on both sides
hist(avgdist_exp, main="Histogram of Expected Distances", xlab="Expected Distance",
     xlim=c(min_value - extension, max_value + extension))
# Create a histogram of the expected distances with an extended x-axis
hist(avgdist_exp, main="Histogram of Expected Distances", xlab="Expected Distance",
     xlim=c(min(avgdist_exp)-1, max(avgdist_exp)))
# Add a vertical line to the histogram at the specified value
abline(v = 2.4798, col = "red", lwd = 1)

############################
#####Perch Height stats#####
############################

#perch height and peak frequency#
perch<-data.frame(
  Species=c("HCD","BD","LD","GRT"),
  Freq=c(3738.158,3297.116,3055.074,2098.720),
  perch=c(31.5,24,13,12.42))
correlation_coeff <- cor(perch$Freq, perch$perch)
correlation_coeff
cor_test_result <- cor.test(perch$Freq, perch$perch)
cor_test_result
??favstats
favstats(~perch|Species,data=perch)

#Chi-square test for perch heights#

data <- data.frame(
  Species = c("BD", "GRT", "LD", "HCD"),
  Low_Canopy = c(0, 17.85, 11.4, 0),
  Mid_Canopy = c(51, 64.3, 80.2, 5.62),
  Upper_Canopy = c(48, 17.85, 8.34, 94.4)
)
# Perform the chi-square test
result <- chisq.test(data[, 2:4])

# Print the result
print(result)

#song perch only#
songperch<-read.csv("CallOnly.csv")
boxplot(Perch.Height~Species,data = songperch)

ggplot(songperch, aes(x = Species, y = Perch.Height, fill = Species)) +
  geom_boxplot(color = "black", width = 0.7, show.legend = TRUE) +
  scale_fill_manual(values = cols, guide = "legend") +
  labs(x = "Species",
       y = "Perch Height")

perchaov<-aov(Perch.Height~Species,data = songperch)
summary(perchaov)

#linear model#
model <- lm(perch ~ Freq, data = perch)
summary(model)
scatter_plot <- ggplot(perch, aes(x = Freq, y = perch)) +
  geom_point() +  # Add points for the scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a linear regression line
  labs(
    x = "Frequency",
    y = "Perch Height",
    title = "Scatter Plot of Frequency vs. Perch Height"
  )
print(scatter_plot)

#height category
library(dplyr)
cat<-read.csv("height_cat.csv")

# Ensure that Height_Category is a factor with levels ordered appropriately
cat$Height_Category <- factor(cat$Height_Category, levels = c("Upper canopy", "Mid canopy", "Low canopy"))

# Calculate proportions for each combination of Species and Height_Category
df_proportions <- cat %>%
  group_by(Species, Height_Category) %>%
  summarize(count = n()) %>%
  mutate(proportion = count / sum(count))

# Plot the bar graph
ggplot(df_proportions, aes(x = Species, y = proportion, fill = Height_Category)) +
  geom_bar(stat = "identity", position = "stack", reverse = TRUE) +
  labs(x = "Species",
       y = "Proportion") +
  scale_fill_manual(values = c("Low canopy" = "#a6611a", "Mid canopy" = "#dfc27d", "Upper canopy" = "#008837")) +  # Choose colors as needed
  theme_minimal()
