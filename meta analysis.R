#in the name of God 
#Analysis of collected data for resistance genes in MTB 
#set working directory and file import 

setwd("D:/meskini/systematic/analysis")
first_table <- read.csv("Data Extraction 3.csv", header = TRUE, na.strings = "?")
View(first_table)
class(first_table)
nrow(first_table)
ncol(first_table)
dim(first_table)

#packages calling (each must hired before need (not at first of run))
install.packages(c("pecman", "ggplot2", "metafor","tidyverse", "meta", "maps"))
library(pecman)
library(ggplot2)
library(metafor)
library(tidyverse)
library(dplyr)
library(meta)
library(MASS)


#tables preparation 
#table 1: genes distribution (forest proportion)
TP1 <- first_table[,c("Gene", "sample.size..Number.of.people.",  "mutation")]
View(TP1) 
#proportion analysis 
library(readxl)
 Proportion_of_each_gene <- read_excel("D:/meskini/systematic/analysis/data/Proportion of each gene.xlsx")
 View(Proportion_of_each_gene)
 
 library(meta)
 proportion_analysis <- data.frame(  
   genename <- Proportion_of_each_gene$Gene,  
   event <- Proportion_of_each_gene$Event,  
   n <- Proportion_of_each_gene$Total  
 )  
 meta_analysis <- metaprop(  
   event <- proportion_analysis$event....Proportion_of_each_gene.Event,  
   n = proportion_analysis$n....Proportion_of_each_gene.Total,  
   studlab = proportion_analysis$genename....Proportion_of_each_gene.Gene,  
   sm = "PLOGIT",   
 )  
 summary(meta_analysis)
 print(meta_analysis$Q)   
 print(meta_analysis$p.value.Q)   
 
 
library(meta)
forest(meta_analysis, pooled.events = TRUE, smlab = "Proportion", col.square = "red", col.diamond = "blue", col.common = "black") 


#each study details 

study_data <- data.frame(
  study = Proportion_of_each_gene$Gene,
  total = Proportion_of_each_gene$Total,
  events = Proportion_of_each_gene$Event
)
view(study_data)
study_data$proportion <- study_data$events / study_data$total
study_data$se <- sqrt((study_data$proportion * (1 - study_data$proportion)) / study_data$total)


meta_analysis <- metaprop(
  event = study_data$events,
  n = study_data$total,
  studlab = study_data$study,
  method = "Inverse",
  sm = "PLOGIT",
  data = study_data
)


print(meta_analysis)
forest(meta_analysis, random = TRUE, pooled.events = TRUE, smlab = "Proportion", col.square = "red", col.diamond = "blue", col.common = "black") 
library(metafor)

study_data$yi <- escalc(measure = "PLO", xi = study_data$events, ni = study_data$total)$yi
study_data$vi <- escalc(measure = "PLO", xi = study_data$events, ni = study_data$total)$vi

res <- rma(yi = study_data$yi, vi = study_data$vi, data = study_data, method = "REML")
rownames(study_data) <- study_data$study
print(res)

loo_results <- leave1out(res)

print(loo_results)




#prevalence of method utilization 
library(readxl)
pie_chart_for_methods <- read_excel("D:/meskini/systematic/analysis/data/pie chart for methods.xlsx")
View(pie_chart_for_methods)

library(plotrix)  

slices <- pie_chart_for_methods$Prevalence
labels <- pie_chart_for_methods$Labels
colors <- pie_chart_for_methods$color

pie3D(slices, explode = 0.1, main = "Method of Sequencing", col = colors, theta = 1, srt = 45)  
legend("bottomright", legend = labels, fill = colors, cex = 0.2) 


#percent 
data <- pie_chart_for_methods$Prevalence
total <- 96
pie_chart_for_methods$percent <- (data/total)*100
View(pie_chart_for_methods)

write.table(pie_chart_for_methods)
#world map distribution 
my.data <- read.table()
install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel","ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library(maps)
library(ggplot2)
setwd("D:/meskini/systematic/analysis/data/world_map")
dir()
x <- read.csv("ethA (Ethionamide and Protionamide).csv")
View(x)


data <- data.frame(
  country = x$Country,
  prevalence = x$Proportion
)

world_map <- map_data("world")
View(world_map)
merged_data <- merge(world_map, data, by.x = "region", by.y = "country", all.x = TRUE)

ggplot() +
    geom_map(data = merged_data, map = world_map,
                         aes(map_id = region, fill = prevalence),
                        color = "white") +
    expand_limits(x = world_map$long, y = world_map$lat) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Prevalence") +
   labs(title = "") +
    theme_minimal()


#time calculation 
# Load the ggplot2 package
library(ggplot2)

setwd("D:/meskini/systematic/analysis/data/time")
dir()
xtime <- read.csv("tlyA.csv", na.strings = "")
View(xtime)
data <- data.frame(
  time = xtime$Year, 
  value =  xtime$Proportion
)


ggplot(data, aes(x = time, y = value)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") + 
  labs(title = "",
       x = "Time",
       y = "Proportion") +
  theme_minimal() 




#publication bias caluculating
install.packages("effectsize")
install.packages("datawizard")



library(effectsize) # For calculating effect size  
library(meta) # For meta-analysis which contains Begg's test  

# Sample data  
data <- c(50, 52, 48, 49, 51, 54, 53, 50, 49)  

# Calculate mean and standard deviation  
mean_data <- mean(data)  
sd_data <- sd(data)  
n <- length(data)  

# Calculate effect size (Cohen's d)  
effect_size <- (mean_data - 0) / sd_data  

# Calculate standard error  
standard_error <- sd_data / sqrt(n)  

# Print results  
cat("Effect Size (Cohen's d):", effect_size, "\n")  
cat("Standard Error:", standard_error, "\n")  


 
study_data <- data.frame(  
  study = 1:8,  
  m = mean_data,  
  sd = sd_data,  
  n = n  
)  


###############################
setwd("D:/meskini/systematic/analysis/data/publication bias")
dir()
bias <- read.csv("gyrA.csv")
view(bias)

study_data <- data.frame(
  study = bias$Gene,
  total = bias$Total,
  events = bias$Mutation
)
study_data$proportion <- study_data$events / study_data$total
study_data$sd <- sd(study_data$proportion)
study_data$sei <- study_data$sd/sqrt(length(study_data$proportion))
study_data$se <- sqrt((study_data$proportion * (1 - study_data$proportion)) / study_data$total)
view(study_data)
meta_analysis <- metacont(n, m, sd, data = study_data, sm = "SMD")  
begg_test <- metabias(meta_analysis)  


sd(x$Efficacy)

print(begg_test)  


res <- rma(yi = study_data$proportion, sei = study_data$se)  
view(res)


egger_test <- regtest(res, model = "lm")  
print(egger_test)  

funnel(res)

##########################################
setwd("D:/meskini/systematic/analysis/data/publication bias")
dir()
library(readxl)
pub <- read_excel("D:/meskini/systematic/analysis/data/publication bias/rrs (other).xlsx")

View(pub)

#hemogenity 
pub$proportion <- pub$Mutation/pub$Total
pub$sd <- sd(pub$proportion)
pub$se <- pub$sd / sqrt (length(pub$proportion))
pub$mean <- mean(pub$proportion)  
pub$median <- median(pub$proportion)
pub$mode <- mode(pub$proportion)
pub$max <- max(pub$proportion)
pub$min <- min(pub$proportion)
View(pub)


#heterogenity 







######  publication bias 
library(readxl)
PB <- read_excel("D:/meskini/systematic/analysis/data/publication bias/ahpC (isoniazid).xlsx")
View(PB)

library(metafor)
PB$Proportion <- PB$Mutation/PB$Total
PB$yi <- PB$Proportion
PB$vi <- escalc(measure = "PLO", xi = PB$Mutation, ni = PB$Total)$vi
PB$sd <- sd(PB$Proportion)
PB$mean <- mean(PB$Proportion)

library(meta)
library(metafor)
meta_analysis <- rma(yi = PB$yi, sei = PB$vi, method = "REML")
egger_test <- regtest(meta_analysis)

begg_test <- metabias(meta_analysis, method.bias = "Begg")

print(egger_test)


print(begg_test)


funnel(meta_analysis, xlab = "Proportion", ylab = "Standard Error",
       main = "Funnel Plot")





######Drug 

library(readxl)
 DR <- read_excel("D:/meskini/systematic/analysis/data/results/test assesment/Pyrazinamide.xlsx")
 View(DR)
DR$Var <- DR$SE^2*length(DR$Proportion)


library(metafor)
res <- rma(yi = DR$Proportion, sei = DR$SE, method = "REML")


DR$Q <- res$QE
DR$p_value_Q <- res$QEp



DR$df <- length(DR$Gene) - 1

# Calculate I^2
DR$I2 <- ((DR$Q - df) / DR$Q) * 100



DR$tau2 <- res$tau2

DR$mean_sample <- mean(DR$Mutation)
DR$mean_Total <- mean(DR$Total)
DR$SD_TOTAL <- sd(DR$Total)
DR$leng <- length(DR$Proportion)
sample_mean <- DR$mean_sample
population_mean <- DR$mean_Total
population_sd <- DR$SD_TOTAL
sample_size <- length(DR$Proportion)

# Calculate the z-value
DR$z_value <- (sample_mean - population_mean) / (population_sd / sqrt(sample_size))


DR$totall <- sum(DR$Total)
DR$Mutations <- sum(DR$Mutation)
DR$totprop <- DR$Mutations/DR$totall

t.test(DR$Proportion)
res




#publication bias 
library(readxl)
 Gene_studies <- read_excel("D:/meskini/systematic/analysis/data/publication bias/RPLC.xlsx")
 View(Gene_studies)

 number_of_successes <- Gene_studies$Mutation
 total_trials <- Gene_studies$Total  
 
Gene_studies$proportion <- number_of_successes/total_trials
 
 # Calculate standard error  
 Gene_studies$standard_error <- sqrt((Gene_studies$proportion * (1 - Gene_studies$proportion)) / total_trials)  
 
 # Define critical value (for 95% confidence, use 1.96)  
 z_score <- 1.96  
 

 
 library(meta)
 library(metafor)
 meta_analysis <- rma(yi = Gene_studies$proportion, sei = Gene_studies$standard_error, method = "REML")
 egger_test <- regtest(meta_analysis)
 
View(Gene_studies)
print(egger_test)
 
funnel(meta_analysis, xlab = "Proportion", ylab = "Standard Error",
        main = "Funnel Plot")
begg_test <- ranktest(meta_analysis)  
begg_test
 
## statistics for each gene
res <- rma(yi = Gene_studies$proportion, sei = Gene_studies$standard_error, method = "REML")
res

##### test fore correctness

margin_of_error = z_score * ((Gene_studies$proportion * (1 - Gene_studies$proportion) / total_trials) ** 0.5)

# Calculate the lower and upper limits of the confidence interval
lower_limit = Gene_studies$proportion - margin_of_error
upper_limit = Gene_studies$proportion + margin_of_error

 
 # Z-score for 95% confidence interval
 z_value <- 1.96
 
 # Calculate standard errors from confidence intervals
 standard_errors_ci <- (upper_limit - lower_limit) / (2 * z_score)
 
 # Print standard errors
 print("Standard Errors from Confidence Intervals:")
 print(standard_errors_ci)
 
 
 
 
 ##### continent analysis 
 
 library(readxl)
 pub <- read_excel("D:/meskini/systematic/analysis/data/continent/tlyA/africa.xlsx")
 
View(pub)

 pub$proportion <- pub$Prevalence/pub$Total
 pub$sd <- sd(pub$proportion)
 pub$se <- pub$sd / sqrt (length(pub$proportion))
 pub$mean <- mean(pub$proportion)  
 pub$median <- median(pub$proportion)
 pub$mode <- mode(pub$proportion)
 pub$max <- max(pub$proportion)
 pub$min <- min(pub$proportion)
 View(pub)
 
 
 
 ########################### income analysis 
 
 
 library(readxl)
 pub <- read_excel("D:/meskini/systematic/analysis/data/income/low.xlsx")
 
 View(pub)
 
 pub$proportion <- pub$positive/pub$total
 pub$sd <- sd(pub$proportion)
 pub$se <- pub$sd / sqrt (length(pub$proportion))
 pub$mean <- mean(pub$proportion)  
 pub$median <- median(pub$proportion)
 pub$mode <- mode(pub$proportion)
 pub$max <- max(pub$proportion)
 pub$min <- min(pub$proportion)
 View(pub)
 