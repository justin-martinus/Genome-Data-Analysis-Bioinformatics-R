#Setup####
#install.packages("BiocManager")
#BiocManager::install(version = "3.17")
#BiocManager::install(c("Biobase", "GEOquery"))

##Library####
library(data.table)
library(Biobase)
library(GEOquery)
library(readr)


##Import####

expdtgeo <- fread("C:/Users/ASUS/Downloads/archive (10)/Breast_GSE45827.csv")

str(expdtgeo)
head(expdtgeo)
#View(expdtgeo)
expdtgeo <- data.frame(expdtgeo)
dim(expdtgeo)

##Transposing and Dividing####

col1 <- colnames(expdtgeo[, -c(1, 2)])

label <- expdtgeo[, c(1, 2)]

label <- t(label)

#View(label)

expdtgeo <- expdtgeo[, -c(1, 2)]

expdtgeo <- as.data.table(expdtgeo)

expdtgeo <- transpose(expdtgeo)

expdtgeo <- data.frame(expdtgeo)

colnames(expdtgeo) <- label[1,]

rownames(expdtgeo) <- col1

#View(expdtgeo)

#Sampling and Filtering####

##Sample 50% of the genes from the expression matrix####
set.seed(2106724832)

row1 <- sample(1:nrow(expdtgeo), 0.5*nrow(expdtgeo))

expdtgeo <- expdtgeo[row1,]

#View(label)

str(expdtgeo)

##Filtering Genes####

library(genefilter)

filt1 <- pOverA(0.25, log2(100))
filt2 <- function(x) (IQR(x) > 0.5)
filtfun <- filterfun(filt1, filt2)
expdtgeo_filtered <- genefilter(expdtgeo, filtfun)
sum(expdtgeo_filtered)

expdtgeo_filtered <- expdtgeo[expdtgeo_filtered,]

str(expdtgeo)

dim(expdtgeo_filtered)

#View(expdtgeo)

par(mfrow=c(1,2))
hist(as.matrix(expdtgeo), main = "original")
hist(as.matrix(expdtgeo_filtered), main="filtered")

#View(label)
str(expdtgeo_filtered)

par(mfrow = c(1, 1))

#Limma Modelling for Comparing 4 Tumor Types####

##Exclude samples that are not part of the 4 categories####
table(label[2,])

library(limma)

cancer <- c("basal", "HER", "luminal_A", "luminal_B")
expdtgeo2 <- expdtgeo[, which(label[2, ] %in% cancer)]

expdtgeo2_filtered <- expdtgeo_filtered[, which(label[2, ] %in% cancer)]

label2 <- label[, which(label[2, ] %in% cancer)]

#View(expdtgeo2)

##Fitting and Showing Top 50 Results####

group <- label2[2, ]
design <- model.matrix(~ group, data=expdtgeo2_filtered)
design

fit <- eBayes(lmFit(expdtgeo2_filtered, design=design))
fit

topResult <- topTable(fit, coef = 2, number = 50)
topResult
library(kableExtra)
knitr::kable(head(topResult, 20), format = "latex", booktabs = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "HOLD_position"))


##Heatmap of each Tumor Type####
# Extract selected genes names
selected <- rownames(expdtgeo2_filtered) %in% rownames(topResult)
selected

# Extract the expression of the selected genes
expdtgeosel <- expdtgeo2_filtered[selected, ]
#View(expdtgeosel)

heatmap(as.matrix(expdtgeosel))

# Heatmap dari kategori grup
expdtgeosel2 <- expdtgeosel
colnames(expdtgeosel2) <- label2[2, ]

# Define the regular expression pattern to match numbers and dots
pattern <- ".*\\b(basal|HER|luminal_A|luminal_B)\\b.*"

# Use gsub() to match and replace the column names
colnames(expdtgeosel2) <- gsub(pattern, "\\1", colnames(expdtgeosel2))
heatmap(as.matrix(expdtgeosel2))

#View(label2)

#View(expdtgeosel2)

##Boxplot for the top 4 genes####
par(mfrow = c(2, 2))
for(i in 1:4){
  temp_data <- data.frame(y = as.vector(t(expdtgeosel2[i, ])), type=label2[2, ])
  temp_data$type <- factor(temp_data$type)
  boxplot(temp_data$y ~ temp_data$type, xlab="Tipe", ylab="EkspresiGen")
}

#str(temp_data)
#View(temp_data)

#colnames(expdtgeosel2)

#View(expdtgeosel2)

par(mfrow = c(1, 1))



#Cancer vs Normal####
##Changing Labels####

table(label[2, ])

expdtgeo3_filtered <- expdtgeo_filtered
label3 <- label
label3[2, ] <- ifelse(label3[2, ] != "normal", "kanker", label3[2, ])

table(label3[2, ])

dim(expdtgeo3_filtered)

##Fitting and Showing Top 50 Results####
group <- label3[2, ]
design <- model.matrix(~ group, data=expdtgeo3_filtered)
design

fit <- eBayes(lmFit(expdtgeo3_filtered, design=design))
fit

topResult2 <- topTable(fit, coef = 2, number = 50)
topResult2
library(kableExtra)
knitr::kable(head(topResult2, 20), format = "latex", booktabs = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "HOLD_position"))

# Extract selected genes names
selected <- rownames(expdtgeo3_filtered) %in% rownames(topResult2)
selected

# Extract the expression of the selected genes
expdtgeosel <- expdtgeo3_filtered[selected, ]
#View(expdtgeosel)

heatmap(as.matrix(expdtgeosel))

##Heatmap of each group####
expdtgeosel2 <- expdtgeosel
colnames(expdtgeosel2) <- label3[2, ]

# Use gsub() to match and replace the column names
colnames(expdtgeosel2) <- gsub(pattern, "\\1", colnames(expdtgeosel2))
heatmap(as.matrix(expdtgeosel2))

#View(label2)

#View(expdtgeosel2)

##Boxplot for the top 4 genes####
par(mfrow = c(2, 2))
for(i in 1:4){
  temp_data <- data.frame(y = as.vector(t(expdtgeosel2[i, ])), type=label3[2, ])
  temp_data$type <- factor(temp_data$type)
  boxplot(temp_data$y ~ temp_data$type, xlab="Tipe", ylab="EkspresiGen")
}



#T-test for Normal vs Cancer####
##The t-test####
vargrp <- label3[2, ]
table(vargrp)

dim(expdtgeo3_filtered)

group <- ifelse(vargrp=="normal", 0, 1)
pval <- apply(as.matrix(expdtgeo3_filtered), 1, function(x) t.test(group, x)$p.value)
dift <- apply(as.matrix(expdtgeo3_filtered), 1, function(x) {
  # Select the columns for the first group (1:86 and 94:151)
  group1 <- c(x[1:86], x[94:151])
  
  # Select the columns for the second group (87:93)
  group2 <- x[87:93]
  
  # Perform the t-test and calculate the difference of means
  diff <- diff(t.test(group1, group2)$estimate)
  return(diff)
})

##Volcano Plot####
par(mfrow = c(1, 1))
library(RColorBrewer)
hist(pval)
sum(pval<0.05)
smoothScatter(dift, -log10(pval))


#Genes' Function Exploration for the 4 Types of Tumor####

library("annotate")
library("hgu133plus2.db")
#View(topResult)

gene_names <- rownames(topResult)
gene_names <- gsub("X", "", gene_names)

GeneSelected <- select(hgu133plus2.db, gene_names,
                       c("SYMBOL", "ENTREZID", "GENENAME"))

GeneSelected <- select(hgu133plus2.db, gene_names, c("SYMBOL", "ENTREZID", "GENENAME", "GO"))
knitr::kable(head(GeneSelected), format = "latex", booktabs = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "HOLD_position"))

# Gene ontology for the top genes
library(GO.db)
GOselected <- select(GO.db, GeneSelected$GO, c("TERM", "GOID"))
# Combine the result
finalres <- cbind(GeneSelected, GOselected)
knitr::kable(head(finalres), format = "latex", booktabs = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "HOLD_position"))

#write.csv2(finalres, file = "C:/Users/ASUS/Downloads/GEOres1.csv")

#Genes' Function Exploration for Healthy vs Cancer####

library("annotate")
library("hgu133plus2.db")
#View(topResult)

gene_names <- rownames(topResult2)
gene_names <- gsub("X", "", gene_names)

GeneSelected <- select(hgu133plus2.db, gene_names,
                       c("SYMBOL", "ENTREZID", "GENENAME"))

GeneSelected <- select(hgu133plus2.db, gene_names, c("SYMBOL", "ENTREZID", "GENENAME", "GO"))
knitr::kable(head(GeneSelected), format = "latex", booktabs = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "HOLD_position"))

# Gene ontology for the top genes
library(GO.db)
GOselected <- select(GO.db, GeneSelected$GO, c("TERM", "GOID"))
# Combine the result
finalres <- cbind(GeneSelected, GOselected)
knitr::kable(head(finalres), format = "latex", booktabs = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "HOLD_position"))


#write.csv2(finalres, file = "C:/Users/ASUS/Downloads/GEOres2.csv")
