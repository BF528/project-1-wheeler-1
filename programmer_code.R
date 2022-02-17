
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))
library(affy)
library(affyPLM)
library(sva)

#Part 3
Data <- ReadAffy(celfile.path='/projectnb/bf528/users/wheeler_2022/project_1/samples')
eset <- rma(Data, background=FALSE)

#Part 4
normalized <-fitPLM(Data,normalize=TRUE,background=TRUE)
nuse_data <- NUSE(normalized,type = 'stats')
rle_data <- RLE(normalized,type = 'stats')
nuse_medians <-nuse_data[1,]
rle_medians <-rle_data[1,]
jpeg(file='nuse_medians.jpeg')
hist(nuse_medians,c ='coral',xlab = 'Median Nuse', main = 'Median of Nuse Values')
dev.off()
jpeg(file='rle_medians.jpeg')
hist(rle_medians,c = 'deeppink4',xlab = 'Median RLE', main = 'Median of RLE Values')
dev.off()

#Part 5
combat_data <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
norm_batch <- as.numeric(factor(combat_data$normalizationcombatbatch)) 
norm_mod <-  as.numeric(factor(combat_data$normalizationcombatmod))
ComBat(eset,batch = norm_batch,mod = norm_mod)
df <- exprs(eset)
write.csv(df,"/projectnb/bf528/users/wheeler_2022/project_1/expression_data.csv", row.names = TRUE)

#Part 6
transposed_scaled_data <- scale(t(df))
scaled_data <- t(transposed_scaled_data)
principal_components <-prcomp(scaled_data,center = FALSE,scale. =FALSE)
values1 <- principal_components$rotation[,1]
values2 <- principal_components$rotation[,2]
jpeg(file='pca_plot.jpeg')
plot(values1,values2, xlab = 'PC1 14.2%',ylab = 'PC2 9.56%', main = 'PCA Plot')
dev.off()

#Part 7
principal_components$importance
summary(principal_components)

