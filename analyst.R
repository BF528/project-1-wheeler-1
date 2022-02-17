library(tidyverse)



expression_data <- read.csv("/projectnb/bf528/users/wheeler_2022/project_1/expression_data.csv",sep = ",",header = TRUE,row.names = 1)






##Noise Filtering and Dimensionality Reduction
#4.1:-For each gene(rows), at least 20% of the gene-expression values must be > log2(15))

threshold_value <- log2(15)
filter1_expressed <- rowSums(expression_data > threshold_value) >= (0.2*ncol(expression_data))
filter1_set<- expression_data[filter1_expressed, ]
# Number of genes left 
print(c("Number of genes left after filter 1:",nrow(filter1_set)))



#4.2:- variance significantly different from the median variance of all probe sets using a threshold of 𝑝<0.01
#row variance
row_var<- apply(filter1_set, 1, var)
#median
med_var <- median(row_var)
dof<- ncol(expression_data) - 1

# chi sq test 
filter2_chi <- filter1_set[((dof*row_var)/med_var)>qchisq(((1 - 0.99)/2),dof,lower.tail = FALSE),]

# Number of genes left
print(c("Number of genes left after filter2:",nrow(filter2_chi )))


#4.3:- Have a coefficient of variation > 0.186.

filter3_var <-subset(filter2_chi, apply(filter2_chi, 1, function(x) sd(x)/mean(x)) > 0.186)
# Number of genes left
print(c("Number of genes left after filter3:", nrow(filter3_var)))

#4.4:-Write out a different file containing the gene expression matrix for genes passing all three of the filters from 4.1, 4.2, and 4.3.
write.csv(filter1_set, sep=',', file="data4.1.csv")
write.csv(filter2_chi, sep=',', file="data4.2.csv")
write.csv(filter3_var, sep=',', file="data4.3.csv")



###############################################################################
##Hierarchical clustering & subtype discovery


#5.1:- hierarchical clustering
h_cluster<-hclust(dist(t(filter3_var)))
plot(h_cluster,main = "HCluster Dendrogram")


#5.2:- Cut the dendrogram such that the samples are divided into 2 clusters
clusters<-cutree(h_cluster, 2)



#Number of samples are in each cluster
print(c("Number of samples in cluster1:", sum(clusters==1)))
print(c("Number of samples in cluster2:", sum(clusters==2)))

#5.3:- heatmap of the gene-expression of each gene across all samples
proj_metadata<- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
proj_metadata<- subset(proj_metadata, select=c(geo_accession, cit.coloncancermolecularsubtype))


color<-ifelse(proj_metadata$cit.coloncancermolecularsubtype == 'C3','red', 'blue')
print(color)
heatmap(as.matrix(filter3_var),ColSideColors = color, main = 'Gene Expression/Samples')
legend(x = "bottomright", inset=c(-0.70,0),
       legend=c('C3', 'OTHER'), title='Subtype', fill=c('red', 'blue'))
# Save heatmap as a file
png('/projectnb/bf528/users/wheeler_2022/project_1/heatmap_result.png')


#5.4:- Using the expression matrix from Part 4.4 and the cluster memberships from Part 5.2, identify genes differentially expressed between the two clusters using a Welch t-test (results in a ranked list). Write out a dataframe containing the probeset ID, t-statistic, p-value, and adjusted p-value (i.e. FDR, see the p.adjust function) to a comma separated file for each comparison. How many genes are differentially expressed at adjusted p<0.05 between the clusters for both lists?
#separate data into cluster1 and cluster2 for t-test  
cluster1<-filter3_var[, clusters==1]
cluster2<-filter3_var[, clusters==2]
plot(h_cluster, hang = -1, cex = 0.8, xlab = "GEO Accession")
rect.hclust(h_cluster, k = 2, border = "red")



#5.5:- Welch t-test 

welch_ttest <-  apply(as.matrix(filter3_var),MARGIN=1,function(x) t.test(x=x[clusters==1],y=x[clusters==2]))

p_value <- c(sapply(welch_ttest,function(x) x$p.value))
t_statistic_val <- c(sapply(welch_ttest,function(x) x$statistic))
pValue_adjusted <- c(p.adjust(p_value ,method = "fdr"))
Probeset_ID = c(row.names(filter3_var))
clustered_data <- data.frame(Probeset_ID, t_statistic_val, p_value, pValue_adjusted)



#differentially expressed at adjusted 𝑝<0.05 between the clusters for both lists
diff_exp<- length(clustered_data$pValue_adjusted < 0.05)
print(diff_exp)
write.csv(clustered_data,sep = ",",row.names = FALSE,file = "data5.5.csv")

#5.6:- the t-test analysis described in 5.4 on the expression matrix from 4.5 
ttest_bio <- apply(as.matrix(filter2_chi),MARGIN=1,function(x) t.test(x=x[clusters==1],y=x[clusters==2]))
pvals_bio <- sapply(ttest_bio,function(x) x$p.value)
tstats_bio <- sapply(ttest_bio,function(x) x$statistic)
p_adjusted_bio <- p.adjust(pvals_bio,"fdr")
biologist_data <- data.frame("Probeset_ID" = c(row.names(filter2_chi)),
                             tstats_bio,pvals_bio,p_adjusted_bio)


write.csv(biologist_data ,row.names = F,file = "data5.6.csv",sep=",")

