library("AnnotationDbi")
library("hgu133plus2.db")
library("GSEABase")


# Repeating the clustering task here,
#   since, though unmentioned, we need that to calculate log2 fold change.

exp_df = read.csv("expression_data.csv", header=TRUE, fill=TRUE, row.names=1)

row_var = apply(exp_df, 1, var)
med_var = median(row_var)
dof = ncol(exp_df) - 1
exp_chi = exp_df[(dof * row_var / med_var) > qchisq(0.01, dof, lower.tail = F),]

rm('row_var', 'med_var', 'dof', 'exp_df')

clustering = hclust(dist(t(exp_chi)))
clusters = cutree(clustering, 2)
cluster1 = exp_chi[, clusters==1]
cluster2 = exp_chi[, clusters==2]
l2fc = apply(log2(cluster2), 1, mean) - apply(log2(cluster1), 1, mean)

rm('exp_chi', 'clustering', 'clusters', 'cluster1', 'cluster2')


# Part 1, loading differential expression results

dif_df = read.csv("data5.6.csv", header=TRUE, fill=TRUE)
mapping = select(hgu133plus2.db, keys=as.character(dif_df$Probeset_ID), columns="SYMBOL", keyset="PROBEID")
gs_df = merge(x=mapping, y=dif_df, by.x="PROBEID", by.y="Probeset_ID")

# Part 2 doesn't happeen in R


# Part 3, assessing top and bottom differentially expressed genes.
gs_df$l2fc = l2fc[as.character(gs_df$PROBEID)]

rm('dif_df', 'mapping', 'l2fc')

gs_df = gs_df[order(-gs_df$l2fc, -gs_df$p_adjusted_bio),]
gs_df = gs_df[!duplicated(gs_df$SYMBOL),]

print(head(gs_df, 10))
print(tail(gs_df, 10))

topk = head(gs_df, 1000)
botk = tail(gs_df, 1000)


# Part 4, GSEA, unexplained, but I'll figure it out
hallcol = getGmt('h.all.v7.5.1.symbols.gmt')
keggcol = getGmt('c2.cp.kegg.v7.5.1.symbols.gmt')
gocol = getGmt('c5.go.v7.5.1.symbols.gmt')

print(paste('Hallmark Collection Size: ', length(hallcol)))
print(paste('KEGG Collection Size: ', length(keggcol)))
print(paste('GO Collection Size: ', length(gocol)))

# Part 5, Fisher Test

N = nrow(gs_df)

fisher_sum = function(expressed, gene_set) {
  k = length(intersect(gene_set@geneIds, expressed))
  q = length(gene_set@geneIds)-k
  ft = fisher.test(matrix(c(k, 1000-k, q, N - 1000 - q), nrow=2), alternative='greater')
  return(c(unname(ft$estimate), unname(ft$p)))
}

gatherem = function(expressed, gene_col) {
  dfo = data.frame(geneset=character(),
                   hyperg_stat=double(),
                   p_nom=double(),
                   stringsAsFactors=FALSE) 
  for (gene_set_name in names(gene_col)) {
    ft = fisher_sum(as.character(expressed$SYMBOL), gene_col[[gene_set_name]])
    nextRow = nrow(dfo) +1
    dfo[nextRow,1] = gene_set_name
    dfo[nextRow,2:3] = ft
  }
  dfo$p_adj = p.adjust(dfo$p_nom, 'fdr')
  return(dfo)
}

HallTopDF = gatherem(topk, hallcol)
KeggTopDF = gatherem(topk, keggcol)
GoTopDF = gatherem(topk, gocol)

print(paste('Hallmark Significantly Enriched Up: ', sum(HallTopDF$p_adj < 0.05)))
print(head(HallTopDF[order(HallTopDF$p_nom),], 50))

print(paste('KEGG Significantly Enriched Up: ', sum(KeggTopDF$p_adj < 0.05)))
print(head(KeggTopDF[order(KeggTopDF$p_nom),], 3))

print(paste('GO Significantly Enriched Up: ', sum(GoTopDF$p_adj < 0.05)))
print(head(GoTopDF[order(GoTopDF$p_nom),], 3))


HallBotDF = gatherem(botk, hallcol)
KeggBotDF = gatherem(botk, keggcol)
GoBotDF = gatherem(botk, gocol)

print(paste('Hallmark Significantly Enriched Down: ', sum(HallBotDF$p_adj < 0.05)))
print(head(HallBotDF[order(HallBotDF$p_nom),], 3))

print(paste('KEGG Significantly Enriched Down: ', sum(KeggBotDF$p_adj < 0.05)))
print(head(KeggBotDF[order(KeggBotDF$p_nom),], 3))

print(paste('GO Significantly Enriched Down: ', sum(GoBotDF$p_adj < 0.05)))
print(head(GoBotDF[order(GoBotDF$p_nom),], 3))