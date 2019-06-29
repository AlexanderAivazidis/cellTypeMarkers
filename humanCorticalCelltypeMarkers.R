require(mfishtools)
require(matrixStats)
require(ComplexHeatmap)

panelSize = 250

data = read.delim('data/Allen/human_MTG_2018-06-14_exon+intron_cpm-matrix.csv', sep = ',')
data = data[,2:dim(data)[2]]
rowdata = read.delim('data/Allen/human_MTG_2018-06-14_genes-rows.csv', sep = ',')
coldata = read.delim('data/Allen/human_MTG_2018-06-14_samples-columns.csv', sep = ',')
rownames(data) = rowdata[,1]
rSums = rowSums(data)
data = data[rSums > 10000,]
data = as.matrix(data)
data = log(data+1, 2)

specific_type = coldata[,'cluster']
names(specific_type) = coldata[,1]
data = data[,specific_type != 'no class']
specific_type = specific_type[specific_type != 'no class']

exprThresh = 1
medianExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMedians(data[,x])))
propExpr   = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(data[,x]>exprThresh)))
rownames(medianExpr) <- rownames(propExpr) <- genes <- rownames(data)
nonZeroMedian = (!rowSums(medianExpr) == 0)
medianExpr = medianExpr[nonZeroMedian,]
propExpr = propExpr[nonZeroMedian,]

runGenes <- filterPanelGenes(
  summaryExpr = 2^medianExpr-1,  # medians (could also try means); We enter linear values to match the linear limits below
  propExpr    = propExpr,    # proportions
  startingGenes  = c(),  # Starting genes 
  numBinaryGenes = 1000,      # Number of binary genes 
  minOn     = 10,   # Minimum required expression in highest expressing cell type
  maxOn     = 10^8,  # Maximum allowed expression
  fractionOnClusters = 0.5,  # Max fraction of on clusters 
  excludeFamilies = c("LOC","Fam","RIK","RPS","RPL","\\-","Gm","Rnf","BC0")) # Avoid LOC markers, in this case

corDist         <- function(x) return(as.dist(1-cor(x)))
clusterDistance <- as.matrix(corDist(medianExpr))

# fishPanel = c('RORB', 'CUX1', 'FOXP2')
# runGenes = c(runGenes, fishPanel[2])
# save(fishPanel, file = 'humanCorticalCellMarkerGenes1000.RData')
load('humanCorticalCellMarkerGenes1000.RData')
for (pS in 1:panelSize){
  print(pS)
  
  fishPanel <- buildMappingBasedMarkerPanel(
    mapDat        = data[runGenes,],                # Data for optimization
    medianDat     = medianExpr[runGenes,],            # Median expression levels of relevant genes in relevant clusters
    clustersF     = as.character(specific_type),                         # Vector of cluster assignments
    panelSize     = pS,                               # Final panel size
    currentPanel  = fishPanel,                        # Starting gene panel
    subSamp       = 50,                          # Maximum number of cells per cluster to include in analysis (20-50 is usually best)
    panelMin      = 1,
    optim = 'clusterDistance',
    clusterDistance = clusterDistance,                # Cluster distance matrix (potentiall multiplied by weight matrix)
    percentSubset = 100                               # Only consider a certain percent of genes each iteration to speed up calculations (in most cases this is not recommeded)
  )
  save(fishPanel, file = 'humanCorticalCellMarkerGenes1000.RData')
}

## Plot of median expression of top50 genes
load('data/Allen/cellTaxonomy_humanMTG.RData')
order = cell_taxonomy[[4]][cell_taxonomy[[3]]]
order = order[order != 'no class']
ht1 = Heatmap(medianExpr[fishPanel[1:50], order], name = "log2(counts+1)",
              column_title = 'Top 50 human MTG cell type markers',
              cluster_columns = FALSE, cluster_rows = FALSE)
pdf(file = 'top50humanMTGcellMarkers.pdf', width = 21, height = 14)
ht1
dev.off()

write.csv(fishPanel, file = 'top100humanMTGcellMarkers.csv', quote = F, row.names = F, col.names = F)






