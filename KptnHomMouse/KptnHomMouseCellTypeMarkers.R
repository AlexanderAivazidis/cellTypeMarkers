### Kptn Mouse Cell type marker selection

require(ComplexHeatmap)

counts = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_counts.rds')
coldata = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_coldata.rds')
specific_type = as.character(coldata[,2])
names(specific_type) = colnames(counts)
meanExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(counts[,x])))

catalog<-fread("/home/jovyan/cellTypeMarkers/RNAscope Probe List December2018_Leica.csv")

catalog.species<-subset(catalog, `Species Name` %in%  'Mus musculus')
catalog.species<-subset(catalog.species, `Product Type` %in%  'RNAscope LS 2.5 Target Probes')

cellGroupsUp = c('Neuron.Gad1Gad2.Drd1-Nefm', 'Neuron.Gad1Gad2.Adora2a-Nefm', 'Neuron.Gad1Gad2.Drd1-Fos')
cellGroupsDown = c('Neuron.Gad1Gad2.Pvalb-Rgs12','Neuron.Gad1Gad2.Pnoc')

genes = c('Drd1', 'Nefm', 'Fos', 'Rgs12', 'Pnoc', 'Adora2a', 'Pvalb')
genes2 = c('Drd1', 'Fos', 'Pvalb', 'Rgs12', 'Pnoc')

probes.selected<-subset(catalog.species, `Gene Symbol` %in% genes)

allGenes = sort(catalog.species$`Gene Symbol`)

ht1 = Heatmap(meanExpr[genes2,c(cellGroupsUp, cellGroupsDown, colnames(meanExpr)[!colnames(meanExpr) %in% c(cellGroupsDown, cellGroupsUp)])],
              name = "Mean counts",
              column_title = 'Cell type markers',
              cluster_columns = FALSE, cluster_rows = FALSE)

### Chosen probes: for marking Neuron.Gad1Gad2.Drd1-Fos and Neuron.Gad1Gad2.Pvalb-Rgs12

# 517708      RNAscope® 2.5 LS Probe - Mm-Rgs12-XHs RNAscope LS 2.5 Target Probes       1 Mus musculus        Rgs12
# 406498-C2   RNAscope® 2.5 LS Probe - Mm-Drd1a-C2 RNAscope LS 2.5 Target Probes        2 Mus musculus        Drd1
# 316928-C3   RNAscope® 2.5 LS Probe - Mm-Fos-C3 RNAscope LS 2.5 Target Probes          3 Mus musculus        Fos
# 421938-C4   RNAscope® 2.5 LS Probe - Mm-Pvalb-C4 RNAscope LS 2.5 Target Probes        4 Mus musculus        Pvalb


