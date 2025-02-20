#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("clusterProfiler")
#BiocManager::install("DOSE")
#BiocManager::install("pathview")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("GOfuncR")
#install.packages("tidyverse")
#install.packages("tidyr")
#install.packages("pathfindR")
#install.packages("igraph")

#Requied libraries
library(ggplot2)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(tidyr)
library(clusterProfiler)
library(DOSE)
library(pathview)
library(EnsDb.Hsapiens.v86)
library(pathfindR)
library(GOfuncR)
library(igraph)
library(ReactomePA)
library(msigdbr)


#Load the R input file generated by the python script

# args <- commandArgs(trailingOnly = TRUE)
# csv_file <- args[1]
# comparisons <- read.csv(file=csv_file)

comparisons <- read.csv(file="/Users/juansola/Desktop/240930_TCGA-LIHC_mapping.csv",header = TRUE)

for (x in 1:nrow(comparisons)){

name=comparisons[x,1]
name
name=strsplit(name, split = "/")
name=strsplit(name[[1]][length(name[[1]])], split = "_RAW.csv")


#change X with the complete path to your raw count file
raw_count_file=comparisons[x,1]

#change X with the complete path to your condition file
condition_file=comparisons[x,2]

#specify where to save the output
output_path = comparisons[x,3]

#name of the output files
name=name[[1]][1]

#cutoff parameters
pvalue_cut=0.05
folchange=1


#Load the conditions
condition <- read.csv(file=condition_file, row.names ="samples")
#Load the read raw count 
raw_count <- as.matrix(read.csv(file=raw_count_file,check.names = FALSE, row.names=1))
raw_count=t(raw_count)
raw_count=raw_count[,rownames(condition)]




#run Diferrentail expression analysis 
dea <- DESeqDataSetFromMatrix(countData =raw_count , colData = condition,design = ~ condition)
dea$condition <- relevel(dea$condition, ref = "Control")
dea <- DESeq(dea)
res <- results(dea)
resOrdered <- res[order(res$padj),]
resOrdered$gene <- row.names(resOrdered)
rownames(resOrdered) <- NULL
sublog <-resOrdered[!is.na(resOrdered$log2FoldChange) & ((resOrdered$log2FoldChange >= folchange) | (resOrdered$log2FoldChange <= -folchange)),]
resOrdered$gene <- sapply(strsplit(resOrdered$gene, "[.]"), `[`, 1) # delet decimal in ensmbl code
resSig <- subset(sublog, padj < pvalue_cut) #select genes with a p-value less than 0.05
resSig$gene <- sapply(strsplit(resSig$gene, "[.]"), `[`, 1) # delet decimal in ensmbl code

resultsNames(dea)

#save DEA results
deseqname=paste(name,"deseq",sep="_")
deseq_folder=paste(output_path,deseqname,sep="/")
dir.create(paste(output_path,deseqname,sep="/"))
dea_name=paste(name,"p-value-005_log1.csv",sep="_")
write.csv(as.data.frame(resSig), file=paste(deseq_folder,dea_name,sep="/"),row.names = FALSE)
dea_name=paste(name,"not_cutoff.csv",sep="_")
write.csv(as.data.frame(resOrdered), file=paste(deseq_folder,dea_name,sep ="/"),row.names = FALSE)


#Volcano Plot

volplotdata <- read.csv(file=paste(deseq_folder,dea_name,sep ="/"),)
volplotdata <- volplotdata[!duplicated(volplotdata$gene), ]
rownames(volplotdata) <- volplotdata$gene
volplotdata$Expression <- "NO"
volplotdata$Expression[volplotdata$log2FoldChange > folchange & volplotdata$padj < pvalue_cut] <- "UP"
volplotdata$Expression[volplotdata$log2FoldChange < -folchange & volplotdata$padj < pvalue_cut] <- "Down"


ggplot(data=volplotdata, aes(x=log2FoldChange, y=-log10(padj), col=Expression)) + geom_point() + theme_minimal() + geom_vline(xintercept=c(-folchange, folchange), col="red") +
  geom_hline(yintercept=-log10(pvalue_cut), col="red")+scale_color_manual(values=c("red", "black", "blue"))+ ggtitle("Deseq result") + theme(plot.title = element_text(hjust = 0.5))

volcano_name=paste(name,"volcano_plot.pdf",sep="_")
ggsave(paste(deseq_folder,volcano_name,sep="/"))



######################################################
######################################################
######################################################
######################################################




#functional enrichment


go_name=paste(name,"p-value-005_log1.csv",sep="_")
fun_enrich_cutoff <- read.csv(paste(deseq_folder,go_name,sep ="/"), row.names = 'gene')


cutoff_genes <- fun_enrich_cutoff %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


annotations_cutoff <- AnnotationDbi::select(org.Hs.eg.db, # database
                                         keys = cutoff_genes$gene,  # data to use for retrieval
                                         columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                         keytype = "ENSEMBL")



# Determine the indices for the non-duplicated genes
non_duplicates_cutoff <- which(duplicated(annotations_cutoff$ENSEMBL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_cutoff <- annotations_cutoff[non_duplicates_cutoff, ]

#delete NA
annotations_cutoff <- na.omit(annotations_cutoff)

## Run GO enrichment analysis with a 0.05 p-value cutoff 
colnames(cutoff_genes)[1] <- "ENSEMBL"
annotations_cutoff<- merge(annotations_cutoff,cutoff_genes,by="ENSEMBL")

positives_cutoff <- annotations_cutoff[annotations_cutoff$log2FoldChange >= folchange,]
up_genes=paste(name,'positive_genes.csv',sep ="_")
write.csv(positives_cutoff, file=paste(deseq_folder,up_genes,sep ="/"),row.names = FALSE)

negatives_cutoff <- annotations_cutoff[annotations_cutoff$log2FoldChange <= -folchange,]
down_genes=paste(name,'negative_genes.csv',sep ="_")
write.csv(negatives_cutoff, file=paste(deseq_folder,down_genes,sep ="/"),row.names = FALSE)

ego_cutoff_genes <- enrichGO(gene = annotations_cutoff$ENSEMBL,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)


ego_cutoff_genes_positives <- enrichGO(gene = positives_cutoff$ENSEMBL,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)

ego_cutoff_genes_negatives <- enrichGO(gene = negatives_cutoff$ENSEMBL,
                                       keyType = "ENSEMBL",
                                       OrgDb = org.Hs.eg.db, 
                                       ont = "BP", 
                                       pAdjustMethod = "BH", 
                                       qvalueCutoff = 0.05, 
                                       readable = TRUE)

go_results <- data.frame(ego_cutoff_genes)
go_results_positives <- data.frame(ego_cutoff_genes_positives)
go_results_negatives <- data.frame(ego_cutoff_genes_negatives)

go_name=paste(name,"go_results",sep="_")
GO_folder=paste(output_path,go_name,sep="/")
dir.create(paste(output_path,go_name,sep="/"))
GO_name=paste(name,"GO_positives.csv",sep="_")
write.csv(as.data.frame(go_results_positives), file=paste(GO_folder,GO_name,sep="/"),row.names = FALSE)
GO_name=paste(name,"GO_negatives.csv",sep="_")
write.csv(as.data.frame(go_results_negatives), file=paste(GO_folder,GO_name,sep="/"),row.names = FALSE)
GO_name=paste(name,"GO_merge_positives_and_negatives.csv",sep="_")
write.csv(as.data.frame(go_results), file=paste(GO_folder,GO_name,sep="/"),row.names = FALSE)


#Create dotplot for each GO enrichment

if (nrow(go_results_positives) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_go_up <-dotplot(ego_cutoff_genes_positives, showCategory=20, font.size = 14)
  plot_go_p=paste(name,"PLOT_positives.png",sep="_")
  png(filename=paste(GO_folder,plot_go_p,sep="/"), width=400, height=700)
  plot(plot_go_up)
  dev.off()
}


if (nrow (go_results_negatives) < 1) {
  print ("there are not GO categories")
  
}else{
  plot_go_down <-dotplot(ego_cutoff_genes_negatives, showCategory=20, font.size = 14)
  plot_go_d=paste(name,"PLOT_negatives.png",sep="_")
  png(filename=paste(GO_folder,plot_go_d,sep="/"), width=400, height=700)
  plot(plot_go_down)
  dev.off() 
}



######################################################
######################################################
######################################################
######################################################


#pathway enrichment with KEGG



pathway_name=paste(name,"pathway_enrichment",sep="_")
path_folder=paste(output_path,pathway_name,sep="/")
dir.create(path_folder)

pathfinder_positive <- positives_cutoff[, c('SYMBOL', 'log2FoldChange', 'padj')]


output_pathfindr_positive <- tryCatch(
  {
    run_pathfindR(pathfinder_positive, p_val_threshold = 0.05)
  },
  error = function(e) {
    output_pathfindr_positive <- data.frame(matrix(ncol = 0, nrow = 0))
  },
  warning = function(w) {
    output_pathfindr_positive <- data.frame(matrix(ncol = 0, nrow = 0))
  }
)



unlink("pathfindR_Results",recursive=TRUE)
pathway_name=paste(name,"path_positives.csv",sep="_")
write.csv(output_pathfindr_positive, file=paste(path_folder,pathway_name,sep="/"),row.names = FALSE)

if (nrow(output_pathfindr_positive) ==0){
  print ("there are not metabolic categories")
}else{
  plot_path_up <- enrichment_chart(output_pathfindr_positive[1:10, ])
  plot_pathway_p=paste(name,"plot_pathway_positives.png",sep="_")
  png(filename=paste(path_folder,plot_pathway_p,sep="/"), width=400, height=700)
  plot(plot_path_up)
  dev.off()
}


pathfinder_negatives <- negatives_cutoff[, c('SYMBOL', 'log2FoldChange', 'padj')]

output_pathfindr_negatives <- tryCatch(
  {
    run_pathfindR(pathfinder_negatives, p_val_threshold = 0.05)
  },
  error = function(e) {
    output_pathfindr_negatives <- data.frame(matrix(ncol = 0, nrow = 0))
  },
  warning = function(w) {
    output_pathfindr_negatives <- data.frame(matrix(ncol = 0, nrow = 0))
  }
)



unlink("pathfindR_Results",recursive=TRUE)
pathway_name=paste(name,"path_negatives.csv",sep="_")
write.csv(output_pathfindr_negatives, file=paste(path_folder,pathway_name,sep="/"),row.names = FALSE)

if (nrow(output_pathfindr_negatives)==0 ){
  print ("there are not metabolic categories")
}else{
  plot_path_down <- enrichment_chart(output_pathfindr_negatives[1:10, ])
  plot_pathway_d=paste(name,"plot_pathway_negatives.png",sep="_")
  png(filename=paste(path_folder,plot_pathway_d,sep="/"), width=400, height=700)
  plot(plot_path_down)
  dev.off()
}




######################################################
######################################################
######################################################
######################################################


#pathway enrichment with Reactome

reactome_pathway_name=paste(name,"reactome_pathway_enrichment",sep="_")
reactome_path_folder=paste(output_path,reactome_pathway_name,sep="/")
dir.create(reactome_path_folder)


positive_entrez_id <- positives_cutoff$ENTREZID[!is.na(positives_cutoff$ENTREZID)]

positive_enrich_reactome <- tryCatch(
  {
    enrichPathway(gene=positive_entrez_id,organism = "human",minGSSize = 10, maxGSSize = 500,)
  },
  error = function(e) {
    positive_enrich_reactome <- data.frame(matrix(ncol = 0, nrow = 0))
  },
  warning = function(w) {
    positive_enrich_reactome <- data.frame(matrix(ncol = 0, nrow = 0))
  }
)


file_name=paste(name,"reactome_positives.csv",sep="_")
write.csv(DataFrame(positive_enrich_reactome), file=paste(reactome_path_folder,file_name,sep="/"),row.names = FALSE)


if (nrow(positive_enrich_reactome)==0 ){
  print ("there are not metabolic categories")
}else{
  plot_path_up <- barplot(positive_enrich_reactome, showCategory=20)
  plot_pathway_up=paste(name,"plot_pathway_positives.png",sep="_")
  png(filename=paste(reactome_path_folder,plot_pathway_up,sep="/"), width=400, height=700)
  plot(plot_path_up)
  dev.off()
}






negative_entrez_id <- negatives_cutoff$ENTREZID[!is.na(negatives_cutoff$ENTREZID)]

negative_enrich_reactome <- tryCatch(
  {
    enrichPathway(gene=negative_entrez_id,organism = "human",minGSSize = 10, maxGSSize = 500,)
  },
  error = function(e) {
    negative_enrich_reactome <- data.frame(matrix(ncol = 0, nrow = 0))
  },
  warning = function(w) {
    negative_enrich_reactome <- data.frame(matrix(ncol = 0, nrow = 0))
  }
)


file_name=paste(name,"reactome_negatives.csv",sep="_")
write.csv(DataFrame(negative_enrich_reactome), file=paste(reactome_path_folder,file_name,sep="/"),row.names = FALSE)


if (nrow(negative_enrich_reactome)==0 ){
  print ("there are not metabolic categories")
}else{
  plot_path_down <- barplot(negative_enrich_reactome, showCategory=20)
  plot_pathway_up=paste(name,"plot_pathway_negatives.png",sep="_")
  png(filename=paste(reactome_path_folder,plot_pathway_up,sep="/"), width=400, height=700)
  plot(plot_path_down)
  dev.off()
}








#Gene Set Enrichment Analysis


dea_name=paste(name,"not_cutoff.csv",sep="_")
gsea_input <- read.csv(paste(deseq_folder,dea_name,sep ="/"))
gsea_input2 <- gsea_input[!duplicated(gsea_input$gene), ]

notfilter_deg <- gsea_input2 %>%
  rownames_to_column() %>% 
  as_tibble()


notfilter_deg <- AnnotationDbi::select(org.Hs.eg.db, # database
                                            keys = notfilter_deg$gene,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")



# Determine the indices for the non-duplicated genes
non_duplicates_notfilter_deg <- which(duplicated(notfilter_deg$ENSEMBL) == FALSE)

# Return only the non-duplicated genes using indices
notfilter_deg <- notfilter_deg[non_duplicates_notfilter_deg, ]

#delete NA
notfilter_deg <- na.omit(notfilter_deg)

## Run gseGO enrichment analysis with a 0.05 p-value cutoff 

colnames(gsea_input2)[7] <- "ENSEMBL"
notfilter_deg<- merge(notfilter_deg,gsea_input2,by="ENSEMBL")

foldchanges <- notfilter_deg$log2FoldChange
names(foldchanges) <- notfilter_deg$ENSEMBL

foldchanges <- sort(foldchanges, decreasing = TRUE)

gseaenrichment_result <- gseGO(geneList     = foldchanges,
                               OrgDb        = org.Hs.eg.db,
                               ont          = "BP",
                               keyType = "ENSEMBL",
                               minGSSize    = 100,
                               maxGSSize    = 500,
                               pvalueCutoff = 0.05,
                               verbose      = FALSE)
require(DOSE)

gesaplot <- dotplot(gseaenrichment_result, showCategory=10, split=".sign",font.size = 8) + facet_grid(.~.sign)

gsea_final_result <- gseaenrichment_result@result

gsea_name=paste(name,"gsea_enrichment",sep="_")
gsea_folder=paste(output_path,gsea_name,sep="/")
dir.create(gsea_folder)
gsea_name=paste(name,"gsea.csv",sep="_")
write.csv(gsea_final_result, file=paste(gsea_folder,gsea_name,sep="/"),row.names = FALSE)

plot_gsea_name=paste(name,"plot_gsea.png",sep="_")
png(filename=paste(gsea_folder,plot_gsea_name,sep="/"), width=400, height=700)
plot(gesaplot)

dev.off()


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################


#GSEA analysis using whole MSigDB database

collection <- c('H','C1','C2','C3','C4','C5','C6','C7','C8')

geneList_MSigDB <- resOrdered$log2FoldChange
names(geneList_MSigDB) <- resOrdered$gene
geneList_MSigDB <- sort(geneList_MSigDB, decreasing = TRUE)


GSEA_MSigDB_name=paste(name,"GSEA-MSigDB_enrichment",sep="_")
GSEA_MSigDB_path_folder=paste(output_path,GSEA_MSigDB_name,sep="")
dir.create(GSEA_MSigDB_path_folder)


for (cole in collection){
all_gene_sets_msigdb <- msigdbr(species = "Homo sapiens",category = cole)

gsea_results_MSigDB <- GSEA(geneList = geneList_MSigDB,pvalueCutoff = 0.05, TERM2GENE = all_gene_sets_msigdb[, c("gs_name", "ensembl_gene")])

enrich_dataframe <- as.data.frame(gsea_results_MSigDB)

msigdb_dataframe <- as.data.frame(all_gene_sets_msigdb)
msigdb_dataframe <- msigdb_dataframe[,c("gs_name","gs_description")]

enrich_dataframe <- merge(x=enrich_dataframe,y=msigdb_dataframe,  by.x= "ID",  by.y= "gs_name")

enrich_dataframe <- enrich_dataframe[!duplicated(enrich_dataframe$ID), ]

csv_file_name= paste('collection_',cole,".csv",sep="")

write.csv(enrich_dataframe,paste(GSEA_MSigDB_path_folder,csv_file_name,sep="/"))


  if (nrow(enrich_dataframe) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_msigdb  <- dotplot(gsea_results_MSigDB, showCategory=10, split=".sign",font.size = 8) + facet_grid(.~.sign)
  #Create name to save the file
  plot_go_p=paste('collection_',cole,".png",sep="")
  #opne the png file
  png(filename=paste(GSEA_MSigDB_path_folder,plot_go_p,sep="/"), width=400, height=700)
  plot(plot_msigdb)
  dev.off()
}

}








#Network Analysis

Humannet_library <- read.csv(file="/Users/juansola/Documents/TC3R/research/092023_Sima_experiment22/Report_RNAseq_SIMA_experiment/240422_Comparison_to_ TCGA-pacreas/240422_Comparison_TCGA_pacreas_inputs/240510_network/HumanNet-XC.tsv",sep = '\t', header=FALSE) 
Humannet_library
colnames(Humannet_library) <- c("gene_1","gene_2","score")


#create a list with our gene of interest
my_entrez <- as.list(annotations_cutoff$ENTREZID)


#create a new dataframe with the genes that show interaction in HUmannet database

df_network <- data.frame()
for (i in 1:nrow(Humannet_library)){
    current_row <- Humannet_library[i, ]
    id_1 <- current_row[,1]
    id_2 <- current_row[,2]
    if ((id_1 %in% my_entrez) & (id_2 %in% my_entrez)){
        df_network <- rbind(df_network, current_row)
    }
}


#anotate the new dataframe with the symbol of each gene
symbol_1 <- AnnotationDbi::select(org.Hs.eg.db, # database
                                         keys = as.character(df_network$gene_1),  # data to use for retrieval
                                         columns = c("SYMBOL"), # information to retreive for given data
                                         keytype = "ENTREZID")
symbol_2 <- AnnotationDbi::select(org.Hs.eg.db, # database
                                         keys = as.character(df_network$gene_2),  # data to use for retrieval
                                         columns = c("SYMBOL"), # information to retreive for given data
                                         keytype = "ENTREZID")



name_in_network=strsplit(name, split = "_")
mage_in_network = name_in_network[[1]][5]


#create folder to save th network results
net_name=paste(name,"network_interaction",sep="_")
dir.create(paste(output_path,net_name,sep="/"))
network_folder=paste(output_path,net_name,sep="/")



#Select symbosl infromation and generate the complete network
df_network["symbol_1"] <- symbol_1$SYMBOL
df_network["symbol_2"] <- symbol_2$SYMBOL

#Save the data used to create the network
connected_symbosl=paste(name,"nextwork_data.csv",sep="_")
write.csv(as.data.frame(df_network), file=paste(network_folder,connected_symbosl,sep="/"),row.names = FALSE)

pdf_network=paste(name,"networks_plots.pdf",sep="_")
pdf(paste(network_folder,pdf_network,sep="/"))


#generate the complete network
graph_node <- df_network[c("symbol_1","symbol_2")]
weights <- c(df_network$score)

g <- graph_from_data_frame(graph_node, directed=FALSE)
E(g)$weight <- weights # Use directed=TRUE for directed graphs

plot(g)
title(main = "General Network")


#most connected node network
degrees <- degree(g)  # Calculate degrees for all nodes
write.csv(as.data.frame(degrees), file=paste(network_folder,"node_degree.csv",sep="/"),row.names = TRUE)


most_connected_node <- which.max(degrees)  # Find the node with the highest degree

high_connected <- c(V(g)[most_connected_node])

# Step 2: Select neighbors of the most connected node
neighbors_of_most_connected <- neighbors(g, high_connected)
# Step 3: Create a sub-graph
# Include both the most connected node and its neighbors
nodes_to_include <- c(neighbors_of_most_connected,high_connected)

# Create the subgraph based on the nodes
sub_g <- induced_subgraph(g, vids = nodes_to_include)

# Step 4: Plot the sub-graph (optional)
plot(sub_g)
title(main = "Most Connected Gene Network")


# MAGE network

network_mage <- c(V(g)[name ==mage_in_network ])

# Step 2: Select neighbors of the most connected node
neighbors_of_mage <- neighbors(g, network_mage)
# Step 3: Create a sub-graph
# Include both the most connected node and its neighbors
nodes_to_include <- c(neighbors_of_mage,network_mage)

# Create the subgraph based on the nodes
sub_g_mage <- induced_subgraph(g, vids = nodes_to_include)

# Step 4: Plot the sub-graph (optional)
plot(sub_g_mage)
title(main = paste(mage_in_network,"Linkages",sep=" "))

dev.off()



#general statistis of the network
file_net_sta <- file(paste(network_folder,"network_statistics.txt",sep="/"), "w")
#Number of edges
cat("NUMBER OF EDGES\n", file = file_net_sta)
cat(gsize(g), file = file_net_sta)
cat("\n", file = file_net_sta)


#Number of nodes
cat("NUMBER OF NODES\n", file = file_net_sta)
cat(gorder(g), file = file_net_sta)
cat("\n", file = file_net_sta)

#articulation points
cat("ARTICULATION POINTS\n", file = file_net_sta)
cat(articulation_points(g)$name, file = file_net_sta)
cat("\n", file = file_net_sta)





#Sumary table

sumarry <- matrix(c(table(condition$condition)[2],table(condition$condition)[1],nrow(positives_cutoff),
                    nrow(negatives_cutoff),nrow(go_results_positives),nrow(go_results_negatives),
                    nrow(output_pathfindr_positive),nrow(output_pathfindr_negatives),resultsNames(dea)[2]),ncol=1,byrow=TRUE)
colnames(sumarry) <- c("Quantity")
rownames(sumarry) <- c("Positive Tumors","Negative Tumors","Up-regulated","Down-regulated","Up-regulated GO-enrichment",
                       "Down-regulated GO-enrichment","Up-regulated Pathway-enrichment","Down-regulated Pathway-enrichment","Comparison")
sumarry <- as.table(sumarry)
summary_name=paste(name,"summary.csv",sep="_")
write.csv(sumarry, file=paste(output_path,summary_name,sep="/"))

nrow(go_results_negatives)

}
