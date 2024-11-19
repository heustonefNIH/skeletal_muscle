# Global parameters -------------------------------------------------------

rnaProject <- "preserved_muscle"
percent.mt.filtered <- 10
cum.var.thresh <- 80
resolution <- 0.5
do.doubletFinder <- TRUE
min.cells <- 3
min.features <- 200
doublet.var.thresh <- 90
predicted.doubletRate <- 0.05


# Directories -------------------------------------------------------------
rna.dir <- ""
path_to_data <- "./cellranger_files"
sourceable.functions <- list.files(path = "./RFunctions", pattern = "*.R$", full.names = TRUE)
metadata.location <- ""

# Load libraries ----------------------------------------------------------

library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggplot2)

##load local functions
invisible(sapply(sourceable.functions, source))

try(setwd(rna.dir), silent = TRUE)
writeLines(capture.output(sessionInfo()), paste0(rnaProject, "_sessionInfo.txt"))

# Load data ---------------------------------------------------------------

## load data list
sc.data <- sapply(list.dirs(
	path = path_to_data, recursive = FALSE, full.names = TRUE), 
	basename, 
	USE.NAMES = TRUE)

# Load metadata and select only relevant columns
metadata <- read.table(file = paste0(metadata.location, "Metadata_file.txt"), header = TRUE, sep = "\t", row.names = 1)[,c(1:3)]
head(metadata)


# Read cellbender data into seurat objects --------------------------------

object.list <- list()
for(i in 1:length(sc.data)){
	object.item <- Read10X_h5(paste0(names(sc.data)[i], "/outs/cb-seurat_feature_bc_matrix_filtered.h5")) # select cellbender files
	object.item <- CreateSeuratObject(
		object.item, 
		project = rnaProject, 
		min.cells = min.cells, 
		min.features = min.features)
	object.item$orig.ident <- sc.data[[i]]
	object.item <- AssignMetadata(metadata.df = metadata, seurat.object = object.item)
	object.item <- PercentageFeatureSet(object.item, pattern = "MT-", col.name = "percent.mt")
	print(paste("adding", sc.data[[i]], "to list"))
	object.list <- c(object.list, object.item)
}

names(object.list) <- sc.data
seurat.raw <- merge(object.list[[1]], y = object.list[2:length(object.list)], add.cell.ids = names(object.list))
saveRDS(seurat.raw, file = paste0(rnaProject, "-rawMergedSeurat.Object.RDS"))


# Reporting metrics -------------------------------------------------------

head(seurat.raw@meta.data)
raw.meta <- as.data.frame(seurat.raw@meta.data)
raw.meta %>% 
	group_by(orig.ident) %>% 
	reframe(umi.median = median(nCount_RNA),
					umi.qs25 = quantile(nCount_RNA, 0.25), 
					umi.qs75 = quantile(nCount_RNA, 0.75), 
					feat.median = median(nFeature_RNA),
					feat.qs25 = quantile(nFeature_RNA, 0.25), 
					feat.qs75 = quantile(nFeature_RNA, 0.75), 
					mt.median = median(percent.mt),
					mt.qs25 = quantile(percent.mt, 0.25), 
					mt.qs75 = quantile(percent.mt, 0.75), 
	)

raw.meta %>% 
	group_by(Prep) %>% 
	reframe(umi.median = median(nCount_RNA),
					umi.qs25 = quantile(nCount_RNA, 0.25), 
					umi.qs75 = quantile(nCount_RNA, 0.75), 
					feat.median = median(nFeature_RNA),
					feat.qs25 = quantile(nFeature_RNA, 0.25), 
					feat.qs75 = quantile(nFeature_RNA, 0.75), 
					mt.median = median(percent.mt),
					mt.qs25 = quantile(percent.mt, 0.25), 
					mt.qs75 = quantile(percent.mt, 0.75), 
	)

raw.meta %>% 
	group_by(Load) %>% 
	reframe(umi.median = median(nCount_RNA),
					umi.qs25 = quantile(nCount_RNA, 0.25), 
					umi.qs75 = quantile(nCount_RNA, 0.75), 
					feat.median = median(nFeature_RNA),
					feat.qs25 = quantile(nFeature_RNA, 0.25), 
					feat.qs75 = quantile(nFeature_RNA, 0.75), 
					mt.median = median(percent.mt),
					mt.qs25 = quantile(percent.mt, 0.25), 
					mt.qs75 = quantile(percent.mt, 0.75), 
	)

summarise(raw.meta)

# Subset object -----------------------------------------------------------

seurat.object <- subset(seurat.raw, 
												subset = nFeature_RNA >= min.features &
													nFeature_RNA <= 2500 &
													percent.mt <= percent.mt.filtered)

# QC for filtered seurat object ----------------------------------------------------------------------

seurat.meta <- as.data.frame(seurat.object@meta.data)
seurat.meta %>% 
	group_by(orig.ident) %>% 
	reframe(umi.median = median(nCount_RNA),
					umi.qs25 = quantile(nCount_RNA, 0.25), 
					umi.qs75 = quantile(nCount_RNA, 0.75), 
					feat.median = median(nFeature_RNA),
					feat.qs25 = quantile(nFeature_RNA, 0.25), 
					feat.qs75 = quantile(nFeature_RNA, 0.75), 
					mt.median = median(percent.mt),
					mt.qs25 = quantile(percent.mt, 0.25), 
					mt.qs75 = quantile(percent.mt, 0.75), 
	)

seurat.meta %>% 
	group_by(Prep) %>% 
	reframe(umi.median = median(nCount_RNA),
					umi.qs25 = quantile(nCount_RNA, 0.25), 
					umi.qs75 = quantile(nCount_RNA, 0.75), 
					feat.median = median(nFeature_RNA),
					feat.qs25 = quantile(nFeature_RNA, 0.25), 
					feat.qs75 = quantile(nFeature_RNA, 0.75), 
					mt.median = median(percent.mt),
					mt.qs25 = quantile(percent.mt, 0.25), 
					mt.qs75 = quantile(percent.mt, 0.75), 
	)

seurat.meta %>% 
	group_by(Load) %>% 
	reframe(umi.median = median(nCount_RNA),
					umi.qs25 = quantile(nCount_RNA, 0.25), 
					umi.qs75 = quantile(nCount_RNA, 0.75), 
					feat.median = median(nFeature_RNA),
					feat.qs25 = quantile(nFeature_RNA, 0.25), 
					feat.qs75 = quantile(nFeature_RNA, 0.75), 
					mt.median = median(percent.mt),
					mt.qs25 = quantile(percent.mt, 0.25), 
					mt.qs75 = quantile(percent.mt, 0.75), 
	)


# Normalize and scale data ----------------------------------------------------------

##normalize & scale data
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)
seurat.object <- ScaleData(seurat.object, features = all.genes, vars.to.regress = c("nFeature_RNA"))



# Variable gene tracking --------------------------------------------------

hvg.table <- HVFInfo(seurat.object)
write.table(hvg.table, file = paste0(rnaProject, "-HVGtable.txt"),
						quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

# Linear dimensional reduction --------------------------------------------

seurat.object <- RunPCA(seurat.object, features = top(object = seurat.object))

##determine dimensionality
tot.var <- percent.variance(seurat.object@reductions$pca@stdev, plot.var = FALSE, return.val = TRUE)
paste0("Num pcs for 80% variance: ", length(which(cumsum(tot.var) <= 80)))
paste0("Num pcs for 85% variance: ", length(which(cumsum(tot.var) <= 85)))
paste0("Num pcs for 90% variance: ", length(which(cumsum(tot.var) <= 90)))
paste0("Num pcs for 95% variance: ", length(which(cumsum(tot.var) <= 95)))

cluster.dims <- 0
if(cum.var.thresh > 0){
	cluster.dims <- length(which(cumsum(tot.var) <= cum.var.thresh))
}


# DoubletFinder -----------------------------------------------------------

seurat.object <- runDoubletFinder(seurat.object, sctransformed = FALSE, predicted.doubletRate = predicted.doubletRate)
seurat.object <- subset(seurat.object, subset = DF.classifications == "Singlet")
saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, ".RDS"))

# Louvain cluster ---------------------------------------------------------

##cluster cells
seurat.object <- FindNeighbors(seurat.object, dims = 1:cluster.dims)
seurat.object <- FindClusters(seurat.object, resolution = resolution)
seurat.object <- RunUMAP(seurat.object, dims = 1:cluster.dims)
saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, "-", as.character(cluster.dims), "dims.RDS"))

#distribution tables
cluster.distribution.by.orig.ident <- table(seurat.object$seurat_clusters,seurat.object$orig.ident)
cluster.distribution.by.orig.ident <- round(prop.table(cluster.distribution.by.orig.ident, 2) * 100, 2)

cluster.distribution.by.load <- table(seurat.object$seurat_clusters,seurat.object$Load)
cluster.distribution.by.load <- round(prop.table(cluster.distribution.by.load, 2) * 100, 2)

cluster.distribution.by.prep <- table(seurat.object$seurat_clusters,seurat.object$Prep)
cluster.distribution.by.prep <- round(prop.table(cluster.distribution.by.prep, 2) * 100, 2)

# Tables ------------------------------------------------

#markers based on cluster ID
Idents(seurat.object) <- "seurat_clusters"
all.markers <- FindAllMarkers(seurat.object, only.pos = FALSE)
markers.all <- all.markers %>%
	group_by(cluster) %>%
	select(-c(pct.1, pct.2, p_val)) %>%
	select(c(gene, avg_log2FC, p_val_adj, cluster)) %>%
	arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
write.table(markers.all, file = paste0(rnaProject, "_SupTable2.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#for webgestalt
pos.markers <- FindAllMarkers(seurat.object, only.pos = TRUE)


#markers based on cells vs nuclei
Idents(seurat.object) <- "Load"
load.markers <- FindAllMarkers(seurat.object, only.pos = FALSE)
markers.load <- load.markers %>%
	group_by(cluster) %>%
	arrange(desc(avg_log2FC), .by_group = TRUE) %>%
	select(-c(pct.1, pct.2, p_val))
write.table(markers.load, file = paste0(rnaProject, "_SupTable3.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

#markers based on sorted vs filtered
Idents(seurat.object) <- "Prep"
prep.markers <- FindAllMarkers(seurat.object, only.pos = FALSE)
markers.prep <- prep.markers %>%
	group_by(cluster) %>%
	arrange(desc(avg_log2FC), .by_group = TRUE) %>%
	select(-c(pct.1, pct.2, p_val)) 
write.table(markers.prep, file = paste0(rnaProject, "_SupTable4.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# Web Gestalt Markers --------------------------------------------

library(WebGestaltR)

# reformat to include only significant markers
markers.pos<- pos.markers %>%
	filter(p_val_adj < 0.05) %>%
	select(-p_val_adj, -pct.1, -pct.2, -p_val) %>%
	relocate(gene, avg_log2FC)
markers.pos <- data.frame(markers.pos)
rownames(markers.pos) <- 1:nrow(markers.pos)

#Select databases to test
enrich.databases <- c(
	"pathway_KEGG",
	"pathway_Reactome",
	"pathway_Wikipathway",
	"geneontology_Biological_Process",
	"geneontology_Molecular_Function_noRedundant"
)

#GSEA for DEGs between clusters
clusters <- unique(seurat.object$seurat_clusters)
for(i in clusters){
	group.subset <- i
	marker.subset <- pos.markers %>%
		filter(cluster == group.subset, avg_log2FC >= 0.5) %>%
		select(-cluster)
	enrichResult <- WebGestaltR(
		enrichMethod="GSEA", organism="hsapiens",
		enrichDatabase=enrich.databases, 
		interestGene = marker.subset,
		interestGeneType="genesymbol", 
		sigMethod = "fdr", minNum=10, fdrThr = 1, gseaPlotFormat = "png",
		isOutput = TRUE, saveRawGseaResult = FALSE, 
		projectName = paste0("posmarkers-", group.subset, "-pwayandgo-", as.character(as.integer(Sys.time())))
	)
	
	marker.subset <- all.markers %>%
		filter(cluster == group.subset, avg_log2FC >= 0.5 | avg_log2FC <= -0.5) %>%
		select(-cluster)
	enrichResult <- WebGestaltR(
		enrichMethod="GSEA", organism="hsapiens",
		enrichDatabase=enrich.databases, interestGene = marker.subset,
		interestGeneType="genesymbol", sigMethod = "fdr", minNum=10, fdrThr = 1, gseaPlotFormat = "png",
		isOutput = TRUE, saveRawGseaResult = FALSE, 
		projectName = paste0("allmarkers-", group.subset, "-pwayandgo-", as.character(as.integer(Sys.time())))
	)
}

#ORA for HVG

hvg.table <- HVFInfo(seurat.object)
hvg.

max(hvg.table$mean)
dim(hvg.table[hvg.table$variance.standardized>2,])
write.table(hvg.table, file = paste0(rnaProject, "-HVGtable.txt"),
						quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)


for(i in enrich.databases){
	enrichResult <- WebGestaltR(
		enrichMethod="ORA", organism="hsapiens",
		enrichDatabase=i, interestGene = beta.markers.trim$gene,
		interestGeneType="genesymbol", 
		referenceGeneType="genesymbol", 
		sigMethod = "fdr", 
		minNum=10, fdrThr = 1, gseaPlotFormat = "png",
		isOutput = TRUE, saveRawGseaResult = FALSE, 
		projectName = paste0("Beta1vsBeta2-", i, "-", as.character(as.integer(Sys.time())))
	)
	enrichResult.list <- c(enrichResult.list, enrichResult)
}

ora.data <- read.table("/Users/heustonef/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/scRNA_ATAC-Seq/huMuscle/DataAnalysis/AllProtect_GEX_markersCellsvsNuclei.txt", 
											 sep = "\t", header = TRUE)
head(ora.data)
markers.ora <- ora.data %>%
	arrange(desc(abs(avg_log2FC)), .by_group = TRUE) %>%
	filter(variance.standardized >=2.5) 

markers.ora <- data.frame(markers.ora)
rownames(markers.ora) <- 1:nrow(markers.ora)






# miloR - Differential abundance ------------------------------------------

# See `MiloR.R` for milo data
