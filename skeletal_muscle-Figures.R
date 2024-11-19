#Read in files



orig.ident.levels <- c("FC", "FN", "SC", "SN")
prep.levels <- c("Filtered", "Sorted")
load.levels <- c("Cells", "Nuclei")

seurat.raw$orig.ident <- factor(x = seurat.raw$orig.ident, levels = orig.ident.levels)
seurat.object$orig.ident <- factor(x = seurat.object$orig.ident, levels = orig.ident.levels)
seurat.object$Prep <- factor(x = seurat.object$Prep, levels = prep.levels)
seurat.raw$Prep <- factor(x = seurat.raw$Prep, levels = prep.levels)
seurat.object$Load <- factor(x = seurat.object$Load, levels = load.levels)
seurat.raw$Load <- factor(x = seurat.raw$Load, levels = load.levels)

# Figure 2 ----------------------------------------------------------------

png(filename = paste0(rnaProject, "-Figure2b.png"), height = 800, width = 800)
VlnPlot(seurat.object, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") + 
	NoLegend() + 
	geom_boxplot(
		width = 0.07,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values=c("white", "gray", "white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-Figure2c.png"), height = 800, width = 800)
VlnPlot(seurat.object, features = "percent.mt", pt.size = 0, group.by = "orig.ident") + 
	NoLegend() + 
	geom_boxplot(
		width = 0.07,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values=c("white", "gray", "white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-Figure 2e.png"), height = 800, width = 800)
DimPlot(seurat.object, reduction = "umap", cols = color.palette, label = F, label.size = 7, repel = T, pt.size = 1, shuffle = TRUE) + NoLegend()
dev.off()

#Generate barplot for Figure 2d
ora.data <- read.table("enrichment_results_wg_result1710432063.txt", 
											 sep = "\t", row.names = 1, header = TRUE)
ora.table <- ora.data %>%
	filter(FDR <= 0.05) %>%
	mutate(neglog10FDR = -log10(FDR)) %>%
	select(c(description, enrichmentRatio, neglog10FDR, FDR)) %>%
	arrange(desc(enrichmentRatio))

hvg.cc.barplot <- ggplot(ora.table, aes(x=reorder(description,enrichmentRatio), y=enrichmentRatio, fill=neglog10FDR)) + 
	geom_bar(stat='identity') + 
	coord_flip() +
	labs(y = "enrichment Ratio", x = "cell compartment") +
	scale_fill_continuous(low="lightgrey", high="blue") +
	theme_minimal_grid()

png(filename = paste0(rnaProject, "-Figure 2d.png"), height = 300, width = 800)
plot(hvg.cc.barplot)
dev.off()

# Figure 3 ----------------------------------------------------------------

library(ComplexHeatmap)
library(RColorBrewer)


figmarkers <- c(
	'APOC1', 'APOE', 'MYF5', #satellite
	'MYOD1', 'PAX7', #myoblast and/or progenitor
	'RGS5', 'ACTA2', 'NOTCH3', 'NDUFA4L2', 'HIGD1B', #pericyte
	'LUM', 'DCN', # LUM+ FAP
	'FBN1', 'PCOLCE2', #PRG4+ FAP
	'APOD', 'GPX3', #adipocyte
	'COL3A1', 'COL1A1', 'PDGFRA', #fibroblast
	'PECAM1', 'CDH5', 'KDR', 'ESAM', 'DARC', 'FABP4', 'IFI27', 'VWF',#endothelial
	'PERGL', 'MYH11', 'PLN', #smooth muscle
	'TPM3', 'TNNC1', 'TNNI1', 'TNNT1', 'MYL2', 'MYL3', 'MYL6B', 'MYH2', #type IIA fiber
	'MYH1', #Type IIX fiber
	'TPM1', 'TNNC2', 'TNNI2', 'TNNT3', 'MYL1', 'MYLPF', 'MYBPC2', 'ENO3',  #sarcomeric
	'ATP2A1', 'SLN','ATP2A2',  #calcium transport
	'MYH7', #type I fiber
	'HBB', 'HBA1', 'HBA2' #erythroid
)
figmarkers <- figmarkers[figmarkers %in% Features(seurat.object)]

#Fig3A
heatmap.subset <- subset(seurat.object, downsample = 400)
heatmap.data <- heatmap.subset[["RNA"]]@scale.data[figmarkers, ] %>% as.matrix()
heatmap.data <- scale(t(heatmap.data))
cluster_anno<- heatmap.subset@meta.data$seurat_clusters
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00")) # set color scale based on 10% and 95% quantiles
ht.palette <- color.palette[1:12]
color.vector <- c(7, 6, 8, 5, 2, 1, 0, 3, 11, 4, 9, 10) +1 # +1 because of the stupid 0 shift
ht.palette <- ht.palette[color.vector]

png(filename = paste0(rnaProject, "-Figure3a.png"), height = 700, width = 1600, bg = "transparent")
ht <- draw(Heatmap(heatmap.data, name = "Expression",  
									 row_split = factor(cluster_anno),
									 cluster_rows = TRUE,
									 show_row_dend = FALSE,
									 cluster_row_slices = TRUE,
									 row_title_gp = gpar(fontsize = 10),
									 row_gap = unit(1, "mm"),
									 cluster_columns = TRUE,
									 show_column_dend = TRUE,
									 col = col_fun,
									 column_names_gp = gpar(fontsize = 10),
									 row_title_rot = 90,
									 left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = ht.palette))),
									 show_row_names = FALSE, use_raster = TRUE))
dev.off()

ht.col <- column_order(ht)
figmarkers[ht.col]

#Fig 3b
png(filename = paste0(rnaProject, "-Figure3b.png"), height = 700, width = 1600)
DotPlot(seurat.object, features = figmarkers[ht.col], col.min = 0, dot.min = 0.01, dot.scale = 9) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# Figure 4 ----------------------------------------------------------------

#umaps
png(filename = paste0(rnaProject, "-Figure4a.png"), height = 800, width = 1500)
DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", split.by = "Load", label = F, pt.size = 1.5, cols = color.palette, shuffle = TRUE) + NoLegend()
dev.off()

png(filename = paste0(rnaProject, "-Figure4e.png"), height = 800, width = 1500)
DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", split.by = "Prep", label = F, pt.size = 1.5, cols = color.palette, shuffle = TRUE) + NoLegend()
dev.off()


png(filename = paste0(rnaProject, "-Figure4i.png"), height = 800, width = 1500)
DimPlot(seurat.subset, reduction = "umap", group.by = "seurat_clusters", split.by = "Prep", label = F, pt.size = 1.5, cols = subset.palette) + NoLegend()
dev.off()

#violin plots
png(filename = paste0(rnaProject, "-Figure4b.png"), height = 800, width = 600)
VlnPlot(seurat.object, features = "nFeature_RNA", pt.size = 0, group.by = "Load") + 
	geom_boxplot(
		width = 0.04,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values = c("white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-Figure4c.png"), height = 800, width = 600)
VlnPlot(seurat.object, features = "percent.mt", pt.size = 0, group.by = "Load") + 
	geom_boxplot(
		width = 0.1,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values = c("white", "gray"))
dev.off()


png(filename = paste0(rnaProject, "-Figure4f.png"), height = 800, width = 600)
VlnPlot(seurat.object, features = "nFeature_RNA", pt.size = 0, group.by = "Prep") + 
	geom_boxplot(
		width = 0.1,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values = c("white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-Figure4g.png"), height = 800, width = 600)
VlnPlot(seurat.object, features = "percent.mt", pt.size = 0, group.by = "Prep") + 
	geom_boxplot(
		width = 0.1,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values = c("white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-Figure4j.png"), height = 800, width = 600)
VlnPlot(seurat.object, group.by = "Prep", features = "nFeature_RNA", pt.size = 0, cols = c("white", "gray")) + 
	NoLegend() +
	stat_summary(fun.y=median, geom="crossbar", size=2, color="black", width = 0.5, linewidth = 0.5) + 
	theme(rect = element_rect(fill = "transparent"))
dev.off()

png(filename = paste0(rnaProject, "-Figure4h.png"), height = 800, width = 600)
VlnPlot(seurat.object, group.by = "Prep", features = "percent.mt", pt.size = 0, cols = c("white", "gray")) + 
	NoLegend() +
	stat_summary(fun.y=median, geom="crossbar", size=2, color="black", width = 0.5, linewidth = 0.5) + 
	theme(rect = element_rect(fill = "transparent"))
dev.off()

#See `MiloR.R` for beeswarm plots

# Extended Data Figure 1 --------------------------------------------------

png(filename = paste0(rnaProject, "-EData1a.png"), height = 800, width = 800)
VlnPlot(seurat.raw, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident") + 
	NoLegend() + 
	geom_boxplot(
		width = 0.07,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values=c("white", "gray", "white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-EData1b.png"), height = 800, width = 800)
VlnPlot(seurat.raw, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") + 
	NoLegend() + 
	geom_boxplot(
		width = 0.07,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	geom_hline(yintercept = 200, linetype = "dashed", color = "gray40", linewidth = 1.1) +
	geom_hline(yintercept = 2500, linetype = "dashed", color = "gray40", linewidth = 1.1) +
	scale_fill_manual(values=c("white", "gray", "white", "gray"))
dev.off()


png(filename = paste0(rnaProject, "-EData1c.png"), height = 800, width = 800)
VlnPlot(seurat.raw, features = "percent.mt", pt.size = 0, group.by = "orig.ident") + 
	NoLegend() + 
	geom_boxplot(
		width = 0.07,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	geom_hline(yintercept = 10, linetype = "dashed", color = "gray40", linewidth = 1.1) +
	scale_fill_manual(values=c("white", "gray", "white", "gray"))
dev.off()

#Extended data figure 1d
for(i in levels(seurat.raw$orig.ident)){
	scatterp <- FeatureScatter(subset(seurat.raw, subset = orig.ident == i), feature1 = "nCount_RNA", feature2 = "nFeature_RNA", split.by = "orig.ident", ncol = 2) + 
		geom_hex(bins = 100) + 
		NoLegend()
	png(filename = paste0(rnaProject, "-", i, "-EData1d.png"), height = 800, width = 800)
	plot(scatterp)
	dev.off()
}

png(filename = paste0(rnaProject, "-EData1e.png"), height = 800, width = 800)
FeatureScatter(seurat.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
														group.by = "orig.ident", pt.size = 4, cols = color.palette) + 
	geom_hex(bins = 80) 
dev.off()



# Extended Data Figure 2 --------------------------------------------------

png(filename = paste0(rnaProject, "-EData2a.png"), height = 800, width = 800)
VlnPlot(seurat.object, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident") + 
	NoLegend() + 
	geom_boxplot(
		width = 0.07,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values=c("white", "gray", "white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-EData2b.png"), height = 800, width = 800)
VlnPlot(seurat.object, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") + 
	NoLegend() + 
	geom_boxplot(
		width = 0.07,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values=c("white", "gray", "white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-EData2c.png"), height = 800, width = 800)
VlnPlot(seurat.object, features = "percent.mt", pt.size = 0, group.by = "orig.ident") + 
	NoLegend() + 
	geom_boxplot(
		width = 0.07,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values=c("white", "gray", "white", "gray"))
dev.off()



#Extended data figure 2d
for(i in levels(seurat.object$orig.ident)){
	scatterp <- FeatureScatter(subset(seurat.object, subset = orig.ident == i), feature1 = "nCount_RNA", feature2 = "nFeature_RNA", split.by = "orig.ident", ncol = 1) + 
		geom_hex(bins = 50) + 
		NoLegend()
	png(filename = paste0(rnaProject, "-", i, "-EData2d.png"), height = 800, width = 800)
	plot(scatterp)
	dev.off()
}

png(filename = paste0(rnaProject, "-EData2e.png"), height = 800, width = 800)
FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
														group.by = "orig.ident", pt.size = 4, cols = color.palette) + 
	geom_hex(bins = 80) 
dev.off()




# Extended Data Figure 3 --------------------------------------------------

png(filename = paste0(rnaProject, "-EDataFigure3a.png"), height = 800, width = 1200)
VlnPlot(seurat.object, features = "nFeature_RNA", pt.size = 0, split.by = "Prep", split.plot = TRUE) + 
	geom_boxplot(
		width = 0.2,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values = c("white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-EDataFigure3b.png"), height = 800, width = 1200)
VlnPlot(seurat.object, features = "nCount_RNA", pt.size = 0, split.by = "Prep", split.plot = TRUE) + 
	geom_boxplot(
		width = 0.2,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values = c("white", "gray"))
dev.off()

png(filename = paste0(rnaProject, "-EDataFigure3c.png"), height = 800, width = 1200)
VlnPlot(seurat.object, features = "percent.mt", pt.size = 0, split.by = "Prep", split.plot = TRUE) + 
	geom_boxplot(
		width = 0.2,
		notch = FALSE,
		notchwidth = 0.1,
		outlier.shape = NA,
		coef = 0) +
	scale_fill_manual(values = c("white", "gray"))
dev.off()

#Generate barplot for Extended Data Figure 3d
gsea.data <- read.table("UnsortedVsSorted-combined.txt", sep = "\t", row.names = 1, header = TRUE)

gsea.table <- gsea.data %>%
	filter(FDR <= 0.05) %>%
	mutate(FDR=replace(FDR, FDR==0, (min(gsea.data$FDR[gsea.data$FDR > 0])/10))) %>%
	mutate(neglog10FDR = -log10(FDR)) %>%
	select(c(description, normalizedEnrichmentScore, neglog10FDR, FDR)) %>%
	arrange(desc(normalizedEnrichmentScore))

head(gsea.table)
hvg.cc.barplot <- ggplot(gsea.table, aes(x=reorder(description,normalizedEnrichmentScore), y=normalizedEnrichmentScore, fill=neglog10FDR)) + 
	geom_bar(stat='identity') + 
	coord_flip() +
	labs(y = "normalized Enrichment ratio", x = "GSEA term") +
	scale_fill_continuous(low="lightgrey", high="blue") +
	theme_minimal_grid()
plot(hvg.cc.barplot)

png(filename = paste0(rnaProject, "-UnsortedvsSorted-GSEA_barplot.png"), height = 600, width = 800)
plot(hvg.cc.barplot)
dev.off()


