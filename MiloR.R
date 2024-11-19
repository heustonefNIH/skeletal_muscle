#MiloR for DA testing

# Global parameters -------------------------------------------------------

rnaProject <- "preserved_muscle"
rna.dir <- ""
try(setwd(rna.dir), silent = TRUE)

# Libraries ---------------------------------------------------------------

library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

# Load data ---------------------------------------------------------------

seurat.object <- readRDS(file = paste0(rnaProject, "-", as.character(cluster.dims), "dims.RDS"))


# Create milo object ------------------------------------------------------

sc.milo <- Milo(as.SingleCellExperiment(seurat.object))
sc.milo

#add graph from Seurat's NN graph
graph(sc.milo) <- graph(buildFromAdjacency(seurat.object@graphs$RNA_nn, k = 29))

# d = n_dim for KNN graph building
# k = n_dim for clustering and UMAP vis
sc.milo <- makeNhoods(sc.milo, prop = 0.05, k = 29, d=29, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(sc.milo)
sc.milo <- countCells(sc.milo, meta.data = data.frame(colData(sc.milo)), samples="orig.ident") # make the nhood X sample counts matrix

## Create design matrix
milo_design <- data.frame(colData(sc.milo))[,c("orig.ident", "Prep", "Load")]
milo_design$Prep <- as.factor(milo_design$Prep)
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$orig.ident
milo_design

sc.milo <- calcNhoodDistance(sc.milo, d=29, reduced.dim = "PCA")
da.results.all.prep <- testNhoods(sc.milo, design = ~ Prep, design.df = milo_design, fdr.weighting = "graph-overlap")
da.results.all.prep <- annotateNhoods(sc.milo, da.results.all.prep, coldata_col = "seurat_clusters")
da.results.all.prep$seurat_clusters <- as.factor(as.numeric(da.results.all.prep$seurat_clusters))
head(da.results.all.prep)
da.results.all.prep %>%
	arrange(SpatialFDR) %>%
	head()

da.results.all.prep %>%
	filter(SpatialFDR <= 0.05 & logFC >=2) %>%
	filter(seurat_clusters == 3) %>%
	nrow()


ggplot(da.results.all.prep, aes(PValue)) + geom_histogram(bins=50)
ggplot(da.results.all.prep, aes(logFC, -log10(SpatialFDR))) + 
	geom_point() +
	geom_hline(yintercept = 1)

sc.milo <- buildNhoodGraph(sc.milo)

## Plot single-cell UMAP
umap_pl <- scater::plotReducedDim(sc.milo, dimred = "UMAP", colour_by="Prep", text_by = "seurat_clusters", text_size = 3, point_size=0.5) +
	guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sc.milo, da.results.all.prep, layout="UMAP", alpha=0.1) 

umap_pl + nh_graph_pl +
	patchwork::plot_layout(guides="collect")

# ggplot(da.results.all.prep, aes(Prep_fraction)) + geom_histogram(bins=50)
ggplot(da.results.all.prep, aes(seurat_clusters_fraction)) + geom_histogram(bins=50)
da.results.all.prep$seurat_clusters_called <- ifelse(da.results.all.prep$seurat_clusters < 0.7, "Mixed", da.results$seurat_clusters)
head(da.results.all.prep, n = 20)


#default bswarm plot
plotDAbeeswarm(da.results.all.prep, group.by = "seurat_clusters") 

#modify bswarm plot
da.res <- da.results.all.prep
beeswarm_pos <- ggplot_build(
	da.res %>%
		mutate(is_signif = ifelse(SpatialFDR < 0.1, 1, 0)) %>%
		arrange(seurat_clusters) %>%
		ggplot(aes(seurat_clusters, logFC)) +
		ggbeeswarm::geom_quasirandom()
)

pos_x <- beeswarm_pos$data[[1]]$x
pos_y <- beeswarm_pos$data[[1]]$y

n_groups <- unique(da.res$seurat_clusters) %>% length()

da.prepplot <- da.res %>%
	mutate(is_signif = ifelse(SpatialFDR < 0.1, 1, 0)) %>%
	mutate(logFC_color = ifelse(is_signif==1, logFC, 0)) %>%
	arrange(seurat_clusters) %>%
	mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
	mutate(pos_x = pos_x, pos_y=pos_y)

p <- ggplot(da.prepplot, aes(pos_x, pos_y, color=logFC_color)) +
	scale_colour_gradient2(mid = "gray87") +
	guides(color="none") +
	xlab("seurat_clusters") + ylab("Log Fold Change") +
	scale_x_continuous(
		breaks = seq(1,n_groups),
		labels = setNames(levels(da.res$seurat_clusters), seq(1,n_groups))) +
	geom_point() +
	ylim(c(-9, 9)) +
	coord_flip() +
	theme_bw(base_size=22) +
	theme(strip.text.y =  element_text(angle=0))

plot(p)
png(filename = paste0(rnaProject, "-miloRbeeswarm-Prep.png"), height = 800, width = 800)
plot(p)
dev.off()



# DA results for Load -----------------------------------------------------

da.results.all.load <- testNhoods(sc.milo, design = ~ Load, design.df = milo_design, fdr.weighting = "graph-overlap")
da.results.all.load <- annotateNhoods(sc.milo, da.results.all.load, coldata_col = "seurat_clusters")
da.results.all.load$seurat_clusters <- as.factor(as.numeric(da.results.all.load$seurat_clusters))
head(da.results.all.load)
da.results.all.load %>%
	arrange(SpatialFDR) %>%
	head()
head(da.results.all.load)

dim(da.results.all.load)
da.res <- da.results.all.load

# Get position with ggbeeswarm
beeswarm_pos <- ggplot_build(
	da.res %>%
		mutate(is_signif = ifelse(SpatialFDR < 0.1, 1, 0)) %>%
		arrange(seurat_clusters) %>%
		ggplot(aes(seurat_clusters, logFC)) +
		ggbeeswarm::geom_quasirandom()
)

pos_x <- beeswarm_pos$data[[1]]$x
pos_y <- beeswarm_pos$data[[1]]$y

n_groups <- unique(da.res$seurat_clusters) %>% length()

da.loadplot <- da.res %>%
	mutate(is_signif = ifelse(SpatialFDR < 0.1, 1, 0)) %>%
	mutate(logFC_color = ifelse(is_signif==1, logFC, 0)) %>%
	arrange(seurat_clusters) %>%
	mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
	mutate(pos_x = pos_x, pos_y=pos_y) %>%
	ggplot(aes(pos_x, pos_y, color=logFC_color)) +
	scale_colour_steps2(midpoint = 0, mid = "gray") +
	guides(color="none") +
	xlab("seurat_clusters") + ylab("Log Fold Change") +
	scale_x_continuous(
		breaks = seq(1,n_groups),
		labels = setNames(levels(da.res$seurat_clusters), seq(1,n_groups))
	) +
	ylim(c(-9, 9))+
	geom_point() +
	coord_flip() +
	theme_bw(base_size=22) +
	theme(strip.text.y =  element_text(angle=0))

plot(da.loadplot)

png(filename = paste0(rnaProject, "-miloRbeeswarm-Load.png"), height = 800, width = 800)
plot(da.loadplot)
dev.off()




# Seurat subset -----------------------------------------------------------

seurat.subset <- subset(seurat.object, idents = c("0", "3"), invert = TRUE)
seurat.subset@meta.data$seurat_clusters <- droplevels(seurat.subset@meta.data$seurat_clusters)
levels(seurat.subset$seurat_clusters)


sc.milo <- Milo(as.SingleCellExperiment(seurat.subset))

#try to add graph from Seurat's NN graph
graph(sc.milo) <- graph(buildFromAdjacency(seurat.subset@graphs$RNA_nn, k = 29))
sc.milo <- makeNhoods(sc.milo, prop = 0.05, k = 29, d=29, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(sc.milo)

sc.milo <- countCells(sc.milo, meta.data = data.frame(colData(sc.milo)), samples="orig.ident")


## Convert create design matrix
milo_design <- data.frame(colData(sc.milo))[,c("orig.ident", "Prep", "Load")]
milo_design$Prep <- as.factor(milo_design$Prep)
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$orig.ident
milo_design

sc.milo <- calcNhoodDistance(sc.milo, d=29, reduced.dim = "PCA")
da.results.subset.prep <- testNhoods(sc.milo, design = ~ Prep, design.df = milo_design, fdr.weighting = "graph-overlap")
da.results.subset.prep <- annotateNhoods(sc.milo, da.results.subset.prep, coldata_col = "Prep")
da.results.subset.prep %>%
	arrange(SpatialFDR) %>%
	head()
da.results.subset.prep <- annotateNhoods(sc.milo, da.results.subset.prep, coldata_col = "seurat_clusters")
da.results.subset.prep$seurat_clusters <- as.factor(as.numeric(da.results.subset.prep$seurat_clusters))

sc.milo <- buildNhoodGraph(sc.milo)

## Plot single-cell UMAP
umap_pl <- scater::plotReducedDim(sc.milo, dimred = "UMAP", colour_by="Prep", text_by = "seurat_clusters", text_size = 3, point_size=0.5) +
	guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sc.milo, da.results.subset.prep, layout="UMAP", alpha=0.1) 

umap_pl + nh_graph_pl +
	patchwork::plot_layout(guides="collect")

ggplot(da.results.subset.prep, aes(Prep_fraction)) + geom_histogram(bins=50)
ggplot(da.results.subset.prep, aes(seurat_clusters_fraction)) + geom_histogram(bins=50)
da.results.subset.prep$seurat_clusters_called <- ifelse(da.results.subset.prep$seurat_clusters < 0.7, "Mixed", da.results.subset.prep$seurat_clusters)
head(da.results.subset.prep, n = 20)

plotDAbeeswarm(da.results.subset.prep, group.by = "seurat_clusters")

# fancy plot
da.res <- da.results.subset.prep

# Get position with ggbeeswarm
beeswarm_pos <- ggplot_build(
	da.res %>%
		mutate(is_signif = ifelse(SpatialFDR < 0.1, 1, 0)) %>%
		arrange(seurat_clusters) %>%
		ggplot(aes(seurat_clusters, logFC)) +
		ggbeeswarm::geom_quasirandom()
)

pos_x <- beeswarm_pos$data[[1]]$x
pos_y <- beeswarm_pos$data[[1]]$y

n_groups <- unique(da.res$seurat_clusters) %>% length()

da.subsetplot <- da.res %>%
	mutate(is_signif = ifelse(SpatialFDR < 0.1, 1, 0)) %>%
	mutate(logFC_color = ifelse(is_signif==1, logFC, 0)) %>%
	arrange(seurat_clusters) %>%
	mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
	mutate(pos_x = pos_x, pos_y=pos_y) %>%
	ggplot(aes(pos_x, pos_y, color=logFC_color)) +
	scale_colour_steps2(midpoint = 0, mid = "gray80") +
	guides(color="none") +
	xlab("seurat_clusters") + ylab("Log Fold Change") +
	scale_x_continuous(
		breaks = seq(1,n_groups),
		labels = setNames(levels(da.res$seurat_clusters), seq(1,n_groups))
	) +
	ylim(c(-9, 9))+
	geom_point() +
	coord_flip() +
	theme_bw(base_size=22) +
	theme(strip.text.y =  element_text(angle=0))

plot(da.subsetplot)

png(filename = paste0(rnaProject, "-miloRbeeswarm-subsetPrep.png"), height = 800, width = 800)
plot(da.subsetplot)
dev.off()


