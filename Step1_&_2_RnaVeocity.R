library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(velocyto.R)
library(dplyr)
library(ggpubr)


##  AstroNeurontest - seurat object with defined clusters

# filtering out N-stage7,A-stage5,6,7
small_data <- subset(AstroNeurontest,idents = c("NSC-stage1","NSC-stage2","RG-like","N-stage1",
                                                "N-stage2","N-stage3","N-stage4","N-stage5","N-stage6",
                                                "A-stage1","A-stage2","A-stage3","A-stage4"))


## separating control and TBI cells ##
control <- subset(small_data,cells = 1:1412)
TBI <- subset(small_data,cells = 1413:3314)


#### Adding spliced and unspliced data to seurat object ####

##load loom files
S5_ldat <- read.loom.matrices("S5_10X_05042019.loom") #control
S6_ldat <- read.loom.matrices("S6_10X_05042019.loom") #TBI

##changing cell names according to AstroNeurontest
#control
colnames(S5_ldat$spliced) <- gsub("S5_10X_05042019:","",S5_ldat$spliced%>%colnames())%>%gsub(pattern = "x",replacement = ".1_1")
colnames(S5_ldat$unspliced) <- colnames(S5_ldat$spliced)
colnames(S5_ldat$ambiguous) <- colnames(S5_ldat$spliced)
#TBI
colnames(S6_ldat$spliced) <- gsub("S6_10X_05042019:","",S6_ldat$spliced%>%colnames())%>%gsub(pattern = "x",replacement = ".1_2")
colnames(S6_ldat$unspliced) <- colnames(S6_ldat$spliced)
colnames(S6_ldat$ambiguous) <- colnames(S6_ldat$spliced)

## subsetting cells according to small_data
#control
S5_cells <- colnames(control)
S5_ldat$spliced <- S5_ldat$spliced[,S5_cells]
S5_ldat$unspliced <- S5_ldat$unspliced[,S5_cells]
S5_ldat$ambiguous <- S5_ldat$ambiguous[,S5_cells]
#TBI
S6_cells <- colnames(TBI)
S6_ldat$spliced <- S6_ldat$spliced[,S6_cells]
S6_ldat$unspliced <- S6_ldat$unspliced[,S6_cells]
S6_ldat$ambiguous <- S6_ldat$ambiguous[,S6_cells]

## adding to seurat object
#control
S5_Sdat <- as.Seurat(x = S5_ldat)
control[["spliced"]] <- S5_Sdat[["spliced"]]
control[["unspliced"]] <- S5_Sdat[["unspliced"]]
control[["ambiguous"]] <- S5_Sdat[["ambiguous"]]
#TBI
S6_Sdat <- as.Seurat(x = S6_ldat)
TBI[["spliced"]] <- S6_Sdat[["spliced"]]
TBI[["unspliced"]] <- S6_Sdat[["unspliced"]]
TBI[["ambiguous"]] <- S6_Sdat[["ambiguous"]]




####  RunVelocity  ####
#control
control <- RunVelocity(object = control, deltaT = 1, kCells = 20, 
                       fit.quantile = 0.02,spliced.average = 0.5,ncores = 10) #Step1
cell.colors_5 <- sccore::fac2col(control@active.ident)
control_velo <- show.velocity.on.embedding.cor(emb = Embeddings(object = control, reduction = "umap"),
                                               vel = Tool(object = control,slot = "RunVelocity"), 
                                               n = 300, scale = "sqrt", cell.colors = ac(x = cell.colors_5, alpha = 0.5), 
                                               cex = 0.8, arrow.scale = 1, show.grid.flow = TRUE, 
                                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                                               do.par = FALSE, cell.border.alpha = 0.1) #Step2

#TBI
TBI <- RunVelocity(object = TBI, deltaT = 1, kCells = 20, fit.quantile = 0.02,spliced.average = 0.5,ncores = 10)
cell.colors_6 <- sccore::fac2col(TBI@active.ident)
TBI_velo <- show.velocity.on.embedding.cor(emb = Embeddings(object = TBI, reduction = "umap"),
                                           vel = Tool(object = TBI,slot = "RunVelocity"), 
                                           n = 300, scale = "sqrt", cell.colors = ac(x = cell.colors_6, alpha = 0.5), 
                                           cex = 0.8, arrow.scale = 1, show.grid.flow = TRUE, 
                                           min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                                           do.par = FALSE, cell.border.alpha = 0.1)

