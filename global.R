rm(list=ls())

gc()

library(SeuratObject)
library(shinyjs)
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(ggplot2)
library(DT)
library(Seurat)



create_metadata_UMAP <- function(obj, col){
    if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")){
    col_df <- data.frame(obj@reductions$umap@cell.embeddings, data = obj@meta.data[,col])
    umap <- ggplot(data = col_df) +
        geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(data)), size = 0.01) +
        scale_colour_gradientn(colours = rainbow(7))
    }  else if (col %in% colnames(obj@meta.data)) {
      if (col == "seurat_clusters") {
        mycols <- c("MG1"="darkorange","MG2"="royalblue1", "MG3"="forestgreen", "MG4"="violet","MF"="yellow3","DC1"="greenyellow", "DC2"="darkturquoise", "DC3"="violetred","DC4"="darkcyan","Tcells1"="seagreen1","Tcells2"="green3","NK"="salmon3","NKL"="gray","ILC"="navy","Neutro"="red","Prol"="lightpink1","My1"="gray","My2"="gray","My3"="gray", "Mix"="gray")
        umap <- DimPlot(obj, reduction = "umap",cols=mycols , label = T, label.size = 4, pt.size = 1)
      } else { umap <- DimPlot(obj, pt.size = .1, label = F, label.size = 4, group.by = col, reduction = "umap")
      }
      
    } else {
    umap <- ggplot() +
        theme_void() +
        geom_text(aes(x = 0.5, y = 0.5, label = "Choose a feature"), size = 15, color = "gray73", fontface = "bold") +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    return(umap)
}



create_VlnPlot <- function(obj, gene, clusters) {
  
  if (length(gene) > 1) {
    
    violin_plot <-VlnPlot(obj,features=gene,combine=T, stack=T, cols= NULL, idents = clusters, flip = T,same.y.lims = F) + 
      theme(axis.title = element_text(size = 20,face="bold"),axis.text.x=element_text(size=20,face="bold",angle=90),axis.text.y.left = element_text(size=15, face="bold"),axis.line=element_line(size=1),axis.ticks=element_line(size=1))+
      NoLegend()
    
  } else {
    
    violin_plot <-VlnPlot(obj,features=gene,combine=T, stack=F, cols= NULL, idents = clusters, flip = T,same.y.lims = F, pt.size = 0) + 
      theme(axis.title = element_text(size = 20,face="bold"),axis.text.x=element_text(size=20,face="bold",angle=90),axis.text.y.left = element_text(size=15, face="bold"),axis.line=element_line(size=1),axis.ticks=element_line(size=1))+
      NoLegend()
    
  }
  gc()
  return(violin_plot)
  gc()
  
}




create_featurePlot <- function(obj, gene) {
  
  if (gene %in% rownames(obj)) {
    
    FeaturePlot(obj,features= gene,label=F,label.size=3, pt.size = 1.5, combine = F, cols = c("yellow", "blue"),order=T)
    
    
  } else { ggplot() + 
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Choose a gene"), size = 15, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))}
  
}



