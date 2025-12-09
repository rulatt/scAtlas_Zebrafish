source('global.R')
obj <- readRDS("data/immune_atlas_dataset.rds")
genes_names = as.vector(read.table("data/gene_names.txt")$V1)
gc()

ui <- dashboardPage(
  dashboardHeader(title = "scRNAseq Analysis"),
  dashboardSidebar(
    tags$head(
      tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))
    ),
    sidebarMenu(
      id = 'tab',
      useShinyjs(),
      menuItem("Home Page", tabName = "home", icon = icon("list")),
      menuItem("Run scRNAseq Analysis", tabName = "input", icon = icon("edit"))
    )
  ), 
  dashboardBody(
    tabItems(
      tabItem(tabName = "input",
              tabsetPanel(
                id = 'main_tabs',
                tabPanel("Summary",
                         HTML("<h3>Summary</h3><br>
                              <p>UMAP: Split Uniform Manifold Approximation Projection (UMAP) related to Figure 2B</p>
                              <p>Gene expression: UMAP plots depicting gene expression related to Figures 2C, 3A, 4A, S4</p>
                              <p>Violin Plot: Violin plot depicting gene expression related to Figures 3D, 4D</p>"
                         )),
                tabPanel("UMAP",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'umap'),
                             downloadButton("downloadUMAP", "Download UMAP")
                           ),
                           column(width = 4,selectizeInput("MetadataUMAP","Metadata", choices = NULL)
                           )
                         )
                ),
                
                tabPanel("Gene Expression",
                         fluidRow(
                           
                           column(
                             width = 6,
                             selectizeInput("gene1", "Genes", choices = NULL),
                             plotOutput(outputId = 'featurePlot1'),
                             downloadButton("downloadFeaturePlot1", "Download Feature Plot")
                             
                           ),
                         )
                ),
                tabPanel("Violin Plot",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'ViolinPlot'),
                             downloadButton("downloadViolinPlot", "Download Violin Plot")
                           ),
                           column(
                             width = 4,
                             selectizeInput("gene", 
                                            "Genes", 
                                            choices = NULL
                                            , multiple = FALSE)
                           ),
                           column(
                             width = 4,
                             selectizeInput("clusters", 
                                            "Clusters", 
                                            choices = NULL,
                                            multiple = TRUE)
                           ),
                           column(
                             width = 4,
                             actionButton("runButtonV", 
                                          "Run" 
                             )
                           )
                           
                         )
                ),
              )
      ),
      tabItem(tabName = "home",
              tags$h1(HTML("<u>A single-cell transcriptomic atlas reveals resident dendritic-like cells in the zebrafish brain parenchyma</u>"),
                      HTML("<h3>Mireia Rovira, Giuliano Ferrero, Magali Miserocchi, Alice Montanari, Ruben Lattuca ,Valérie Wittamer</h3><br>
                              <p>Institut de Recherche Interdisciplinaire en Biologie Humaine et Moléculaire (IRIBHM)</p>
                              <p>ULB Institute of Neuroscience (UNI), Université Libre de Bruxelles (ULB), Brussels, Belgium</p>"
                      ))
      )
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 5000 * 1024^2)
  
  
  output$umap <- renderPlot({
    if (!is.null(input$MetadataUMAP)) {
      create_metadata_UMAP(obj, input$MetadataUMAP)
    }
  })
  
  output$featurePlot1 <- renderPlot({
    if (!is.null(input$gene1)) {
      create_featurePlot(obj, input$gene1)
    }
  })
  
  generate_ViolinPlot <- eventReactive(input$runButtonV, {
    gc()
    create_VlnPlot(obj, input$gene, input$clusters)
  })
  
  gc()
  
  output$ViolinPlot <- renderPlot({
    gc()
    generate_ViolinPlot()
  })
  
  gc()
  
  
  output$downloadFeaturePlot1 <- downloadHandler(
    filename = function(){
      paste0(input$gene1, '_feature_plot', '.png')
    },
    content = function(file){
      plot <- create_featurePlot(obj, input$gene1)
      ggsave(filename=file, width = 10, height = 5, type = "cairo")
    }
  )
  
  
  output$downloadUMAP <- downloadHandler(
    filename = function(){
      paste0('UMAP', '.png')
    },
    content = function(file){
      plot <- create_metadata_UMAP(obj, input$MetadataUMAP)
      ggsave(filename=file, width = 10, height = 5, type = "cairo")
    }
  )
  
  output$downloadViolinPlot <- downloadHandler(
    filename = function(){
      paste0('violin_plot', '.png')
    },
    content = function(file){
      plot <- create_VlnPlot(obj, input$gene, input$clusters)
      ggsave(filename=file, width = 10, height = 5, type = "cairo")
    }
  )
  
  output$downloadVolcanoPlot <- downloadHandler(
    filename = function(){
      paste0('volcano_plot', '.png')
    },
    content = function(file){
      plot <- generate_diff_genes()$plotDE
      ggsave(filename=file, width = 10, height = 5, type = "cairo")
    }
  )
  
  
  updateSelectizeInput(session, "MetadataUMAP", choices = c("seurat_clusters", "nCount_RNA","nFeature_RNA"), server = TRUE)
  updateSelectizeInput(session, "gene1", choices = genes_names, server = TRUE)
  updateSelectizeInput(session, "gene", choices = genes_names, server = TRUE)
  
  updateSelectizeInput(session, "clusters", choices = c("MG1","MG2","MG3","MG4","MF","DC1","DC2","DC3","DC4","Tcells1","Tcells2", "NK","ILC","Neutro","Prol") , selected = c("MG1","MG2","MG3","MG4","MF","DC1","DC2","DC3","DC4","Tcells1","Tcells2", "NK","ILC","Neutro","Prol") ,server = TRUE)
  updateSelectizeInput(session, "clusters_1", choices = c("MG1","MG2","MG3","MG4","MF","DC1","DC2","DC3","DC4","Tcells1","Tcells2", "NK","ILC","Neutro","Prol") , server = TRUE)
  updateSelectizeInput(session, "clusters_2", choices = c("MG1","MG2","MG3","MG4","MF","DC1","DC2","DC3","DC4","Tcells1","Tcells2", "NK","ILC","Neutro","Prol") , server = TRUE)
  
  
}


shinyApp(ui, server)