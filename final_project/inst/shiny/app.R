options(shiny.maxRequestSize = 30*1024^2)

#load required libraries
library(plotrix)
library(shiny)
library(ggplot2)
library(colourpicker)
library(tidyverse)
library(rlang)
library(dplyr)
library(shinyjs)
library(DT)
library(reader)
library(abind)
library(pheatmap)
library(broom)
library(DESeq2)
library(scales)
library(RColorBrewer)
library(ggbeeswarm)
library(reshape2)
library(cowplot)
library(psych)

ui <- navbarPage(
  
  title = "Gene Expression Exploration Project",
  id="Gene Expression Exploration Project",
  fluid=TRUE,
  
  #Set up UI for Sample Exploration tabs
  tabPanel("Sample Metadata Exploration",
           
           # Application title
           titlePanel("Sample Metadata Exploration"),
           
           # Place for uploading data
           sidebarLayout(
             sidebarPanel(
               h3("Upload metadata data table"),
               tags$div(tags$p(
                 'Data table should be a matrix of metadata information. Metadata information must be submitted here for use in other tabs.'
               )),
               fileInput(
                 "meta",
                 "Metadata table",
                 multiple = FALSE,
                 accept = ".csv"
               ),
               actionButton(inputId = 'm_submit',label = 'Upload'
               ),
               selectizeInput('metadata_col','Choose variable to plot',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               tags$div(tags$p(
                 "Variable in metadata table. Selected variable must be of numeric values."
               )),
               actionButton(inputId = 'metadata_plot',label = 'Plot'
               )
             ),
             
             # Show a table of the inputted data
             mainPanel(
               tabsetPanel(
                 tabPanel('Input Matrix',
                          h4(strong("Data Preview")),
                          DTOutput('meta_header'),
                 ),
                 tabPanel('Summary Table',
                          DTOutput('summary')
                 ),
                 tabPanel('Summary Plots',
                          plotOutput('metadata_hist')
                 )
               )
             )
           )),
  
  #Set up UI for Counts Exploration tabs
  tabPanel("Expression Counts Exploration",
           
           # Application title
           titlePanel("Expression Counts Exploration"),
           
           # Place for uploading data
           sidebarLayout(
             sidebarPanel(
               h3("Upload counts data table"),
               tags$div(tags$p(
                 'Data table should be a matrix of normalized counts.'
               )),
               fileInput(
                 "counts",
                 "Counts matrix",
                 multiple = FALSE,
                 accept = ".csv"
               ),
               sliderInput(inputId="v_slider", 
                           label="Select the gene expression variance quantile threshold to filter genes by:", 
                           min=0, max=1, value=0.5, step = NULL,round = FALSE,ticks = TRUE),
               sliderInput(inputId="e_slider", 
                           label="Select the gene non-zero count threshold to filter genes by:", 
                           min=0, max=1, value=0.5, step = NULL,round = FALSE,ticks = TRUE),
               selectizeInput('pca_1','Choose principle component to plot',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               selectizeInput('pca_2','Choose principle component to plot',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               selectizeInput('pca_col','Choose variable to color PCA plot',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               selectizeInput('top_PCAs','Choose number of top principle components to plot',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               tags$div(tags$p(
                 "Plots and filter statistics may take some time to load for large datasets."
               )),
               actionButton(inputId = 'c_submit',label = 'Filter & Plot'
               )
             ),
             
             # Show a table of the inputted data
             mainPanel(
               tabsetPanel(
                 tabPanel('Input Matrix',
                          h4(strong("Data Preview")),
                          DTOutput('counts_header'),
                 ),
                 tabPanel('Filter Statistics',
                          tableOutput('filter_stats'),
                 ),
                 tabPanel('Scatter Plot',
                          plotOutput('scattermv'),
                          plotOutput('scattermz'),
                 ),
                 tabPanel('Heatmap',
                          plotOutput('count_hmap'),
                 ),
                 tabPanel('PCA',
                          plotOutput('pca'),
                 ),
                 tabPanel('Beeswarm',
                          plotOutput('bees'),
                 )
               )
             )
           )),
  
  #Set up UI for Differential Expression Exploration tabs
  tabPanel("Differential Expression Exploration",
           
           # Application title
           titlePanel("Differential Expression Exploration"),
           
           # Place for uploading data
           sidebarLayout(
             sidebarPanel(
               h3("Upload differential expression results data table"),
               tags$div(tags$p(
                 'Data table should be a table of differential expression results.'
               )),
               fileInput(
                 "de_res",
                 "Differential Expression Analysis Results Table",
                 multiple = FALSE,
                 accept = ".csv"
               ),
               selectizeInput('volcano_col1','Choose log fold change variable to plot',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               selectizeInput('volcano_col2','Choose p-value variable to plot',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               sliderInput(inputId="de_slider", 
                           label="Select the magnitude of significance value coloring:", 
                           min=-300, max=0, value=-150, step = NULL,round = FALSE,ticks = TRUE),
               actionButton(inputId = 'de_plot',label = 'Plot'
               )
             ),
             
             # Show a table of the inputted data
             mainPanel(
               tabsetPanel(
                 tabPanel('Input Matrix',
                          h4(strong("Data Preview")),
                          DTOutput('de_res_header'),
                 ),
                 tabPanel('Volcano Plot',
                          plotOutput('volcano_p'),
                 )
               )
             )
           )),
  
  #Set up UI for Individual Gene Counts Exploration tabs
  tabPanel("Individual Gene Counts Exploration",
           
           # Application title
           titlePanel("Individual Gene Counts Exploration"),
           
           # Place for uploading data
           sidebarLayout(
             sidebarPanel(
               h3("Upload counts data table"),
               tags$div(tags$p(
                 'Data tables should be a matrix of normalized counts associated with the previously uploaded metadata. Selection drop-downs may take a moment to load for large datasets.'
               )),
               fileInput(
                 "individual_counts",
                 "Counts matrix",
                 multiple = FALSE,
                 accept = ".csv"
               ),
               selectizeInput('individual_counts_row','Choose gene to individually visualize',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               selectizeInput('individual_metadata_col','Choose metadata variable to group gene counts by',
                              multiple=F,choices = c(''),selected = NULL,
                              options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))
               ),
               selectizeInput('plot_type','Choose plot type',
                              multiple=F,choices = c('bar', 'pie'),selected = NULL
               ),
               actionButton(inputId = 'individual_plot',label = 'Plot'
               )
             ),
             
             mainPanel(
               tabsetPanel(
                 tabPanel('Input Data',
                          h4(strong("Data Preview")),
                          DTOutput('individual_counts_header'),
                 ),
                 tabPanel('Chosen Plot',
                          plotOutput('individual_plot'),
                          DTOutput('individual_stats')
                 )
               )
             )
           ))
)
server <- function(input, output, session) { 
  #set up reactive values for user inputs
  reactivevalue=reactiveValues(meta=NULL,
                               meta_location=NULL,
                               metadata=NULL,
                               counts=NULL,
                               counts_location=NULL,
                               counts_matrix=NULL,
                               pca_fit=NULL,
                               de_res=NULL,
                               de_res_location=NULL,
                               de_res_matrix=NULL,
                               individual_counts=NULL,
                               icl=NULL,
                               individual_counts_matrix=NULL)
  
  #function for updating Sample Information Exploration selections
  setupSelections1 = function(){
    updateSelectizeInput(session=session, inputId="metadata_col",
                         choices=colnames(reactivevalue$metadata),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
  }
  
  #function for updating Counts Exploration selections
  setupSelections2 = function(){
    updateSelectizeInput(session=session, inputId="pca_col",
                         choices=colnames(reactivevalue$metadata),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
    updateSelectizeInput(session=session, inputId="pca_1",
                         choices=colnames(reactivevalue$pca_fit$x),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
    updateSelectizeInput(session=session, inputId="pca_2",
                         choices=colnames(reactivevalue$pca_fit$x),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
    updateSelectizeInput(session=session, inputId="top_PCAs",
                         choices=c(1:length(colnames(reactivevalue$pca_fit$x))),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
  }
  
  #function for updating Differential Expression Exploration selections
  setupSelections3 = function(){
    updateSelectizeInput(session=session, inputId="volcano_col1",
                         choices=colnames(reactivevalue$de_res_matrix),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
    updateSelectizeInput(session=session, inputId="volcano_col2",
                         choices=colnames(reactivevalue$de_res_matrix),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
    
  }
  
  #function for updating Individual Gene Counts Exploration selections
  setupSelections4 = function(){
    updateSelectizeInput(session=session, inputId="individual_counts_row",
                         choices=rownames(reactivevalue$individual_counts_matrix),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
    updateSelectizeInput(session=session, inputId="individual_metadata_col",
                         choices=colnames(reactivevalue$metadata),selected=NULL,
                         options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))
    
  }
  #Create metadata table preview  
  observeEvent( input$meta, {
    req(input$meta)
    reactivevalue$meta_location=input$meta$datapath
    reactivevalue$metadata=read.table(reactivevalue$meta_location,header = T,row.names = 1,
                                      sep = get.delim(reactivevalue$meta_location,n = 10,delims = c('\t',',')),check.names = F)
    output$meta_header=renderDT((datatable(reactivevalue$metadata)))
    setupSelections1()
  }
  )
  #Create sample summary table when data table is submitted
  observeEvent(
    eventExpr = input[["m_submit"]],
    handlerExpr = {
      #get vectors for column names, column types, and mean/standard deviation or unique values
      reactivevalue$metadata[is.na(reactivevalue$metadata)] = 0
      CN <- colnames(reactivevalue$metadata)
      Type <- sapply(reactivevalue$metadata, class)
      MDV <- c()
      for (i in 1:length(colnames(reactivevalue$metadata))) {
        if (Type[i]=="numeric") {
          st <- paste("+/-", round(sd(reactivevalue$metadata[,i])))
          MDV <- c(MDV,paste(round(mean(reactivevalue$metadata[,i])), st))
        }
        else if (Type[i]=="integer") {
          st <- paste("+/-", round(sd(reactivevalue$metadata[,i])))
          MDV <- c(MDV,paste(round(mean(reactivevalue$metadata[,i])), st))
        }
        else {
          MDV <- c(MDV,paste(unique(reactivevalue$metadata[,i]), sep=", ", collapse=", "))
        }
      }
      #make summary table
      sum_tab <- abind(CN,Type,MDV,along=2)
      col_names <- c('Column Name','Type','Mean (sd) or Distinct Values')
      colnames(sum_tab) <- col_names
      output$summary <- renderDT((datatable(sum_tab)))
    } 
  )
  
  #Create sample summary plot when plot is pressed
  observeEvent(
    eventExpr = input[["metadata_plot"]],
    handlerExpr = {
      col_to_plot <- which(colnames(reactivevalue$metadata)==input$metadata_col)
      output$metadata_hist <- renderPlot(hist(reactivevalue$metadata[,col_to_plot],main = paste("Histogram of" , input$metadata_col),xlab = input$metadata_col))
    } 
  )
  
  #Create counts table preview  
  observeEvent( input$counts, {
    req(input$counts)
    reactivevalue$counts_location=input$counts$datapath
    reactivevalue$counts_matrix=read.table(reactivevalue$counts_location,header = T,row.names = 1,
                                           sep = get.delim(reactivevalue$counts_location,n = 10,delims = c('\t',',')),check.names = F)
    
    reactivevalue$counts_matrix[is.na(reactivevalue$counts_matrix)] = 0
    
    output$counts_header=renderDT((datatable(reactivevalue$counts_matrix)))
    
    #pca is calculated to get vector of PCs for PC selection update from pca_fit reactive value
    pca_counts_matrix <- as_tibble(t(reactivevalue$counts_matrix[rowSums(reactivevalue$counts_matrix[])>0,]))
    
    reactivevalue$pca_fit <- pca_counts_matrix %>% 
      select(where(is.numeric)) %>% 
      prcomp(scale = TRUE)
    
    setupSelections2()
  }
  )
  
  #Create filter statistics table, filtered data scatter plots, filtered data heatmap 
  #PCA plot, and beeswarm plot when button is pressed
  observeEvent(
    eventExpr = input[["c_submit"]],
    handlerExpr = {
      #get vector of variance per gene
      gene_vars <- apply(
        reactivevalue$counts_matrix,
        MARGIN = 1,var)
      #get vector of median per gene
      gene_meds <- apply(
        reactivevalue$counts_matrix,
        MARGIN = 1,median)
      #get variance threshold from input percentile
      v_threshold <- quantile(gene_vars, input$v_slider)
      
      #get logical vector of genes that pass variance threshold
      to_filter <- apply(
        reactivevalue$counts_matrix,
        MARGIN = 1,var)>=v_threshold
      
      cords <- which(to_filter==TRUE)
      
      #save variance threshold genes for genes passing table
      fills1 <- to_filter
      
      #get logical vector of genes that pass zero count threshold for genes passing table
      fills2 <- apply(
        reactivevalue$counts_matrix,
        MARGIN = 1,
        function(x){length(which(x!=0))/length(x)})>=input$e_slider
      #get names of passed genes
      to_filter <- names(to_filter[cords])
      #filter matrix by passed genes
      filtered_counts_matrix <- reactivevalue$counts_matrix[to_filter,]
      
      #get logical vector of genes that pass zero count threshold for matrix filtering
      to_filter <- apply(
        filtered_counts_matrix,
        MARGIN = 1,
        function(x){length(which(x!=0))/length(x)})>=input$e_slider
      
      cords <- which(to_filter==TRUE)
      
      #make table of filter logical vectors 
      fills <- cbind(fills1,fills2)
      
      #make new logical vector for genes passing both filters
      passed <- c()
      
      for (i in 1:length(rownames(fills))) {
        if (FALSE %in% fills[i,]) {
          passed <- c(passed,FALSE)
        }
        else {
          passed <- c(passed,TRUE)
        }
      }
      
      to_filter <- names(to_filter[cords])
      #filter counts matrix for a second time for genes passing zero count filter
      filtered_counts_matrix <- filtered_counts_matrix[to_filter,]
      
      #make filter statistics table
      Statistic <- c('Number of Samples','Total Number of Genes','Number of Genes Passing','Percent of Genes Passing','Number of Genes Not Passing','Percent of Genes Not Passing')
      
      Value <- round(c(ncol(reactivevalue$counts_matrix),nrow(reactivevalue$counts_matrix),nrow(filtered_counts_matrix),(nrow(filtered_counts_matrix)/nrow(reactivevalue$counts_matrix))*100,nrow(reactivevalue$counts_matrix)-nrow(filtered_counts_matrix),(nrow(reactivevalue$counts_matrix)-nrow(filtered_counts_matrix))/nrow(reactivevalue$counts_matrix)*100))
      
      output$filter_stats <- renderTable(cbind(Statistic,Value))
      
      #make scatter plots of log10 median vs variance 
      gene_meds[gene_meds == 0] <- 1
      gene_meds <- log10(gene_meds)
      
      to_plotmv <- data.frame(abind(gene_vars,gene_meds,passed,along=2,make.names=TRUE))
      
      output$scattermv <- renderPlot(ggplot(to_plotmv, aes(x=gene_meds, y=gene_vars, color=as.character(passed))) +
                                       scale_colour_manual(name = '0:filtered 1:not filtered', values = c("1" = "blue4", "0" ="cyan3")) +
                                       ggtitle("Log10 Median Count  vs Variance") + xlab("Log10 Median Count") + ylab("Variance") +
                                       geom_point())
      #make scatter plots of log10 median vs zero counts
      #the loop makes a vector of columns with zero counts for each gene
      zero_counts <- c()
      for (i in 1:length(rownames(reactivevalue$counts_matrix))) {
        zero_counts <- c(zero_counts,length(which(reactivevalue$counts_matrix[i,]==0)))
      }
      
      to_plotmz <- data.frame(abind(zero_counts,gene_meds,passed,along=2,make.names=TRUE))
      
      output$scattermz <- renderPlot(ggplot(to_plotmz, aes(x=gene_meds, y=zero_counts, color=as.character(passed))) +
                                       scale_colour_manual(name = '0:filtered 1:not filtered', values = c("1" = "blue4", "0" ="cyan3")) +
                                       ggtitle("Log10 Median Count  vs Number of Zeroes") + xlab("Log10 Median Count") + ylab("Zero Counts") +
                                       geom_point())
      
      #make heatmap
      #if statement restricts gene variances to plot to the top 20 most variable
      #unless there are less than 20 genes in the filtered matrix. In that case, 
      #half the most variable genes are plotted
      
      if (length(rownames(filtered_counts_matrix))<20) {
        most_variable_genes <- names(tail(sort(gene_vars),(length(rownames(filtered_counts_matrix))/2)))
      }
      else {
        most_variable_genes <- names(tail(sort(gene_vars),20))
      }
      
      output$count_hmap <- renderPlot(pheatmap(as.matrix(filtered_counts_matrix)[c(most_variable_genes),],main='Heatmap of Filtered Counts in Highly Variable Genes'))
      
      #make PCA plot
      pca_counts_matrix <- as_tibble(t(reactivevalue$counts_matrix))
      
      col_to_plot <- which(colnames(reactivevalue$metadata)==input$pca_col)
      
      groups <- reactivevalue$metadata[,col_to_plot]
      
      #adds user input metadata column for data point coloring
      pca_counts_matrix <- pca_counts_matrix %>%
        mutate(
          groups_col = groups,
        )
      
      pca_fit <- pca_counts_matrix %>% 
        select(where(is.numeric)) %>% 
        prcomp(scale = TRUE)
      
      pca_model <- pca_fit %>%
        augment(pca_counts_matrix)
      
      #gets the principle component numbers from user inputs
      pc_1 <- as.numeric(str_replace_all(input$pca_1, "PC", ""))
      pc_2 <- as.numeric(str_replace_all(input$pca_2, "PC", ""))
      
      eigs <- pca_fit$sdev^2
      
      #the principle component number taken from the user input is added to 1 and the number of 
      #columns in the count matrix to ge the principle component column indexes
      output$pca <- renderPlot(pca_model %>% 
                                 ggplot(aes(as.matrix(pca_model[,ncol(pca_counts_matrix)+1+pc_1]), 
                                            as.matrix(pca_model[,ncol(pca_counts_matrix)+1+pc_2]), 
                                            color = groups_col)) + 
                                 ggtitle("Gene Expression PCA Plot") + xlab(paste(input$pca_1,round(eigs[pc_1] / sum(eigs), digits = 2))) + ylab(paste(input$pca_2,round(eigs[pc_2] / sum(eigs), digits = 2))) +
                                 geom_point(size = 1.5) +
                                 theme_half_open(12) + background_grid())
      
      #make the beeswarm plot
      #the principle component columns are isolated from the augmented matrix for plotting
      for_bees <- as.matrix(pca_model)
      
      for_bees <- for_bees[,-c(1:(nrow(reactivevalue$counts_matrix)+2))]
      
      for_bees <- for_bees[,c(1:input$top_PCAs)]
      
      output$bees <- renderPlot(ggplot(data = melt(as.data.frame(for_bees),id.vars=NULL),aes(x = variable, y = value, color = variable, labels = NULL)) +
                                  geom_beeswarm(cex = 3))
      
      
    } 
  )
  
  #Create differential expression analysis results table
  observeEvent( input$de_res, {
    req(input$de_res)
    reactivevalue$de_res_location=input$de_res$datapath
    reactivevalue$de_res_matrix=read.table(reactivevalue$de_res_location,header = T,row.names = 1,
                                           sep = get.delim(reactivevalue$de_res_location,n = 10,delims = c('\t',',')),check.names = F)
    
    reactivevalue$de_res_matrix[is.na(reactivevalue$de_res_matrix)] = 0
    
    output$de_res_header=renderDT((datatable(reactivevalue$de_res_matrix)))
    
    setupSelections3()
  }
  )
  
  #Create volcano plot when button is pressed
  observeEvent(
    eventExpr = input[["de_plot"]],
    handlerExpr = {
      
      volcano_data <- as.data.frame(reactivevalue$de_res_matrix)
      slider_factor <- (1 * (10^input$de_slider))
      
      #gets the index of the user input columns to plot
      col1 <- which(colnames(volcano_data)==input$volcano_col1)
      col2 <- which(colnames(volcano_data)==input$volcano_col2)
      
      volcano_data <- volcano_data %>%
        mutate(slider_cond = case_when(volcano_data[, col2] < slider_factor ~ "TRUE", 
                                       volcano_data[, col2] >= slider_factor ~ "FALSE", TRUE ~ 'NA'))
      
      output$volcano_p <- renderPlot(ggplot(data=volcano_data, aes(x=volcano_data[,col1], y=-log10(volcano_data[,col2]), color=slider_cond)) + 
                                       geom_point() +
                                       scale_color_manual(values = c('FALSE' = 'red', 'TRUE' = 'orange', 'NA'='black')) + 
                                       xlab(input$volcano_col1) +
                                       ylab(input$volcano_col2) +
                                       labs(title=paste0('Volcano Plot: ', input$volcano_col1, '  VS  -log10(', input$volcano_col2, ')' )) +
                                       theme(legend.position="bottom"))
      
    } 
  )
  
  #Create counts matrix
  observeEvent( input$individual_counts, {
    req(input$individual_counts)
    reactivevalue$icl=input$individual_counts$datapath
    reactivevalue$individual_counts_matrix=read.table(reactivevalue$icl,header = T,row.names = 1,
                                                      sep = get.delim(reactivevalue$icl,n = 10,delims = c('\t',',')),check.names = F)
    
    reactivevalue$individual_counts_matrix[is.na(reactivevalue$individual_counts_matrix)] = 0
    
    
    output$individual_counts_header=renderDT((datatable(reactivevalue$individual_counts_matrix)))
    
    setupSelections4()
  }
  )
  
  #Creates bar and pie chart options as well as summary statistcs table
  observeEvent(
    eventExpr = input[["individual_plot"]],
    handlerExpr = {
      
      #gets column index of user input metadata and subsets this column
      chosen_metadata_idx <- which(colnames(reactivevalue$metadata)==input$individual_metadata_col)
      chosen_metadata <- reactivevalue$metadata[,chosen_metadata_idx]
      
      chosen_data <- reactivevalue$individual_counts_matrix
      
      #changes column names of count matrix to user input metadata
      colnames(chosen_data) <- chosen_metadata
      
      #makes matrix with counts summed by metadata groups
      grouped <- t(rowsum(t(chosen_data), colnames(chosen_data)))
      grouped <- as.data.frame(grouped)
      
      #gets row index of user input gene and gets this row
      chosen_gene_idx <- which(rownames(grouped)==input$individual_counts_row)
      chosen_gene_grouped <- grouped[chosen_gene_idx,]
      
      chosen_gene <- chosen_data[chosen_gene_idx,]
      
      #gets metadata group names
      col_groups <- unique(colnames(chosen_gene))
      
      #gets the variance per metadata group
      group_vars <- c()
      for (i in col_groups) {
        idxs <- which(colnames(chosen_gene)==i)
        
        group_vars <- c(group_vars,var(as.numeric(chosen_gene[,c(idxs)])))
        
      }
      
      # gets the mean count per metadata group
      group_means <- c()
      for (i in col_groups) {
        idxs <- which(colnames(chosen_gene)==i)
        
        group_means <- c(group_means,mean(as.numeric(chosen_gene[,c(idxs)])))
        
      }
      
      #set up data for grouped bar graph
      data <- c(group_means,group_vars)
      
      labs <- c(rep("Group mean count",length(col_groups)),rep("Group variance",length(col_groups)))
      
      col_groups <- rep(c(col_groups),length(col_groups)) 
      
      grouped_bar <- data.frame(labs,col_groups,data)
      
      #plot options
      #makes grouped bar chart if selected by user
      if (input$plot_type=='bar') {
        output$individual_plot <- renderPlot(ggplot(grouped_bar, aes(fill=col_groups, y=data, x=labs)) + 
                                               geom_bar(position="dodge", stat="identity")+
                                               xlab("Count Statistic") +
                                               ylab("Statistic Value") +
                                               labs(title="Grouped Bar Chart of Average Gene Expression Count and Gene Expression Variance in Each Metadata Group"))}
      
      #makes pie chart if selected by user
      else if (input$plot_type=='pie') {
        
        output$individual_plot <- renderPlot(pie3D(as.numeric(chosen_gene_grouped), 
                                                   explode=0.1, 
                                                   theta=pi/3, 
                                                   labels = colnames(chosen_gene_grouped),
                                                   col=c(brewer.pal(length(colnames(chosen_gene_grouped)),name = 'Pastel1')),
                                                   main="Pie Chart of Gene Expression Counts in Each Metadata Group"))
        
      }
      #gets summary statistics by metadata group and makes data table 
      stats <- describeBy(t(chosen_gene), colnames(chosen_gene))
      stats <- datatable(do.call(rbind.data.frame, stats))
      chosen_gene_group <- t(chosen_gene_grouped)
      colnames(chosen_gene_group) <- "count"
      var <- group_vars
      stats <- datatable(abind(var,chosen_gene_group,stats$x$data[-c(1,2)],along = 2,make.names = T))
      output$individual_stats <- renderDT(stats)
    } 
  )
  
  
}

shinyApp(ui = ui, server = server)



