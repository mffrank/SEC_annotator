library(ggplot2)
library(data.table)
library(reshape)
library(plotly)
library(DT)

# library(devtools)
# install_github("PedrioliLab/SECprofiler",ref="peakpicker")
# library('SECprofiler')

# library(Cairo)   # For nicer ggplot2 output when deployed on Linux
source("plot.Features.R")
source("tracesMethods.R")

ann_mode<- "all"
# traces <- readRDS("data/pepTracses.annotated.consec3.SibPepCorr.filtered.all.rda")
traces <- readRDS("data/protTraces_raw_a_cs3_spcf_fdr1pc_top2.rda")

if(traces$trace_type == "peptide"){
  if(ann_mode == "features"){
    # features <- readRDS("data/ProteinFeatures.tab.pcorr08.c05_subs.rds")
    load("data/ProteinFeatures")
    features <- ProteinFeatures
    hyp_names <- unique(features$protein_name)
  } else if(ann_mode == "all"){
    # load("data/pepTracses.annotated.consec3.SibPepCorr.filtered")
    # traces <- pepTracses.annotated.consec3.SibPepCorr.filtered
    hyp_names <- unique(traces$trace_annotation$protein_id)
  }
  
} else {
  # features <- filtered_corum_table
  features <- readRDS("data/corum_table_filtered_subset_current.rda")
  # Filter out those proteins that have no trace
  features <- subset(features, protein_id %in% traces$traces$id)
  hyp_names <- unique(features$complex_name)
  
}

instance <- round(runif(1)*10^7,0)



server <- function(input, output, session) {
  annotations <- reactiveValues()
  annotations$dt <- data.table(FeatureName = numeric(0),
                               leftBoundary = numeric(0),
                               rightBoundary = numeric(0),
                               Confidence = numeric(0),
                               Apex = numeric(0))
  ranges <- reactiveValues(x = NULL, y = NULL)
  apex <- reactiveValues(sel = NULL, bound_left = NULL, bound_right = NULL)
  
  
  ## Observe the Button Events ------------------
  ###################################################
  
  plot.data <- eventReactive(input$select,{
    complexID <- input$select
    ranges$x <- NULL
    apex$sel <- NULL
    apex$bound_right <- NULL
    apex$bound_left <- NULL
    if(traces$trace_type == "peptide"){
      if(ann_mode == "features"){
        plot.peptideFeatures.data(features,traces,
                                  complexID)
      } else if(ann_mode == "all"){
        peptides <- traces$trace_annotation$id[traces$trace_annotation$protein_id == complexID]
        pepTraces_sel <- subset(traces,trace_ids=peptides)
        toLongFormat(pepTraces_sel$traces)
      }
    } else {
      proteins <- features$protein_id[features$complex_name == complexID]
      Traces_sel <- subset(traces,trace_ids = proteins)
      toLongFormat(Traces_sel$traces)
    }
  })
  
  observeEvent(input$prev, {
    pos <- which(input$select == hyp_names)
    if(pos != 1){
      updateSelectInput(session, "select", selected = hyp_names[pos-1])
      apex$sel <- NULL
      apex$bound_right <- NULL
      apex$bound_left <- NULL
      ranges$x <- NULL
    }
  })
  
  observeEvent(input$forw, {
    pos <- which(input$select == hyp_names)
    if(pos != length(hyp_names)){
      updateSelectInput(session, "select", selected = hyp_names[pos+1])
      apex$sel <- NULL
      apex$bound_right <- NULL
      apex$bound_left <- NULL
      ranges$x <- NULL
    }
  })
  

  observeEvent(input$save, {
    if(!is.null(apex$sel) & !is.null(apex$bound_right)){
      newLine <- isolate(data.table(FeatureName = input$select,
                                    Apex = apex$sel,
                                    leftBoundary = apex$bound_left,
                                    rightBoundary = apex$bound_right,
                                    Confidence = "Low"
      ))
      isolate(annotations$dt <- rbind(newLine,annotations$dt))
      apex$sel <- NULL
      apex$bound_right <- NULL
      apex$bound_left <- NULL
    }
    # ranges$x <- NULL
    # pos <- which(input$select == hyp_names)
    # if(pos != length(hyp_names)){
    #   updateSelectInput(session, "select", selected = hyp_names[pos+1])
    # }
    
    })
  
  observeEvent(input$save_sure, {
    if(!is.null(apex$sel) & !is.null(apex$bound_right)){
      newLine <- isolate(data.table(FeatureName = input$select,
                                    Apex = apex$sel,
                                    leftBoundary = apex$bound_left,
                                    rightBoundary = apex$bound_right,
                                    Confidence = "High"
      ))
      isolate(annotations$dt <- rbind(newLine,annotations$dt))
      apex$sel <- NULL
      apex$bound_right <- NULL
      apex$bound_left <- NULL
    }
    # ranges$x <- NULL
    # pos <- which(input$select == hyp_names)
    # if(pos != length(hyp_names)){
    #   updateSelectInput(session, "select", selected = hyp_names[pos+1])
    # }
    
    })
  
  observeEvent(input$unzoom, {
    ranges$x <- NULL
    # ranges$y <- NULL
  })
  
  observeEvent(input$saveall, {
    saveRDS(isolate(annotations$dt), file = paste0("data/Manual_annotation",gsub("-","_",Sys.Date()),"_",instance,".rda"))
  })
  
  observeEvent(input$delete, {
    if(!is.null(input$view_rows_selected)){
      annotations$dt<- annotations$dt[-input$view_rows_selected]
    }
  })
  
  observeEvent(input$keypress[1],{
    if(input$keypress[1] == "3" |input$keypress[1] == "ArrowRight"){
      pos <- which(input$select == hyp_names)
      if(pos != length(hyp_names)){
        updateSelectInput(session, "select", selected = hyp_names[pos+1])
        apex$sel <- NULL
        apex$bound_right <- NULL
        apex$bound_left <- NULL
        ranges$x <- NULL
      } 
    } else if(input$keypress[1] == "ArrowLeft"){
      pos <- which(input$select == hyp_names)
      if(pos != 1){
        updateSelectInput(session, "select", selected = hyp_names[pos-1])
        apex$sel <- NULL
        apex$bound_right <- NULL
        apex$bound_left <- NULL
        ranges$x <- NULL
      }
    } else if(input$keypress[1] == "2"){
      if(!is.null(apex$sel) & !is.null(apex$bound_right)){
        newLine <- isolate(data.table(FeatureName = input$select,
                                      Apex = apex$sel,
                                      leftBoundary = apex$bound_left,
                                      rightBoundary = apex$bound_right,
                                      Confidence = "Low"
        ))
        isolate(annotations$dt <- rbind(newLine,annotations$dt))
        apex$sel <- NULL
        apex$bound_right <- NULL
        apex$bound_left <- NULL
      }
    } else if(input$keypress[1] == "1"){
      if(!is.null(apex$sel) & !is.null(apex$bound_right)){
        newLine <- isolate(data.table(FeatureName = input$select,
                                      Apex = apex$sel,
                                      leftBoundary = apex$bound_left,
                                      rightBoundary = apex$bound_right,
                                      Confidence = "High"))
        isolate(annotations$dt <- rbind(newLine,annotations$dt))
        apex$sel <- NULL
        apex$bound_right <- NULL
        apex$bound_left <- NULL
      }
      
    }
    print(input$keypress)
  })
  ## Observe the Mouse Events in Plots ---------------------------
  ###############################################################
  
  observeEvent(input$plot2_dblclick, {
    brush <- input$plot2_brush
    if (is.null(brush)) {
      ranges$x <- NULL
    }
  })

  observeEvent(input$plot2_brush, {
    brush <- input$plot2_brush
    ranges$x <- c(round(brush$xmin,0), round(brush$xmax,0))
  })
  
  
  observeEvent(input$plot1_click, {
    # if (!is.null(input$plot1_brush)) {
      apex$sel <- round(input$plot1_click$x,digits = 0)
    # }
  })
  
  observeEvent(input$plot1_brush, {
      apex$bound_left <- round(input$plot1_brush$xmin,digits = 0)
      apex$bound_right <- round(input$plot1_brush$xmax,digits = 0)
  })

  ## Output the Plots and Data table -------------------------
  ###############################################
  
  
  output$plot1 <- renderPlot({
    # print(apex_sel())
    if(traces$trace_type == "peptide"){
      if(ann_mode == "features"){
        p <- plot.complexFeatures(plot.data(),ranges,apex_man = apex, plot_peak = F, log = F,
                                  plot_monomer = F, plot_apex = F)
      } else if(ann_mode == "all"){
        p <- plot.traces(plot.data(),ranges,apex_man = apex, plot=FALSE, ledgend = FALSE)
      } 
    } else if (traces$trace_type == "protein"){
      p <- plot.traces(plot.data(),ranges,apex_man = apex, plot=FALSE,
                       ledgend = FALSE, trace_type = "Protein")
    }
    p
  })
  
  output$plot2 <- renderPlot({
    if(ann_mode == "features"){
      p <- plot.complexFeatures(plot.data(), plot_peak = F, log = T,
                                plot_monomer = F, plot_apex = F,
                                plot.title = F) 
    } else if(ann_mode == "all"){
      p <- plot.traces(plot.data(), plot=TRUE,
                       ledgend = FALSE,log = TRUE)
    } 
    p
  })
  
  output$view <-DT::renderDataTable(annotations$dt)
  
  # ## Used for debugging, no longer needed-------------------
  # output$click_info <- renderPrint({
  #   # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
  #   # were a base graphics plot, we'd need those.
  #   # nearPoints(traces, input$plot1_click, addDist = TRUE)
  #   list(is.null(plot.data()),
  #        input$select,
  #        input$plot1_brush$xmin,
  #        input$plot1_brush$xmax,
  #        apex$sel,
  #        session)
  # })
  
  # output$brush_info <- renderPrint({
  #   # brushedPoints(mtcars2, input$plot1_brush)
  #   input$plot1_brush
  # })
  #----------------------------------
  
}
