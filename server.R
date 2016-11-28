library(ggplot2)
library(data.table)
library(reshape)
library(plotly)
library(DT)

# library(Cairo)   # For nicer ggplot2 output when deployed on Linux
source("plot.complexFeatures_peptides.R")
source("tracesMethods.R")

ann_mode<- "all"

if(ann_mode == "features"){
  
  # features <- readRDS("data/ProteinFeatures.tab.pcorr08.c05_subs.rds")
  load("data/ProteinFeatures")
  features <- ProteinFeatures
  hyp_names <- unique(features$protein_name)
  
  pepTraces.annotated.consec3.SibPepCorr.filtered <- readRDS("data/pepTracses.annotated.consec3.SibPepCorr.filtered.all.rda")

  # load("data/pepTracses.annotated.consec3.SibPepCorr.filtered")
  # pepTraces.annotated.consec3.SibPepCorr.filtered <- pepTracses.annotated.consec3.SibPepCorr.filtered
  # ## NOTE: Traces were characters, need to be converted to numeric first!
  # traces.numeric<- subset(pepTraces.annotated.consec3.SibPepCorr.filtered$traces, select=-id) 
  # traces.numeric[, names(traces.numeric) := lapply(.SD, as.numeric)]
  # # sapply(traces.numeric, class)
  # pepTraces.annotated.consec3.SibPepCorr.filtered$traces <- data.table(traces.numeric,
  #                                                                      id = pepTraces.annotated.consec3.SibPepCorr.filtered$traces$id)
} else if(ann_mode == "all"){
  # load("data/pepTracses.annotated.consec3.SibPepCorr.filtered")
  # pepTraces.annotated.consec3.SibPepCorr.filtered <- pepTracses.annotated.consec3.SibPepCorr.filtered
  
  pepTraces.annotated.consec3.SibPepCorr.filtered <- readRDS("data/pepTracses.annotated.consec3.SibPepCorr.filtered.all.rda")
  hyp_names <- unique(pepTraces.annotated.consec3.SibPepCorr.filtered$trace_annotation$protein_id)
  #   pepTraces.annotated.consec3.SibPepCorr.filtered <- pepTracses.annotated.consec3.SibPepCorr.filtered
  # ## NOTE: Traces were characters, need to be converted to numeric first!
  # traces.numeric<- subset(pepTraces.annotated.consec3.SibPepCorr.filtered$traces, select=-id)
  # traces.numeric[, names(traces.numeric) := lapply(.SD, as.numeric)]
  # # sapply(traces.numeric, class)
  # pepTraces.annotated.consec3.SibPepCorr.filtered$traces <- data.table(traces.numeric,
  #                                                                      id = pepTraces.annotated.consec3.SibPepCorr.filtered$traces$id)

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
    if(ann_mode == "features"){
      plot.complexFeatures.data(features,pepTraces.annotated.consec3.SibPepCorr.filtered,
                                complexID)
    } else if(ann_mode == "all"){
      peptides <- pepTraces.annotated.consec3.SibPepCorr.filtered$trace_annotation$id[pepTraces.annotated.consec3.SibPepCorr.filtered$trace_annotation$protein_id == complexID]
      pepTraces_sel <- subset(pepTraces.annotated.consec3.SibPepCorr.filtered,trace_ids=peptides)
      toLongFormat(pepTraces_sel$traces)
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
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  # observeEvent(input$zoom, {
  #   brush <- input$plot2_brush
  #   if (!is.null(brush)) {
  #     ranges$x <- c(brush$xmin, brush$xmax)
  #     # ranges$y <- c(brush$ymin, brush$ymax)
  #   }
  # })
  
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
  # apex_sel <- reactive(input$plot1_click$x)
  
  # output$firstPlot <- reactive({
  #   # return(is.null(output$plot1))
  #   TRUE
  # })
  
  
  ## Output the Plots and Data table -------------------------
  ###############################################
  
  
  output$plot1 <- renderPlot({
    # print(apex_sel())
    if(ann_mode == "features"){
      p <- plot.complexFeatures(plot.data(),ranges,apex_man = apex, plot_peak = F, log = F,
                                plot_monomer = F, plot_apex = F)
    } else if(ann_mode == "all"){
      p <- plot.traces(plot.data(),ranges,apex_man = apex, plot=FALSE, ledgend = FALSE)
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
  #   # nearPoints(pepTraces.annotated.consec3.SibPepCorr.filtered, input$plot1_click, addDist = TRUE)
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
