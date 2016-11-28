

jscode <- '
$(function() {
$(document).keypress(function(e) {
Shiny.onInputChange("keypress", [e.key, Math.random()]);

});
});
'

#Load data
ann_mode <- "all"

if(ann_mode == "features"){
  # features <- readRDS("data/ProteinFeatures.tab.pcorr08.c05_subs.rds")
  load("data/ProteinFeatures")
  features <- ProteinFeatures
  
  hyp_names <- unique(features$protein_name)
  
} else if(ann_mode == "all"){
  # load("data/pepTracses.annotated.consec3.SibPepCorr.filtered")
  pepTraces.annotated.consec3.SibPepCorr.filtered <- readRDS("data/pepTracses.annotated.consec3.SibPepCorr.filtered.all.rda")
  hyp_names <- unique(pepTraces.annotated.consec3.SibPepCorr.filtered$trace_annotation$protein_id)
}


ui <- fluidPage(
  tags$head(tags$script(HTML(jscode))),
  titlePanel("SEC-WATH MS Protein Viewer"),
  fluidRow(
    column(4,
      wellPanel(
        h3("Select Protein"),
        br(),
        selectizeInput("select", label = h5("UniprotID"), 
                    choices = hyp_names,
                    options = list(maxOptions=10000)),
        actionButton("prev","",icon = icon("arrow-left")),
        actionButton("forw","",icon = icon("arrow-right")),
        br(),
        br(),
        fluidRow(
          # actionButton("zoom",label = "Zoom"),
          actionButton("unzoom",label = "Reset Zoom")
        ),
        br(),
        fluidRow(
          actionButton("save_sure", label = "Save Feature",icon = icon("thumbs-o-up")),
          actionButton("save", label = "Save Feature",icon = icon("thumbs-o-down"))
        ),
        actionButton("saveall", label = "Save Data Frame")
      ),
      plotOutput("plot2",
                 dblclick = "plot2_dblclick",
                 brush = brushOpts(id = "plot2_brush",
                                   direction = "x",
                                   resetOnNew = TRUE),
                 height = 200
                 ),
      h4("How to annotate Features"),
      tags$ul(
        tags$li("Use small plot to zoom to desired region"), 
        tags$li("Drag your mouse over the peak boundaries"), 
        tags$li("Double-click the peak apex"),
        tags$li("Press Save Feature and check the table"),
        tags$li("Repeat for all features in a plot, then go to the next plot"),
        tags$li("Delete wrong entries by selecting them and pressing delete"),
        tags$li("Press Save Data Frame often, to save a copy of the table")
      )
      
    ),
    column(8,
      # conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
      #       h4("The data is loading...")
      # ),
           plotOutput("plot1",
                 # Equivalent to: click = clickOpts(id = "plot_click")
                 dblclick = "plot1_click",
                 brush = brushOpts(id = "plot1_brush",
                                   direction = "x",
                                   resetOnNew = TRUE)
      ),
      
      # ## Debugging
      # fluidRow(
      #   column(width = 6,
      #          h6("system vars"),
      #          verbatimTextOutput("click_info")
      #   )),
      fluidRow(
        div(actionButton("delete", "Delete",icon = icon("remove-sign",lib="glyphicon")), style="float:right")
      ),
      # div(style="display:inline-block",actionButton("delete", "Delete"), style="float:right"),
      br(),
      DT::dataTableOutput("view")
    )
  )
)
