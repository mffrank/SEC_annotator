#' Plot the result of the sliding window algorithm.
#' @param sw.result An object of type 'complexFeaturesSW'.
#' @param trace.mat The same matrix that was passend to
#'        \code{findComplexFeatures}.
#' @param protein.names The names of proteins whose traces are contained in
#'        the matrix \code{trace.mat}.
#' @param n.largest Only plot the n largest subgroups. Optional.
#' @export
plot.complexFeatures.data <- function(res,
                                 traces.obj,
                                 proteinID
                                 ) {
  
    features <- subset(res, protein_id == proteinID)
    peptides <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
    traces <- subset(traces.obj,trace_ids = peptides)
    traces.long <- toLongFormat(traces$traces)
    
    # add monomer MW
    monomer = as.numeric(unlist(features$monomer_sec))
    monomer_max = as.numeric(unlist(features$monomer_sec_max))
    
    setkey(traces.long, id)
    
    proteinName = unique(features$protein_name)[1]
    
    group.colors <- c("TRUE" = "green3", "FALSE" = "firebrick1")
    
    n_annotatedSubunits = features$n_subunits_annotated[1]
    n_signalSubunits = features$n_subunits_with_signal[1]
    n_detectedSubunits = features$n_subunits_detected[1]
    complexCompleteness = round(features$completeness[1],digits=2)
    if(!("in_complex" %in% colnames(features))){
      features$in_complex <- features$apex < monomer_max
    }
    res <-list(traces.long = traces.long,
               proteinName = proteinName,
               n_annotatedSubunits = n_annotatedSubunits,
               n_signalSubunits = n_signalSubunits,
               n_detectedSubunits = n_detectedSubunits,
               complexCompleteness = complexCompleteness,
               group.colors = group.colors,
               monomer = monomer,
               monomer_max = monomer_max,
               features = features)
  
    return(res)
    
}

plot.complexFeatures <- function(plot.data,ranges = NULL,
                                 apex_man = NULL,
                                 plot_peak=TRUE,
                                 plot_monomer=TRUE,
                                 log=FALSE,
                                 plot_apex=TRUE,
                                 plot.title =TRUE) {
  list2env(plot.data, environment())
  if(!is.null(ranges$x[1])){
  traces.long.s <- subset(traces.long,fraction %in% ranges$x[1]:ranges$x[2] )
  # print(ranges$x[1]:ranges$x[2])
  # print(traces.long.s)
  }  else {
    traces.long.s <- traces.long
  }
  p <- ggplot(traces.long.s) +
    geom_line(aes_string(x='fraction', y='intensity', color='id')) +
    xlab('fraction') +
    ylab('intensity') +
    guides(color=FALSE) +
    coord_cartesian(xlim = ranges$x)
  
  if(plot.title){
    p <- p + ggtitle(proteinName)
  }
  if (log) {
    p <- p + scale_y_log10('log(intensity)')
  }
  
  if(plot_peak==TRUE){
    p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = 0, ymax = Inf, fill=in_complex),alpha = 0.25) +
      scale_fill_manual(values=group.colors)
    p <- p + geom_vline(data=features,aes(xintercept = left_pp), colour="darkgrey", linetype="dashed")
    p <- p + geom_vline(data=features,aes(xintercept = right_pp), colour="darkgrey",linetype="dashed")
  }
  
  if (plot_monomer==TRUE){
    p <- p + geom_vline(xintercept=monomer,linetype="solid")
    p <- p + geom_vline(xintercept=monomer_max,linetype="dashed")
  }
  if (!is.null(apex_man$bound_right)){
    p <- p + geom_vline(xintercept=apex_man$bound_left,linetype="dashed")
    p <- p + geom_vline(xintercept=apex_man$bound_right,linetype="dashed")
  }
  if(plot_apex==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="darkgrey", linetype="solid")
  }
  if(!is.null(apex_man$sel)){
    p <- p + geom_vline(data=features,aes(xintercept=apex_man$sel), colour="darkgrey", linetype="solid")
  }
  #p <- p + theme(legend.position="none")
  return(p)
}


