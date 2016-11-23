#' Subset a traces.obj by trace_ids or fraction_ids.
#' @param traces.obj An object of type \code{traces.obj}.
#' @param trace_ids A character vector specifying the trace identifiers
#'        for subsetting \code{traces.obj}.
#' @param fraction_ids A numeric vector specifying the fraction identifiers
#'        for subsetting \code{traces.obj}.
#' @return traces.obj An object of type \code{traces.obj}.
#' @export
subset.traces <- function(traces.obj,trace_ids=NULL,fraction_ids=NULL){
    if (!is.null(trace_ids)) {
      traces.obj$traces <- subset(traces.obj$traces,id %in% trace_ids)
      traces.obj$trace_annotation <- subset(traces.obj$trace_annotation ,id %in% trace_ids)
    }
    if (!is.null(fraction_ids)){
      traces.obj$traces <- traces.obj$traces[,c(fraction_ids,"id"),with=FALSE]
      traces.obj$fraction_annotation <- subset(traces.obj$fraction_annotation ,id %in% fraction_ids)
    }
    traces.obj
}

#' Get a matrix of intensity values from a traces object.
#' @param traces.obj An object of type \code{traces.obj}.
#' @return A matrix with intensity values.
#' @export
getIntensityMatrix <- function(traces.obj) {
    ids <- traces.obj$traces$id
    intensity.mat <- as.matrix(sapply(subset(traces.obj$traces,
                                      select=-id),as.numeric))
    rownames(intensity.mat) <- ids
    intensity.mat
}


#' Convert a data.table containing traces from wide format to long format.
#' @param traces.dt A data.table with an id column \code{id} and
#'        columns of continuously numbered fractions.
#' @return A data.table with columns
#'          \itemize{
#'           \item \code{id}
#'           \item \code{fraction}
#'           \item \code{intensity}
#'          }
toLongFormat <- function(traces.dt) {
  traces.dt.long <-
    melt(traces.dt, id.var='id', variable.name='fraction',
         value.name='intensity', variable.factor=FALSE)
  traces.dt.long[, fraction := as.numeric(fraction)]
  setkey(traces.dt.long,id)
  data.table(traces.dt.long)
  traces.dt.long
}

#' Plot a traces.obj.
#' @param traces.obj An object of type \code{traces.obj}.
#' @export
plot.traces <- function(traces.long,
                        ranges = NULL,
                        apex_man = NULL,
                        plot=TRUE,
                        ledgend = TRUE,
                        log = FALSE,
                        trace_type = "Peptide") {
  
  if(!is.null(ranges$x[1])){
  traces.long.s <- subset(traces.long,fraction %in% ranges$x[1]:ranges$x[2] )
  }  else {
    traces.long.s <- traces.long
  }

    pl <- ggplot(traces.long.s) +
      ggtitle(paste(trace_type, 'traces')) +
      xlab('fraction') + ylab('intensity') +
      coord_cartesian(xlim = ranges$x)
      
    pl <- pl + geom_line(aes(x=fraction, y=intensity, color=id))
    if (!ledgend) {
      pl <- pl + theme(legend.position="none")
    }
    if (plot) plot(pl)
    if (log) {
      pl <- pl + scale_y_log10('log(intensity)')
    }
    if(!is.null(apex_man$sel)){
      pl <- pl + geom_vline(xintercept=apex_man$sel, colour="darkgrey", linetype="solid")
    }
    if (!is.null(apex_man$bound_right)){
      pl <- pl + geom_vline(xintercept=apex_man$bound_left,linetype="dashed")
      pl <- pl + geom_vline(xintercept=apex_man$bound_right,linetype="dashed")
    }
    pl
}
