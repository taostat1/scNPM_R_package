.make_sampler_progress <- function(verbose, progress_file) {
  if (!verbose && is.null(progress_file)) return(NULL)
  run_id <- paste(format(Sys.time(), "%Y%m%dT%H%M%OS6"), Sys.getpid(), sep="-")
  clock <- function(s) {
    if (!is.finite(s)) return("estimating")
    s <- max(0, floor(s))
    sprintf("%02d:%02d:%02d", as.integer(s %/% 3600), as.integer(s %/% 60 %% 60), as.integer(s %% 60))
  }
  function(record) {
    if (verbose) message(sprintf(
      "MCMC iteration %d/%d (%.1f%%) | %s | elapsed %s | ETA %s | stable %d/%d | freeze %s | retained %d/%d",
      record$iteration,record$total,100*record$iteration/record$total,record$stage,
      clock(record$elapsed_seconds),clock(record$eta_seconds),record$stable_count,
      record$stable_required,record$freeze_reason,
      record$retained_draws,record$retained_target))
    if (verbose) utils::flush.console()
    if (!is.null(progress_file)) {
      row <- as.data.frame(c(list(run_id=run_id,timestamp=format(Sys.time(), "%Y-%m-%d %H:%M:%S%z")),record),
        stringsAsFactors=FALSE)
      header <- !file.exists(progress_file) || file.info(progress_file)$size == 0
      utils::write.table(row,progress_file,sep=",",row.names=FALSE,col.names=header,
        append=!header,quote=TRUE,na="")
    }
    invisible(NULL)
  }
}
