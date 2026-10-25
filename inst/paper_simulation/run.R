#!/usr/bin/env Rscript
# Command-line entry. Paths are resolved relative to this file, not the working directory.
args<-commandArgs(FALSE)
script<-sub('^--file=','',args[grepl('^--file=',args)])
if(length(script)!=1L)stop('Run with Rscript path/to/run.R; see README.md.')
SIM_HOME<-dirname(normalizePath(script,winslash='/',mustWork=TRUE))
cli<-commandArgs(TRUE)
if(!length(cli)||'--help'%in%cli) {
  cat('Rscript run.R --phase prepare|run|summarize|reference --n both|3000|1500\n',
      '  --output PATH --replicates 1,2,...,10 --methods all|scNPM,... --threads 8\n',
      'Run is serial; independent-prior fits follow their matching reference.\n')
  quit(status=0)
}
if(length(cli)%%2L)stop('Every option needs a value.')
opts<-setNames(as.list(cli[seq(2,length(cli),2)]),cli[seq(1,length(cli),2)])
if(anyDuplicated(names(opts))||any(!names(opts)%in%c('--phase','--n','--output','--replicates','--methods','--threads')))
  stop('Unknown or duplicate option; use --help.')
getopt<-function(key,default)if(is.null(opts[[key]]))default else opts[[key]]
phase<-getopt('--phase','prepare')
if(!phase%in%c('prepare','run','summarize','reference'))stop('Invalid phase.')
sizes<-getopt('--n','both');sizes<-if(sizes=='both')c(3000L,1500L)else as.integer(sizes)
reps<-as.integer(strsplit(getopt('--replicates',paste(1:10,collapse=',')),',',fixed=TRUE)[[1]])
threads<-as.integer(getopt('--threads','8'))
stopifnot(length(sizes)>0,all(sizes%in%c(3000L,1500L)),!anyNA(reps),all(reps%in%1:10),
          !anyDuplicated(reps),length(threads)==1L,!is.na(threads),threads>=1L)
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads)
for(file in c('config','generation','baselines','evaluation','methods','portable','summary'))
  source(file.path(SIM_HOME,'R',paste0(file,'.R')),local=globalenv())
methods<-getopt('--methods','all');methods<-if(methods=='all')METHODS else strsplit(methods,',',fixed=TRUE)[[1]]
stopifnot(all(methods%in%METHODS))
output<-getopt('--output',file.path(getwd(),'paper_simulation_output'))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
output<-normalizePath(output,winslash='/',mustWork=TRUE)
if(phase=='reference') {
  raw<-read.csv(file.path(SIM_HOME,'reference_results','all_replicates.csv'))
  summarize_table(raw,output)
} else if(phase=='summarize') {
  summarize_paper(output)
} else {
  for(n in sizes)prepare_paper(output,n,reps)
  if(phase=='run') {
    required<-c('scNPM',if('BLGGM'%in%methods)'BLGGM',
      if(any(methods%in%c('glasso','mixtureGGM','JGL_fused')))'glasso',
      if('JGL_fused'%in%methods)'JGL',if('AJGM'%in%methods)'jsonlite')
    for(pkg in required)if(!requireNamespace(pkg,quietly=TRUE))stop('Install dependency: ',pkg)
    for(n in sizes)run_paper(output,n,reps,methods,threads)
    summarize_paper(output)
  }
}
