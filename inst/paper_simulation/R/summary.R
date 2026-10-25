summarize_paper <- function(output) {
  rows<-states<-list()
  for(n in c(3000L,1500L))for(r in 1:10)for(method in METHODS) {
    p<-file.path(output,paste0('N',n),'fits',sprintf('%02d_%s',r,method),'result.rds')
    if(!file.exists(p))next
    z<-readRDS(p)
    states[[length(states)+1L]]<-data.frame(N=n,replicate=r,method=method,status=z$status,error=z$error)
    if(z$status!='ok')next
    for(scope in c('overall_micro','differential_directional')) {
      row<-z$evaluation$metrics[z$evaluation$metrics$scope==scope,]
      ap<-if(scope=='overall_micro')z$evaluation$graphs$AP else
        z$evaluation$metrics$AP[grepl('^differential_T',z$evaluation$metrics$scope)]
      stopifnot(length(ap)==if(scope=='overall_micro')6L else 3L)
      row$N<-n;row$method<-method;row$replicate<-r;row$mAP<-mean(ap)
      row$matching_accuracy<-unname(z$evaluation$summary['matching_accuracy'])
      rows[[length(rows)+1L]]<-row
    }
  }
  if(length(states))write.csv(do.call(rbind,states),file.path(output,'results_status.csv'),row.names=FALSE)
  if(!length(rows))stop('No completed fits to summarize.')
  raw<-do.call(rbind,rows)
  write.csv(raw,file.path(output,'all_replicates.csv'),row.names=FALSE)
  summarize_table(raw,output)
}
summarize_table <- function(raw,output) {
  stopifnot(!anyDuplicated(raw[,c('N','scope','method','replicate')]))
  fields<-c('precision','recall','F1','mAP','AP','FDP','TP','FP','FN','selected','matching_accuracy')
  tab<-paired<-list()
  for(n in c(3000L,1500L))for(scope in c('differential_directional','overall_micro'))for(m in METHODS) {
    x<-raw[raw$N==n&raw$scope==scope&raw$method==m,]
    if(!nrow(x))next
    row<-data.frame(N=n,scope=scope,method=m,n=nrow(x))
    for(v in fields) {
      a<-x[[v]];a<-a[is.finite(a)];nn<-length(a)
      row[[paste0(v,'_n')]]<-nn
      row[[paste0(v,'_mean')]]<-if(nn)mean(a)else NA_real_
      row[[paste0(v,'_SE')]]<-if(nn>1)sd(a)/sqrt(nn)else NA_real_
    }
    tab[[length(tab)+1L]]<-row
    if(m=='scNPM')next
    a<-raw[raw$N==n&raw$scope==scope&raw$method=='scNPM',]
    b<-merge(a,x,by='replicate',suffixes=c('_scNPM','_other'))
    if(!nrow(b))next
    for(v in c('precision','recall','F1','mAP')) {
      d<-b[[paste0(v,'_scNPM')]]-b[[paste0(v,'_other')]];d<-d[is.finite(d)]
      paired[[length(paired)+1L]]<-data.frame(N=n,scope=scope,comparison=paste('scNPM -',m),
        metric=v,n=length(d),mean=if(length(d))mean(d)else NA_real_,
        SE=if(length(d)>1)sd(d)/sqrt(length(d))else NA_real_)
    }
  }
  write.csv(do.call(rbind,tab),file.path(output,'summary.csv'),row.names=FALSE)
  if(length(paired))write.csv(do.call(rbind,paired),file.path(output,'paired_differences.csv'),row.names=FALSE)
  invisible(do.call(rbind,tab))
}
