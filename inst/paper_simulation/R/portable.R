# File handling and launch adapters; the statistical functions are in adjacent files.
stage <- function(folder,text) {
  line <- sprintf('[%s] %s',format(Sys.time(),'%Y-%m-%d %H:%M:%S'),text)
  writeLines(line,file.path(folder,'stage.txt'));cat(line,'\n');flush.console()
}
paper_data <- function(r,n) {
  seeds <- replicate_seeds(r,n)
  net <- make_networks(seeds$network)
  dat <- generate_expression(net,generate_original_means(DATA_CFG,seeds$means),seeds$dropout)
  if(n==1500L) dat <- fresh_expression(dat$truth,seeds$expression)
  net <- dat$truth$networks
  candidate <- Reduce(`|`,lapply(c(net$precision$T,net$precision$C),adjacency))
  up <- upper.tri(candidate);set.seed(seeds$ppi)
  ids <- take(which(up&candidate),round(.6*sum(up&candidate)))
  noise <- take(which(up&!candidate),round(length(ids)*.1/.9))
  ppi <- matrix(FALSE,100,100,dimnames=dimnames(candidate))
  ppi[c(ids,noise)]<-TRUE;ppi<-ppi|t(ppi)
  dat$observed$PPI<-ppi;net$ppi<-ppi;net$ppi_eligible<-candidate
  net$ppi_true_ids<-ids;net$ppi_false_ids<-noise
  net$ppi_origin<-'six_network_union';net$ppi_noise_pool<-'six_network_union_complement'
  net$config$ppi_reference<-'union';dat$truth$networks<-net
  dat$seeds<-seeds
  dat
}
prepare_paper <- function(output,n,reps) {
  for(r in reps) {
    folder<-file.path(output,paste0('N',n),'data',sprintf('replicate_%02d',r))
    files<-file.path(folder,c('observed_inputs.rds','evaluation_truth.rds'))
    manifest<-file.path(folder,'manifest.rds')
    if(file.exists(manifest)) {
      m<-readRDS(manifest)
      stopifnot(identical(m$seeds,replicate_seeds(r,n)),m$n==n,
                identical(m$md5,unname(tools::md5sum(files))))
      next
    }
    if(any(file.exists(files)))stop('Incomplete data directory: ',folder)
    dat<-paper_data(r,n);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    saveRDS(dat$observed,files[1]);saveRDS(dat$truth,files[2])
    saveRDS(list(n=n,seeds=dat$seeds,md5=unname(tools::md5sum(files))),manifest)
    write.csv(dat$truth$networks$audit,file.path(folder,'network_audit.csv'),row.names=FALSE)
    write.csv(dat$summary,file.path(folder,'expression_audit.csv'),row.names=FALSE)
    write.csv(as.data.frame(dat$seeds),file.path(folder,'seeds.csv'),row.names=FALSE)
    cat('Generated',n,'cells per condition; replicate',r,'\n')
  }
}
wsl_path <- function(path) {
  x<-normalizePath(path,winslash='/',mustWork=TRUE)
  if(!grepl('^[A-Za-z]:/',x))stop('WSL needs a Windows drive path: ',x)
  paste0('/mnt/',tolower(substr(x,1,1)),substring(x,3))
}
fit_ajgm <- function(obs,cfg,seed) {
  work<-file.path(tempdir(),paste0('ajgm_',Sys.getpid()))
  dir.create(work,recursive=TRUE,showWarnings=FALSE)
  for(d in c('T','C'))write.table(t(obs[[paste0('X_',d)]]),file.path(work,paste0(d,'.csv')),
    sep=',',row.names=FALSE,col.names=FALSE)
  py<-file.path(SIM_HOME,'python','ajgm_adapter.py')
  args<-c('-u',shQuote(py),shQuote(work),seed,cfg$ajgm_maxiter,cfg$ajgm_edge_cut)
  if(nzchar(cfg$wsl_distro)) {
    user<-if(nzchar(cfg$wsl_user))c('-u',cfg$wsl_user)else character()
    args<-c('-d',cfg$wsl_distro,user,'--exec','env',paste0('OMP_NUM_THREADS=',cfg$threads),
      paste0('OPENBLAS_NUM_THREADS=',cfg$threads),cfg$ajgm_python,'-u',shQuote(wsl_path(py)),
      shQuote(wsl_path(work)),seed,cfg$ajgm_maxiter,cfg$ajgm_edge_cut)
    status<-system2('wsl.exe',args)
  } else status<-system2(cfg$ajgm_python,args)
  if(status!=0)stop('AJGM process failed; inspect stdout/stderr.')
  P<-lapply(c('T','C'),function(d)lapply(1:3,function(k)
    as.matrix(read.csv(file.path(work,paste0(d,'_P',k,'.csv')),header=FALSE))))
  pc<-lapply(P,function(z)lapply(z,partial_cor))
  labels<-lapply(c('T','C'),function(d)as.integer(read.csv(file.path(work,paste0(d,'_labels.csv')),header=FALSE)[,1]))
  list(labels=labels,scores=lapply(pc,function(z)lapply(z,abs)),
    masks=lapply(pc,function(z)lapply(z,function(a)abs(a)>=cfg$ajgm_edge_cut)),
    mapping=match_graphs(pc[[1]],pc[[2]]),frozen_at=NA_integer_,
    selection='Prespecified absolute partial correlation cutoff',edge_cut=cfg$ajgm_edge_cut,
    convergence=lapply(c('T','C'),function(d)jsonlite::read_json(file.path(work,paste0(d,'_convergence.json')),simplifyVector=TRUE)))
}
fit_model <- function(obs,method,seed,folder,cfg=CFG) {
  ans<-fit_adapter(obs,method,seed,folder,cfg)
  if(method=='scNPM_independent_fixedNPM') {
    r<-as.integer(substr(basename(folder),1,2))
    ref<-readRDS(file.path(dirname(folder),sprintf('%02d_scNPM',r),'fit_summary.rds'))
    audit<-data.frame(replicate=r,same_mapping=identical(ans$mapping,ref$mapping),
      same_freeze_iteration=identical(ans$frozen_at,ref$frozen_at),
      same_adaptation=identical(ans$adaptation_audit,ref$adaptation_audit))
    write.csv(audit,file.path(folder,'fixed_mapping_audit.csv'),row.names=FALSE)
    if(!all(audit$same_mapping,audit$same_freeze_iteration,audit$same_adaptation))
      stop('Fixed-matching replay differs from the reference fit.')
  }
  ans
}
run_paper <- function(output,n,reps,methods,threads) {
  stopifnot(as.character(packageVersion('scNPM'))=='1.1.1')
  cfg<-CFG;cfg$threads<-threads
  getFromNamespace('set_scNPM_threads','scNPM')(threads)
  Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads)
  root<-file.path(output,paste0('N',n))
  for(r in reps)for(method in METHODS[METHODS%in%methods]) {
    folder<-file.path(root,'fits',sprintf('%02d_%s',r,method))
    result<-file.path(folder,'result.rds')
    if(file.exists(result)&&readRDS(result)$status=='ok')next
    if(method=='scNPM_independent_fixedNPM'&&!file.exists(file.path(root,'fits',sprintf('%02d_scNPM',r),'fit_summary.rds')))
      stop('Run scNPM first for this replicate.')
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    start<-Sys.time();stage(folder,paste('Starting',method,'N',n,'replicate',r))
    seed<-replicate_seeds(r,n)$fit
    saveRDS(cfg,file.path(folder,'settings.rds'))
    writeLines(capture.output(sessionInfo()),file.path(folder,'sessionInfo.txt'))
    z<-tryCatch({
      input<-file.path(root,'data',sprintf('replicate_%02d',r))
      obs<-readRDS(file.path(input,'observed_inputs.rds'))
      f<-fit_model(obs,method,seed,folder,cfg)
      truth<-readRDS(file.path(input,'evaluation_truth.rds'))
      ev<-evaluate_all(f,truth)
      saveRDS(f,file.path(folder,'fit_summary.rds'))
      write.csv(ev$graphs,file.path(folder,'overall_by_network.csv'),row.names=FALSE)
      write.csv(ev$metrics,file.path(folder,'metrics.csv'),row.names=FALSE)
      write.csv(ev$edge_audit,file.path(folder,'differential_edges.csv'),row.names=FALSE)
      list(status='ok',evaluation=list(summary=ev$summary,graphs=ev$graphs,metrics=ev$metrics),error='')
    },error=function(e)list(status='failed',error=conditionMessage(e)))
    z$method<-method;z$replicate<-r;z$seed<-seed;z$threads<-threads
    z$elapsed<-as.numeric(difftime(Sys.time(),start,units='secs'))
    saveRDS(z,result);stage(folder,paste('Finished',z$status,z$error))
    cat(sprintf('%s N=%d replicate=%d method=%s status=%s\n',Sys.time(),n,r,method,z$status),
        file=file.path(output,'progress.log'),append=TRUE)
    if(z$status!='ok')stop(z$error)
    gc(FALSE)
  }
}
