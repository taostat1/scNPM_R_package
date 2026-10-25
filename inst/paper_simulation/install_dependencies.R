options(repos=c(CRAN='https://cloud.r-project.org'))
if(!requireNamespace('remotes',quietly=TRUE))install.packages('remotes')
versions<-c(glasso='1.11',JGL='2.3.2',jsonlite='2.0.0')
for(p in names(versions))
  if(!requireNamespace(p,quietly=TRUE)||as.character(packageVersion(p))!=versions[[p]])
    remotes::install_version(p,version=versions[[p]],upgrade='never')
commit<-'2acabd4c583b49801460865ef7ccb8b944edc7fa'
work<-tempfile('BLGGM_');dir.create(work)
archive<-file.path(work,'source.zip')
download.file(paste0('https://github.com/WgitU/BLGGM/archive/',commit,'.zip'),archive,mode='wb')
unzip(archive,exdir=work)
src<-file.path(work,paste0('BLGGM-',commit))
for(f in c('Makevars','Makevars.win')) {
  path<-file.path(src,'src',f)
  if(file.exists(path))writeLines(gsub('CXX11','CXX14',readLines(path)),path)
}
remotes::install_local(src,upgrade='never')
