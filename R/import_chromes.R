
import_chromes <- function(paths){
  dat<-list()
  for (path in paths){
    files <- list.files(path=path, pattern = "CSV",full.names = T)
    file_names <- gsub(pattern = ".CSV",x = basename(files), replacement = "")
    mydata <- lapply(X=files, FUN=read.csv, row.names = 1, header=T,
                     fileEncoding="utf-16",check.names = F)
    mydata <- lapply(mydata, FUN=as.matrix)
    names(mydata) <- file_names
    dat <- append(dat,mydata)
  } 
  dat
}
