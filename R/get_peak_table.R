# get_peak_table <- function(pks,sigma.r, sigma.t)
# rt <- sapply(pks, function(x){
#   rt <- x[[1]]$rt
#   rt
# }, simplify=T)
# for (i in 1:length(rt)){
# names(rt[[i]]) <- paste(names(rt)[[i]],rt[[i]], sep='_')
# }
# d.rt <- as.matrix(dist(unlist(rt)))
# 
# sp <- sapply(1:length(dat.pr), function(i){
#   sp <- scales::rescale(t(dat.pr[[i]][(rt[[i]]-.09)*100,]))
#   colnames(sp)<-paste(i, colnames(sp),sep='_')
#   sp
# }, simplify=T)
# sp2 <- do.call(cbind,sp)
# 
# c <- cor(sp2,sp2,method = "pearson")
# sigma.r=1
# #sigma.t=45
# sigma.t=1
# S<-exp((-(1-abs(c))^2)/(2*sigma.r^2))*exp(-(mint^2)/2*sigma.t^2)
# #S<-exp((-(1-abs(c))^2))*exp(-(mint^2)/(2*sigma.t^2))
# D<-1-S
# hist(D)
# 
# library(fastcluster)
# library(dynamicTreeCut)
# linkage = "average"
# clust <- fastcluster::hclust(as.dist(D), 
#                     method = linkage)
# hmax<-NULL
# deepSplit<-FALSE
# if (is.null(hmax)) {
#   hmax <- 0.3
# }
# clus <- dynamicTreeCut::cutreeDynamicTree(clust, 
#                                           maxTreeHeight = hmax, deepSplit = deepSplit, minModuleSize = 2)
# sing <- which(clus == 0)
# clus[sing] <- max(clus) + 1:length(sing)
# 
# 
# sapply(pks,function(x)nrow(x[[1]]))
# clus
# 
# cl.rt <- aggregate(unlist(rt), by = list(clus), FUN = "mean")
# 
# 
# new.ts<-rownames(dat.pr[[1]])
# plot(new.ts,dat.pr[[1]][,'280'],type='l')
# points(cl.rt$x, dat.pr[[1]][cl.rt$x*100,'280'],col='red')
# 
# mint <- abs(outer(unlist(rt),unlist(rt), FUN="-"))
# 
# # ff<-function(i,j) mean(diag(A[i,]) %*% ( (diag(B[i,])-diag(B[j,])) %*% (diag(C[i,])-diag(C[j,])) %*% diag(A[j,])))
# # distMatrix<-outer(1:3,1:3,Vectorize(ff))
# 
# 
# mint <- outer(rt,rt, FUN = "-")
# # fn<-function(x,y){x-y}
# # proxy::dist(rt, method = fn)
