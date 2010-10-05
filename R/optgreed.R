optgreed <- function(dist, n.tr, l2, l1names, valid, validvar, validlb, validub, verbose){
p = (length(unique(l1names)) %/% n.tr) + as.integer((length(unique(l1names)) %% n.tr) > 0)
out <- .C("optgreed",
           vec = as.double(dist[row(dist) < col(dist)]),
           n = as.integer(choose(nrow(dist), 2)),
           nrow = as.integer(nrow(dist)),
           ntr = as.integer(n.tr),
           l2 = as.integer(l2),
           l1names = as.integer(as.factor(l1names)),
           valid = as.integer(valid),
           validvar = as.double(validvar),
           validlb = as.double(validlb),
           validub = as.double(validub),
           verbose = as.integer(verbose),
           pairdist = numeric(p),
           result = integer(p * n.tr))
result <- data.frame(matrix(out$result, ncol=(n.tr), byrow=TRUE), out$pairdist)
return(result)
}