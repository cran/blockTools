mahalnaive <- function(x, block.vars, vcov, n.tr, l2, l1names, valid, validvar, validlb, validub, verbose){
vcd <- as.matrix(x[,block.vars])
as.double(vcd)
vcovi <- solve(vcov)
p = (length(unique(l1names)) %/% n.tr) + as.integer((length(unique(l1names)) %% n.tr) > 0)
out <- .C("mahalnaive",
          data = as.double(vcd),
          nrow = as.integer(nrow(x)),
          ncol = as.integer(ncol(vcd)),
          vcovi = as.double(vcovi),
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
