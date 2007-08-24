
data(x100)
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
## out$blocks contains 3 data frames

## To illustrate two-level blocking, with multiple level two units per
##  level one unit:
for(i in (1:nrow(x100))){if(even(i)){x100$id[i] <- x100$id[i-1]}}
rm(i)

out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id", "id2"),
              block.vars = c("b1", "b2"), algorithm="optGreedy",
              distance = "mahalanobis", level.two = TRUE, valid.var =
              "b1", valid.range = c(0,500), verbose = TRUE) 
