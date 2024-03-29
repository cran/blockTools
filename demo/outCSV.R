
data(x100)

# First, block
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
# Second, assign
assg <- assignment(out, seed = 123)

# create three .csv files of blocks
# outCSV(out)

# create three .csv files of assigned blocks
#   (note: overwrites blocked .csv files)
# outCSV(assg)

# create three .csv files with custom file names
# outCSV(assg, file.names = list("file1", "file2", "file3"))
