
optimised_trees <- read.tree('clockstar_example_data/optimised_tr.trees') # optimised_tr.trees are the gene trees with the same topology but different branch lengths
bsd_mat <- bsd.matrix(optimised_trees)
part_data <- partitions(bsd_mat)
write.table(part_data$parts.mat, file = 'partitions_output.csv', sep = ',') # partitions_output.csv is the output file with the partitions for different values of k. This can be opened in excel or other spreadsheet software.

