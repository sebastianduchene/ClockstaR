
read.phylip <- function(file_name){
    require(ape)
    raw_lines <- readLines(file_name)
    raw_seqs <- raw_lines[-grep('^ *[0-9]', raw_lines)]
    list_seqs <- list()
    for(i in 1:length(raw_seqs)){
        line_split <- strsplit(raw_seqs[i], '\t+| +')
        list_seqs[[i]] <- line_split[[1]][2]
        names(list_seqs)[i] <- line_split[[1]][1]
    }
    seq_lens <- sapply(1:length(list_seqs), function(x) nchar(list_seqs[[x]]))

    if(length(unique(seq_lens)) > 1) stop('sequences are not the same length')

    matrix_output <- matrix(NA, length(seq_lens), seq_lens[1])
    for(i in 1:length(list_seqs)){
        matrix_output[i, ] <- strsplit(tolower(list_seqs[[i]]), '')[[1]]
    }
    rownames(matrix_output) <- names(list_seqs)
    matrix_output <- as.DNAbin(matrix_output)
    return(matrix_output)
}
