# make sure the working directory has the PF partition specification and the alignment.
# The partitions are also written to this directory.
## Sebastian Duchene
## Feb 27 2015
# To run. source this function and type:
partition_data_blocks('partition_finder.cfg', 'alignment.phy')

partition_data_blocks <- function(pf_config, phy_alignment){
    require(ape)
    read.phylip <- function(file_name){
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

    partition_lines <- readLines(pf_config)
    sequence_matrix<- read.phylip(phy_alignment)

    pf_arguments <- partition_lines[grep('=', partition_lines)]
    partitions_args<- pf_arguments[-grep('alignment|models|model_selection|branchlengths|search', pf_arguments)]

    parts_output <- list()
    for(i in 1:length(partitions_args)){
        part_split <- strsplit(partitions_args[i], '=')[[1]]
        part_name <- gsub(' ', '', part_split[1])
        part_range <- strsplit(part_split[2], '[-]')[[1]]
        part_start <- as.numeric(part_range[1])
        if(length(grep('.*\\\\.*', part_range[2])) > 0){
            codon_temp <- strsplit(part_range[2], '\\\\')
            part_end <- as.numeric(codon_temp[[1]][1])
            part_by <- as.numeric(gsub(';', '', codon_temp[[1]][2]))
            part_seq <- seq(from = part_start, to = part_end, by = part_by)
        }else{
            part_end <- as.numeric(gsub(';', '', part_range[2]))
            part_seq <- seq(from = part_start, to = part_end)
        }
        parts_output[[i]] <- part_seq
        names(parts_output)[i] <- part_name
    }

    for(i in 1:length(parts_output)){
        part_temp <- sequence_matrix[, parts_output[[i]]]
        print(paste(names(parts_output)[i], 'saved in', getwd(), sep = ' '))
        write.dna(part_temp, file = paste0(names(parts_output)[i], '.fasta'), format = 'fasta', colsep = '', nbcol = -1)
    }
}
