
partition_data_partitionfinder <- function(data_file_fasta, partition_finder_file){
   partition_finder_file <- readLines(partition_finder_file)
   data_file <- read.dna(data_file_fasta, format = 'fasta')

  if(!is.matrix(data_file)) stop('The data have not been read as a matrix. Please check that it is aligned')

  raxml_format_parts <- partition_finder_file[(grep('RaxML[-]style', partition_finder_file) + 1):length(partition_finder_file)]
  raxml_format_parts <- raxml_format_parts[1:(which(raxml_format_parts == '')[1] - 1)]

  partitions_raw <- gsub('^[A-Z]+,', '', raxml_format_parts)
  parts_list  <- list()
  for(i in 1:length(partitions_raw)){
    parts_dat_temp <- gsub(' ', '', strsplit(partitions_raw[i], '=')[[1]])  
    parts_name_temp <- parts_dat_temp[1]
    parts_ranges_temp <- strsplit(parts_dat_temp[2:length(parts_dat_temp)], ',')[[1]]
    parts_ranges_proc <- list()
    parts_list[[i]] <- list()
    for(p in 1:length(parts_ranges_temp)){
      range_raw <- strsplit(parts_ranges_temp[p], '-')[[1]]

      if(grepl('[\\]', parts_ranges_temp[p])){
        end_p_temp <- substr(parts_ranges_temp[p], gregexpr('-', parts_ranges_temp[p])[[1]] + 1, nchar(parts_ranges_temp[p]))
        split_end <- strsplit(end_p_temp, '[\\]')[[1]]
        end_p <- split_end[1]
        by_p <- split_end[2]
        seq_range <- seq(from = as.numeric(range_raw[1]), to = as.numeric(end_p), by = as.numeric(by_p))
      }else{
        end_p <- range_raw[2]
        seq_range <- seq(from = as.numeric(range_raw[1]), to = as.numeric(end_p)) 
      }

      parts_list[[i]] <- c(parts_list[[i]], seq_range)
    }
    names(parts_list)[i] <- parts_name_temp
    parts_list[[i]] <- unlist(parts_list[[i]])
  }

  for(pdat in 1:length(parts_list)){
    cat('saving partition', names(parts_list)[pdat], 'to', paste0(names(parts_list)[pdat], '.fasta\n'))
    write.dna(data_file[, parts_list[[pdat]]], file = paste0(names(parts_list)[pdat], '.fasta'), format = 'fasta', nbcol = -1, colsep = '')
  }

}