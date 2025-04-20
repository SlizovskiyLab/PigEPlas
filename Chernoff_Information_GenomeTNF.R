install.packages(c("Biostrings", "stringr", "data.table"
                   , "pheatmap"))
library(Biostrings)
library(stringr)
library(pheatmap)

setwd("/Users/islizovs/Downloads")
#GCA_018041865.1_PDT000863098.1_genomic.fna >>ETEC Genome 1
#GCA_019078505.1_PDT001079507.1_genomic.fna >>ETEC Genome 2
genome <- readDNAStringSet("GCA_018041865.1_PDT000863098.1_genomic.fna")  # Replace with your FASTA file
seq <- as.character(genome[[1]])  # Assuming single contig

window_size <- 10000
step <- 1000

if (nchar(seq) >= window_size) {
  starts <- seq(1, nchar(seq) - window_size + 1, by = step)
  windows <- substring(seq, starts, starts + window_size - 1)
} else {
  stop("Genome is shorter than the window size. Try a smaller window or a longer genome.")
}

library(stringr)
get_tnf <- function(s) {
  s <- toupper(s)
  trinucs <- substring(s, 1:(nchar(s)-2), 3:(nchar(s)))
  trinucs <- trinucs[!grepl("N", trinucs)]
  freq <- table(factor(trinucs, levels = apply(expand.grid(c("A","C","G","T"), 
                                                           c("A","C","G","T"), 
                                                           c("A","C","G","T")), 1, paste, collapse="")))
  return(as.numeric(freq) / sum(freq))
}

tnf_list <- lapply(windows, get_tnf)
tnf_list

chernoff_info <- function(p, q, lambdas = seq(0.01, 0.99, by = 0.05)) {
  max_info <- -Inf
  for (lambda in lambdas) {
    mixed <- sum(p^lambda * q^(1 - lambda))
    if (mixed > 0) {
      info <- -log(mixed)
      if (info > max_info) max_info <- info
    }
  }
  return(max_info)
}


n <- length(tnf_list)
mat <- matrix(0, nrow = n, ncol = n)

for (i in 1:n) {
  for (j in i:n) {
    val <- chernoff_info(tnf_list[[i]], tnf_list[[j]])
    mat[i, j] <- val
    mat[j, i] <- val  # symmetry
  }
}
range(mat)
library(viridis)
# Clamp small negatives to 0
mat[mat < 0] <- 0

window_step <- 1000  # in bp
genome_positions_kb <- seq(0, by = window_step, length.out = nrow(mat)) / 1000  # in kb

tick_interval <- 50  # in kb
tick_indices <- which(genome_positions_kb %% tick_interval == 0)
tick_labels <- paste0(genome_positions_kb[tick_indices], " kb")

row_labels <- rep("", nrow(mat))
col_labels <- rep("", ncol(mat))

row_labels[tick_indices] <- tick_labels
col_labels[tick_indices] <- tick_labels



pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         labels_row = row_labels,
         labels_col = col_labels,
         color = viridis(100),
         breaks = seq(0, 0.01, length.out = 101),
         fontsize = 8,
         # Add manual tickmarks
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE)
