
af_by_bamstat <- function(df_bam, sample, POS, REF, ALT) {
	df_tmp <- df_bam[df_bam$sample==sample,]
	print(df_tmp$sample[1])
	stopifnot(length(POS)==length(ALT))
	stopifnot(length(POS)==length(REF))
	tmp <- lapply(seq_along(POS), function(i){
		pos_i <- POS[i]
		ref_i <- REF[i]
		alt_i <- ALT[i]
		if(nchar(ref_i)>nchar(alt_i)){ # deletion
			pos_i <- pos_i + 1
			mut <- "deletions"
		} else if(nchar(ref_i)<nchar(alt_i)){ # insertion
			mut <- "insertions"
		} else { # sub
			mut <- toupper(alt_i)
		}
		depth_i <- df_tmp$reads_all[pos_i]
		depth_alt_i <- df_tmp[[mut]][pos_i]
		return(c(depth_alt_i, depth_i))
	})
	tmp <- matrix(unlist(tmp), ncol=2, byrow=T)
	colnames(tmp) <- c("ALT_depth", "Total_depth")
	return(as_tibble(tmp))
}
