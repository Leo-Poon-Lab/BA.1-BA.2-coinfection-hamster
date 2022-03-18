comp_seqs <- function(seq_int, seqs_compare, ignore_nt=c("-", "N", "?"), gap_filter=NA, primer_file='../../../2021/2021-10-11_nCoV_primers/results/primers_20211011.bed', outfile_prefix=NA){
	stopifnot(width(seq_int)==29903)
	pos_ignore <- which(strsplit(as.character(seq_int), "")[[1]] %in% ignore_nt)
	char_seq_int <- strsplit(as.character(seq_int), "")[[1]]
	
	suppressMessages(df_primer <- readr::read_tsv(primer_file, col_names = F))
	df_primer$X2 <- df_primer$X2+1 # 0 based

	# initialize results
	len_rst <- length(seqs_compare)
	quantile_len_rst <- round(quantile(seq_len(len_rst)))

	num_gap <- rep(NA, len_rst)
	mut <- rep(NA, len_rst)
	num_diff <- rep(NA, len_rst)
	num_diff_spike <- rep(NA, len_rst)
	primer_binding_check <- rep(NA, len_rst)
	muts_excluding_primer_binding <- rep(NA, len_rst)
	num_diff_excluding_primer_binding <- rep(NA, len_rst)
	num_diff_spike_excluding_primer_binding <- rep(NA, len_rst)
	# cal results
	invisible(sapply(seq_along(seqs_compare), function(i){
		if(i %in% quantile_len_rst){print(names(quantile_len_rst)[quantile_len_rst==i])}
		seq_t <- seqs_compare[i]
		check <- strsplit(compareStrings(seq_int, seq_t), "")[[1]] %in% c("?", "-", "N")
		idx <- which(check)
		idx <- idx[!idx %in% pos_ignore]
		ref_base <- char_seq_int[idx]
		alt_base <- strsplit(as.character(seq_t), "")[[1]][idx]
		check2 <- !alt_base %in% ignore_nt
		
		pos_diff <- idx[check2]
		num_gap[i] <<- sum(!check2)
		mut_i <<- paste0(ref_base[check2], pos_diff, alt_base[check2])
		mut[i] <<- paste(mut_i, collapse="|")
		num_diff[i] <<- length(pos_diff)
		num_diff_spike[i] <<- sum(pos_diff <= 25384 & pos_diff >= 21563)

		primer_binding_check_t <- sapply(pos_diff, function(y){
			tmp_y <- df_primer$X4[df_primer$X2<=y & df_primer$X3>=y]
			if(length(tmp_y)==0){return(NA)}else{return(paste(tmp_y,collapse=","))}
		})
		check_primer_bind <- !is.na(primer_binding_check_t)
		if(sum(check_primer_bind)>0){
			primer_binding_check[i] <<- paste0(paste0(mut_i[check_primer_bind], ":", primer_binding_check_t[check_primer_bind]), collapse="|")
		} else{
			primer_binding_check[i] <<- NA
		}

		muts_excluding_primer_binding[i] <<- mapply(function(x,y){
			paste0(x[is.na(y)], collapse = "|")
		},mut_i, primer_binding_check_t) %>% .[.!=""] %>% paste0(collapse="|")

		num_diff_excluding_primer_binding[i] <<- sum(is.na(primer_binding_check_t))
		num_diff_spike_excluding_primer_binding[i] <<- sum(pos_diff[is.na(primer_binding_check_t)] <= 25384 & pos_diff[is.na(primer_binding_check_t)] >= 21563)
		
		return(NA)
	}))

	df_diff <- tibble(sample=names(seqs_compare), num_diff=num_diff, num_diff_excluding_primer_binding=num_diff_excluding_primer_binding, num_diff_spike = num_diff_spike,num_diff_spike_excluding_primer_binding=num_diff_spike_excluding_primer_binding, num_gap=num_gap, mutations = mut, primer_binding_site=primer_binding_check, mutations_excluding_primer_binding=muts_excluding_primer_binding) 
	
	if(!is.na(gap_filter)){
		df_diff <- df_diff %>% filter(num_gap<=gap_filter)
	}
	if(any(is.na(df_diff$num_diff_excluding_primer_binding))){
		df_diff$num_diff_excluding_primer_binding[is.na(df_diff$num_diff_excluding_primer_binding)] <- 0
	}
	# df_diff <- df_diff %>% arrange(num_diff_excluding_primer_binding, num_diff)
	
	if(!is.na(outfile_prefix)){
		write_csv(df_diff, paste0(outfile_prefix, ".csv"))
		writeXStringSet(c(seq_int, seqs_compare), paste0(outfile_prefix, ".fasta"))
	}
	
	return(df_diff)
}
