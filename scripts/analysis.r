library(Biostrings)
library(tidyverse)
library(naturalsort)

## analysis consensus sequences to check BA1 and BA2 specific mutations
### specific mutations of VOC588 and VOC589
source("./helper/helper_comp_seq.r")
files_seqs <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus/", "fa$", full.names=T)
check <- grepl("VOC0195-P3-25", files_seqs) | grepl("VOC588-P2_", files_seqs) | grepl("/HL-H", files_seqs) | grepl("/HL-NC", files_seqs)
seq_all_raw <- lapply(files_seqs[check], readDNAStringSet)
seq_all_raw <- do.call(c, seq_all_raw)
names(seq_all_raw) <- sapply(strsplit(names(seq_all_raw), "_"), function(x) {x[2]})
writeXStringSet(seq_all_raw, "../results/seq_all_raw.fasta")

system("mafft --auto --thread -1 --addfragments ../results/seq_all_raw.fasta ../../../2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta > ../results/seq_all_raw_aln.fasta")
seq_all <- readDNAStringSet("../results/seq_all_raw_aln.fasta") 
# seq_all <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus/con_all_combined_aln.fasta")

check <- grepl("HL-H", names(seq_all)) | grepl("HL-NC", names(seq_all))
sample_names_exp <- names(seq_all)[check]
seq_exp <- seq_all[check]
sample_names_exp <- sapply(strsplit(sample_names_exp, "-"), function(x) {x[grepl("^N", x)]})
sample_names_num <- gsub("[NC]","",sample_names_exp)

seqs_filled <- mclapply(c(1:3, 4:6, 10:12, 22:24), function(i) {
	print(i)
	seq_i <- seq_exp[sample_names_num==i]
	pos_tofill <- 21290:24270
	check_i <<- grepl("N\\d", names(seq_i))
	seq_filled <- seq_i[check_i]
	seq_segment <- seq_i[!check_i]
	sapply(pos_tofill, function(x){
		subseq(seq_filled, x, x) <<- subseq(seq_segment, x, x)
		return(NA)
	})
	return(seq_filled)
}, mc.cores=8)
seqs_filled <- do.call(c, seqs_filled)
names(seqs_filled) <- paste0(names(seqs_filled), "_filled")

check_ins_pos <- which(strsplit(as.character(seq_all[1]), "")[[1]] == "-")
seq_all_trim <- c(seq_all[1], seq_all[grepl("VOC", names(seq_all))], seqs_filled)
sapply(check_ins_pos-(0:((length(check_ins_pos)-1))), function(pos) {
	subseq(seq_all_trim, pos, pos) <<- ""
})
stopifnot(!any(strsplit(as.character(seq_all_trim[1]), "")[[1]] == "-"))

seq_ref <- seq_all_trim[grepl("MN90", names(seq_all_trim))] # Ref
seq_voc195 <- seq_all_trim[grepl("VOC0195-P3-25", names(seq_all_trim))] # BA1
seq_voc588 <- seq_all_trim[grepl("VOC588-P2", names(seq_all_trim))][1] # BA2

seq_h16 <- seq_all_trim[grepl("HL-H", names(seq_all_trim)) | grepl("HL-NC", names(seq_all_trim))]
seq_h16 <- seq_h16[naturalorder(names(seq_h16))]

comp_seqs(seq_voc588, c(seq_h16, seq_voc195, seq_ref), outfile_prefix="../results/voc588_comp")
comp_seqs(seq_voc195, c(seq_h16, seq_voc588, seq_ref), outfile_prefix="../results/voc195_comp")
# comp_seqs(seq_voc195, seq_voc588, outfile_prefix="../results/ref_comp")

## analysis VCFs
samples <- c(gsub("_NA", "", c(names(seq_voc195), names(seq_voc588))), names(seq_h16))
samples <- sort(samples)

# files_tsv <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/variant_caller/ivar/")
# files_tsv_full <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/variant_caller/ivar/", full.names = T)

# files_int <- files_tsv_full[files_tsv %in% paste0("ivar_", samples, ".tsv")]
source("./helper/convert_to_vcf.r")
# convert_to_vcf(files_int)
# files_vcf <- list.files("../results/vcf/", full.names = T)
# annotate_snpeff(files_vcf)
# files_snpeff_csv <- list.files("../results/vcf/", "csv$", full.names = T)
# df_vcf <- read_snpeff(files_snpeff_csv)


file_bam <- list.files("../data/", "bam$", full.names = T)
# file_bam <- file_bam[grepl("NC", file_bam)]
mclapply(file_bam, function(x){
	system(paste0("./helper/call_vcf_from_bam.sh ", x))
}, mc.cores=8)
files_vcf <- list.files("../data/", "\\.normvcf$", full.names = T)
# files_vcf <- files_vcf[grepl("NC", files_vcf)]
annotate_snpeff(files_vcf, mc.cores=8)
files_snpeff_csv <- list.files("../data", "normvcf\\.csv$", full.names = T)
df_vcf_ori <- read_snpeff(files_snpeff_csv)

files_depth <- list.files("../data", "cov$", full.names = T)
files_depth_sim <- list.files("../data", "cov$")
sample_name <- gsub("-trimmed.masked.bam.cov", "", files_depth_sim)
df_depth <- lapply(seq_along(sample_name), function(i) {
	tmp <- read_tsv(files_depth[i], col_names=F)
	tmp$sample <- sample_name[i]
	return(tmp)
})
df_depth <- bind_rows(df_depth)
names(df_depth)[2:3] <- c("POS", "depth")

### plot composition of interested sites 
df_vcf <- df_vcf_ori
df_vcf$mutation <- paste0(df_vcf$gene, ": ", gsub("p\\.", "", df_vcf$mut_aa))
df_vcf %>% filter(POS==25000) %>% .$mutation
df_vcf <- df_vcf %>% filter(!is.na(mut_aa))
# df_vcf <- df_vcf %>% filter(!grepl("del", mut_aa))
# df_vcf <- df_vcf %>% filter(!grepl("ins", mut_aa))
# df_vcf$mutation[is.na(df_vcf$mut_aa)] <- paste0("Nt_", df_vcf$POS, ": ", df_vcf$mut_cdna
# )[is.na(df_vcf$mut_aa)]
# df_vcf$mutation[is.na(df_vcf$mut_aa)] <- gsub("del\\D+", "del", df_vcf$mutation[is.na(df_vcf$mut_aa)])

# df_vcf$gene[is.na(df_vcf$mut_aa)] <- NA
df_vcf$sample <- gsub("-trimmed", "", df_vcf$sample)
unique(df_vcf$sample)
df_vcf_voc195 <- df_vcf %>% filter(sample == "VOC0195-P3-25-S3-iseq")
df_vcf_voc588 <- df_vcf %>% filter(sample == "VOC588-P2")
df_vcf_hl <- df_vcf %>% filter(!sample %in% c("VOC0195-P3-25-S3-iseq", "VOC588-P2"))

mut_BA1BA2 <- c(unique(df_vcf_voc195$mutation), unique(df_vcf_voc588$mutation))
mut_dup <- mut_BA1BA2[duplicated(mut_BA1BA2)]

muts_int <- table(df_vcf$mutation)[table(df_vcf$mutation)>1]
muts_int <- names(muts_int)
muts_int <- muts_int[(muts_int %in% mut_BA1BA2)|grepl("ins", muts_int)|grepl("del", muts_int)]
muts_int <- muts_int[!muts_int %in% mut_dup]

df_vcf_plot <- df_vcf %>% filter(mutation %in% muts_int)
df_vcf_plot$mutation_fct <- reorder(df_vcf_plot$mutation, df_vcf_plot$POS)
levels(df_vcf_plot$mutation_fct) <- gsub(".+: ", "", levels(df_vcf_plot$mutation_fct))
df_vcf_plot$gene_fct <- reorder(df_vcf_plot$gene, df_vcf_plot$POS)
df_vcf_plot$gene_fct <- fct_rev(df_vcf_plot$gene_fct)
df_vcf_plot$sample_fct <- naturalfactor(df_vcf_plot$sample)
levels(df_vcf_plot$sample_fct) <- sapply(levels(df_vcf_plot$sample_fct), function(x){
	tmp <- strsplit(x, "-")[[1]]
	tmp[grepl("VOC", tmp) | grepl("^N+", tmp)]
})
df_vcf_plot$group <- sapply(as.character(df_vcf_plot$sample_fct), function(x){
	if(grepl("VOC", x)){return("VOC")}
	x <- gsub("[NC]", "", x)
	if(x %in% 1:3){return("BA.1 infected")}
	if(x %in% 4:6){return("BA.2 infected")}
	if(x %in% 10:12){return("Direct contact 1")}
	if(x %in% 22:24){return("Direct contact 2")}
})

df_vcf_plot$DP <- sapply(df_vcf_plot$mut, function(x){
	tmp <- strsplit(x, ";")[[1]]
	depth <- tmp[grepl("DP=", tmp)]
	depth <- as.numeric(gsub("DP=","",depth))
}, USE.NAMES=F)
df_vcf_plot$ALT_FREQ <- sapply(df_vcf_plot$mut, function(x){
	tmp <- strsplit(x, ";")[[1]]
	depth <- tmp[grepl("DP=", tmp)]
	depth <- as.numeric(gsub("DP=","",depth))
	AD <- tmp[grepl("AD=", tmp)]
	AD <- as.numeric(strsplit(AD, ",")[[1]][2])
	AD/depth
}, USE.NAMES=F)
df_vcf_plot <- df_vcf_plot %>% filter(ALT_FREQ>0.1)
# df_vcf_plot$ALT_FREQ[df_vcf_plot$DP<=10]

df_vcf_plot_gap <- df_vcf_plot %>% select(POS,mutation_fct,gene_fct) %>% unique()
df_vcf_plot_gap <- lapply(unique(df_vcf_plot$sample), function(x) {
	df_tmp <- df_vcf_plot_gap
	df_tmp$sample <- x
	df_tmp
})
df_vcf_plot_gap <- bind_rows(df_vcf_plot_gap)
df_vcf_plot_gap <- left_join(df_vcf_plot_gap, df_depth, c("sample", "POS"))
df_vcf_plot_gap <- df_vcf_plot_gap %>% filter(depth<=10)
df_vcf_plot_gap$sample_fct <- naturalfactor(df_vcf_plot_gap$sample)
levels(df_vcf_plot_gap$sample_fct) <- sapply(levels(df_vcf_plot_gap$sample_fct), function(x){
	tmp <- strsplit(x, "-")[[1]]
	tmp[grepl("VOC", tmp) | grepl("^N+", tmp)]
})
df_vcf_plot_gap$group <- sapply(as.character(df_vcf_plot_gap$sample_fct), function(x){
	if(grepl("VOC", x)){return("VOC")}
	x <- gsub("[NC]", "", x)
	if(x %in% 1:3){return("BA.1 infected")}
	if(x %in% 4:6){return("BA.2 infected")}
	if(x %in% 10:12){return("Direct contact 1")}
	if(x %in% 22:24){return("Direct contact 2")}
})

p_out <- ggplot(df_vcf_plot)+
	geom_tile(aes(x=sample_fct, y=mutation_fct, fill=ALT_FREQ))+
	geom_tile(aes(x=sample_fct, y=mutation_fct), fill="grey", data=df_vcf_plot_gap)+
	scale_fill_gradient2(name= "Allele frequency", low = "blue", mid="yellow",
	high = "dark red", midpoint = 0.5)+
	facet_grid(rows=vars(gene_fct), cols=vars(group), scales="free", space="free")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	theme_classic()+
	ylab("SNVs")+
	xlab("Sample")
ggsave("../results/alt_freq.pdf", height=12, width=15, plot=p_out)

df_vcf_plot$mutation_fct2 <- paste0(df_vcf_plot$POS, ": ", df_vcf_plot$mutation_fct)
df_vcf_plot$mutation_fct2 <- naturalfactor(df_vcf_plot$mutation_fct2)
df_vcf_plot_gap$mutation_fct2 <- paste0(df_vcf_plot_gap$POS, ": ", df_vcf_plot_gap$mutation_fct)
df_vcf_plot_gap$mutation_fct2 <- naturalfactor(df_vcf_plot_gap$mutation_fct2)
p_out2 <- ggplot(df_vcf_plot)+
	geom_tile(aes(x=sample_fct, y=mutation_fct2, fill=ALT_FREQ))+
	geom_tile(aes(x=sample_fct, y=mutation_fct2), fill="grey", data=df_vcf_plot_gap)+
	scale_fill_gradient2(name= "Allele frequency", low = "blue", mid="yellow",
	high = "dark red", midpoint = 0.5)+
	facet_grid(rows=vars(gene_fct), cols=vars(group), scales="free", space="free")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	ylab("SNVs")+
	xlab("Sample")+
	theme_classic()
ggsave("../results/alt_freq_nt.pdf", height=12, width=15, plot=p_out2)

### generate the table for comparision
df_vcf_ba2 <- df_vcf_plot %>% filter(grepl("VOC588", sample))
df_vcf_ba1 <- df_vcf_plot %>% filter(grepl("195", sample))

files_bamstat <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", full.names=T)
check <- sapply(paste0(samples, ".tsv"), function(x) {
	tmp <- grep(x, files_bamstat)
	if(length(tmp)==0){return(NA)}else{return(tmp)}
})
df_bamstat <- mclapply(files_bamstat[check[!is.na(check)]], read_tsv, mc.cores=8)
df_bamstat <- bind_rows(df_bamstat)
df_bamstat$sample <- rep(samples[!is.na(check)], each=29903)

source("./helper/af_by_bamstat.r")

df_vcf_int <- df_vcf_ba2
df_af <- lapply(naturalsort(samples), function(x) {
	tmp <- af_by_bamstat(df_bamstat, x, df_vcf_int$POS, df_vcf_int$REF, df_vcf_int$ALT)
	tmp$sample <- x
	tmp$POS <- df_vcf_int$POS
	tmp$REF <- df_vcf_int$REF
	tmp$ALT <- df_vcf_int$ALT
	return(tmp)
})
df_af <- bind_rows(df_af)
df_af$Alt_freq <- df_af$ALT_depth/df_af$Total_depth
df_af <- left_join(df_af, df_vcf_int %>% select(POS:ALT, mutation))
df_af$sample <- sapply(df_af$sample, function(x){
	tmp <- strsplit(x, "-")[[1]]
	tmp[grepl("VOC", tmp) | grepl("^N+", tmp)]
})
df_af_af <- df_af %>% filter(Total_depth>=100) %>% select(sample:ALT, mutation, Alt_freq) %>% pivot_wider(names_from=sample, values_from=Alt_freq)
df_af_ad <- df_af %>% select(sample:ALT, mutation, ALT_depth) %>% pivot_wider(names_from=sample, values_from=ALT_depth)
df_af_td <- df_af %>% select(sample:ALT, mutation, Total_depth) %>% pivot_wider(names_from=sample, values_from=Total_depth)
writexl::write_xlsx(df_af_ad, "../results/df_ad_ba2.xlsx")
writexl::write_xlsx(df_af_af, "../results/df_af_ba2_100.xlsx")
writexl::write_xlsx(df_af_td, "../results/df_td_ba2.xlsx")

df_vcf_int <- df_vcf_ba1
df_af <- lapply(naturalsort(samples), function(x) {
	tmp <- af_by_bamstat(df_bamstat, x, df_vcf_int$POS, df_vcf_int$REF, df_vcf_int$ALT)
	tmp$sample <- x
	tmp$POS <- df_vcf_int$POS
	tmp$REF <- df_vcf_int$REF
	tmp$ALT <- df_vcf_int$ALT
	return(tmp)
})
df_af <- bind_rows(df_af)
df_af$Alt_freq <- df_af$ALT_depth/df_af$Total_depth
df_af <- left_join(df_af, df_vcf_int %>% select(POS:ALT, mutation))
df_af$sample <- sapply(df_af$sample, function(x){
	tmp <- strsplit(x, "-")[[1]]
	tmp[grepl("VOC", tmp) | grepl("^N+", tmp)]
})
df_af_af <- df_af %>% filter(Total_depth>=100) %>% select(sample:ALT, mutation, Alt_freq) %>% pivot_wider(names_from=sample, values_from=Alt_freq)
df_af_ad <- df_af %>% select(sample:ALT, mutation, ALT_depth) %>% pivot_wider(names_from=sample, values_from=ALT_depth)
df_af_td <- df_af %>% select(sample:ALT, mutation, Total_depth) %>% pivot_wider(names_from=sample, values_from=Total_depth)
writexl::write_xlsx(df_af_af, "../results/df_af_ba1_100.xlsx")
writexl::write_xlsx(df_af_ad, "../results/df_ad_ba1.xlsx")
writexl::write_xlsx(df_af_td, "../results/df_td_ba1.xlsx")

## analysis whether there is any reassortment
