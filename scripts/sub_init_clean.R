################################################################################
# MAKEFILE: clean rules

clean_str = ''
if(bool_bis_samp) {
	clean_str = paste(clean_str, 'clean_bisulfite_align_tmp')
}

if(bool_pull_samp) {
	clean_str = paste(clean_str, 'clean_pulldown_sample_tmp')
}

if(bool_bis_samp || bool_pull_samp) {
	clean_str = paste(clean_str, 'clean_sample_classification_tmp')
}

if(bool_bis_comp) {
	clean_str = paste(clean_str, 'clean_bisulfite_compare_tmp')
}

if(bool_pull_comp) {
	clean_str = paste(clean_str, 'clean_pulldown_compare_tmp')
}

if(bool_bis_comp || bool_pull_comp) {
	clean_str = paste(clean_str, 'clean_compare_classification_tmp')
}

make_clean_tmp_rule = sprintf('.PHONY : clean_tmp
clean_tmp : %s', clean_str)

cat(make_clean_tmp_rule, file = file_make, sep = '\n', append = TRUE)
