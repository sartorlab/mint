make_clean_tmp_rule = '.PHONY : clean_tmp
clean_tmp :
	clean_bisulfite_align_tmp clean_bisulfite_compare_tmp clean_pulldown_sample_tmp clean_pulldown_compare_tmp clean_sample_classification_tmp clean_compare_classification_tmp'

cat(make_clean_tmp_rule, file = file_make, sep = '\n', append = TRUE)
