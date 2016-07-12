library(annotatr)
library(readr)
library(ggplot2)
library(regioneR)
library(optparse)

################################################################################
# Deal with command line options

option_list = list(
  make_option('--file', type='character', help='[Required] Tab-delimited file with genomic locations and possibly associated data.'),
  make_option('--genome', type='character', help='[Required] The shortname for the genome used, e.g. hg19, mm9, rn4.'),
  make_option('--annot_type', type='character', help='[Required] One of bismark, simple_bisulfite, macs2, simple_pulldown, sample_class, methylSig, PePr, or compare_class. Indicates what type of data is being annotated.'),
  make_option('--group1', type='character', help='[Required] A character indicating the name of group1 or NULL.'),
  make_option('--group0', type='character', help='[Required] A character indicating the name of group2 or NULL.')
)

file = opt$file
genome = opt$genome
annot_type = opt$annot_type

# If the group names are not NULL, use them
if(!is.null(opt$group1)) {
	group1 = opt$group1
}
if(!is.null(opt$group0)) {
	group0 = opt$group0
}

# Deal with the annot_type
if(annot_type == 'bismark') {
	# gunzip -c ./bisulfite/bismark/IDH2mut_2_mc_hmc_bisulfite_trimmed_bismark_bt2.bismark.cov.gz | head
	# chr21	9437517	9437517	86.241610738255	257	41
	# chr21	9437519	9437519	75.503355704698	225	73
	# chr21	9437531	9437531	58.9225589225589	175	122
	# chr21	9437533	9437533	86.5771812080537	258	40
	col_names = c('chr','start','perc_meth','coverage')
	col_types = 'ci-di-'
	stranded = FALSE
	prefix = gsub('_bt2.bismark.cov.gz','', basename(file))
} else if (annot_type == 'simple_bisulfite') {
	# head ./classifications/simple/IDH2mut_1_mc_hmc_bisulfite_simple_classification.bed
	# chr1	249239097	249239098	mc_hmc_high	1000	.	249239097	249239098	102,0,102
	# chr1	249239100	249239101	mc_hmc_high	1000	.	249239100	249239101	102,0,102
	# chr1	249239104	249239105	mc_hmc_high	1000	.	249239104	249239105	102,0,102
	# chr1	249239117	249239118	mc_hmc_high	1000	.	249239117	249239118	102,0,102
	col_names = c('chr','start','end','class')
	col_types = 'ciic-----'
	stranded = FALSE
	prefix =
} else if (annot_type == 'macs2') {
	# head ./pulldown/macs2_peaks/IDH2mut_2_hmc_pulldown_macs2_peaks.narrowPeak
	# chr21	9909565	9910066	IDH2mut_2_hmc_pulldown_macs2_peak_1	126	.	3.04226	16.46495	12.64028	213
	# chr21	15471800	15472211	IDH2mut_2_hmc_pulldown_macs2_peak_2	146	.	7.33008	18.53385	14.66105	186
	# chr21	15855465	15855886	IDH2mut_2_hmc_pulldown_macs2_peak_3	77	.	4.05129	11.47467	7.78926	178
	# chr21	15865175	15865647	IDH2mut_2_hmc_pulldown_macs2_peak_4	322	.	6.33401	36.40632	32.21310	240
	col_names = c('chr','start','end','name','fold','pval')
	col_types = 'ciic--dd--'
	stranded = FALSE
	prefix = gsub('_peaks.narrowPeak','', basename(file))
} else if (annot_type == 'simple_pulldown') {
	# head ./classifications/simple/IDH2mut_1_hmc_pulldown_simple_classification.bed
	# chr21	9944160	9944341	hmc_low	1000	.	9944160	9944341	102,102,255
	# chr21	9999539	9999701	hmc_low	1000	.	9999539	9999701	102,102,255
	# chr21	10132638	10133209	hmc_low	1000	.	10132638	10133209	102,102,255
	# chr21	10135234	10135504	hmc_low	1000	.	10135234	10135504	102,102,255
	col_names = c('chr','start','end','class')
	col_types = 'ciic----'
	stranded = FALSE
	prefix =
} else if (annot_type == 'sample_class') {
	# head ./classifications/sample/NBM_1_sample_classification.bed
	# chr1	10468	10472	unclassifiable	1000	.	10468	10472	192,192,192
	# chr1	10483	10485	unclassifiable	1000	.	10483	10485	192,192,192
	# chr1	10488	10490	unclassifiable	1000	.	10488	10490	192,192,192
	# chr1	10492	10494	unclassifiable	1000	.	10492	10494	192,192,192
	col_names = c('chr','start','end','class')
	col_types = 'ciic----'
	stranded = FALSE
	prefix = gsub('.bed','', basename(file))
} else if (annot_type == 'methylSig') {
	# head ./bisulfite/methylsig_calls/IDH2mut_v_NBM_mc_hmc_bisulfite_DMR_methylSig.txt
	# chr	start	end	strand	pvalue	qvalue	meth.diff	logLikRatio	theta	df	mu1	mu0
	# chr21	9437472	9437521	*	0.000228579898882106	0.0519394743163992	28.6754246208103	158.6873713159	1e+06	4	81.4241341438466	52.7487095230364
	# chr21	9437522	9437571	*	0.0208553163921381	0.283830375900465	22.738738618692	13.6819351620397	49.7403416134298	4	85.8751984050635	63.1364597863715
	# chr21	9438322	9438371	*	0.777292550560113	0.922829646963547	-2.62511686393549	0.0915541817157646	23.6729453104579	4	75.0320958189817	77.6572126829172
	col_names = c('chr','start','end','pval','meth_diff','mu1','mu0')
	col_types = 'cii-d-d---dd'
	stranded = FALSE
	prefix = gsub('.txt','', basename(file))
} else if (annot_type == 'PePr') {
	# head ./pulldown/pepr_peaks/INTpos_v_INTneg_hmc_pulldown__PePr_chip1_peaks.bed
	# chr1	70760670	70761490	INTpos	245.943822054	.	6.5588413905	7.4943442631e-16	1.41424309192e-09
	# chr14	68731580	68732810	INTneg	333.178100172	.	8.8852092139	2.11278442242e-15	2.6579961895e-09
	# chr19	49923650	49924470	INTpos	167.462258567	.	4.4658913717	8.00430994576e-15	7.5523886073e-09
	# chr10	32465030	32466260	INTneg	117.736116164	.	3.13979227207	1.16899565361e-14	8.82395561008e-09
	col_names = c('chr','start','end','name','fold','pval')
	col_types = 'ciic--dd-'
	stranded = FALSE
	prefix =
} else if (annot_type == 'compare_class') {
	# head ./classifications/comparison/IDH2mut_v_NBM_compare_classification.bed
	# chr1	0	10000	unclassifiable	1000	.	0	10000	192,192,192
	# chr1	10000	10162	no_DM	1000	.	10000	10162	0,0,0
	# chr1	10162	10185	unclassifiable	1000	.	10162	10185	192,192,192
	# chr1	10185	10276	no_DM	1000	.	10185	10276	0,0,0
	col_names = c('chr','start','end','class')
	col_types = 'ciic----'
	stranded = FALSE
	prefix = gsub('.bed','', basename(file))
} else {
	stop('annot_type is invalid. Must be one of bismark, simple_bisulfite, macs2, simple_pulldown, sample_class, methylSig, PePr, or compare_class.')
}
