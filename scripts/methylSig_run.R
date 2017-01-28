library(optparse)
library(methylSig)

option_list = list(
	make_option('--project', type='character'),
	make_option('--cytfiles', type='character'),
	make_option('--sampleids', type='character'),
	make_option('--assembly', type='character'),
	make_option('--pipeline', type='character'),
	make_option('--context', type='character'),
	make_option('--resolution', type='character'),
	make_option('--treatment', type='character'),
	make_option('--destranded', type='logical'),
	make_option('--maxcount', type='integer'),
	make_option('--mincount', type='integer'),
	make_option('--filterSNPs', type='logical'),
	make_option('--ncores', type='integer', default=1),
	make_option('--quiet', type='logical', default=FALSE),
	make_option('--dmtype', type='character'),
	make_option('--winsize.tile', type='integer'),
	make_option('--dispersion', type='character'),
	make_option('--local.disp', type='logical'),
	make_option('--winsize.disp', type='integer'),
	make_option('--local.meth', type='logical'),
	make_option('--winsize.meth', type='integer'),
	make_option('--minpergroup', type='character'),
	make_option('--T.approx', type='logical'),
	make_option('--outprefix', type='character')
)

opt = parse_args(OptionParser(option_list=option_list))

# Setup data directory string
project = opt$project
prefix = opt$outprefix

# Parse comma-separated arguments
cyt_files = unlist(strsplit(opt$cytfiles, ','))
sample_id = unlist(strsplit(opt$sampleids, ','))
treatment = as.integer(unlist(strsplit(opt$treatment, ',')))
min_per_group = as.integer(unlist(strsplit(opt$minpergroup, ',')))

if(any(!file.exists(cyt_files))){
  stop(sprintf('Some input files do not exist: %s', str(file.exists(cyt_files))))
}

meth = methylSigReadData(
	fileList = cyt_files,
	sample.ids = sample_id,
	assembly = opt$assembly,
	pipeline = opt$pipeline,
	header= FALSE,
	context = opt$context,
	resolution = opt$resolution,
	treatment = treatment,
	destranded = opt$destranded,
	maxCount = opt$maxcount,
	minCount = opt$mincount,
	filterSNPs = opt$filterSNPs,
	num.cores = opt$ncores,
	quiet= opt$quiet)

if(opt$dmtype == 'DMR') {
	message('Doing tiled analysis')
	meth_tiled = methylSigTile(meth, win.size = opt$winsize.tile)

	diff_meth = methylSigCalc(
		meth_tiled,
		group = c('Treatment' = max(treatment),'Control' = min(treatment)),
		dispersion = opt$dispersion,
		local.disp = opt$local.disp,
		winsize.disp = opt$winsize.disp,
		local.meth = opt$local.meth,
		winsize.meth = opt$winsize.meth,
		min.per.group = min_per_group,
		T.approx = opt$T.approx,
		num.cores= opt$ncores)

	write.methylSigDiff(diff_meth, file=sprintf('bisulfite/methylsig_calls/%s.txt', prefix), row.names=F,quote=F,sep='\t')
} else if (opt$dmtype == 'DMC' ){
	message('Doing CpG analysis')
	diff_meth = methylSigCalc(
		meth,
		group = c('Treatment' = max(treatment),'Control' = min(treatment)),
		dispersion = opt$dispersion,
		local.disp = opt$local.disp,
		winsize.disp = opt$winsize.disp,
		local.meth = opt$local.meth,
		winsize.meth = opt$winsize.meth,
		min.per.group = min_per_group,
		T.approx = opt$T.approx,
		num.cores= opt$ncores)

	write.methylSigDiff(diff_meth, file=sprintf('bisulfite/methylsig_calls/%s.txt', prefix),
		row.names=F, col.names=F, quote=F, sep='\t')

} else {
	stop('Error in methylSig run. Invalid OPT_DM_TYPE in config.mk. Must be DMC for CpG resolution or DMR for regions of winsize.tile resolution.')
}

save.image(file=sprintf('RData/%s_analysis.RData', prefix))
