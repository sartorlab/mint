library(optparse)
library(methylSig)

# Example call
#    Rscript ~/latte/Methylation/Methylation_Code/process_errbs_comparison-wise_run_methylSig.R  --cytfiles   --sampleids   --assembly   --pipeline  --context  --resolution  --treatment  --destranded  --maxcount  --mincount  --filterSNPs  --ncores  --quiet  --tile  --dispersion  --minpergroup

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
    make_option('--tile', type='logical'),
    make_option('--dispersion', type='character'),
    make_option('--minpergroup', type='character'),
    make_option('--comparison', type='character')
)

opt = parse_args(OptionParser(option_list=option_list))

# Setup data directory string
project = opt$project
comparison = opt$comparison

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

save(meth, file=sprintf('projects/%s/bisulfite/methylsig_calls/%s_raw_data.RData', project, comparison))

if(opt$tile) {
    message('Doing tiled analysis')
    meth_tiled = methylSigTile(meth, win.size=100)

    diff_meth = methylSigCalc(
        meth_tiled,
        group = c('Treatment'=1,'Control'=0),
        dispersion = opt$dispersion,
        min.per.group = min_per_group,
        num.cores= opt$ncores
        )

    write.methylSigDiff(diff_meth, file=sprintf('project/%s/bisulfite/methylsig_calls/%s.txt', project, comparison), row.names=F,quote=F,sep='\t')
} else {
    message('Doing CpG analysis')
    diff_meth = methylSigCalc(
        meth,
        group = c('Treatment'=1,'Control'=0),
        dispersion = 'both',
        min.per.group = min_per_group,
        num.cores = opt$ncores
        )

    write.methylSigDiff(diff_meth, file=sprintf('projects/%s/bisulfite/methylsig_calls/%s.txt', project, comparison), row.names=F,quote=F,sep='\t')

}
