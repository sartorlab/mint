library(optparse)

option_list = list(
    make_option('--project', type='character'),
    make_option('--comparison', type='character'),
    make_option('--inbed', type='character'),
    make_option('--outbed', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

project=opt$project
comparison=opt$comparison
inbed=opt$inbed
outbed=opt$outbed

# Set working directory
setwd(sprintf('~/latte/mint/%s', project))

classification = read.table(inbed, header=F, sep='\t', comment.char='', quote='', stringsAsFactors=F)
colnames(classification) = c('chr','start','end','ms_class','p_class','code')
scheme = data.frame(
    name = c(
        'hyper5mC_5hmC', 'hyper5mC_hypo5hmC', 'hyper5mC', 'hyper5mC',
        'hypo5mC_hyper5hmC', 'hypo5mC_5hmC', 'hypo5mC', 'hypo5mC',
        'hyper5hmC', 'hypo5hmC', 'noDM', 'noDM',
        'hyper5hmC', 'hypo5hmC', 'noDM', 'unclassifiable'),
    color = c(
        'dark violet', 'green', 'red', 'red',
        'yellow', 'light violet', 'orange', 'orange',
        'dark blue', 'light blue', 'black', 'black',
        'dark blue', 'light blue', 'black', 'gray'),
    rgb = c(
        '102,0,204', '0,255,0', '255,0,0', '255,0,0',
        '255,255,0', '204,153,255', '255,102,102', '255,102,102',
        '0,0,255', '102,102,255', '0,0,0', '0,0,0',
        '0,0,255', '102,102,255', '0,0,0', '192,192,192'),
    code = c(
        6, 12, 18, 24,
        10, 20, 30, 45,
        14, 28, 42, 56,
        22, 44, 66, 88),
    stringsAsFactors=F
)
merge_scheme = scheme[,c('code','name','rgb')]

classification$strand = '.'
classification$thickStart = classification$start
classification$thickEnd = classification$end
classification$score = 1000

prelim_classification = merge(classification, merge_scheme, by='code', sort=F)
final_classification = prelim_classification[, c('chr','start','end','name','score','strand','thickStart','thickEnd','rgb')]
write.table(final_classification, file=outbed, sep='\t', col.names=F, row.names=F, quote=F)
