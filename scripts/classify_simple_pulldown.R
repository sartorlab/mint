library(dplyr)
library(optparse)

option_list = list(
    make_option('--project', type='character'),
    make_option('--humanID', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))
project = opt$project
humanID = opt$humanID

# Set working directory
setwd(sprintf('~/latte/mint/%s', project))

# Setup file paths
pulldown_file = sprintf('./analysis/macs_peaks/%s_pulldown_macs2_peaks.narrowPeak', humanID)
class_bed_file = sprintf('./analysis/classification_simple/%s_pulldown_classification_simple.bed', humanID)
class_bb_file = sprintf('./analysis/summary/%s_hub/hg19/%s_pulldown_classification_simple.bb', project, humanID)

# Interpret antibody and corresponding coding
if(grepl('_hmc', humanID)) {
  antibody = 'hmc'
} else if (grepl('_mc', humanID)) {
  antibody = 'mc'
}
if(antibody == 'mc') {
    colors = c('255,102,102','255,0,0','102,0,0')
    classes = c('5mc_low','5mc_med','5mc_high')
} else if (antibody == 'hmc') {
    colors = c('102,102,255','0,0,255','0,0,102')
    classes = c('5hmc_low','5hmc_med','5hmc_high')
} else {
    stop('Antibody not recognized. Should be mc or hmc')
}

classification = data.frame(
    code = c(1,2,3),
    class = classes,
    color = colors,
    stringsAsFactors=F)

peaks = read.table(pulldown_file,
    header=F,
    sep='\t',
    quote='',
    comment.char='',
    col.names=c('chrom','start','end','name','score','strand','fold','pval','qval','summit'),
    stringsAsFactors=F)

peaks[which(peaks$score > 1000),'score'] = 1000

# Rather than doing 'tertiles' where the groups have equal number of members,
# split the range of fold change into three equal parts
peaks$fold = ceiling(peaks$fold * 10)

# Determine the 1%-tile and use the smallest fold change in the top 1%-tile
# as the maximum of the range of the fold changes.
fold_percentiles = ntile(peaks$fold, 100)
minimum_fold = min(peaks$fold)
maximum_fold = min(peaks$fold[which(fold_percentiles==100)])
fold_range = c(minimum_fold, maximum_fold)
fold_thirds = seq(floor(min(fold_range)), ceiling(max(fold_range)), by=(ceiling(max(fold_range))-floor(min(fold_range)))/3 )

# Find the indices of peaks in the first, second, and third thirds
first_third = which(peaks$fold < fold_thirds[2])
second_third = which(peaks$fold >= fold_thirds[2] & peaks$fold < fold_thirds[3])
third_third = which(peaks$fold >= fold_thirds[3])

peaks$group_fold = 1
peaks$group_fold[second_third] = 2
peaks$group_fold[third_third] = 3

merged = merge(peaks, classification, by.x='group_fold', by.y='code')
merged$thickStart = merged$start
merged$thickEnd = merged$end

merged_final = merged[,c('chrom','start','end','class','fold','strand','thickStart','thickEnd','color')]
merged_final = merged_final[order(merged_final$chrom, merged_final$start),]

write.table(merged_final, file=class_bed_file, sep='\t', row.names=F, col.names=F, quote=F)

# Convert to bigBed
message('Converting BED to bigBED...')
command = sprintf('bedToBigBed %s ~/latte/Homo_sapiens/chromInfo_hg19.txt %s', class_bed_file, class_bb_file)
system(command)
