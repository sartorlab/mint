library(dplyr)
library(readr)
library(optparse)

option_list = list(
    make_option('--project', type='character'),
    make_option('--inFile', type='character'),
	make_option('--outFile', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))
project = opt$project
inFile = opt$inFile		# Can be macs2 peaks or bismark bedGraphs
outFile = opt$outFile	# Is the bed file with *_mark_pulldown_simple_classification.bed

# Set working directory
# Working directory should be projects/project

# Interpret platform
if(grepl('pulldown', inFile)) {
	platform = 'pulldown'

	# Interpret macs2 or PePr peaks
	if(grepl('PePr', inFile)) {
		peak_type = 'PePr'
	} else if (grepl('macs2', inFile)) {
		peak_type = 'macs2'
	}
} else if (grepl('bisulfite', inFile)) {
	platform = 'bisulfite'
}

# Interpret mark
# This order matters very much!
if(grepl('mc_hmc', inFile)) {
	mark = 'mc_hmc'
} else if (grepl('hmc', inFile)) {
	mark = 'hmc'
} else if (grepl('mc', inFile)) {
	mark = 'mc'
}

# Set the colors and classes according to platform and mark
# and then construct the classification table for joining
if(platform == 'pulldown') {

	if(peak_type == 'macs2') {
		column_names = c('chr','start','end','name','score','strand','fold','pval','qval','summit')
		skip = 0

		if(mark == 'mc') {
			colors = c('255,102,102','255,0,0','102,0,0')
			classes = c('mc_low','mc_med','mc_high')
		} else if (mark == 'hmc') {
			colors = c('102,102,255','0,0,255','0,0,102')
			classes = c('hmc_low','hmc_med','hmc_high')
		}
	} else if (peak_type == 'PePr') {
		column_names = c('chr','start','end','DM_status','fold','strand','pval')
		skip = 0

		if(mark == 'mc') {
			colors = c('255,102,102','255,0,0','102,0,0')
			classes = c('diff_mc_weak','diff_mc_mod','diff_mc_strong')
		} else if (mark == 'hmc') {
			colors = c('102,102,255','0,0,255','0,0,102')
			classes = c('diff_hmc_weak','diff_hmc_mod','diff_hmc_strong')
		}
	}

	classification = data.frame(
		code = c(1,2,3),
		class = classes,
		color = colors,
		stringsAsFactors=F)

} else if (platform == 'bisulfite') {
	if(mark == 'mc_hmc') {
		colors = c('0,0,0', '255,153,255', '255,0,255', '102,0,102')
		classes = c('mc_hmc_none', 'mc_hmc_low', 'mc_hmc_med', 'mc_hmc_high')
	} else if (mark == 'mc') {
		colors = c('0,0,0', '255,102,102', '255,0,0', '102,0,0')
		classes = c('mc_none', 'mc_low', 'mc_med', 'mc_high')
	} else if (mark == 'hmc') {
		colors = c('0,0,0', '102,102,255', '0,0,255', '0,0,102')
		classes = c('hmc_none', 'hmc_low', 'hmc_med', 'hmc_high')
	}

	classification = data.frame(
	    code = c(1,2,3,4),
	    class = classes,
	    color = colors,
	    stringsAsFactors=F)

	column_names = c('chr','start','end','perc_meth')
	skip = 1
}

# Read data
peaks = readr::read_tsv(inFile, col_names = column_names, skip = skip)

# Create the code column in the appropriate manner
if(platform == 'pulldown') {

	if(peak_type == 'macs2') {
		peaks[which(peaks$score > 1000),'score'] = 1000
	}

	# Rather than doing 'tertiles' where the groups have equal number of members,
	# split the range of fold change into three equal parts
	peaks$fold = ceiling(peaks$fold * 10)

	# Determine the 1%-tile and use the smallest fold change in the top 1%-tile
	# as the maximum of the range of the fold changes.
	fold_percentiles = dplyr::ntile(peaks$fold, 100)
	minimum_fold = min(peaks$fold)
	maximum_fold = min(peaks$fold[which(fold_percentiles==100)])
	fold_range = c(minimum_fold, maximum_fold)
	fold_thirds = seq(floor(min(fold_range)), ceiling(max(fold_range)), by=(ceiling(max(fold_range))-floor(min(fold_range)))/3 )

	# Find the indices of peaks in the first, second, and third thirds
	first_third = which(peaks$fold < fold_thirds[2])
	second_third = which(peaks$fold >= fold_thirds[2] & peaks$fold < fold_thirds[3])
	third_third = which(peaks$fold >= fold_thirds[3])

	# Encode the peaks
	peaks$code = 1
	peaks$code[second_third] = 2
	peaks$code[third_third] = 3

} else if (platform == 'bisulfite') {

	peaks$code = sapply(peaks$perc_meth, function(pm) {
		if(pm < 5) {
			return(1)
		} else if (pm < 33 && pm >= 5) {
			return(2)
		} else if (pm >= 33 && pm < 66) {
			return(3)
		} else if (pm >= 66) {
			return(4)
		}
	}, USE.NAMES = FALSE)

}

merged = dplyr::inner_join(x = peaks, y = classification, by='code')
merged$thickStart = merged$start
merged$thickEnd = merged$end
merged$score = 1000
merged$strand = '.'

merged_final = merged[,c('chr','start','end','class','score','strand','thickStart','thickEnd','color')]
merged_final = merged_final[order(merged_final$chr, merged_final$start),]

readr::write_tsv(merged_final, path=outFile, col_names=F)
