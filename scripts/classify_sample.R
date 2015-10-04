library(optparse)
library(GenomicRanges)
library(methylSig)

option_list = list(
    make_option('--project', type='character'),
    make_option('--bisulfiteID', type='character'),
    make_option('--pulldownID', type='character'),
    make_option('--pulldownInputID', type='character'),
    make_option('--humanID', type='character'),
    make_option('--winsize', type='integer')
)
opt = parse_args(OptionParser(option_list=option_list))

project = opt$project
bisulfiteID = opt$bisulfiteID
pulldownID = opt$pulldownID
pulldownInputID = opt$pulldownInputID
humanID = opt$humanID
winsize = opt$winsize

# Setup directory paths
basedir = '~/latte/mint'
referencedir = sprintf('%s/data/reference', basedir)
macsdir = sprintf('%s/analysis/%s/macs_peaks', basedir, project)
extractordir = sprintf('%s/analysis/%s/bismark_extractor_calls', basedir, project)
coveragesdir = sprintf('%s/analysis/%s/pulldown_coverages', basedir, project)
sampleclassifydir = sprintf('%s/analysis/%s/classification_sample', basedir, project)
hubdir = sprintf('%s/analysis/%s/summary/ucsc_trackhub/hg19', basedir, project)

# Setup file paths
pulldown_file = sprintf('%s/%s_pulldown_macs2_peaks.narrowPeak', macsdir, pulldownID)
pulldown_input_file = sprintf('%s/%s_pulldown_zero.bdg', coveragesdir, pulldownInputID)
cov_file = sprintf('%s/%s_trim.fastq.gz_bismark.bismark.cov', extractordir, bisulfiteID)
cyt_file = sprintf('%s/%s_trim.fastq.gz_bismark.CpG_report.txt', extractordir, bisulfiteID)
class_bed_file = sprintf('%s/%s_sample_classification_regions.bed', sampleclassifydir, humanID)
class_bb_file = sprintf('%s/%s_sample_classification_regions.bb', hubdir, humanID)

chr_lengths = read.table(sprintf('%s/chromInfo_hg19.txt',referencedir), header=F, sep='\t', stringsAsFactors=F)

classification = data.frame(
    code = c(0, 12, 14, 18, 21, 24, 28, 30, 35),
    class = c('5hmC','none_or_5mC','unclassifiable','none','none','5mC_(low)','5mC_or_5hmC_(low)','5mC','5mC_or_5hmC'),
    color = c('0,0,255','0,255,255','192,192,192','0,0,0','0,0,0','255,0,255','255,255,0','255,0,0','0,255,0'),
    color_name = c('blue','cyan','gray','black','black','magenta','yellow','red','green'),
    stringsAsFactors=F)

message(sprintf('Processing sample: %s',sample))

    # File checks
        if(!file.exists(pulldown_file)) {
            stop(sprintf('%s does not exist!',pulldown_file))
        }
        if(!file.exists(pulldown_input_file)) {
            stop(sprintf('%s does not exist!',pulldown_input_file))
        }
        if(!file.exists(cov_file)) {
            stop(sprintf('%s does not exist!',cov_file))
        }
        if(!file.exists(cyt_file)) {
            stop(sprintf('%s does not exist!',cyt_file))
        }

# Peaks from MACS2
# Presumably 0-based because they are narrowPeaks which are a kind of BED
    message('Reading pulldown narrowPeak...')
    pulldown = read.table(pulldown_file, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
    colnames(pulldown) = c('chrom','start','end','name','score','strand','fold','log10p','log10q','summit')

    pulldown_gr = GRanges(
        seqnames=pulldown$chrom,
        ranges=IRanges(start=pulldown$start, end=pulldown$end),
        fold=pulldown$fold,
        pval=pulldown$log10p,
        qval=pulldown$log10q)
    seqlengths(pulldown_gr) = chr_lengths[match(names(seqlengths(pulldown_gr)), chr_lengths[,1]), 2]
    genome(pulldown_gr) = 'hg19'
    pulldown_gr = trim(pulldown_gr)
    rm(pulldown)

# Input bedGraph for pulldown with only 0 coverage
# 0-based
# bedtools genomecov -bga -ibam $bamFile -g ~/latte/Methylation/Data/chromInfo_hg19.txt | grep -w '0$' > $bedFile
    message('Reading pulldown zero input regions...')
    pulldown_input = read.table(pulldown_input_file, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
    colnames(pulldown_input) = c('chrom','start','end','coverage')

    pulldown_input_gr = GRanges(
        seqnames=pulldown_input$chrom,
        ranges=IRanges(start=pulldown_input$start, end=pulldown_input$end))
    pulldown_input_gr = restrict(pulldown_input_gr, start=1, )
    seqlengths(pulldown_input_gr) = chr_lengths[match(names(seqlengths(pulldown_input_gr)), chr_lengths[,1]), 2]
    genome(pulldown_input_gr) = 'hg19'
    pulldown_input_gr = trim(pulldown_input_gr)
    rm(pulldown_input)

# Methylation calls from bismark_methylation_extractor
# 1-based?
    message('Reading bismark_methylation_extractor into methylSig...')
    bisulfite = readBismarkData(
        bismarkCovFiles=cov_file,
        cytosineCovFiles=cyt_file,
        sample.ids=sample,
        assembly='hg19',
        pipeline='bismark',
        context='CpG',
        resolution='base',
        treatment=c(1),
        destranded=TRUE,
        maxCount=500,
        minCount=10,
        filterSNPs=TRUE,
        num.cores=1,
        quiet=FALSE)

    # Tile the input
    # Each region will be classified and output as a BigBed
    bisulfite_tiled = methylSigTile(bisulfite, win.size=winsize)
    rm(bisulfite)

    bisulfite_gr = GRanges(
        seqnames=bisulfite_tiled@data.chr,
        ranges=IRanges(start=bisulfite_tiled@data.start, end=bisulfite_tiled@data.end),
        strand=bisulfite_tiled@data.strand,
        bisulfite_coverage=as.integer(bisulfite_tiled@data.coverage),
        perc_meth=round(100*(bisulfite_tiled@data.numCs / bisulfite_tiled@data.coverage),2))
    names(bisulfite_gr@elementMetadata) = c('bisulfite_coverage','perc_meth')
    seqlengths(bisulfite_gr) = chr_lengths[match(names(seqlengths(bisulfite_gr)), chr_lengths[,1]), 2]
    genome(bisulfite_gr) = 'hg19'
    bisulfite_gr = trim(bisulfite_gr)
    rm(bisulfite_tiled)

    # Add columns to include pulldown data
    bisulfite_gr$pulldown_peak = 0
    bisulfite_gr$pulldown_input = 1

# Overlaps
    message('Computing overlaps between pulldown peaks and tiled bisulfite data...')
    pulldown_bisulfite = findOverlaps(pulldown_gr, bisulfite_gr, minoverlap=40)
    pulldown_input_bisulfite = findOverlaps(pulldown_input_gr, bisulfite_gr, minoverlap=90)

# Give attributes to cpgs_gr
    message('Overlaying pulldown peak/input and bisulfite methylation information onto CpG index...')
    bisulfite_gr[pulldown_bisulfite@subjectHits]$pulldown_peak = 1
    bisulfite_gr[pulldown_input_bisulfite@subjectHits]$pulldown_input = 0

    # Consider methylation rates below 3% to be unmethylated
    bisulfite_gr[which(bisulfite_gr$perc_meth < 4)]$perc_meth = 0

# Apply coding to cpgs_gr
    message('Numerically coding regions for classification...')
    bisulfite_gr$bisulfite_code = 2
    bisulfite_gr[which(bisulfite_gr$perc_meth == 0)]$bisulfite_code = 3
    bisulfite_gr[which(bisulfite_gr$perc_meth > 0 & bisulfite_gr$perc_meth < 50)]$bisulfite_code = 4
    bisulfite_gr[which(bisulfite_gr$perc_meth != 999 & bisulfite_gr$perc_meth >= 50)]$bisulfite_code = 5

    bisulfite_gr$pulldown_code = 0
    bisulfite_gr[which(bisulfite_gr$pulldown_peak == 0 & bisulfite_gr$pulldown_input == 1)]$pulldown_code = 6
    bisulfite_gr[which(bisulfite_gr$pulldown_peak == 0 & bisulfite_gr$pulldown_input == 0)]$pulldown_code = 7

    bisulfite_gr$final_code = bisulfite_gr$bisulfite_code * bisulfite_gr$pulldown_code

# Add columns to GRanges for BED file
# Subtract one from thickStart because BEDs are supposed to be 0-based
# and the
    message('Annotating GRanges object for BED conversion...')
    bisulfite_gr$thickStart = as.integer(start(bisulfite_gr) - 1)
    bisulfite_gr$thickEnd = as.integer(end(bisulfite_gr))
    bisulfite_gr$rgb = classification[match(bisulfite_gr$final_code, classification$code),'color']
    bisulfite_gr$class = classification[match(bisulfite_gr$final_code, classification$code),'class']

# Convert to data.frame, tidy up, and write as .bed
# Note that the width column should really just be 1000 for the .bed
    message('Building BED file...')
    cpgs_df = as.data.frame(bisulfite_gr)[,c('seqnames','start','end','final_code','width','strand','thickStart','thickEnd','rgb')]
    cpgs_df$start = as.integer(cpgs_df$start - 1)
    cpgs_df$end = as.integer(cpgs_df$end)
    cpgs_df$width = 1000
    cpgs_df$strand = '.'

    message('Writing BED file...')
    write.table(cpgs_df, file=class_bed_file, row.names=F, col.names=F, quote=F, sep='\t')

# Convert to bigBed
    message('Converting BED to bigBED...')
    command = sprintf('bedToBigBed %s ~/latte/Methylation/Data/chromInfo_hg19.txt %s', class_bed_file, class_bb_file)
    system(command)

    # Add new track to the custom track file
    # cat(sprintf('track type=bigBed name=%s_class description=%s_region_classification db=hg19 itemRgb=on bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/%s',opt$sample, opt$sample, class_bb_filename),
    # file=sprintf('~/latte/Methylation/Data/GSE52945/%s_ucsc_customtracks.txt',opt$sample), sep='\n', append=T)
