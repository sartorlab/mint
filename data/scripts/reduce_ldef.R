# Locus definitions from chipenrich.data are not always as reduced as possible.
# For example, the following lines:
#
# chr1	10002	14318	10251
# chr1	14319	17257	653635
# chr1	17258	17905	653635
# chr1	17906	18909	653635
# chr1	18910	22329	653635
# chr1	22330	27134	653635
# chr1	27135	29664	653635
# chr1	29665	33020	653635
# chr1	33021	41081	654835
#
# indicate that a single line would suffice to describe geneid 653635. The
# purpose of this script is, where possible, to reduce such geneids to the
# minimal number of lines possible. bedtools is not a good option for this
# job because it ignores the geneid (4th) column.

library(optparse)

option_list = list(
    make_option('--inpath', type='character'),
    make_option('--outpath', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

inpath = opt$inpath
outpath = opt$outpath

if(!file.exists(inpath)) {
  stop(sprintf('The does not exist: %s', inpath))
}

library(GenomicRanges)

message('Reading data and converting to GRanges...')
data = read.table(inpath, sep='\t', header=F, comment.char='', quote='', stringsAsFactors=F)
gr_data = GRanges(seqnames=data$V2, ranges=IRanges(start=data$V3, end=data$V4), strand='*', geneid=data$V1)

# Group GRanges by the geneid
message('Splitting GRanges by geneid...')
gr_split_data = split(gr_data, gr_data$geneid)
# Within each geneid group, reduce
message('Reducing geneid chunks...')
gr_split_reduce_data = mclapply(gr_split_data, reduce, mc.cores=2)
# Call to reduce destroys geneid metadata so we have to add it back
message('Annotating geneid chunks...')
gr_split_annotated_reduce_data = lapply(names(gr_split_reduce_data), function(id) { l=gr_split_reduce_data[[id]]; l$geneid = id; return(l) })

# Use Reduce (different from reduce) to stitch together the split pieces
message('Recombining geneid chunks...')
gr_final_data = Reduce(c, gr_split_annotated_reduce_data)
# Ranges are ordered by geneid, so order by chrom + first position
message('Sorting...')
gr_final_data_sort = sort(gr_final_data)

message('Forming output data.frame...')
gr_final_data_df = as.data.frame(gr_final_data_sort)
gr_final_data_df = gr_final_data_df[,c('seqnames','start','end','geneid')]
colnames(gr_final_data_df) = c('chrom','start','end','geneid')

message(sprintf('Writing reduced locus definition to %s', outpath))
write.table(gr_final_data_df, file=outpath, sep='\t', col.names=F, row.names=F, quote=F)
