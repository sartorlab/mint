################################################################################
# MAKEFILE: bisulfite_align rules

if(bool_bis_samp) {

make_var_bis_align_prefix = sprintf('
################################################################################
# Workflow for bisulfite_align
BISULFITE_ALIGN_PREFIXES := %s', paste(bisulfite_samples$fullHumanID, collapse=' '))

make_var_bis_align = 'BISULFITE_ALIGN_PREREQS := 	$(patsubst %,$(DIR_TRACK)/%_simple_classification.bb,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_CLASS_SIMPLE)/%_simple_classification.bed,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_TRACK)/%_trimmed.fq.gz_bismark_bt2.bw,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_SUM_FIGURES)/%_bismark_counts.png,$(BISULFITE_ALIGN_PREFIXES))\\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bam,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES))'

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_bis_align = '
.PHONY : bisulfite_align
bisulfite_align : $(BISULFITE_ALIGN_PREREQS)

# Rule for UCSC bigBed track of simple classifiation
$(DIR_TRACK)/%_bisulfite_simple_classification.bb : $(DIR_CLASS_SIMPLE)/%_bisulfite_simple_classification.bed
	bedToBigBed $< $(CHROM_PATH) $@

# Simple classification for percent methylation
$(DIR_CLASS_SIMPLE)/%_bisulfite_simple_classification.bed : $(DIR_BIS_BISMARK)/%_bisulfite_trimmed.fq.gz_bismark_bt2.bedGraph.gz
	Rscript ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@

# Rules for UCSC bigWig track
$(DIR_TRACK)/%_trimmed.fq.gz_bismark_bt2.bw : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph
	bedGraphToBigWig $< $(CHROM_PATH) $@

.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph.gz
	gunzip -c $< | awk \'NR > 1 {print $$0}\' | sort -T . -k1,1 -k2,2n > $@

# Rule for annotatr of extractor results
$(DIR_SUM_FIGURES)/%_bismark_counts.png : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt
	Rscript ../../scripts/annotatr_bis_align.R --file $< --genome $(GENOME)

# Rule for methylSig input
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -f ../../scripts/to_methylSig.awk $< | sort -T . -k2,2 -k3,3n > $@

# Rule for annotatr input
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -f ../../scripts/to_annotatr.awk $< | sort -T . -k1,1 -k2,2n > $@

# Rule for bismark methylation extractor
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph.gz $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bam
	cd $(DIR_BIS_BISMARK); \\
	bismark_methylation_extractor $(OPTS_EXTRACTOR) $(<F)

# Rule for bismark alignment
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bam : $(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz $(DIR_BIS_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip
	bismark $(OPTS_BISMARK) --output_dir $(@D) --temp_dir $(@D) $<
	samtools sort $@ $(patsubst %.bam,%,$@)
	samtools index $@

# Rule for FastQC on trimmed
$(DIR_BIS_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip : $(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz
	fastqc $(OPTS_FASTQC) --outdir $(@D) $<

# Rule for trim_galore
$(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz : $(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip
	trim_galore $(OPTS_TRIMGALORE_BISULFITE) --output_dir $(@D) $(DIR_BIS_RAW_FASTQS)/$*.fastq.gz

# Rule for FastQC on raw
$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip :
	fastqc $(OPTS_FASTQC) --outdir $(@D) $(DIR_BIS_RAW_FASTQS)/$*.fastq.gz
'

cat(make_var_bis_align_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_var_bis_align, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_bis_align, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
bisulfite_align_q = c(
	'#!/bin/bash',
	'#### Begin PBS preamble',
	'#PBS -N bis_align',
	'#PBS -l procs=10,mem=80gb,walltime=6:00:00',
	'#PBS -A sartor_lab',
	'#PBS -q first',
	'#PBS -M rcavalca@umich.edu',
	'#PBS -m abe',
	'#PBS -j oe',
	'#PBS -V',
	'#### End PBS preamble',
	'# Put your job commands after this line',
	sprintf('cd ~/latte/mint/projects/%s/',project),
	'make -j 2 bisulfite_align')
cat(bisulfite_align_q, file=sprintf('projects/%s/bisulfite_align.q', project), sep='\n')

for(i in 1:nrow(bisulfite_samples)) {
	# trackDb.txt entry for Bismark methylation calls
	trackEntry = c(
	  sprintf('track %s_pct_meth', bisulfite_samples[i,'fullHumanID']),
	  sprintf('parent %s_sample', bisulfite_samples[i,'humanID']),
	  sprintf('bigDataUrl %s_trimmed.fq.gz_bismark_bt2.bw', bisulfite_samples[i,'fullHumanID']),
	  sprintf('shortLabel %s_pct_meth', bisulfite_samples[i,'fullHumanID']),
	  sprintf('longLabel %s_percent_methylation', bisulfite_samples[i,'fullHumanID']),
	  'visibility full',
	  'viewLimits 0:100',
	  'type bigWig',
	  'priority 1.4',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

	# trackDb.txt entry for bisulfite simple classification results
	trackEntry = c(
	  sprintf('track %s_simple_class', bisulfite_samples[i,'fullHumanID']),
	  sprintf('parent %s_sample', bisulfite_samples[i,'humanID']),
	  sprintf('bigDataUrl %s_simple_classification.bb', bisulfite_samples[i,'fullHumanID']),
	  sprintf('shortLabel %s_simp_class', bisulfite_samples[i,'fullHumanID']),
	  sprintf('longLabel %s_simple_classification', bisulfite_samples[i,'fullHumanID']),
	  'visibility pack',
	  'itemRgb on',
	  'type bigBed 9 .',
	  'priority 1.2',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}
