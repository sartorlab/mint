library(optparse)

option_list = list(
    make_option('--project', type='character'),
    make_option('--comparison', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

project = opt$project
comparison = opt$comparison

# Directories
basedir = '~/latte/mint'
analysisdir = sprintf('%s/analysis/%s', basedir, project)
scriptdir = sprintf('%s/scripts', basedir)
projectscriptdir = sprintf('%s/scripts/%s', basedir, project)

# Expect a tab-delimited text file (${PROJECT}_annotation.txt) to be placed in data/${PROJECT} with the columns:
# 1. projectID - Same name used for project_init.sh
# 2. sampleID - SRA, GEO, Sequencing Core, etc ID (usually alpha-numeric)
# 3. humanID - Corresponding human readable name to more quickly interpret file
# 4. pulldown - 0/1 indicating data is from pulldown (1) or not (0)
# 5. bisulfite - 0/1 indicating data is bisulfite converted (1) or not (0)
# 6. mc - 0/1 indicating data uses a 5mc antibody or, if input == 1, is 5mc matched input
# 7. hmc - 0/1 indicating data uses a 5hmc antibody or, if input == 1, is 5hmc matched input
#    NOTE: These two can be 1 when pulldown == 0 once TAB-seq / oxBS-seq are reproducible
#    and an appropriate pipeline will be run for pure bisulfite conversion data.
# 8. input - 0/1 indicating data is input control for pulldown, if (input == 1 && 5mc == 0 && 5hmc == 0)
#    then the input is not matched to antibody and is used (possibly) twice, otherwise, if
#    (input == 1 && ((5mc == 1 && 5hmc == 0) || (5mc == 0 && 5hmc == 1))) then the input is antibody matched.
# 9. group - 0/1 indicates data grouping for comparison in methylSig and/ or PePr
#    if all are 0 then the assumption is a sample-wise analysis only.
annotfile = sprintf('%s/data/%s/%s_annotation.txt', basedir, project, project)
annotation = read.table(annotfile, header=T, sep='\t', quote='', comment.char='', stringsAsFactors=F)
colnames(annotation) = c('projectID','sampleID','humanID','pulldown','bisulfite','mc','hmc','input','group')

# Create fullHumanID column to simplify file tracking
annotation$fullHumanID = ''
for(i in 1:nrow(annotation)) {
  if(annotation[i,'pulldown'] == 1){
    if(annotation[i,'mc'] == 1 && annotation[i,'input'] == 0) {
      annotation[i,'fullHumanID'] = sprintf('%s_%s', annotation[i,'humanID'], 'mc')
    } else if (annotation[i,'mc'] == 1 && annotation[i,'input'] == 1) {
      annotation[i,'fullHumanID'] = sprintf('%s_%s', annotation[i,'humanID'], 'mc_input')
    } else if (annotation[i,'hmc'] == 1 && annotation[i,'input'] == 0) {
      annotation[i,'fullHumanID'] = sprintf('%s_%s', annotation[i,'humanID'], 'hmc')
    } else if (annotation[i,'hmc'] == 1 && annotation[i,'input'] == 1) {
      annotation[i,'fullHumanID'] = sprintf('%s_%s', annotation[i,'humanID'], 'hmc_input')
    } else if (annotation[i,'input'] == 1) {
      annotation[i,'fullHumanID'] = sprintf('%s_%s', annotation[i,'humanID'], 'input')
    }
  }

  if(annotation[i,'bisulfite'] == 1) {
    if(annotation[i,'mc'] == 1 && annotation[i,'hmc'] == 0) {
      annotation[i,'fullHumanID'] = sprintf('%s_%s', annotation[i,'humanID'], 'mc')
    } else if (annotation[i,'mc'] == 0 && annotation[i,'hmc'] == 1) {
      annotation[i,'fullHumanID'] = sprintf('%s_%s', annotation[i,'humanID'], 'hmc')
    } else if (annotation[i,'mc'] == 1 && annotation[i,'hmc'] == 1) {
      annotation[i,'fullHumanID'] = sprintf('%s_%s', annotation[i,'humanID'], 'mc_hmc')
    }
  }
}

##################################
# Pull out subtables
bisulfite = subset(annotation, annotation$bisulfite == 1)
pulldown = subset(annotation, annotation$pulldown == 1)

# Control checks for script creation
boolBis = nrow(bisulfite) > 0
boolPull = nrow(pulldown) > 0
boolComp = comparison != ''

# Bisulfite scripts
  if(boolBis) {

    if(0 %in% bisulfite$mc || 0 %in% bisulfite$hmc) {
      # If a bisulfite experiment is just mc or just hmc then indicate that
      # we do not support oxBS-seq and/or TAB-seq right now
      message('oxBS-seq and TAB-seq are not yet supported in the pipeline.')
    } else {
      # Bisulfite experiment measures mc + hmc (RRBS, WGBS, etc.), proceed
      message('Creating bisulfite alignment scripts.')

      bisAlignScript = sprintf('%s/%s_bisulfite_alignment.sh', projectscriptdir, project)
      for(i in 1:nrow(bisulfite)) {
        command = sprintf('sh %s/process_bisulfite_align.sh -project %s -sampleID %s -humanID %s',
          scriptdir,
          project,
          bisulfite[i,'sampleID'],
          bisulfite[i,'fullHumanID'])
        cat(command, file=bisAlignScript, sep='\n', append=T)
      }

      if(boolComp) {
        message('Creating bisulfite comparison scripts.')

        cytfiles = paste(paste(analysisdir, '/bismark_extractor_calls/',bisulfite$fullHumanID, '_trim', '.fastq.gz_bismark.CpG_report_for_methylSig.txt', sep=''), collapse=',')
        samples = paste(bisulfite$fullHumanID, collapse=',')
        treatment = paste(bisulfite$group, collapse=',')
        destrand = T
        max = 500
        min = 5
        filter = T
        tile = T
        cores = 2

        bisCompareScript = sprintf('%s/%s_bisulfite_comparison.sh', projectscriptdir, project)
        command = sprintf('sh %s/process_bisulfite_comparison.sh -project %s -cyt %s -samples %s -treatment %s -destrand %s -max %s -min %s -filter %s -tile %s -cores %s -comparison %s',
        scriptdir,
        project,
        cytfiles,
        samples,
        treatment,
        destrand,
        max,
        min,
        filter,
        tile,
        cores,
        comparison)
        cat(command, file=bisCompareScript, sep='\n', append=T)

      } else {
        message('No groups annotated. Creating simple classification scripts.')

        bisSimpClassScript = sprintf('%s/%s_bisulfite_classification_simple.sh', projectscriptdir, project)
        for(i in 1:nrow(bisulfite)) {
          command = sprintf('sh %s/classify_simple_bisulfite.sh -project %s -humanID %s',
            scriptdir,
            project,
            bisulfite[i,'fullHumanID'])
          cat(command, file=bisSimpClassScript, sep='\n', append=T)
        }
      }

    }

  } else {
    message('No bisulfite experiments annotated.')
  }

##################################
# Pulldown scripts
  # If pulldown, then create pulldown alignment scripts
  if(boolPull) {
    message('Creating pulldown alignment scripts.')

    pullAlignScript = sprintf('%s/%s_pulldown_alignment.sh', projectscriptdir, project)
    for(i in 1:nrow(pulldown)) {

      command = sprintf('sh %s/process_pulldown_align.sh -project %s -sampleID %s -humanID %s',
        scriptdir,
        project,
        pulldown[i,'sampleID'],
        pulldown[i,'fullHumanID'])
      cat(command, file=pullAlignScript, sep='\n', append=T)
    }

    if(boolComp) {
      message('Creating pulldown comparison scripts.')

      # Split pulldown into input and not input, then split the notinputs by hmc status.
      # This will catch the instance where there is only pulldown and no bisulfite
      inputs = subset(pulldown, input == 1)
      inputsList = split(inputs, inputs$hmc)
      notinputs = subset(pulldown, input == 0)

      compList = split(notinputs, notinputs$hmc)
      for(i in names(compList)) {
        # i gives the mc/hmc context
        if(i == '0') {
          compContext = 'mc'
        } else {
          compContext = 'hmc'
        }

        subpull = compList[[i]]

        # This detects whether we need to share the inputs, or whether they are matched
        if(!is.null(inputsList[[i]])) {
          subpull = rbind(subpull, inputsList[[i]])
        } else {
          subpull = rbind(subpull, inputs)
        }

        input1 = subset(subpull, input == 1 & group == 1)
        input2 = subset(subpull, input == 1 & group == 0)
        chip1 = subset(subpull, input == 0 & group == 1)
        chip2 = subset(subpull, input == 0 & group == 0)

        input1files = paste(paste(analysisdir, '/bowtie2_bams/', input1$fullHumanID, '_pulldown_aligned.bam', sep=''), collapse=',')
        input2files = paste(paste(analysisdir, '/bowtie2_bams/', input2$fullHumanID, '_pulldown_aligned.bam', sep=''), collapse=',')
        chip1files = paste(paste(analysisdir, '/bowtie2_bams/', chip1$fullHumanID, '_pulldown_aligned.bam', sep=''), collapse=',')
        chip2files = paste(paste(analysisdir, '/bowtie2_bams/', chip2$fullHumanID, '_pulldown_aligned.bam', sep=''), collapse=',')
        pulldownComparison = sprintf('%s_%s', comparison, compContext)

        pullCompScript = sprintf('%s/%s_%s_pulldown_comparison.sh', projectscriptdir, project, compContext)
        command = sprintf('sh %s/process_pulldown_comparison.sh -project %s -input1 %s -input2 %s -chip1 %s -chip2 %s -comparison %s',
          scriptdir,
          project,
          input1files,
          input2files,
          chip1files,
          chip2files,
          pulldownComparison)
        cat(command, file=pullCompScript, sep='\n', append=T)
      }

    } else {
      message('No groups annotated. Creating pulldown sample analysis and simple classification scripts.')

      # The humanIDs should be the same for matched samples, with differentiating
      # information provided by the rest of the annotation table.
      sampList = split(pulldown, pulldown$humanID)
      pullSampleScript = sprintf('%s/%s_pulldown_sample.sh', projectscriptdir, project)
      for(hID in names(sampList)) {

        subSampList = sampList[[hID]]
        chipID = subset(subSampList, subSampList$input == 0)$fullHumanID
        inputID = subset(subSampList, subSampList$input == 1)$fullHumanID

        chipFile = sprintf('%s_pulldown_aligned.bam', chipID)
        inputFile = sprintf('%s_pulldown_aligned.bam', inputID)

        command = sprintf('sh %s/process_pulldown_sample.sh -project %s -chip %s -input %s -name %s',
          scriptdir,
          project,
          chipFile,
          inputFile,
          hID)
        cat(command, file=pullSampleScript, sep='\n', append=T)

        pullSimpClassScript = sprintf('%s/%s_pulldown_classification_simple.sh', projectscriptdir, project)
        command = sprintf('Rscript %s/classify_simple_pulldown.R --project %s --humanID %s',
          scriptdir,
          project,
          chipID)
        cat(command, file=pullSimpClassScript, sep='\n', append=T)

      }

    }

  } else {
    message('No pulldown experiments annotated.')
  }

##################################
# Non-simple classification scripts
  if(boolBis && boolPull && !boolComp) {
    # Prepare sample classifier scripts in hybrid case
    message('Creating sample classifier scripts for hybrid experimental setup.')

      sampList = split(annotation, annotation$humanID)
      winsize=100

      for(hID in names(sampList)) {
        subSampList = sampList[[hID]]
        bisulfiteID = subset(subSampList, subSampList$bisulfite == 1)$fullHumanID
        pulldownID = subset(subSampList, subSampList$pulldown == 1 & subSampList$input == 0)$fullHumanID
        pulldownInputID = subset(subSampList, subSampList$pulldown == 1 & subSampList$input == 1)$fullHumanID

        sampleClassScript = sprintf('%s/%s_classification_sample.sh', projectscriptdir, project)
        command = sprintf('Rscript %s/classify_sample.R --project %s --bisulfiteID %s --pulldownID %s --pulldownInputID %s --humanID %s --winsize %s',
          scriptdir,
          project,
          bisulfiteID,
          pulldownID,
          pulldownInputID,
          hID,
          winsize)
        cat(command, file=sampleClassScript, sep='\n', append=T)
      }

  } else if (boolBis && boolPull && boolComp) {
    # Prepare comparison classifier scripts in hybrid case
    # Note we know comparison refers to bisulfite comparison name and pulldownComparison
    # is non-empty and has the correct context
    message('Creating comparison classifier scripts for hybrid experimental setup.')

      compList = split(subset(pulldown, input == 0), subset(pulldown, input == 0)$group)

      # Paths to pulldown coverage bedgraphs
      exp1cov = paste(paste(sprintf('%s/analysis/%s/pulldown_coverages/', basedir, project), compList[['1']]$fullHumanID, '_pulldown_coverage.bdg', sep=''), collapse=',')
      exp2cov = paste(paste(sprintf('%s/analysis/%s/pulldown_coverages/', basedir, project), compList[['0']]$fullHumanID, '_pulldown_coverage.bdg', sep=''), collapse=',')

      # Paths to comparison_prepare script results
      methylfiles = paste(apply(expand.grid(sprintf('%s/analysis/%s/methylsig_calls/', basedir, project), comparison, '_methylSig_', c('DM_up','DM_down','noDM_signal','noDM_nosignal'), '.bed', stringsAsFactors=F),1,paste,collapse=''), collapse=',')
      hydroxyfiles = paste(apply(expand.grid(sprintf('%s/analysis/%s/pepr_peaks/', basedir, project), comparison, '_hmc_PePr_', c('up_peaks','down_peaks','noDM_signal','noDM_nosignal'), '.bed', stringsAsFactors=F),1,paste,collapse=''), collapse=',')

      compClassScript = sprintf('%s/%s_classification_comparison.sh', projectscriptdir, project)
      command1 = sprintf('sh %s/classify_comparison_prepare_bisulfite.sh -project %s -comparison %s',
        scriptdir,
        project,
        comparison)
      command2 = sprintf('sh %s/classify_comparison_prepare_pulldown.sh -project %s -comparison %s -exp1cov %s -exp2cov %s',
        scriptdir,
        project,
        pulldownComparison,
        exp1cov,
        exp2cov)
      command3 = sprintf('sh %s/classify_comparison.sh -project %s -comparison %s -methylfiles %s -hydroxyfiles %s',
        scriptdir,
        project,
        comparison,
        methylfiles,
        hydroxyfiles)
      cat(command1, file=compClassScript, sep='\n', append=T)
      cat(command2, file=compClassScript, sep='\n', append=T)
      cat(command3, file=compClassScript, sep='\n', append=T)

  } else if (!boolBis && boolPull && !boolComp) {
    # Prepare sample classifier scripts in pulldown case
    message('WARNING: sample classifier scripts for pulldown experimental setup is not yet setup.')

  } else if (!boolBis && boolPull && boolComp) {
    # Prepare comparison classifier scripts in pulldown case
    message('Creating comparison classifier scripts for pulldown experimental setup.')

      compList = split(subset(pulldown, input == 0), subset(pulldown, input == 0)$group)

      # Methyl is with respect to hmc
      # The index of compList is with respect to group
      compClassScript = sprintf('%s/%s_classification_comparison.sh', projectscriptdir, project)
      for(methyl in c(0,1)) {

        subCompList1 = subset(compList[['1']], hmc == methyl)
        subCompList2 = subset(compList[['0']], hmc == methyl)

        if(methyl == 0) {
          pulldownComparison = sprintf('%s_mc', comparison)
        } else {
          pulldownComparison = sprintf('%s_hmc', comparison)
        }

        # Paths to pulldown coverage bedgraphs
        exp1cov = paste(paste(sprintf('%s/analysis/%s/pulldown_coverages/', basedir, project), subCompList1$fullHumanID, '_pulldown_coverage.bdg', sep=''), collapse=',')
        exp2cov = paste(paste(sprintf('%s/analysis/%s/pulldown_coverages/', basedir, project), subCompList2$fullHumanID, '_pulldown_coverage.bdg', sep=''), collapse=',')

        command1 = sprintf('sh %s/classify_comparison_prepare_pulldown.sh -project %s -comparison %s -exp1cov %s -exp2cov %ss',
          scriptdir,
          project,
          pulldownComparison,
          exp1cov,
          exp2cov)

        cat(command1, file=compClassScript, sep='\n', append=T)

      }
      # Paths to comparison_prepare script results
      methylfiles = paste(apply(expand.grid(sprintf('%s/analysis/%s/pepr_peaks/', basedir, project), comparison, '_mc_PePr_', c('up_peaks','down_peaks','noDM_signal','noDM_nosignal'), '.bed', stringsAsFactors=F),1,paste,collapse=''), collapse=',')
      hydroxyfiles = paste(apply(expand.grid(sprintf('%s/analysis/%s/pepr_peaks/', basedir, project), comparison, '_hmc_PePr_', c('up_peaks','down_peaks','noDM_signal','noDM_nosignal'), '.bed', stringsAsFactors=F),1,paste,collapse=''), collapse=',')

      command2 = sprintf('sh %s/classify_comparison.sh -project %s -comparison %s -methylfiles %s -hydroxyfiles %s',
        scriptdir,
        project,
        comparison,
        methylfiles,
        hydroxyfiles)
      cat(command2, file=compClassScript, sep='\n', append=T)
  }
