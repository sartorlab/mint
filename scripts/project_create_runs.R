# The logic of this script is as follows:
#
# - Read annotation file and determine if setup is hybrid/pulldown and comparison/sample.
# - IF BISULFITE
#   - Create alignment scripts
#   - Create tracks for sample methylation rates
#   - IF COMPARISON
#     - Create methylSig scripts
#     - Create tracks for methylSig differential methylation rates
#   - ELSE (SAMPLE)
#     - Create simple classification scripts
#     - Create tracks for simple classification of sample methylation rates
# - IF PULLDOWN
#   - Create alignment scripts (both IP and input)
#   - Create tracks for pulldown coverage (both IP and input)
#   - IF COMPARISON
#     - Create PePr scripts (for both 5mC and 5hmC if !BISULFITE)
#     - Create tracks for PePr regions of differential methylation
#   - ELSE (SAMPLE)
#     - Create MACS2 scripts (for both 5mC and 5hmC if !BISULFITE)
#     - Create simple classification scripts (for both 5mC and 5hmC if !BISULFITE)
#     - Create tracks for MACS2 peaks and simple classification of sample qualitative methylation rates
# - IF BISULFITE && PULLDOWN && !COMPARISON
#   - Create hybrid sample classification scripts
#   - Create tracks for hybrid sample classifications
# - ELSE IF BISULFITE && PULLDOWN && COMPARISON
#   - Create hybrid comparison preparation scripts and hybrid comparison classification scripts
#   - Create tracks for hybrid comparison classifications
# - ELSE IF !BISULFITE && PULLDOWN && !COMPARISON
#   - Create pulldown sample classification scripts
#   - Create tracks for pulldown sample classifications
# - ELSE IF !BISULFITE && PULLDOWN && COMPARISON
#   - Create pulldown comparison preparation scripts (both 5mC and 5hmC) and pulldown comparison classification scripts
#   - Create tracks for pulldown comparison classifications

#######################################################################################################
# Argument parsing

library(optparse)

option_list = list(
    make_option('--project', type='character'),
    make_option('--comparison', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

project = opt$project
comparison = opt$comparison

#######################################################################################################
# Preamble

  # Directories
  basedir = '~/latte/mint'

  scriptdir = sprintf('%s/scripts', basedir)
  projectscriptdir = sprintf('%s/%s/scripts', basedir, project)
  analysisdir = sprintf('%s/%s/analysis', basedir, project)

  rawfastqdir = './data/raw_fastqs'
  rawfastqcsdir = './analysis/raw_fastqcs'
  trimfastqsdir= './analysis/trim_fastqs'
  trimfastqcsdir= './analysis/trim_fastqcs'
  bowtie2bamdir = './analysis/bowtie2_bams'
  pulldowncovdir = './analysis/pulldown_coverages'
  bismarkbamdir = './analysis/bismark_bams'
  extractordir = './analysis/bismark_extractor_calls'
  macsdir = './analysis/macs_peaks'
  peprdir = './analysis/pepr_peaks'
  methylsigdir = './analysis/methylsig_calls'
  classsimpledir = './analysis/classification_simple'
  classsampledir = './analysis/classification_sample'
  classcomparedir = './analysis/classification_comparison'
  hubdir = sprintf('../%s/analysis/summary/%s_hub', project, project)

  # The trackDb.txt file is used no matter the situation
  # Make sure we start fresh each time
  hubtrackdbfile = sprintf('../%s/analysis/summary/%s_hub/hg19/trackDb.txt', project, project)
  cat('', file=hubtrackdbfile)

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
  annotfile = sprintf('%s/%s/data/%s_annotation.txt', basedir, project, project)
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
  boolComp = comparison != 'none'

#######################################################################################################
# Workflow-agnostic track hub files
# hubdir/genomes.txt
# hubdir/hub.txt
# hubdir/project_name_hub.html

# genomes.txt
  genomes = c(
    'genome hg19',
    'trackDb hg19/trackDb.txt')
  cat(genomes, file=sprintf('%s/genomes.txt', hubdir), sep='\n')
# hub.txt
  hub = c(
    sprintf('hub %s_hub', project),
    sprintf('shortLabel %s', project),
    sprintf('longLabel %s', project),
    'genomesFile genomes.txt',
    'email rcavalca@umich.edu',
    sprintf('descriptionUrl %s_hub.html', project))
  cat(hub, file=sprintf('%s/hub.txt', hubdir), sep='\n')
# project_name_hub.html
  cat(' ', file=sprintf('%s/%s_hub.html',hubdir, project))

priority = 1
# superTracks
if(boolComp) {
  trackEntry = c(
    sprintf('track %s_group_comparison', comparison),
    'superTrack on show',
    sprintf('group %s_group_comparison', comparison),
    sprintf('shortLabel %s comparison', comparison),
    sprintf('longLabel %s comparison', comparison),
    sprintf('priority %s', priority),
    ' ')
  cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
  priority = priority + 1
}

for(sample in sort(unique(annotation$humanID))) {
  trackEntry = c(
    sprintf('track %s_sample', sample),
    'superTrack on show',
    sprintf('group %s_sample', sample),
    sprintf('shortLabel %s_sample', sample),
    sprintf('longLabel Sample-level tracks for %s', sample),
    sprintf('priority %s', priority),
    ' ')
  cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
  priority = priority + 1
}

#######################################################################################################
# Bisulfite scripts
  if(boolBis) {

    if(0 %in% bisulfite$mc || 0 %in% bisulfite$hmc) {
      # If a bisulfite experiment is just mc or just hmc then indicate that
      # we do not support oxBS-seq and/or TAB-seq right now
      message('oxBS-seq and TAB-seq are not yet supported in the pipeline.')
    } else {
      #######################################################################################################
      # Bisulfite experiment measures mc + hmc (RRBS, WGBS, etc.), proceed
      message('Creating bisulfite alignment scripts.')

      bisAlignScript = sprintf('%s/%s_bisulfite_alignment.sh', projectscriptdir, project)
      cat('', file=bisAlignScript)
      bisAlignVisScript = sprintf('%s/%s_visualize_bisulfite_alignment.sh', projectscriptdir, project)
      cat('', file=bisAlignVisScript)
      for(i in 1:nrow(bisulfite)) {
        command = sprintf('bash %s/process_bisulfite_align.sh -project %s -sampleID %s -humanID %s',
          scriptdir,
          project,
          bisulfite[i,'sampleID'],
          bisulfite[i,'fullHumanID'])
        cat(command, file=bisAlignScript, sep='\n', append=T)

        command = sprintf('bash %s/visualize_align_prepare_bisulfite.sh -project %s -humanID %s',
          scriptdir,
          project,
          bisulfite[i,'fullHumanID'])
        cat(command, file=bisAlignVisScript, sep='\n', append=T)
        command = sprintf('Rscript %s/visualize_align_bisulfite.R --project %s --humanID %s',
          scriptdir,
          project,
          bisulfite[i,'fullHumanID'])
        cat(command, file=bisAlignVisScript, sep='\n', append=T)

        # INSERT RSCRIPT COMMAND FOR ACTUAL VISUALIZATION HERE

        # trackDb.txt entry for Bismark methylation calls
        trackEntry = c(
          sprintf('track %s_pct_meth', bisulfite[i,'fullHumanID']),
          sprintf('parent %s_sample', bisulfite[i,'humanID']),
          sprintf('bigDataUrl %s_trim.fastq.gz_bismark.bw', bisulfite[i,'fullHumanID']),
          sprintf('shortLabel %s_pct_meth', bisulfite[i,'fullHumanID']),
          sprintf('longLabel %s_percent_methylation', bisulfite[i,'fullHumanID']),
          'visibility full',
          'viewLimits 0:100',
          'type bigWig',
          'priority 1.4',
          ' ')
        cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

      }

      if(boolComp) {
        #######################################################################################################
        message('Creating bisulfite comparison scripts.')

        # Files are in ./analysis/bismark_extractor_calls/
        cytfiles = paste(paste(extractordir, '/', bisulfite$fullHumanID, '_trim', '.fastq.gz_bismark.CpG_report_for_methylSig.txt', sep=''), collapse=',')
        samples = paste(bisulfite$fullHumanID, collapse=',')
        treatment = paste(bisulfite$group, collapse=',')
        destrand = T
        max = 500
        min = 5
        filter = T
        tile = T
        cores = 2

        bisCompareScript = sprintf('%s/%s_bisulfite_comparison.sh', projectscriptdir, project)
        cat('', file=bisCompareScript)
        command = sprintf('bash %s/process_bisulfite_comparison.sh -project %s -cyt %s -samples %s -treatment %s -destrand %s -max %s -min %s -filter %s -tile %s -cores %s -comparison %s',
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

        bisCompareVisScript = sprintf('%s/%s_visualize_bisulfite_comparison.sh', projectscriptdir, project)
        cat('', file=bisCompareVisScript)
        command = sprintf('bash %s/visualize_comparison_prepare_bisulfite.sh -project %s -comparison %s',
          scriptdir,
          project,
          comparison)
        cat(command, file=bisCompareVisScript, sep='\n', append=T)
        command = sprintf('Rscript %s/visualize_comparison_bisulfite.R --project %s --comparison %s',
          scriptdir,
          project,
          comparison)
        cat(command, file=bisCompareVisScript, sep='\n', append=T)

        # INSERT RSCRIPT COMMAND FOR ACTUAL VISUALIZATION HERE

        # trackDb.txt entry for methylSig result
        trackEntry = c(
          sprintf('track %s_DM_mSig', comparison),
          sprintf('parent %s_group_comparison', comparison),
          sprintf('bigDataUrl %s_methylSig.bw', comparison),
          sprintf('shortLabel %s_DM_mSig', comparison),
          sprintf('longLabel %s_DM_mSig_regions', comparison),
          'visibility full',
          'autoScale on',
          'alwaysZero on',
          'at y=0.0 on',
          'type bigWig',
          'priority 1.2',
          ' ')
        cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

      } else {
        #######################################################################################################
        message('No groups annotated. Creating simple classification scripts.')

        bisSimpClassScript = sprintf('%s/%s_bisulfite_classification_simple.sh', projectscriptdir, project)
        cat('', file=bisSimpClassScript)
        bisSimpClassVisScript = sprintf('%s/%s_visualize_bisulfite_classification_simple.sh', projectscriptdir, project)
        cat('', file=bisSimpClassVisScript)
        for(i in 1:nrow(bisulfite)) {
          command = sprintf('bash %s/classify_simple_bisulfite.sh -project %s -humanID %s',
            scriptdir,
            project,
            bisulfite[i,'fullHumanID'])
          cat(command, file=bisSimpClassScript, sep='\n', append=T)

          command = sprintf('bash %s/visualize_classification_prepare.sh -project %s -classBED %s',
            scriptdir,
            project,
            sprintf('./analysis/classification_simple/%s_bisulfite_classification_simple.bed', bisulfite[i,'fullHumanID']))
          cat(command, file=bisSimpClassVisScript, sep='\n', append=T)
          command = sprintf('Rscript %s/visualize_classification.R --project %s --classAnnot %s',
            scriptdir,
            project,
            sprintf('./analysis/summary/tables/%s_bisulfite_classification_simple_annot.bed', bisulfite[i,'fullHumanID']))
          cat(command, file=bisSimpClassVisScript, sep='\n', append=T)

          # INSERT RSCRIPT COMMAND FOR ACTUAL VISUALIZATION HERE

          # trackDb.txt entry for bisulfite simple classification results
          trackEntry = c(
            sprintf('track %s_simple_class', bisulfite[i,'fullHumanID']),
            sprintf('parent %s_sample', bisulfite[i,'humanID']),
            sprintf('bigDataUrl %s_bisulfite_classification_simple.bb', bisulfite[i,'fullHumanID']),
            sprintf('shortLabel %s_bis_simp_class', bisulfite[i,'fullHumanID']),
            sprintf('longLabel %s_bisulfite_simple_classification', bisulfite[i,'fullHumanID']),
            'visibility pack',
            'itemRgb on',
            'type bigBed 9 .',
            'priority 1.2',
            ' ')
          cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
        }
      }

    }

  } else {
    message('No bisulfite experiments annotated.')
  }

#######################################################################################################
# Pulldown scripts
  # If pulldown, then create pulldown alignment scripts
  if(boolPull) {
    #######################################################################################################
    message('Creating pulldown alignment scripts.')

    pullAlignScript = sprintf('%s/%s_pulldown_alignment.sh', projectscriptdir, project)
    cat('', file=pullAlignScript)
    pullAlignVisScript = sprintf('%s/%s_visualize_pulldown_alignment.sh', projectscriptdir, project)
    cat('', file=pullAlignVisScript)
    for(i in 1:nrow(pulldown)) {

      command = sprintf('bash %s/process_pulldown_align.sh -project %s -sampleID %s -humanID %s',
        scriptdir,
        project,
        pulldown[i,'sampleID'],
        pulldown[i,'fullHumanID'])
      cat(command, file=pullAlignScript, sep='\n', append=T)

      command = sprintf('bash %s/visualize_align_prepare_pulldown.sh -project %s -humanID %s',
        scriptdir,
        project,
        pulldown[i,'fullHumanID'])
      cat(command, file=pullAlignVisScript, sep='\n', append=T)
      command = sprintf('Rscript %s/visualize_align_pulldown.R --project %s --humanID %s',
        scriptdir,
        project,
        pulldown[i,'fullHumanID'])
      cat(command, file=pullAlignVisScript, sep='\n', append=T)

      # INSERT RSCRIPT COMMAND FOR ACTUAL VISUALIZATION HERE

      # trackDb.txt entry for chip/input pulldown coverages
      trackEntry = c(
        sprintf('track %s_cov', pulldown[i,'fullHumanID']),
        sprintf('parent %s_sample', pulldown[i,'humanID']),
        sprintf('bigDataUrl %s_pulldown_coverage.bw', pulldown[i,'fullHumanID']),
        sprintf('shortLabel %s_cov', pulldown[i,'fullHumanID']),
        sprintf('longLabel %s_coverage', pulldown[i,'fullHumanID']),
        'visibility full',
        'autoScale on',
        'type bigWig',
        'priority 1.6',
        ' ')
      cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

    }

    if(boolComp) {
      #######################################################################################################
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

        # Files are in ../bowtie2_bams/ relative to analysis/pepr_peaks/ (execution context)
        input1files = paste(paste('../bowtie2_bams/', input1$fullHumanID, '_pulldown_aligned.bam', sep=''), collapse=',')
        input2files = paste(paste('../bowtie2_bams/', input2$fullHumanID, '_pulldown_aligned.bam', sep=''), collapse=',')
        chip1files = paste(paste('../bowtie2_bams/', chip1$fullHumanID, '_pulldown_aligned.bam', sep=''), collapse=',')
        chip2files = paste(paste('../bowtie2_bams/', chip2$fullHumanID, '_pulldown_aligned.bam', sep=''), collapse=',')
        pulldownComparison = sprintf('%s_%s', comparison, compContext)

        # This is okay in the for loop because it creates a fresh contextual pulldown comparison
        pullCompScript = sprintf('%s/%s_%s_pulldown_comparison.sh', projectscriptdir, project, compContext)
        cat('', file=pullCompScript)
        command = sprintf('bash %s/process_pulldown_comparison.sh -project %s -input1 %s -input2 %s -chip1 %s -chip2 %s -comparison %s',
          scriptdir,
          project,
          input1files,
          input2files,
          chip1files,
          chip2files,
          pulldownComparison)
        cat(command, file=pullCompScript, sep='\n', append=T)

        pullCompareVisScript = sprintf('%s/%s_visualize_pulldown_comparison.sh', projectscriptdir, project)
        cat('', file=pullCompareVisScript)
        command = sprintf('bash %s/visualize_comparison_prepare_pulldown.sh -project %s -comparison %s',
          scriptdir,
          project,
          pulldownComparison)
        cat(command, file=pullCompareVisScript, sep='\n', append=T)
        command = sprintf('Rscript %s/visualize_comparison_pulldown.R --project %s --comparison %s',
          scriptdir,
          project,
          pulldownComparison)
        cat(command, file=pullCompareVisScript, sep='\n', append=T)

        # trackDb.txt entry for PePr output
        trackEntry = c(
          sprintf('track %s_DM_PePr', comparison),
          sprintf('parent %s_group_comparison', comparison),
          sprintf('bigDataUrl %s_%s_PePr_peaks.bb', comparison, compContext),
          sprintf('shortLabel %s_DM_PePr', comparison),
          sprintf('longLabel %s_DM_PePr_peaks', comparison),
          'visibility pack',
          'itemRgb on',
          'type bigBed 9 .',
          'priority 1.3',
          ' ')
        cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

      }

    } else {
      #######################################################################################################
      message('No groups annotated. Creating pulldown sample analysis and simple classification scripts.')

      # The humanIDs should be the same for matched samples, with differentiating
      # information provided by the rest of the annotation table.
      inputs = subset(pulldown, input == 1)
      inputsList = split(inputs, inputs$hmc)
      notinputs = subset(pulldown, input == 0)

      sampList = split(notinputs, notinputs$hmc)

      pullSampleScript = sprintf('%s/%s_pulldown_sample.sh', projectscriptdir, project)
      cat('', file=pullSampleScript)
      pullSimpClassScript = sprintf('%s/%s_pulldown_classification_simple.sh', projectscriptdir, project)
      cat('', file=pullSimpClassScript)

      pullSampleVisScript = sprintf('%s/%s_visualize_pulldown_sample.sh', projectscriptdir, project)
      cat('', file=pullSampleVisScript)
      pullSimpClassVisScript = sprintf('%s/%s_visualize_pulldown_classification_simple.sh', projectscriptdir, project)
      cat('', file=pullSimpClassVisScript)

      for(i in names(sampList)) {
        # i gives the mc/hmc context
        if(i == '0') {
          compContext = 'mc'
        } else {
          compContext = 'hmc'
        }

        subPull = sampList[[i]]

        # This detects whether we need to share the inputs, or whether they are matched
        if(!is.null(inputsList[[i]])) {
          subPull = rbind(subPull, inputsList[[i]])
        } else {
          subPull = rbind(subPull, inputs)
        }

        subPull = split(subPull, subPull$humanID)

        for(hID in names(subPull)) {
          subSampList = subPull[[hID]]

          chipID = subset(subSampList, subSampList$input == 0)$fullHumanID
          inputID = subset(subSampList, subSampList$input == 1)$fullHumanID

          # Files are in ./analysis/bowtie2_bams/
          chipFile = chipID
          inputFile = inputID

          command = sprintf('bash %s/process_pulldown_sample.sh -project %s -chip %s -input %s -name %s',
            scriptdir,
            project,
            chipFile,
            inputFile,
            chipID)
          cat(command, file=pullSampleScript, sep='\n', append=T)

          command = sprintf('bash %s/visualize_sample_prepare_pulldown.sh -project %s -humanID %s',
            scriptdir,
            project,
            chipID)
          cat(command, file=pullSampleVisScript, sep='\n', append=T)
          command = sprintf('Rscript %s/visualize_sample_pulldown.R --project %s --humanID %s',
            scriptdir,
            project,
            chipID)
          cat(command, file=pullSampleVisScript, sep='\n', append=T)

          # trackDb.txt entry for MACS2 output
          trackEntry = c(
            sprintf('track %s_peaks', chipID),
            sprintf('parent %s_sample', hID),
            sprintf('bigDataUrl %s_pulldown_macs2_peaks.bb', chipID),
            sprintf('shortLabel %s_peaks', chipID),
            sprintf('longLabel %s_MACS_peaks', chipID),
            'visibility dense',
            'type bigBed',
            'priority 1.5',
            ' ')
          cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

          command = sprintf('Rscript %s/classify_simple_pulldown.R --project %s --humanID %s',
            scriptdir,
            project,
            chipID)
          cat(command, file=pullSimpClassScript, sep='\n', append=T)

          command = sprintf('bash %s/visualize_classification_prepare.sh -project %s -classBED %s',
            scriptdir,
            project,
            sprintf('./analysis/classification_simple/%s_pulldown_classification_simple.bed', chipID))
          cat(command, file=pullSimpClassVisScript, sep='\n', append=T)
          command = sprintf('Rscript %s/visualize_classification.R --project %s --classAnnot %s',
            scriptdir,
            project,
            sprintf('./analysis/summary/tables/%s_pulldown_classification_simple_annot.bed', chipID))
          cat(command, file=pullSimpClassVisScript, sep='\n', append=T)

          # trackDb.txt entry for pulldown simple classification output
          trackEntry = c(
            sprintf('track %s_simple_class', chipID),
            sprintf('parent %s_sample', hID),
            sprintf('bigDataUrl %s_pulldown_classification_simple.bb', chipID),
            sprintf('shortLabel %s_pull_simp_class', chipID),
            sprintf('longLabel %s_pulldown_simple_classification', chipID),
            'visibility pack',
            'itemRgb on',
            'type bigBed 9 .',
            'priority 1.2',
            ' ')
          cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
        }
      }

    }

  } else {
    message('No pulldown experiments annotated.')
  }

#######################################################################################################
# Non-simple CLASSIFICATION scripts
  if(boolBis && boolPull && !boolComp) {
    #######################################################################################################
    # Prepare SAMPLE classifier scripts in HYBRID case
    message('Creating sample classifier scripts for hybrid experimental setup.')

      sampList = split(annotation, annotation$humanID)
      winsize=100

      sampleClassScript = sprintf('%s/%s_classification_sample_hybrid.sh', projectscriptdir, project)
      cat('', file=sampleClassScript)
      sampleClassVisScript = sprintf('%s/%s_visualize_classification_sample_hybrid.sh', projectscriptdir, project)
      cat('', file=sampleClassVisScript)

      for(hID in names(sampList)) {
        subSampList = sampList[[hID]]
        bisulfiteID = subset(subSampList, subSampList$bisulfite == 1)$fullHumanID
        pulldownID = subset(subSampList, subSampList$pulldown == 1 & subSampList$input == 0)$fullHumanID
        pulldownInputID = subset(subSampList, subSampList$pulldown == 1 & subSampList$input == 1)$fullHumanID

        # Note, the appropriate file paths are established in the classify_sample_hybrid.R script
        command = sprintf('Rscript %s/classify_sample_hybrid.R --project %s --bisulfiteID %s --pulldownID %s --pulldownInputID %s --humanID %s --winsize %s',
          scriptdir,
          project,
          bisulfiteID,
          pulldownID,
          pulldownInputID,
          hID,
          winsize)
        cat(command, file=sampleClassScript, sep='\n', append=T)

        command = sprintf('bash %s/visualize_classification_prepare.sh -project %s -classBED %s',
          scriptdir,
          project,
          sprintf('./analysis/classification_sample/%s_hybrid_classification_sample.bed', hID))
        cat(command, file=sampleClassVisScript, sep='\n', append=T)
        command = sprintf('Rscript %s/visualize_classification.R --project %s --classAnnot %s',
          scriptdir,
          project,
          sprintf('./analysis/summary/tables/%s_hybrid_classification_sample_annot.bed', hID))
        cat(command, file=sampleClassVisScript, sep='\n', append=T)

        # trackDb.txt entry for hybrid sample classification
        trackEntry = c(
          sprintf('track %s_hybrid_sample_classification', hID),
          sprintf('parent %s_sample', hID),
          sprintf('bigDataUrl %s_hybrid_classification_sample.bb', hID),
          sprintf('shortLabel %s_hybrid_sample_class', hID),
          sprintf('longLabel %s_hybrid_sample_classification', hID),
          'visibility pack',
          'itemRgb on',
          'type bigBed 9 .',
          'priority 1.1',
          ' ')
        cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
      }

  } else if (boolBis && boolPull && boolComp) {
    #######################################################################################################
    # Prepare COMPARISON classifier scripts in HYBRID case
    # Note we know comparison refers to bisulfite comparison name and pulldownComparison
    # is non-empty and has the correct context
    message('Creating comparison classifier scripts for hybrid experimental setup.')

      compList = split(subset(pulldown, input == 0), subset(pulldown, input == 0)$group)

      # Files are in ./analysis/pulldown_coverages/
      exp1cov = paste(paste(pulldowncovdir, '/', compList[['1']]$fullHumanID, '_pulldown_coverage.bdg', sep=''), collapse=',')
      exp2cov = paste(paste(pulldowncovdir, '/', compList[['0']]$fullHumanID, '_pulldown_coverage.bdg', sep=''), collapse=',')

      # Paths to comparison_prepare script results
      methylfiles = paste(apply(expand.grid(methylsigdir, '/', comparison, '_methylSig_', c('DM_up','DM_down','noDM_signal','noDM_nosignal'), '.bed', stringsAsFactors=F),1,paste,collapse=''), collapse=',')
      hydroxyfiles = paste(apply(expand.grid(peprdir, '/', comparison, '_hmc_PePr_', c('up_peaks','down_peaks','noDM_signal','noDM_nosignal'), '.bed', stringsAsFactors=F),1,paste,collapse=''), collapse=',')

      compClassScript = sprintf('%s/%s_classification_comparison.sh', projectscriptdir, project)
      cat('', file=compClassScript)
      command1 = sprintf('bash %s/classify_comparison_prepare_bisulfite.sh -project %s -comparison %s',
        scriptdir,
        project,
        comparison)
      command2 = sprintf('bash %s/classify_comparison_prepare_pulldown.sh -project %s -comparison %s -exp1cov %s -exp2cov %s',
        scriptdir,
        project,
        pulldownComparison,
        exp1cov,
        exp2cov)
      command3 = sprintf('bash %s/classify_comparison.sh -project %s -comparison %s -methylfiles %s -hydroxyfiles %s',
        scriptdir,
        project,
        comparison,
        methylfiles,
        hydroxyfiles)
      cat(command1, file=compClassScript, sep='\n', append=T)
      cat(command2, file=compClassScript, sep='\n', append=T)
      cat(command3, file=compClassScript, sep='\n', append=T)

      compareClassVisScript = sprintf('%s/%s_visualize_classification_comparison.sh', projectscriptdir, project)
      cat('', file=compareClassVisScript)
      command = sprintf('bash %s/visualize_classification_prepare.sh -project %s -classBED %s',
        scriptdir,
        project,
        sprintf('./analysis/classification_comparison/%s_classification_comparison.bed', comparison))
      cat(command, file=compareClassVisScript, sep='\n', append=T)
      command = sprintf('Rscript %s/visualize_classification.R --project %s --classAnnot %s',
        scriptdir,
        project,
        sprintf('./analysis/summary/tables/%s_classification_comparison_annot.bed', comparison))
      cat(command, file=compareClassVisScript, sep='\n', append=T)

      # trackDb.txt entry for hybrid comparison classification
      trackEntry = c(
        sprintf('track %s_DM_classification', comparison),
        sprintf('parent %s_group_comparison', comparison),
        sprintf('bigDataUrl %s_classification_comparison.bb', comparison),
        sprintf('shortLabel %s_hybrid_class_comp', comparison),
        sprintf('longLabel %s_hybrid_classification_comparison', comparison),
        'visibility pack',
        'itemRgb on',
        'type bigBed 9 .',
        'priority 1.1',
        ' ')
      cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

  } else if (!boolBis && boolPull && !boolComp) {
    #######################################################################################################
    # Prepare SAMPLE classifier scripts in PULLDOWN case
    message('Creating sample classifier scripts for pulldown experimental setup.')

      sampList = split(annotation, annotation$humanID)

      sampleClassScript = sprintf('%s/%s_classification_sample_pulldown.sh', projectscriptdir, project)
      cat('', file=sampleClassScript)
      sampleClassVisScript = sprintf('%s/%s_visualize_classification_sample_pulldown.sh', projectscriptdir, project)
      cat('', file=sampleClassVisScript)
      for(hID in names(sampList)) {
        subSampList = sampList[[hID]]
        mcPeaksID = subset(subSampList, subSampList$pulldown == 1 & subSampList$mc == 1)$fullHumanID
        hmcPeaksID = subset(subSampList, subSampList$pulldown == 1 & subSampList$hmc == 1)$fullHumanID
        if(nrow(subSampList) == 3) {
          # Share input
          mcZeroID = subset(subSampList, subSampList$pulldown == 1 & subSampList$input == 1)$fullHumanID
          hmcZeroID = subset(subSampList, subSampList$pulldown == 1 & subSampList$input == 1)$fullHumanID
        } else {
          # Matched input
          mcZeroID = subset(subSampList, subSampList$pulldown == 1 & subSampList$mc == 1 & subSampList$input == 1)$fullHumanID
          hmcZeroID = subset(subSampList, subSampList$pulldown == 1 & subSampList$hmc == 1 & subSampList$input == 1)$fullHumanID
        }

        # Note, the appropriate file paths are established in the classify_sample_pulldown.sh script
        command = sprintf('bash %s/classify_sample_pulldown.sh -project %s -mcPeaksID %s --mcZeroID %s -hmcPeaksID %s -hmcZeroID %s -humanID %s',
          scriptdir,
          project,
          mcPeaksID,
          mcZeroID,
          hmcPeaksID,
          hmcZeroID,
          hID)
        cat(command, file=sampleClassScript, sep='\n', append=T)

        command = sprintf('bash %s/visualize_classification_prepare.sh -project %s -classBED %s',
          scriptdir,
          project,
          sprintf('./analysis/classification_sample/%s_pulldown_classification_sample.bed', hID))
        cat(command, file=sampleClassVisScript, sep='\n', append=T)
        command = sprintf('Rscript %s/visualize_classification.R --project %s --classAnnot %s',
          scriptdir,
          project,
          sprintf('./analysis/summary/tables/%s_pulldown_classification_sample_annot.bed', hID))
        cat(command, file=sampleClassVisScript, sep='\n', append=T)

        # trackDb.txt entry for pulldown sample classification
        trackEntry = c(
          sprintf('track %s_pulldown_sample_classification', hID),
          sprintf('parent %s_sample', hID),
          sprintf('bigDataUrl %s_pulldown_classification_sample.bb', hID),
          sprintf('shortLabel %s_pull_class_sample', hID),
          sprintf('longLabel %s_pulldown_classification_sample', hID),
          'visibility pack',
          'itemRgb on',
          'type bigBed 9 .',
          'priority 1.1',
          ' ')
        cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
      }

  } else if (!boolBis && boolPull && boolComp) {
    #######################################################################################################
    # Prepare COMPARISON classifier scripts in PULLDOWN case
    message('Creating comparison classifier scripts for pulldown experimental setup.')

      compList = split(subset(pulldown, input == 0), subset(pulldown, input == 0)$group)

      # Methyl is with respect to hmc
      # The index of compList is with respect to group
      compClassScript = sprintf('%s/%s_classification_comparison_pulldown.sh', projectscriptdir, project)
      cat('', file=compClassScript)
      compClassVisScript = sprintf('%s/%s_visualize_classification_comparison_pulldown.sh', projectscriptdir, project)
      cat('', file=compClassVisScript)

      for(methyl in c(0,1)) {

        subCompList1 = subset(compList[['1']], hmc == methyl)
        subCompList2 = subset(compList[['0']], hmc == methyl)

        if(methyl == 0) {
          pulldownComparison = sprintf('%s_mc', comparison)
        } else {
          pulldownComparison = sprintf('%s_hmc', comparison)
        }

        # Paths to pulldown coverage bedgraphs
        exp1cov = paste(paste(pulldowncovdir, '/', subCompList1$fullHumanID, '_pulldown_coverage.bdg', sep=''), collapse=',')
        exp2cov = paste(paste(pulldowncovdir, '/', subCompList2$fullHumanID, '_pulldown_coverage.bdg', sep=''), collapse=',')

        command1 = sprintf('bash %s/classify_comparison_prepare_pulldown.sh -project %s -comparison %s -exp1cov %s -exp2cov %s',
          scriptdir,
          project,
          pulldownComparison,
          exp1cov,
          exp2cov)
        cat(command1, file=compClassScript, sep='\n', append=T)

      }
      # Paths to comparison_prepare script results
      methylfiles = paste(apply(expand.grid(peprdir, '/', comparison, '_mc_PePr_', c('up_peaks','down_peaks','noDM_signal','noDM_nosignal'), '.bed', stringsAsFactors=F),1,paste,collapse=''), collapse=',')
      hydroxyfiles = paste(apply(expand.grid(peprdir, '/', comparison, '_hmc_PePr_', c('up_peaks','down_peaks','noDM_signal','noDM_nosignal'), '.bed', stringsAsFactors=F),1,paste,collapse=''), collapse=',')

      command2 = sprintf('bash %s/classify_comparison.sh -project %s -comparison %s -methylfiles %s -hydroxyfiles %s',
        scriptdir,
        project,
        comparison,
        methylfiles,
        hydroxyfiles)
      cat(command2, file=compClassScript, sep='\n', append=T)

      command = sprintf('bash %s/visualize_classification_prepare.sh -project %s -classBED %s',
        scriptdir,
        project,
        sprintf('./analysis/classification_comparison/%s_classification_comparison.bed', comparison))
      cat(command, file=compClassVisScript, sep='\n', append=T)
      command = sprintf('Rscript %s/visualize_classification.R --project %s --classAnnot %s',
        scriptdir,
        project,
        sprintf('./analysis/summary/tables/%s_classification_comparison_annot.bed', comparison))
      cat(command, file=compClassVisScript, sep='\n', append=T)

      # trackDb.txt entry for pulldown comparison classification
      trackEntry = c(
        sprintf('track %s_DM_pulldown_classification', comparison),
        sprintf('parent %s_group_comparison', comparison),
        sprintf('bigDataUrl %s_classification_comparison.bb', comparison),
        sprintf('shortLabel %s_pull_class_comp', comparison),
        sprintf('longLabel %s_pulldown_classification_comparison', comparison),
        'visibility pack',
        'itemRgb on',
        'type bigBed 9 .',
        'priority 1.1',
        ' ')
      cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
  }
