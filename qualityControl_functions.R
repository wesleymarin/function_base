
## Check to make sure java is accessible
java <- system2('which', c('java'), stdout=T, stderr=T)
check.system2_output(java, 'java not found')

## This function marks duplicate reads
picard.mark_duplicates <- function(java_command, currentSample, resultsDirectory){
  
  ## Initialize an output path for the duplicate BAM file
  currentSample[['duplicatesBamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.dupMarked.bam'))
  ## Initialize an output path for the duplicate stats file
  currentSample[['duplicatesStatsPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.dupMetrics.tsv'))
  
  ## Building up the run command
  optionsCommand <- c('-jar','/home/LAB_PROJECTS/tools/picard.jar','MarkDuplicates',paste0('I=',currentSample$sortedBamPath),
                      paste0('O=',currentSample$duplicatesBamPath),paste0('M=',currentSample$duplicatesStatsPath))
  
  cat('\n\n',java_command, optionsCommand)
  output.markDups <- system2(java_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.markDups, 'picard tools mark duplicates failed')
  
  ## Print the output
  cat('\n',paste0(output.markDups), collapse='\n')
  
  ## Check to make sure the duplicate stats file exists
  currentSample[['duplicatesStatsPath']] <- normalizePath(currentSample$duplicatesStatsPath, mustWork=T)
  
  cat('\n\nSuccessfully marked duplicates ',currentSample$duplicatesStatsPath)
  
  return(currentSample)
}

## This function generates HsMetrics
picard.CollectHsMetrics <- function(java_command, currentSample, resultsDirectory, intervalList){
  
  ## Initialize an output path for the HsMetrics file
  currentSample[['HsMetricsPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.HsMetrics.tsv'))
  
  ## Building up the run command
  optionsCommand <- c('-jar','/home/LAB_PROJECTS/tools/picard.jar','CollectHsMetrics',paste0('I=',currentSample$sortedBamPath),
                      paste0('O=',currentSample$HsMetricsPath),'R=/home/LAB_PROJECTS/reference/fasta/hg38/allChr.fa',
                      paste0('BAIT_INTERVALS=/home/LAB_PROJECTS/reference/interval_list/',intervalList),
                      paste0('TARGET_INTERVALS=/home/LAB_PROJECTS/reference/interval_list/',intervalList)
  )
  
  cat('\n\n',java_command, optionsCommand)
  output.HsMetrics <- system2(java_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.HsMetrics, 'picard tools CollectHsMetrics failed')
  
  ## Print the output
  cat('\n',paste0(output.HsMetrics), collapse='\n')
  
  ## Check to make sure the metrics file exists
  currentSample[['HsMetricsPath']] <- normalizePath(currentSample$HsMetricsPath, mustWork=T)
  
  cat('\n\nSuccessfully ran CollectHsMetrics ',currentSample$HsMetricsPath)
  
  return(currentSample)
}

## This function reads picard stats and stores them in depthDF
collect.picard_stats <- function(currentSample, resultsDirectory, depthDF){
  cat('\nCollecting picard stats.')
  
  ## This two files need to exist
  dupMetricsPath <- currentSample$duplicatesStatsPath
  hsMetricsPath <- currentSample$HsMetricsPath
  
  ## Read in the duplication metrics
  dupMetricsDF <- read.delim(dupMetricsPath,skip=6,sep='\t',header=T,check.names=F,stringsAsFactors=F)
  
  ## Record the one duplication metric we are interested in
  depthDF[currentSample$name,'PERCENT_DUPLICATION'] <- dupMetricsDF[1,'PERCENT_DUPLICATION']
  
  ## Read in the hsMetrics
  hsMetricsDF <- read.delim(hsMetricsPath,skip=6,sep='\t',header=T,check.names=F,stringsAsFactors=F)
  
  ## Find all columns that are also in depthDF (these are the columns we want to record)
  grabCols <- intersect(colnames(hsMetricsDF), colnames(depthDF))
  
  ## Record the metrics we are interested in
  depthDF[currentSample$name,grabCols] <- hsMetricsDF[1,grabCols]
  
  cat('\nFinished with picard stats.')
  return(depthDF)
}
