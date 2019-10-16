

## This function runs a bowtie2 alignment looking for exact matches to all KIR references in subsetKirReference
bowtie2.c4_alignment <- function(bowtie2_command, reference_index, threads, current_sample, resultsDirectory){
  
  ## Intitialize an output path for the SAM file
  current_sample[['samPath']] <- file.path(resultsDirectory,paste0(current_sample$name,'.sam'))
  
  ## Building up the run command
  optionsCommand <- c(paste0('-x ',reference_index),
                      '-5 0', '-3 6', '-N 0', '--end-to-end', paste0('-p ',threads), '--score-min "L,0,-0.10"',
                      '-I 75', '-X 1000',
                      paste0('-1 ',current_sample$fastq1path),
                      paste0('-2 ',current_sample$fastq2path),
                      '--no-unal','-a','--mp 2,2', '--rdg 1,1', '--rfg 1,1',
                      paste0('-S ',current_sample$samPath))
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2_command,optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 failed, retrying alignment...')
    output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  }
  
  check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
  
  ## Print the bowtie2 output
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  ## Check to make sure the SAM file actually exists
  current_sample[['samPath']] <- normalizePath(current_sample$samPath, mustWork=T)
  
  cat('\n\nSuccessfully aligned',current_sample$name,'to',reference_index)
  
  return(current_sample)
}

## This function returns a list of dataframes of allele sequences found in the reference fasta
read.c4_dataframe_from_reference_fasta <- function(fasta_path, referenceKeyList){
  
  ## Make sure the fasta can be found
  fasta_path <- normalizePath(fasta_path, mustWork=T)
  
  ## Initialize a list to store the sequence strings from the file
  alleleSeqList <- list()
  
  ## Read in the fasta file and store the allele names and sequences
  for(currentLine in readLines(fasta_path)){
    alleleNameBool <- grepl('>',currentLine,fixed=T)
    
    if(alleleNameBool){
      alleleName <- strsplit(currentLine, '_',fixed=TRUE)[[1]][3]
      alleleName <- strsplit(alleleName,' ',fixed=T)[[1]][1]
      alleleName <- referenceKeyList[[alleleName]]
    }else{
      alleleSeq <- currentLine
      
      if(alleleName %in% c('C4A','C4B')){
        alleleSeqList[[alleleName]] <- alleleSeq
      }
    }
  }
  
  ## Initialize the output list
  output.alleleSeqDFList <- list()
  
  cat('\n\tProcessing C4 sequence...')
  
  largestAlleleSize <- max(sapply(alleleSeqList,nchar))
  minAlleleSize <- min(sapply(alleleSeqList,nchar))
  
  ## Make sure the largest and smallest allele are the same sizes
  if(largestAlleleSize != minAlleleSize){
    stop(fasta_path,' has alleles of different sizes. Cannot transform into dataframe.')
  }
  
  ## Initialize a dataframe for storing the allele sequence strings
  alleleSeqDF <- data.frame(matrix(0, length(alleleSeqList), largestAlleleSize),row.names=names(alleleSeqList),check.names=F,stringsAsFactors=F)
  
  ## For each allele for the current locus, input the sequence into the dataframe
  for(currentAlleleName in names(alleleSeqList)){
    alleleSeqDF[currentAlleleName,] <- strsplit(alleleSeqList[[currentAlleleName]],'')[[1]]
  }
  
  return(alleleSeqDF)
}

## This function builds a list of deletion indices for each kir allele
build.c4_inverse_deletion_index_list <- function(c4AlleleDF){
  cat('\nBuilding up a deletion index for faster lookup during read assignment.')
  
  ## Initialize a list for storing deletion position conversions
  deletionIndexList <- list()
  for(currentAllele in rownames(c4AlleleDF)){
    deletionIndexList[[currentAllele]] <- which(c4AlleleDF[currentAllele,] == '.')
  }
  return(deletionIndexList)
}

## This function builds a list of deletion indices for each kir allele
build.c4_deletion_index_list <- function(c4AlleleDF){
  cat('\nBuilding up a deletion index for faster lookup during read assignment.')
  
  ## Initialize a list for storing deletion position conversions
  deletionIndexList <- list()
  for(currentAllele in rownames(c4AlleleDF)){
    deletionIndexList[[currentAllele]] <- which(c4AlleleDF[currentAllele,] != '.')
  }
  return(deletionIndexList)
}

## This function initializes a list of dataframes for storing the assembled nucleotides for each kir locus
build.c4_nuc_frame <- function(c4AlleleDF){
  refSeqDF <- data.frame(matrix(0,5,ncol(c4AlleleDF)),check.names=F,stringsAsFactors=F)
  colnames(refSeqDF) <- as.character(colnames(c4AlleleDF))
  rownames(refSeqDF) <- c('1','2','3','4','5')
  return(refSeqDF)
}

## This function counts how many reads map to a unique locus or allele
c4.count_read_matches <- function(currentSample, samTable, alignedLocusVect, maxReadThreshold){
  
  ## Pull out the unique read names
  uniqueReadNames <- unique(samTable$read_name)
  
  ## Randomize the read name order
  set.seed(001) # just to make it reproducible
  randomUniqueReadNames <- sample(uniqueReadNames)
  
  ## Check if there are more reads than the threshold and take some out if so
  if(length(randomUniqueReadNames) > maxReadThreshold){
    randomUniqueReadNames <- randomUniqueReadNames[1:maxReadThreshold]
  }
  
  ## Initialize the list for storing reference matches
  uniqueLocusMatchList <- as.list(alignedLocusVect)
  names(uniqueLocusMatchList) <- alignedLocusVect
  uniqueLocusMatchList[alignedLocusVect] <- 0
  
  
  ## Initialize variables to check on the progress of the next for loop
  i = 1
  max_i = length(randomUniqueReadNames)
  checkAmount = ceiling(max_i/10)
  j=0
  
  ## Find the reference matches for each unique read name
  for(currentReadName in randomUniqueReadNames){
    
    ## Pull out the lines of the SAM file that match the current read name
    samSubsetTable <- samTable[read_name == currentReadName]
    
    if('C4ins' %in% samSubsetTable$locus){
      uniqueLocusMatchList['C4ins'] = uniqueLocusMatchList['C4ins'][[1]] + 1
      uniqueLocusMatchList['C4'] = uniqueLocusMatchList['C4'][[1]]+1
      
      ## Will display the percent completion every 10%
      i = i+1
      if(i%%checkAmount == 0){
        j = j+10
        cat(paste0(j,'% ', collapse = ''))
      }
      
      next
    }
    
    ## Find the best alignment score for this read
    maxAlignmentScore <- max(samSubsetTable$alignment_score)
    
    ## Pull out the current read name alignments that have the best alignment score
    samSubsetTable <- samSubsetTable[alignment_score == maxAlignmentScore]
    
    ## Pull out the unique locus names
    matchedLocusList <- samSubsetTable$locus
    
    ###### This section will count all locus matches (as opposed to only unique locus matches)
    #for(matchedLocus in matchedLocusList){
    #  uniqueLocusMatchList[matchedLocus] = uniqueLocusMatchList[matchedLocus][[1]] + 1
    #}
    #if('KIR2DL5A' %in% matchedLocusList & 'KIR2DL5B' %in% matchedLocusList){
    #  uniqueLocusMatchList['KIR2DL5'] = uniqueLocusMatchList['KIR2DL5'][[1]] + 1
    #}
    ###### /s
    
    ###### This section will count only unique locus matches (as opposed to all locus matches)
    ## If there is only 1 unique locus, then add 1 to the unique match count for that locus
    
    if(length(matchedLocusList) == 1){
      uniqueLocusMatchList[matchedLocusList] = uniqueLocusMatchList[matchedLocusList][[1]] + 1
      
      ## If there is only a single matching reference allele, then iterate the count of that allele
      #if(length(matchedAlleleList) == 1){
      #  uniqueAlleleMatchList[matchedAlleleList] = uniqueAlleleMatchList[matchedAlleleList][[1]] + 1
      #}
    }else if('C4A' %in% matchedLocusList & 'C4B' %in% matchedLocusList){
      uniqueLocusMatchList['C4'] = uniqueLocusMatchList['C4'][[1]] + 1
    }
    ###### /s
    
    ## Will display the percent completion every 10%
    i = i+1
    if(i%%checkAmount == 0){
      j = j+10
      cat(paste0(j,'% ', collapse = ''))
    }
  }
  cat('\n',i)
  cat("\n\nFinished counting!")
  
  ## Cutting the function short for now to test what a good threshold would be
  return(list(locusMatches = uniqueLocusMatchList))
  #locusRatio <- sapply(uniqueLocusMatchList, function(x) x/uniqueLocusMatchList$KIR3DL3)
}

## This functions counts how many reads match CFF probes
c4.count_cff_matches <- function(currentSample,probeDF,samTable){
  cat('\n\nCounting CFF probes for the current sample..')
  
  ## Initialize the list for cff probe match counts
  cffProbeMatchList <- as.list(probeDF$Name)
  names(cffProbeMatchList) <- probeDF$Name
  cffProbeMatchList[probeDF$Name] <- 0
  
  for(probeName in probeDF$Name){
    cffProbeMatchList[probeName] <- length(grep(probeDF[probeName,'Sequence'], samTable$read_seq, fixed=T)) + as.integer(cffProbeMatchList[[probeName]])
  }
  cat('\tFinished.')
  return(cffProbeMatchList)
}

## Format reads for snp phasing (fill in deletion positions with '.')
c4.read_formatter <- function(readName, readSeq, startPos, endPos, refAlleleName, cigarString){
  
  cigarOperationVect <- c('M','I','D')
  
  cigarCheckVect <- c(cigarOperationVect,'1','2','3','4','5','6','7','8','9','0')
  
  ## Format the start position
  startPos <- as.numeric(startPos)
  
  ## Format the end position
  endPos <- as.numeric(endPos)
  
  ## Testing if there deletions found in the reference for the current read
  if(nchar(cigarString) > (nchar(nchar(readSeq)) + 1)){
    
    ## Split the cigar string up into a vector
    cigarVect <- strsplit(cigarString,'')[[1]]
    
    ## Make sure there are no weird codes that we havent accounted for
    if(!all(cigarVect %in% cigarCheckVect)){
      stop('c4.read_formatter found cigar operators that are not coded for ',readName)
    }
    
    ## Determine which position of the vector denote cigar operators
    operationPosVect <- which(cigarVect %in% cigarOperationVect)
    
    ## Initialize a list for storing the different operations and their positions
    operationList <- list()
    prevPos <- 1
    countInt <- 1
    
    ## Convert the cigar vect into the operation list
    for(operationPos in operationPosVect){
      operationList[[countInt]] <- cigarVect[prevPos:operationPos]
      prevPos <- operationPos+1
      countInt <- countInt+1
    }
    
    ## Initialize a start position for CIGAR operations
    prevPos <- 1
    subStrVect <- c()
    for(operationVect in operationList){
      ## Pull out the end position of the current operation
      operationEndPos <- as.integer(paste0(operationVect[1:(length(operationVect)-1)], collapse='')) + prevPos - 1
      
      ## Pull out the operation type
      operationTypeChar <- operationVect[length(operationVect)]
      
      ## If the end position is not an integer something went horribly wrong
      if(is.na(operationEndPos) | !(operationTypeChar %in% cigarOperationVect)){
        stop('c4.read_formatter something weird happened for ',readName)
      }
      
      if(operationTypeChar == 'M'){
        ## Add the good matches to the read string vect
        subStrVect <- c(subStrVect, substr(readSeq,prevPos,operationEndPos))
        prevPos <- operationEndPos + 1
      }else if(operationTypeChar == 'I'){
        
        ## If the insertion is in this specific position we want to keep it
        if((startPos + prevPos) == 14796 & refAlleleName == 'C4B'){
          subStrVect <- c(subStrVect, substr(readSeq,prevPos,operationEndPos))
        }
        
        ## Otherwise we ignore it
        prevPos <- operationEndPos + 1
      }else if(operationTypeChar == 'D'){
        ## Create a string of '.' of the current operation length
        delStr <- paste0(replicate((operationEndPos - prevPos) + 1, '.'), collapse='')
        
        ## Add the del string to the read string vect
        subStrVect <- c(subStrVect, delStr)
      }
    }
    
    ## Collapse the subStrVect into the read sequence
    readSeq <- paste0(subStrVect,collapse='')
  }else if((endPos-startPos+1) > nchar(readSeq)){
    
    delIndexList <- which(startPos:endPos %in% currentDelIndex[[refAlleleName]])
    
    ## If the index list is empty, something went wrong
    if(length(delIndexList) == 0){
      stop('Something weird happened with the deletion index.')
    }
    
    ## Split apart the read sequence and add in deletion symbols
    subStringList <- c()
    for(delIndex in delIndexList){
      preString <- substr(readSeq, 1, delIndex-1)
      postString <- substr(readSeq, delIndex, nchar(readSeq))
      
      readSeq <- paste0(preString, '.', postString)
    }
  }
  
  ## Create a full index that cooresponds to the sequence nucleotides and where they go
  fullIndex <- startPos:(nchar(readSeq)+startPos-1)
  
  ## Turn the sequence string into a list
  seqList <- strsplit(readSeq,'')[[1]]
  
  if(length(seqList) != length(fullIndex)){
    cat('\n',startPos,endPos,refAlleleName)
  }
  
  names(seqList) <- fullIndex
  
  seqListList <- list(seqList)
  return(seqListList)
}

