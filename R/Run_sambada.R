#' @title Run sambada on parallel cores
#' @description Read samBada's input file to retrieve necessary information (number of individuals etc...), split the dataset using SamBada's Supervision tool, run sambada on the splitted dataset and merge all using Supervision. For this function you need SamBada to be installed on your computer; if this is not already the case, you can do this with downloadSambada() - for Mac users, please read the details in downloadSambada's documentation. This function produces the following output files: \code{outputFile}-Out-0.csv to \code{outputFile}-Out-\code{dimMax}.csv as well as outputFile-storey.csv (\code{outputFile} and \code{dimMax} are parameters of the function). See sambada's documentation for more information. In case you have to specify several words in one parameter, you can either specify them in one string and separate them with a space or add a vector string
#' @author Solange Duruz, Sylvie Stucki
#' @param genoFile The name of the file in the current directory of genetic information, compliant with samBada's format (use prepareGeno to transform it)
#' @param envFile  The name of the file in the current directory of environmental information (use \code{link{createEnv}} to create it and \code{link{prepareEnv}} to reduce the correlated dataset and check order)
#' @param idGeno Name of the column in the \code{genoFile} corresponding to the id of the animals
#' @param idEnv Name of the column in the \code{envFile} corresponding to the id of the animals
#' @param outputFile char Base name(s) for the results file(s). Output files will be created from the base name with suffixes (e.g. -Out-) 
#' @param dimMax Maximum number of environmental variables included in the logistic models. Use 1 for univariate models, 2 for univariate and bivariates models
#' @param cores Number of cores to use. If NULL, the #cores-1 will be used where #cores corresponds to all cores available on your computer.
#' @param wordDelim char Word delimiter of input file(s). Default ' ', 
#' @param saveType composed of three words 1) one of 'end' or 'real' to save the result during the analysis or at the end (allows sorting of result) 2) one of 'all' or 'best' to save all models or only significant models 3) If 'best' specify the threshold of significance (before applying Bonferroni's correction). Default 'END BEST 0.05', 
#' @param populationVar one of 'first' or 'last'. This option indicates whether any explanatory variables represent the population structure. If present, the said population variables must be gathered in the input file, either on the left or on the right side of the group of environmental variables. Default null.
#' @param spatial composed of 5 words 1) Column name (or number) for longitude 2) Column name (or number) for latitude 3) one of 'spherical' or 'cartesian': to indicate the type of coordinate 4) one of 'distance', 'gaussian', bisquare' or 'nearest': type of weighting scheme (see sambadoc) 5) Number bandwidth of weighting function: Units are in [m] for spherical coordinates; for cartesian coordinates, units match those of the samples' positions (see sambadoc)
#' @param autoCorr composed of 3 words. 1) one of global, local or both: to indicate the type of spatial autocorrelation to compute. 2) one of env, mark or both: to indicate the variables on which to compute the analysis 3) integer The number of permutation to compute the pseudo p-value. Ex 'global both 999'
#' @param shapeFile one of yes or no. With this option, the LISA are saved as a shapefile (in addition to the usual output)
#' @param colSupEnv char or vector of char Name(s) of the column(s) in the environmental data to be excluded from the analysis. Default NULL 
#' @param colSupMark char or vector of char Name(s) of the column(s) in the molecular data to be excluded from the analysis. Default NULL 
#' @param subsetVarEnv char or vector of char Name(s) of the column(s) in the environmental data to be included in the analysis while the other columns are set as inactive. Default NULL 
#' @param subsetVarMark char or vector of char Name(s) of the column(s) in the molecular data to be included in the analysis while the other columns are set as inactive. Default NULL 
#' @param headers logical Presence or absence of variable names in input files Default TRUE
#' @param directory char The directory where binaries of sambada are saved. This parameter is not necessary if directory path is permanently stored in the PATH environmental variable or if a function invoking sambada executable (prepareGeno or sambadaParallel) has been already run in the R active session.
#' @param keepAllFiles logical If TRUE, all parameter files and split \code{genoFile} and log-files are not removed. Default FALSE
#' @examples
#' # Example with data from the package
#' # You first need to download sambada with downloadSambada(tempdir())
#' # Example without population structure, using only one core
#' sambadaParallel(genoFile=system.file("extdata", "uganda-subset-mol.csv", package = "R.SamBada"), 
#'      envFile=system.file("extdata", "uganda-subset-env-export.csv", package = "R.SamBada"), 
#'      idGeno='ID_indiv', idEnv='short_name', dimMax=1, cores=1, saveType='END ALL', 
#'      outputFile=file.path(tempdir(),'uganda-subset-mol')) 
#' \donttest{
#' # Example with population structure, using multiple core
#' sambadaParallel(genoFile=system.file("extdata", "uganda-subset-mol.csv", package = "R.SamBada"), 
#'      envFile=system.file("extdata", "uganda-subset-env-export.csv", package = "R.SamBada"), 
#'      idGeno='ID_indiv', idEnv='short_name', dimMax=2, cores=2, saveType='END ALL', 
#'      populationVar='LAST', outputFile=file.path(tempdir(),'uganda-subset-mol'))
#' }
#' @export
sambadaParallel = function(genoFile, envFile, idGeno, idEnv, outputFile, dimMax=1, cores=NULL, wordDelim=' ', saveType='END BEST 0.05', populationVar=NULL, spatial=NULL, autoCorr=NULL, shapeFile=NULL, colSupEnv=NULL, colSupMark=NULL, subsetVarEnv=NULL, subsetVarMark=NULL, headers=TRUE, directory=NULL, keepAllFiles=FALSE){


  ### Check inputs
  required_string=c('genoFile', 'envFile', 'idGeno', 'idEnv')
  for(i in 1:length(required_string)){
    if(typeof(eval(parse(text = required_string[i])))!='character') stop(paste(required_string[i],"argument supposed to be character"))
  }

  
  if (!file.exists(genoFile)) stop("genoFile not found.")
  if (!file.exists(envFile)) stop("envFile not found.")
  if(typeof(dimMax)!='double') stop("dimMax supposed to be integer")
  if(dimMax%%1!=0) stop("dimMax supposed to be integer")
  if(dimMax>2) warning('dimMax is >2. Are you sure this is what you want?')
  if(typeof(wordDelim)!='character') stop("wordDelim supposed to be a single character")
  if(nchar(wordDelim)!=1) stop("wordDelim supposed to be a single character")
  if(!is.null(directory)){
    if(typeof(directory)!='character') stop("directory supposed to be a character")
  }

  
  add_opt=c('saveType', 'spatial', 'autoCorr', 'shapeFile','outputFile', 'colSupEnv', 'colSupMark', 'subsetVarEnv', 'subsetVarMark')
  for(i in 1:length(add_opt)){
    if(!is.null(eval(parse(text = add_opt[i])))){
      if(typeof(eval(parse(text = add_opt[i])))!='character') stop(paste(add_opt[i],"argument supposed to be character or NULL"))
      assign(add_opt[i],unlist(strsplit(eval(parse(text = add_opt[i])), split=" "))) #to build character vector if all in one string (if already vector, won't change)
    }
  }

  
  #check saveType
  saveType=toupper(saveType)
  if(!(length(saveType)==2 | length(saveType)==3)) stop("saveType should be composed of 2 or 3 words (either in vector or delimited by a space)")
  if(!(saveType[1] %in% c('END','REAL'))) stop("The first word of saveType should be either END or REAL")
  if(!(saveType[2] %in% c('ALL','SIGNIF','BEST'))) stop("The second word of saveType should be either ALL, SIGNIF or BEST")
  if(saveType[2] %in% c('SIGNIF','BEST')){
    if(!(grepl("[[:digit:]\\.-]",saveType[3]))) stop("The third word of saveType should be a number")
    #if(as.numeric(saveType[3])>1 | as.numeric(saveType[3])<0) stop("The third word of saveType should be a number between 0 and 1")
  }

  #check populationVar
  if(!is.null(populationVar)){
    if(!(populationVar %in% c('FIRST','LAST','NONE'))) stop("populationVar should be either FIRST,LAST or NONE")
  }

  #check spatial (first two columns checked later)
  if(!is.null(spatial)){
    spatial[3:5]=toupper(spatial[3:5])
    if(length(spatial)!=5) stop("spatial should be composed of 5 words (either in vector or delimited by a space)")
    if(!(spatial[3] %in% c('SPHERICAL','CARTESIAN'))) stop("The third word of spatial should be either SPHERICAL or CARTESIAN")
    if(!(spatial[4] %in% c('DISTANCE','GAUSSIAN','BISQUARE','NEAREST'))) stop("The fourth word of spatial should be either DISTANCE,GAUSSIAN,BISQUARE or NEAREST")
    if(!(grepl("[[:digit:]\\.-]",spatial[5]))) stop("The fifth word of spatial should be a number")
    if(spatial[4]=='NEAREST' & as.numeric(spatial[5])%%1!=0 ) stop("If the fourth word of spatial is NEAREST, then the fifth word an integer representing the number of neighbours")
  }
  
  #Check autoCorr
  if(!is.null(autoCorr)){
    if(is.null(spatial)) stop("Parameter spatial should be not null when autoCorr is not null")
    if(length(autoCorr)!=3) stop("autoCorr should be composed of 3 words (either in vector or delimited by a space)")
    if(!(autoCorr[1] %in% c('GLOBAL','LOCAL','BOTH'))) stop("The first word of autoCorr should be either GLOBAL, LOCAL or BOTH")
    if(!(autoCorr[2] %in% c('ENV','MARK','BOTH'))) stop("The second word of autoCorr should be either ENV, MARK or BOTH")
    if(!(grepl("[[:digit:]\\.-]",autoCorr[3]))) stop("The third word of autoCorr should be an integer representing the number of permutation")
    if(as.numeric(autoCorr[3])%%1!=0) stop("The third word of autoCorr should be an integer representing the number of permutation")
  }
  
  #Check shapeFile
  if(!is.null(shapeFile)){
    shapeFile=toupper(shapeFile)
    if(!(shapeFile %in% c('YES','NO'))) stop("shapeFile should be either YES or NO")    
  }

  #Check if sambada is installed
  #if(Sys.info()['sysname']=='Windows'){
  #  tryCatch(suppressWarnings(system2('sambada', intern=TRUE, show.output.on.console=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)), error=function(e){stop("sambada is not available. You should first download sambada and either put the binary folder to the path environmental variable or specify the path in the directory input argument")})
  #} else {
  tryCatch(suppressWarnings(system2('sambada', wait=TRUE, stdout=FALSE, stderr=FALSE)), error=function(e){stop("sambada is not available. You should first download sambada and either put the binary folder to the path environmental variable or specify the path in the directory input argument")})
  #}
  file.remove(list=list.files(pattern ='Sambada-error-log*'))
  
  ### Load required libraries if cores=NULL or cores>1
  if(is.null(cores)){
    loadParallel=TRUE
  } else if (cores>1) {
    loadParallel=TRUE
  } else {
    loadParallel=FALSE
  }
  if(loadParallel==TRUE){
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package \"doParallel\" needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package \"foreach\" needed for this function to work. Please install it.", call. = FALSE)
    }
    requireNamespace('doParallel')
  }
  
  
  
  ### Open genofile to count number of mark and indiv
  print("Setting up parameter file for Supervision") 
  
  # Count num mark and check that idGeno, colSupMark and subsetVarMark present in genoFile
  con <- file(genoFile,"r")
  first_line <- readLines(con,n=1)
  #Check column present
  if(regexpr(idGeno, first_line)[1]==-1) stop("idGeno not present in header line of genoFile")
  if(!is.null(colSupMark)){ #Check that all colSupMark in header of envFile
    if(sum(colSupMark %in% first_line)!=length(colSupMark)) stop('not all colSupMark in header line of genoFile')
  }
  if(!is.null(subsetVarMark)){ #Check that all subsetVarMark in header of envFile
    if(sum(subsetVarMark %in% first_line)!=length(subsetVarMark)) stop('not all subsetVarMark in header line of genoFile')
  }
  #Count num mark
  nummark=lengths(gregexpr(wordDelim,trimws(first_line)))
  rm(first_line)
  close(con)
  
  # numVarEnv
  env=read.csv(envFile, sep=wordDelim) 
  numVarEnv=ncol(env)
  if(!(idEnv %in% colnames(env))) stop('idEnv not in envFile')
  if(!is.null(colSupEnv)){ #Check that all colSupEnv in header of envFile
    if(sum(colSupEnv %in% colnames(env))!=length(colSupEnv)) stop('not all colSupEnv in header line of envFile')
  }
  if(!is.null(subsetVarEnv)){ #Check that all subsetVarEnv in header of envFile
    if(sum(subsetVarEnv %in% colnames(env))!=length(subsetVarEnv)) stop('not all subsetVarEnv in header line of envFile')
  }
  if(!is.null(spatial)){ #Check that all spatial in header of envFile
    if(sum(spatial[1:2] %in% colnames(env))!=length(spatial)) stop('The first two words of spatial parameter should be in header line of envFile')
  }
  

  
  # Count numIndiv from genoFile
  con <- file(genoFile,"rb")
  numIndiv = 0L
  while (length(chunk <- readBin(con, "raw", 200000)) > 0) {
    numIndiv = numIndiv + sum(chunk == as.raw(10L))
  }
  numIndiv=numIndiv - 1
  close(con)

  
  if(nrow(env)!=numIndiv){
    stop('length of envfile and genofile do not match!')
  }
  
  genoFileShort=sub('.csv','',sub('-recode','',genoFile)) #without .csv and without -recode if from prepare_geno
  
  working_dir=getwd()
  if(keepAllFiles==FALSE){
    active_dir=tempdir()
  } else {
    active_dir=getwd()
  }
  
  ### Check order ID. If not same order use prepare_env. Could also use the -env from supervision
  
  con <- file(genoFile,"r")
  #Read file line by line. Since big file, more rapid than read.table
  if(exists('IDgeno2')){
    rm(IDgeno2)
  }
  j=1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0){
    if(j==1){ #header line
      col_pos=regexpr(idGeno, line) #position beginning of Id column in header line (already checked that present)
      if(col_pos==1){
        col_pos2=0
      } else{
        col_pos2=lengths(gregexpr(wordDelim, substr(line, 1, col_pos))) #number of columns before ID
      }
    } else{ # line after header
      if(col_pos2==0){ #Id in first column
        length_col=regexpr(wordDelim,line)[1]-1
        ID=substr(line, 1, length_col)       
      } else{ #id in subsequent column
        col_delim=gregexpr(wordDelim,line)[[1]] #stores position of column delimiter
        ID=substr(line, col_delim[col_pos2]+1, col_delim[col_pos2+1]-1) 
      }
      
      if(j==2){ #First line after header, not use rbind
        IDgeno2=ID
      } else{
        IDgeno2=rbind(IDgeno2, ID) #If already exists, rbind IDgeno2
      } 
    } 
    j=j+1
  }
  rm(line)
  close(con)
  
  #If Id from geno and env file do not match => error. Sometimes one of the vector is transformed into char. Set both to char
  if(!identical(as.character(as.vector(IDgeno2)), as.character(as.vector(env[,idEnv])))){
    stop('IDs or ID order between the envfile and genofile do not match. Please double check your idGeno and idEnv input or use prepare_env function to reorder your envfile')
  }
  
  ### Preparing samBada's parameter file
  
  #Cleaning supervisions mess on.exit
  # if(keepAllFiles==FALSE){ #write to temp dir
  #   #file.remove(paste0(genoFileShort,'_paramSupervision.txt'), showWarnings = FALSE)
  #   on.exit(unlink(paste0(genoFileShort,'_paramSupervision.txt')))
  #   on.exit(unlink(paste0(genoFileShort,'_param',0:(cores-1),'.txt')))
  #   on.exit(unlink(paste0(genoFileShort,'_param.txt')))
  #   on.exit(unlink(paste0(genoFileShort,'-mark-',0:(cores-1),'-',(0:(cores-1))*sizeBlock,'.csv'))) 
  #   on.exit(unlink(paste0(genoFileShort,'-mark-',0:(cores-1),'-',(0:(cores-1))*sizeBlock,'-log.csv'))) 
  #   for(dimI in 0:dimMax){
  #     on.exit(unlink(paste0(genoFileShort,'-mark-',0:(cores-1),'-',(0:(cores-1))*sizeBlock,'-Out-',dimI,'.csv'))) 
  #   }
  # }
  
  # Sambada base parameter (nummark will be changed if supervision is used)
  params=c()
  params["HEADERS"]="Yes"
  if(wordDelim!=' '){
    params["WORDDELIM"]=wordDelim
  }
  params["NUMVARENV"]=numVarEnv
  params["NUMMARK"]=nummark+1
  params["NUMINDIV"]=numIndiv
  params["IDINDIV"]=paste(idEnv, idGeno)
  params["STOREY"]="YES"
  
  add_opt=c('dimMax', 'saveType', 'populationVar', 'spatial', 'autoCorr', 'shapeFile', 'colSupEnv', 'colSupMark', 'subsetVarEnv', 'subsetVarMark')
  for(i in 1:length(add_opt)){
    if(!is.null(eval(parse(text = add_opt[i])))){
      params[toupper(add_opt[i])]=paste(eval(parse(text = add_opt[i])), collapse=' ')
    }
  }
  
  ### Detect number of clusters for parallel computing
  
  #Define the number of cores to be used and load libraries for parallel computing
  if(is.null(cores)){
    cores=parallel::detectCores()
    cores=cores[1]-1
  }
  
  setwd(active_dir)
  if(keepAllFiles==FALSE){
    if(basename(genoFile)==genoFile | regexpr('.',genoFile, fixed=TRUE)[1]==1){ #If relative path
      file.copy(file.path(working_dir, genoFile), file.path(active_dir, genoFile), overwrite = TRUE)
    } else {
      file.copy(genoFile, file.path(active_dir, basename(genoFile)), overwrite = TRUE)
    }
  }
  on.exit(setwd(working_dir))
  
  # If cores=1 run directly sambada
  if(cores==1 | !is.null(autoCorr)){
    #write sambada parameterfile
    print("Number of cores used: 1, running sambada without splitting marker file")
    paramFile=paste0(basename(genoFileShort),'-param.txt') #name of paramFile: genoFile (without -recode and without .csv) + -param.txt
    write_sambada_parameter(paramFile, params, 'sambada')
    if(keepAllFiles==TRUE){
      system2('sambada',args=c(paramFile, paste0('"',envFile,'"'), paste0('"',genoFile,'"')))
    } else {
      if (basename(envFile)==envFile | regexpr('.',envFile, fixed=TRUE)[1]==1){
        system2('sambada',args=c(paramFile, paste0('"',file.path(working_dir,envFile),'"'), paste0('"',basename(genoFile),'"')))
      } else {
        system2('sambada',args=c(paramFile, paste0('"',envFile,'"'), paste0('"',basename(genoFile),'"')))
      }
    }
    #Copy files to right location
    if(keepAllFiles==FALSE & dirname(outputFile)!=file.path(dirname(tempdir()),basename(tempdir()))){
      if(basename(outputFile)==outputFile | regexpr('.',outputFile, fixed=TRUE)[1]==1){
        file.copy(paste0(basename(genoFileShort),'-storey.csv'), paste0(file.path(working_dir,outputFile),'-storey.csv'), overwrite=TRUE)
        for(dim in 0:dimMax){
          file.copy(file.path(active_dir,(paste0(basename(genoFileShort),'-Out-',dim,'.csv'))), paste0(file.path(working_dir,outputFile),'-Out-',dim,'.csv'), overwrite = TRUE)
        }
      } else{
        file.copy(paste0(basename(genoFileShort),'-storey.csv'), paste0(outputFile,'-storey.csv'), overwrite=TRUE)
        for(dim in 0:dimMax){
          file.copy(file.path(active_dir,(paste0(basename(genoFileShort),'-Out-',dim,'.csv'))), paste0(outputFile,'-Out-',dim,'.csv'), overwrite = TRUE)
        }
      }
    }
    
    print(paste("Output in ",paste0(outputFile,'-Out-',0:dimMax,'.csv', collapse=' and ')))
    
    return(NA)
  }
  
  ### Split datafile with supervision
  
  #Write parameter file for supervision
  fileSupervision=paste0(basename(genoFileShort),'_paramSupervision.txt')
  paramsSuper=c()
  paramsSuper["genofile"]=basename(genoFile)
  paramsSuper["paramFile"]="null"
  paramsSuper["numEnv"]=1
  paramsSuper["numMark"]=nummark
  paramsSuper["numLigne"]=numIndiv+1 #For the header line
  sizeBlock=ceiling(nummark/cores)
  paramsSuper["sizeBlock"]=sizeBlock
  write_sambada_parameter(fileSupervision, paramsSuper, 'supervision')
  
  #Run supervision
  print("Running Supervision to divide the input genomic file")
  #a=supervision(nomFichier=fileSupervision)
  #a=supervision(nomFichier=fileSupervision, numBlock = 0, blockSize = 0, maxDimension = 0, selScore = "", scoreThreshold = 0, sortScore = "", wordDelim = ' ')
  if(!is.null(directory)){
    changePath(directory)
  }
  a=system2('supervision', args=c(fileSupervision), wait=TRUE)
  ### Run sambada
  
  #Prepare Sambada's parameter file for each supervision split (file 1:n-1 are the same)
  params["NUMMARK"]=paste(sizeBlock, nummark,sep=" ") #Sabada need both the number of marker of the block + total number of markers
  params=params[names(params)!="COLSUPENV" & names(params)!="IDINDIV"] #Supervision has deleted the id from molecular file, and we can't have ID for env file only
  params["COLSUPENV"]=paste(c(colSupEnv,idEnv), collapse=" ")
  
  for(i in 0:(cores-1)){
    fileParam=paste0(basename(genoFileShort),'_param',i,'.txt')
    if(i==(cores-1)){
      SizeLastBlock=nummark-((cores-1)*sizeBlock) #Last bock, number of markers might differ
      params["NUMMARK"]=paste(SizeLastBlock, nummark,sep=" ")
    }
    write_sambada_parameter(fileParam, params, 'sambada')
  }
  
  #Build cluster
  print("Running sambada on parallel cores")
  cl = parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`
  #Close cluster on exit
  on.exit(tryCatch({ parallel::stopCluster(cl)}, error=function(e){})) #If already closed, do nothing
  
  ### Run sambada on parallel 
  # use foreach function
  ##finalMatrix = foreach(i=0:(cores-1), .combine=cbind, .packages='R.SamBada') %dopar% {tempMatrix = sambada(paste0(genoFileShort,'_param',i,'.txt'),envFile,paste0(genoFileShort,'-mark-',i,'-',i*sizeBlock,'.csv'))}
  
  if(basename(envFile)==envFile | regexpr('.',envFile, fixed=TRUE)[1]==1) {
    envFile2=file.path(working_dir,envFile)
  } else {
    envFile2=envFile
  }
  finalMatrix = foreach::foreach(i=0:(cores-1), .combine=cbind, .packages='base') %dopar% {tempMatrix = system2('sambada',args=c(paste0(basename(genoFileShort),'_param',i,'.txt'),paste0('"',envFile2,'"'),paste0(basename(genoFileShort),'-mark-',i,'-',i*sizeBlock,'.csv')))}
  
  #Close cluser
  parallel::stopCluster(cl)
  
  #Supervision merge
  print("Running supervision to merge the files")
  #supervision( base-name.txt, numBlock, blockSize, maxDimension, selScore, scoreThreshold, sortScore, wordDelim)
  #system(paste('supervision',basename(genoFile), cores, sizeBlock, dimMax,"Both", 0.001, "Wald", ' '))
  system2('supervision',args=c(basename(genoFile), cores, sizeBlock, dimMax))
  

  for(i in 0:(cores-1)){
    #add up histograms
    fileStorey=paste0(basename(genoFileShort),'-mark-',i,'-',i*sizeBlock,'-storey.csv')
    storey=read.table(fileStorey)
    if(i==0){
      storeyTot=storey
    } else {
      storeyTot=rbind(storeyTot[1:2,2:ncol(storeyTot)],storeyTot[3:nrow(storeyTot),2:ncol(storeyTot)]+storey[3:nrow(storeyTot),2:ncol(storeyTot)])
      storeyTot=cbind(storey[,1],storeyTot, row.names=NULL)
    }
    # if(keepAllFiles==FALSE){
    #   #Cleaning supervisions mess
    #   file.remove(fileStorey)
    #   file.remove(paste0(genoFileShort,'-mark-',i,'-',i*sizeBlock,'.csv'))
    #   file.remove(paste0(genoFileShort,'-mark-',i,'-',i*sizeBlock,'-log.csv'))
    #   file.remove(paste0(genoFileShort,'_param',i,'.txt'))
    #   for(dim in 0:dimMax){
    #     file.remove(paste0(genoFileShort,'-mark-',i,'-',i*sizeBlock,'-Out-',dim,'.csv'))
    #   }
    # }
    
  }
  # if(keepAllFiles==FALSE){
  #   file.remove(paste0(genoFileShort,'_paramSupervision.txt'))
  # }
  if(basename(outputFile)==outputFile | regexpr('.',outputFile, fixed=TRUE)[1]==1){  
    write.table(storeyTot, file.path(working_dir,paste0(outputFile,'-storey.csv')), row.names=FALSE, col.names=FALSE)
  } else {
    write.table(storeyTot, paste0(outputFile,'-storey.csv'), row.names=FALSE, col.names=FALSE)
  }
  
  for(dim in 0:dimMax){
    if(keepAllFiles==FALSE){
      # }
      if(basename(outputFile)==outputFile | regexpr('.',outputFile, fixed=TRUE)[1]==1){  
        file.copy(file.path(active_dir,(paste0(basename(genoFileShort),'-res-',dim,'.csv'))), file.path(working_dir,paste0(outputFile,'-Out-',dim,'.csv')), overwrite = TRUE)
      } else {
        file.copy(file.path(active_dir,(paste0(basename(genoFileShort),'-res-',dim,'.csv'))), paste0(outputFile,'-Out-',dim,'.csv'), overwrite = TRUE)
      }
    } else {
      file.rename(paste0(basename(genoFileShort),'-res-',dim,'.csv'), paste0(outputFile,'-Out-',dim,'.csv'))
    }

  }
  print(paste("Output in ",paste0(outputFile,'-Out-',0:dimMax,'.csv', collapse=' and ')))
  setwd(working_dir)
  
}

write_sambada_parameter = function(paramFile, paramMatrix, programType){
  #For sambada, write the corresponding category at beginning of lines
  if(programType=='sambada'){
    write(paste(names(paramMatrix[1]),paramMatrix[1]), paramFile, append=FALSE)
    for(i in 2:length(paramMatrix)){
      write(paste(names(paramMatrix[i]),paramMatrix[i]), paramFile, append=TRUE, sep="\n")
    }    
    # For supervision, just write the value of the parameter (without category)
  } else if(programType=='supervision'){
    write(paramMatrix[1], paramFile, append=FALSE)
    for(i in 2:length(paramMatrix)){
      write(paramMatrix[i], paramFile, append=TRUE, sep="\n")
    }  
  }
}


