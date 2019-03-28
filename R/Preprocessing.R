#' @title Prepare genomic input
#' @description Writes a new genomic file that sambada can work with after having applied the selected genomic filtering options. The output file has the same name as the input file but with a .csv extension
#' @author Solange Duruz, Oliver Selmoni
#' @param fileName char Name of the input file (must be in active directory). Can be .gds, .ped, .bed, .vcf. If different from .gds, a gds file (SNPrelate specific format) will be created unless no filtering options are chosen
#' @param outputFile char Name of the output file. Must be a .csv
#' @param saveGDS logical If true (and if the input file extension is different from GDS) the GDS file will be saved. We recommend to set this parameter to TRUE to save time in subsequent functions that rely on GDS file
#' @param mafThresh double A number between 0 and 1 specifying the Major Allele Frequency (MAF) filtering (if null no filtering on MAF will be computed)
#' @param missingnessThresh double A number between 0 and 1 specifying the missing rate filtering (if null no filtering on missing rate will be computed)
#' @param ldThresh double A number between 0 and 1 specifying the linkage disequilibrium (LD) rate filtering (if null no filtering on LD will be computed)
#' @param mgfThresh double A number between 0 and 1 specifying the Major Genotype Frequency (MGF) rate filtering (if null no filtering on MGF will be computed). NB: sambada computations rely on genotypes
#' @param directory char The directory where binaries of sambada are saved. This parameter is not necessary if directory path is permanently stored in the PATH environmental variable or if a function invoking sambada executable (\code{prepareGeno} or \code{sambadaParallel}) has been already run in the R active session.
#' @param interactiveChecks logical If TRUE, plots will show up showing distribution of allele frequency etc... and the user can interactively change the chosen threshold for \code{mafThresh}, \code{missingnessThresh}, \code{mgfThresh} (optional, default value=FALSE)
#' @param verbose logical Turn on verbose mode
#' @return None
#' @examples
#' \dontrun{
#' #With ped input file
#' prepareGeno('myPlinkFile.ped','mySambadaFile.csv',TRUE, mafThresh=0.05, 
#'      missingnessThresh=0.05, mgfThresh=0.8,interactiveChecks=TRUE)
#'
#' #With gds input file
#' prepareGeno('myGDSFile.gds','mySambadaFile.csv',FALSE, mafThresh=0.05, 
#'      missingnessThresh=0.05,mgfThresh=0.8,interactiveChecks=FALSE)
#' }
#' @export
prepareGeno=function(fileName,outputFile,saveGDS,mafThresh=NULL, missingnessThresh=NULL,ldThresh=NULL,mgfThresh=NULL, directory=NULL, interactiveChecks=FALSE, verbose=FALSE){

  ### Check inputs ###
  
  if(typeof(fileName)!='character') stop("fileName argument supposed to be character")
  if (!file.exists(fileName)) stop("Input file not found.")
  if(typeof(outputFile)!='character') stop("outputFile argument supposed to be character")
  extensionO=substr(outputFile,gregexpr("\\.", outputFile)[[1]][length(gregexpr("\\.", outputFile)[[1]])]+1, nchar(outputFile))
  if(extensionO!='csv') stop("outputFile must have a .csv extension")
  if(typeof(saveGDS)!='logical') stop('saveGDS argument supposed to be logical')
  
  if(!is.null(mafThresh)){
    if(typeof(mafThresh)!='double') stop("mafThresh argument supposed to be decimal number")
    if(mafThresh>1 | mafThresh<0) stop("mafThresh argument supposed to be between 0 and 1")    
  }
  if(!is.null(missingnessThresh)){
    if(typeof(missingnessThresh)!='double') stop("missingnessThresh argument supposed to be decimal number")
    if(missingnessThresh>1 | missingnessThresh<0) stop("missingnessThresh argument supposed to be between 0 and 1")
  }

  if(!is.null(mgfThresh)){
    if(typeof(mgfThresh)!='double') stop("mgfThresh argument supposed to be decimal number")
    if(mgfThresh>1 | mgfThresh<0) stop("mgfThresh argument supposed to be between 0 and 1")
  }
  if(!is.null(ldThresh)){
    if(typeof(ldThresh)!='double') stop("ldThresh argument supposed to be decimal number")
    if(ldThresh>1 | ldThresh<0) stop("ldThresh argument supposed to be between 0 and 1")
  }
  if(typeof(interactiveChecks)!='logical') stop('interactiveChecks argument supposed to be logical')
  if(typeof(verbose)!='logical') stop('verbose argument supposed to be logical')
  
  # Get file extension
  filename_short=substr(fileName,1,gregexpr("\\.", fileName)[[1]][length(gregexpr("\\.", fileName)[[1]])]-1)
  extension=substr(fileName,gregexpr("\\.", fileName)[[1]][length(gregexpr("\\.", fileName)[[1]])]+1, nchar(fileName))
  
  if(!is.null(directory)){
    changePath(directory)
  }

  ### If ped file with no filters => no need to code to GDS ###
  
  if(extension=='ped' & is.null(mafThresh) & is.null(missingnessThresh) & is.null(mgfThresh) & is.null(ldThresh)){
    
    if (!file.exists(paste(filename_short,'.map',sep=''))) stop(".map input file not found. Same name as .ped mandatory")
    numSNP=nrow(read.table(paste(filename_short,'.map',sep=''), colClasses=c("character", rep("NULL",3)),sep="\t"))
    #numIndiv=nrow(read.table(filename, colClasses=c(rep("character",2), rep("NULL",4+numSNP*2)),sep=" "))
    #Calculate numIndiv by counting number of new line in pef file
    f = file(fileName, open="rb")
    numIndiv = 0L
    while (length(chunk <- readBin(f, "raw", 20000)) > 0) {
      numIndiv = numIndiv + sum(chunk == as.raw(10L))
    }
    close(f)
    
    #recodePlink(numIndiv,numSNP,filename_short,paste(filename_short,'.csv'))
    system(paste('recode-plink',numIndiv,numSNP,filename_short,outputFile))
    return(NA)
  }
  
  ### Checks that SNPRelate is loaded
  
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("Package \"SNPRelate\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("Package \"gdsfmt\" needed for this function to work. Please install it.", call. = FALSE)
  }
  
  ### Creation and opening of GDS file ###
  
  tmp=tempdir()
  current=getwd()
  
  if(verbose==TRUE){
    print('Creating and opening GDS file')
    print('===================================')
  }
  #Creates GDS file
  if(saveGDS==FALSE){
    outputDir=tmp
  } else {
    outputDir=getwd()
  }
  if(interactiveChecks==TRUE & extension!='gds' & file.exists(paste0(filename_short,'.gds'))){
    useGDS = readline(prompt="A .gds file has been found with the same name. Would you like to use it (press y) or rewrite a new file (press any other key): ")
    if(useGDS=='y'){
      gds_file=paste0(filename_short,'.gds')
    } else{
      gds_file=createGDSfile(fileName,outputDir)
    }
  } else{
    gds_file=createGDSfile(fileName,outputDir)
  }
  
  #Open gds file
  gds_obj=SNPRelate::snpgdsOpen(gds_file)
  on.exit(SNPRelate::snpgdsClose(gds_obj))
  
  ### Filtering ###
  MGF <- NULL
  Rcpp::cppFunction("
    NumericVector MGF(NumericMatrix xx, double maxMGFAllowed){
      int xrow = xx.nrow() ;
      int xcol = xx.ncol();
      int aa;
      int Aa;
      int AA;
      int sum_max;
      int mgf;
      NumericVector yy(xcol);
      NumericVector yybool(xcol);
      NumericVector yysnpid(xcol);
      int k=0;
      for(int i = 0; i < xcol; i++) {
      aa=0;
      Aa=0;
      AA=0;
      for(int j = 0; j < xrow; j++){
      if(xx(j,i)==0){
      aa++;
      }
      else if(xx(j,i)==1){
      Aa++;
      }
      else if(xx(j,i)==2){
      AA++;
      }
      }
      if(aa>=Aa){
      sum_max=aa;
      if(AA>aa){
      sum_max=AA;
      }
      }
      else{
      if(AA>=Aa){
      sum_max=AA;
      }
      else{
      sum_max=Aa;
      }
      }
      if(aa+AA+Aa>0){
      mgf=(sum_max*100)/(aa+Aa+AA);
      }
      else{
      mgf=0;
      }
      yy(i)=mgf;
      if(mgf>maxMGFAllowed*100){
      yybool(i)=0;
      
      }
      else{
      yybool(i)=1;
      yysnpid(k)=i+1;
      k++;
      }
      
      }
      
      //NumericVector yysnpid2(k);
      //yysnpid2(yysnpid.begin() , yysnpid.begin() + k);
      if(maxMGFAllowed>=0){
      return yysnpid;
      }
      else{
      return yy;
      }
  }")
  
  if(verbose==TRUE){
    print('Filtering using SNPRelate in process')
    print('===================================')  
  }
  # Interactive histograms of filtering
  if(interactiveChecks==TRUE){ #If true, shows histogram of pruning values and asks the user to check chosen thresholds
    snpSummary=SNPRelate::snpgdsSNPRateFreq(gds_obj, with.snp.id = TRUE)
    # Creates MAF histograms
    if(is.null(mafThresh)==FALSE){
      hist(snpSummary$MinorFreq, breaks=100, xlab='Minor allele Frequency', main='Histogram of minor allele frequency', xlim=c(0,max(mafThresh,max(snpSummary$MinorFreq)))) 
      abline(v=mafThresh,col="red")
      mafThresh2 = readline(prompt="Would you like to change your maf threshold? (press n if no, or enter a new threshold): ")
      if(grepl("[[:digit:]\\.-]",mafThresh2)){
        if(as.numeric(mafThresh2)>1 | as.numeric(mafThresh2)<0) stop("mafThresh argument supposed to be between 0 and 1")
        mafThresh=as.numeric(mafThresh2)
      }      
    }
    # Creates Missingness histograms
    if(is.null(missingnessThresh)==FALSE){
      hist(snpSummary$MissingRate, breaks=100, xlab='Missingness',main='Histogram of missingness', xlim=c(0,max(missingnessThresh,max(snpSummary$MissingRate)))) 
      abline(v=missingnessThresh,col="red")
      missingnessThresh2 = readline(prompt="Would you like to change your missingness threshold? (press n if no, or enter a new threshold): ")
      if(grepl("[[:digit:]\\.-]",missingnessThresh2)){
        if(as.numeric(missingnessThresh2)>1 | as.numeric(missingnessThresh2)<0) stop("missingnessThresh argument supposed to be between 0 and 1")
        missingnessThresh=as.numeric(missingnessThresh2)
      }
    }
    # Number of pruned SNP by LD
    if(is.null(ldThresh)==FALSE){
      if(length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.id")))>1000000){
        compute.LD = readline(prompt="Would you like to compute the graph of discarded SNP vs LD (takes a long time given the number of SNPs) ? (press y for yes, any other key for no): ")
      } else {
        compute.LD = 'y'
      }
      if(compute.LD=='y'){
        #ldMatrix = SNPRelate::snpgdsLDMat(gds_obj, slide=-1, method="composite")
        #image(t(ldMatrix$LD^2), col=gray(8:0 / 8))
        snp_pruned=vector()
        i=1
        for(ld in seq(0,1,0.05)){
          snp_pruned[i]=length(unlist(SNPRelate::snpgdsLDpruning(gds_obj, ld.threshold=ld, verbose=FALSE, remove.monosnp=FALSE)))
          i=i+1
        }
        barplot(snp_pruned, names.arg=seq(0,1,0.05), space=0, col="white",xlab='LD',ylab='Number of SNPs pruned',main='Histogram of number of pruned SNP by LD')
        abline(v=ldThresh*length(seq(0,1,0.05))+1-ldThresh, col="red")
        ldThresh2 = readline(prompt="Would you like to change your LD threshold? (press n if no, or enter a new threshold): ")
        if(grepl("[[:digit:]\\.-]",ldThresh2)){
          if(as.numeric(ldThresh2)>1 | as.numeric(ldThresh2)<0) stop("ldThresh argument supposed to be between 0 and 1")
          ldThresh=as.numeric(ldThresh2)
        }
      }
    }
    # Major Genotype Frequency
    if(is.null(mgfThresh)==FALSE){
      geno=SNPRelate::snpgdsGetGeno(gds_obj)
      mgf_freq=MGF(geno,-1)
      hist(mgf_freq/100, breaks=100, main='Histogram of MGF')
      abline(v=mgfThresh, col='red')
      mgfThresh2 = readline(prompt="Would you like to change your MGF threshold? (press n if no, or enter a new threshold): ")
      if(grepl("[[:digit:]\\.-]",mgfThresh2)){
        if(as.numeric(mgfThresh2)>1 | as.numeric(mgfThresh2)<0) stop("mgfThresh argument supposed to be between 0 and 1")
        mgfThresh=as.numeric(mgfThresh2)
      }
    }
  }
  
  #LD, MAF and Missing rate pruning
  if(is.null(mafThresh)){
    mafThresh=NaN
  } 
  if(is.null(missingnessThresh)){
    missingnessThresh=NaN
  } 
  if(is.null(ldThresh)){ #Manual filter or use snpgdsSelectSNP
    #if(!exists('snpSummary')){
    #  snpSummary=SNPRelate::snpgdsSNPRateFreq(gds_obj, with.snp.id=TRUE)
    #}
    #maf_filtered=snpSummary$snp.id[snpSummary$MinorFreq>mafThresh] #List of snps with maf > than threshold (pass the filter)
    #miss_filtered=snpSummary$snp.id[snpSummary$MissingRate<missingnessThresh] #List of snps with missingness < than threshold (pass the filter)
    #snp_filtered=maf_filtered[maf_filtered %in% miss_filtered] #List of snps that pass both tests
    snp_filtered=SNPRelate::snpgdsSelectSNP(gds_obj, autosome.only=FALSE, maf=mafThresh, missing.rate=missingnessThresh, verbose=FALSE)
  } else {
    gds_pruned=SNPRelate::snpgdsLDpruning(gds_obj, maf=mafThresh, missing.rate=missingnessThresh, ld.threshold=ldThresh, verbose=FALSE)
    snp_filtered=unlist(gds_pruned)
  }
  #Filter out insertion
  list_snptype=SNPRelate::snpgdsSNPList(gds_obj)
  snp_filtered2=list_snptype$snp.id[!(list_snptype$allele %in% c('A/C','C/A','A/G','G/A','A/T','T/A','C/G','G/C','C/T','T/C','G/T','T/G','A/0','C/0','G/0','T/0'))]
  snp_filtered=snp_filtered[!(snp_filtered %in% snp_filtered2)]
  
  
  #Major genotype frequency pruning
  if(!is.null(mgfThresh)){
    if(interactiveChecks==FALSE){ #If true, geno already calculated
      geno=SNPRelate::snpgdsGetGeno(gds_obj)
    }
    snpMGF=MGF(geno,mgfThresh)
    snpMGF=snpMGF[snpMGF>0] #End full of zero
    #snp present in general and MGF filters
    list_snp=snpMGF[snpMGF %in% snp_filtered] 
  }
  else{
    list_snp=snp_filtered
  }
  
  if(interactiveChecks==TRUE & file.exists(outputFile)){
    cont=readline(prompt=paste0(outputFile," already exists and will be replaced. Do you want to continue? (press x to exit, any other key to continue): "))
    if(cont=='x'){
      print('Function ended on user input')
      return(NA)
    }
  } 
  
  #Write ped file 
  SNPRelate::snpgdsGDS2PED(gds_obj, ped.fn=file.path(tmp,paste(filename_short,'_filtered',sep='')), snp.id=list_snp, verbose=FALSE)

  if(verbose==TRUE){
    print(paste0('Filtering finished: ',length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.id")))-length(list_snp),' deleted out of ',length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.id"))),' SNPs'))
    print('=====================================')
  }  
    
  ### Run recodePlink ###
  
  if(verbose==TRUE){
    print('Running Recodeplink to comply with sambada input file')
    print('=====================================')
  }
  #!Change: take pruned number of snps length(list_snp)?
  numSNP=length(list_snp)
  numIndiv=length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "sample.id")))
  #test2::recodePlink(numIndiv,numSNP,paste(filename_short,'_filtered',sep=''),paste0(filename_short,'.csv'))
  system(paste('recode-plink',numIndiv,numSNP,file.path(tmp,paste(filename_short,'_filtered',sep='')),outputFile))
}

#' @title Set the location of samples through a local web-application with interactive map
#' @description Helps the user defining the location of samples by opening a local web page. If the html fails to open, one must open georeftool.html manually in any web browser: the file is located within the extdata folder of the package. Once opened, the user must upload a file with at least one column corresponding to sample IDs. He can then specify the name of the column corresponding to lat/long if present. For samples without location, he must select the individuals on the list shown and click on a point of the map. The location of the map will be assigned to the chosen samples. When finished, the new file can be downloaded.
#' @author Oliver Selmoni, Solange Duruz
#' @examples 
#' \dontrun{
#' setLocation()
#' }
#' @export
setLocation=function(){
  html=system.file("extdata", "georeftool.html", package = "R.SamBada")
  utils::browseURL(html)
}

#' @title Create env file from raster file(s) and/or global database present in the raster r package
#' @description Create env file as an input for SamBada (it is recommended to run prepare_env function before running samBada) raster file(s) and/or global database present in the raster r package
#' @author Solange Duruz
#' @param locationFileName char Name of the file containing location of individuals. Must be in the active directory. Supported extension are .csv, .shp. All columns present in this file will also be present in the output file
#' @param outputFile char Name of the output file. Must have a .csv extension.
#' @param x char Name of the x (or longitude if not projected coordinate system) column in the \code{locationFileName}. Required if \code{locationFileName} extension is .csv
#' @param y char Name of the y (or latitude if not projected coordinate system) column in the \code{locationFileName}. Required if \code{locationFileName} extension is .csv
#' @param separator char The separator used to separate columns in your \code{locationFileName}
#' @param locationProj integer Coordinate system EPSG code of the \code{locationFileName}. If \code{locationFileName} is already georeferenced, this argument will be skipped. Required if \code{locationFileName} extension is csv.
#' @param worldclim logical If TRUE worldclim bio, tmin, tmax and prec variables will be downloaded at a resolution of 0.5 minutes of degree (the finest resolution). Rely rgdal and gdalUtils R package to merge the tiles. The downloaded tiles will be stored in the (new) wc0.5 directory of the active directory
#' @param srtm logical If TRUE the SRTM (altitude) variables will be downloaded at a resolution ... Rely rgdal and gdalUtils R package to merge the tiles. The downloaded tiles will be stored in the (new) wc0.5 directory of the active directory
#' @param saveDownload logical If TRUE (and if wordclim or srtm is TRUE), the tiles downloaded from global databases will be saved in a non-temporary directory. We recommend setting this parameter to true so that rasters can be used later (post-processing). If wordclim and srtm are FALSE, either value (TRUE/FALSE) will have no effect
#' @param rasterName char or list Name or list of name of raster files to import. Supported format are the one of raster package. If \code{directory} is TRUE then the path to the directory. Can be set to null if worldclim or srtm are set to TRUE.
#' @param rasterProj integer or list of integer Coordinate system EPSG code of the rasterlayer. If rasterlayer is already georeferenced, this argument will be skipped. If \code{rasterName} is a list, can be either a single number if all projections are the same or a list of projection for all files if different. If \code{directory} is TRUE, can only contain one number (all projections must be equal or rasters must be georeferenced)
#' @param directory logical If true, all .tif, .gtiff, .img, .sdat, . present in \code{rasterName} will be loaded
#' @param interactiveChecks logical If TRUE, shows loaded rasters and point locations
#' @param verbose logical If TRUE, indication on process will be shown
#' @return None 
#' @examples
#' \dontrun{
#' #Own raster + worldclim download
#' createEnv(rasterName=c('prec.tif','tmin.sdat'),locationFileName='MyFile.shp',
#'       outputFile='MyFile-env.csv', rasterProj=c(4326,21781), worldclim=TRUE,
#'       saveDownload=TRUE,interactiveChecks=TRUE)
#'
#' #Worldclim download only
#' createEnv(locationFileName='MyFile.csv',outputFile='MyFile-env.csv',
#'       x='Longitude',y='Latitude',locationProj=4326, 
#'       worldclim=TRUE,saveDownload=TRUE,interactiveChecks=FALSE)
#' }
#' @export
createEnv=function(locationFileName,outputFile, x=NULL,y=NULL,locationProj=NULL, separator=',', worldclim=TRUE, srtm=FALSE, saveDownload, rasterName=NULL, rasterProj=NULL,directory=FALSE, interactiveChecks, verbose=TRUE){
  
  ### Load required library
  
  #raster: always needed
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Package \"raster\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package \"sp\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if(worldclim==TRUE | srtm==TRUE){
    if (!requireNamespace("rgdal", quietly = TRUE)) {
      stop("Package \"rgdal\" needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("gdalUtils", quietly = TRUE)) {
      stop("Package \"gdalUtils\" needed for this function to work. Please install it.", call. = FALSE)
    }
  }
  
  ### Check inputs ###
  
  if(!is.null(rasterName)){
    if(typeof(rasterName)!='character') stop("rasterName is supposed to be character or vector of char")
  }
  if(!is.null(rasterProj)){
    if(typeof(rasterProj)!='double' & sum(rasterProj%%1)!=0) stop("rasterProj is supposed to be integer or vector of integers")
    if(length(rasterProj)>1 & length(rasterName)!=length(rasterProj)) stop('The size of rasterProj is supposed to be 1 or equals size of rasterName')
  }
  if(!is.null(directory)){
    if(typeof(directory)!='logical') stop("directory is supposed to be logical")
  }
  if(!is.null(locationFileName)){
    if(typeof(locationFileName)!='character') stop("locationFileName is supposed to be a character")
  }
  if(!is.null(x)){
    if(typeof(x)!='character') stop("x is supposed to be a character")
  }
  if(!is.null(y)){
    if(typeof(y)!='character') stop("y is supposed to be a character")
  }
  if(!is.null(locationProj)){
    if(typeof(locationProj)!='double' & locationProj%%1!=0) stop("locationProj is supposed to be an integer")
    if(length(locationProj)>1) stop('rasterProj is not supposed to be a vector')
  }
  if(typeof(outputFile)!='character') stop("outputFile argument supposed to be character")
  extensionO=substr(outputFile,gregexpr("\\.", outputFile)[[1]][length(gregexpr("\\.", outputFile)[[1]])]+1, nchar(outputFile))
  if(extensionO!='csv') stop("outputFile must have a .csv extension")
  
  if(typeof(saveDownload)!='logical') stop("saveDownload is supposed to be logical")
  if(!is.null(worldclim)){
    if(typeof(worldclim)!='logical') stop("worldclim is supposed to be logical")
  }
  if(!is.null(srtm)){
    if(typeof(srtm)!='logical') stop("srtm is supposed to be logical")
  }
  if(!is.null(interactiveChecks)){
    if(typeof(interactiveChecks)!='logical') stop("interactiveChecks is supposed to be logical")
  }
  if(is.null(rasterName) & worldclim==FALSE & srtm==FALSE) stop('Either provide a rasterName or set worldclim or srtm to true')
  if(directory==TRUE & is.null(rasterName)) stop('The name of the directory must be indicated in the rasterName argument if directory=TRUE')
  if(directory==FALSE & !is.null(rasterName)){
    for(i in 1:length(rasterName)){
      if(!file.exists(rasterName[i])) stop(paste0(rasterName[i],' not found'))
    }
  }
  if(!file.exists(locationFileName))stop("locationFileName not found")

  locationextension=substr(locationFileName,gregexpr("\\.", locationFileName)[[1]][length(gregexpr("\\.", locationFileName)[[1]])]+1, nchar(locationFileName))
  if(locationextension=='csv'){
    if(is.null(x))stop("x must be provided if locationFileName is a .CSV")
    if(is.null(y))stop("y must be provided if locationFileName is a .CSV")
    if(is.null(locationProj))stop("locationProj must be provided if locationFileName is a .CSV")
  }
  
  #Save active directory (will be changed)
  working_dir=getwd()
  if(saveDownload==FALSE){
    active_dir=tempdir()
  } else {
    active_dir=getwd()
  }
  
  
  ### Open location file and set projection ###
  
  if(verbose==TRUE){
    print('loading location file')
  }
  
  
  locations=openEnvData(locationFileName, separator = separator)
  #Create data to be exported (original columns of locationfile for now)
  data=locations
  #Transform locations into a spatial object.! If shapefile
  if(locationextension=='csv'){
    sp::coordinates(locations) = c(x,y)
  }
  
  #If location file is not georeferenced and locationProj is given, assign projection
  if(is.na(raster::projection(locations)) & !is.null(locationProj)){
    raster::projection(locations) = sp::CRS(paste0('+init=epsg:',locationProj))
  }
  
  ### Rasters provided by the user
  
  if(!is.null(rasterName)){
    
    if(verbose==TRUE){
      print('loading rasters provided in rasterName')
    }
    
    if(directory==TRUE){
      
      #Find all raster files in the directory
      setwd(rasterName)
      files=list.files(pattern = "\\.tif$")
      files=c(files, list.files(pattern = "\\.gtiff$")) 
      files=c(files, list.files(pattern = "\\.img$")) 
      files=c(files, list.files(pattern = "\\.sdat$")) 
      files=c(files, list.files(pattern = "\\.ascii$")) 
      files=c(files, list.files(pattern = "\\.grd$")) 
      files=c(files, list.files(pattern = "\\.bil$")) 
      rasterName=files
    
    }
    
    for(i in 1:length(rasterName)){
      
      rasterName2=rasterName[i]
      rasterName2_short=substr(rasterName2,1,gregexpr("\\.", rasterName2)[[1]][length(gregexpr("\\.", rasterName2)[[1]])]-1)
      if(!file.exists(rasterName2))stop(paste(rasterName2,"not found"))
      #Load raster using raster library
      raster=raster::raster(rasterName2)
      
      #Set projection if not georeferenced
      if(is.na(raster::projection(raster)) & is.null(rasterProj)) stop(paste0("Rasters must be georeferenced or rasterProj must be provided! Check file ",rasterName2))
      if(is.na(raster::projection(raster)) & !is.null(rasterProj)){
        if(length(rasterProj)>1){ #one proj per raster
          rasterProj2=rasterProj[i]
        }
        else{ # one proj for all rasters
          rasterProj2=rasterProj
        }
        raster::projection(raster) = sp::CRS(paste0('+init=epsg:',rasterProj2))
      }
      
      #Plot raster
      if(interactiveChecks==TRUE){
        plot(raster, y=1, main="First band of ") #raster::?
        plot(sp::spTransform(locations, raster::projection(raster)), add=TRUE)  
        proceed=readline("Would you like to proceed? (Press x to exit or any other key to continue) ")
        if(proceed=='x'){
          print('Function ended on user input')
          return(NA)
        }
      }
      
      #Extract data at the location of points (ad-hoc function)
      extractdataraster=extractVar(raster, locations)
      #Save column names
      colnames_data=colnames(data)
      data = data.frame(data, extractdataraster)
      #Specify name of columns
      colnames(data)=c(colnames_data, rasterName2_short)
    }
  }
  
  ### Worldclim ###
  
  if(worldclim==TRUE){
    
    if(verbose==TRUE){
      print('Downloading variables from Worldclim database')
    }
    
    setwd(active_dir)
    
    #Load wordclim
    locations2=sp::spTransform(locations,sp::CRS('+init=epsg:4326'))
    coord=sp::coordinates(locations2)
    
    for(i in 1:nrow(locations)){
      latt=coord[i,y]
      long=coord[i,x]
      
      #Download worldclim data (if already downloaded, will not download it twice). Stored in directory wc0.5 of active directory
      raster=raster::getData("worldclim",var="bio",path=active_dir,res=0.5, lon=long, lat=latt)
      raster=raster::getData("worldclim",var="tmin",path=active_dir,res=0.5, lon=long, lat=latt)
      raster=raster::getData("worldclim",var="tmax",path=active_dir,res=0.5, lon=long, lat=latt)
      raster=raster::getData("worldclim",var="prec",path=active_dir,res=0.5, lon=long, lat=latt)
    }
    
    setwd(paste0(active_dir,'/wc0.5'))
    
    for(i in 1:19){ #Loop on number of bio
      
      #Search all raster files of biox
      files=list.files(pattern = paste0('bio',i,'_[0-9]*\\.bil'))
      #merge them using library gdalUtils
      gdalUtils::mosaic_rasters(files, paste0('bio',i,'.tif'), gdalwarp_index=NULL,gdalwarp_params=list(t_srs='+proj=longlat +datum=WGS84',s_srs='+proj=longlat +datum=WGS84'),verbose=FALSE)
      
    }
    var= c('tmin','tmax','prec')
    for(v in 1:length(var)){
      for(i in 1:12){ #Loop on number of bio
        
        #Search all raster files of biox
        files=list.files(pattern = paste0(var[v],i,'_[0-9]*\\.bil'))
        #merge them using library gdalUtils
        gdalUtils::mosaic_rasters(files, paste0(var[v],i,'.tif'), gdalwarp_index=NULL,gdalwarp_params=list(t_srs='+proj=longlat +datum=WGS84',s_srs='+proj=longlat +datum=WGS84'),verbose=FALSE)
        
      }
    }
    
    file_vector=paste0('bio',1:19,'.tif')
    for(v in 1:length(var)){
      file_vector=c(file_vector,paste0(var[v],1:12,'.tif'))
    }
    file_vector_short=paste0('bio',1:19)
    for(v in 1:length(var)){
      file_vector_short=c(file_vector_short,paste0(var[v],1:12))
    }
    
    #Load raster while stacking them
    raster=raster::stack(file_vector)
    raster::projection(raster)=sp::CRS("+init=epsg:4326")
    if(interactiveChecks==TRUE){
      raster::plot(raster, y=1,main='bio1')
      raster::plot(locations, add=TRUE)
      biovar=readline("Would you like to continue? Press x to stop the process, any other letter to continue: ")
      if(biovar=='x'){
        print('Function ended on user input')
        return(NA)
      }
    }
    
    #Extract WC data at point location
    colname_data=colnames(data)
    extractdatawc=extractVar(raster, locations)
    data = data.frame(data, extractdatawc)
    colnames(data)=c(colname_data,file_vector_short)
  }
  
  ### SRTM
  
  if(srtm==TRUE){
    
    if(verbose==TRUE){
      print('Downloading SRTM data')
    }
    
    dir.create(file.path(active_dir, 'srtm'), showWarnings = FALSE)
    setwd(file.path(active_dir, 'srtm'))
    
    
    #Transform to WGS84
    locations2=sp::spTransform(locations, sp::CRS("+init=epsg:4326"))
    coord=sp::coordinates(locations2)
    
    #To find the name of the tile (getData does it automatically but bugs to check if the file is already downloaded)
    rs = raster::raster(nrows=24, ncols=72, xmn=-180, xmx=180, ymn=-60, ymx=60 )
    rowTiles=c()
    colTiles=c()
    
    for(i in 1:nrow(locations)){
      latt=coord[i,2]
      long=coord[i,1]
      
      #Get row and column tile and store them
      
      rowTile = raster::rowFromY(rs, latt)
      if(rowTile<10){ #If <10, add 0 in file name
        rowTile=paste0('0',rowTile)
      }
      rowTiles=c(rowTiles, rowTile)
      colTile = raster::colFromX(rs, long)
      if(colTile<10){
        colTile=paste0('0',colTile)
      }      
      colTiles=c(colTiles, colTile)
      
      #Download srtm data 
      
      if(!file.exists(paste0('srtm_',colTile,'_',rowTile,'.tif'))){ #Check if file already downloaded. Should be done by getData but bugs
        #Download srtm data . Stored in active directory
        tryCatch({raster=raster::getData("SRTM",path=active_dir, lon=long, lat=latt)}, error=function(e){})        
      }

    }
    
    # Merge all files
    files=unique(paste0('srtm_',colTiles,'_',rowTiles,'.tif'))
    gdalUtils::mosaic_rasters(files,'srtm.tif', verbose=FALSE)
    
    #Load raster
    raster=raster::raster('srtm.tif')
    
    #Assign projection
    raster::projection(raster)=sp::CRS("+init=epsg:4326")
    
    #If interactive check, plot data
    if(interactiveChecks==TRUE){
      raster::plot(raster,main='srtm')
      raster::plot(locations, add=TRUE)
      srtm=readline("Would you like to continue?. Press x to stop the process, any other letter to continue: ")
      if(srtm=='x'){
        print('Function ended on user input')
        return(NA)
      }
    }
    
    #Extractdata
    colname_data=colnames(data)
    extractdatasrtm=extractVar(raster, locations)
    data = data.frame(data, extractdatasrtm)
    colnames(data)=c(colname_data,'srtm')
  }
  
  #Export data
  setwd(working_dir)
  locationfilename_short=substr(locationFileName,1,gregexpr("\\.", locationFileName)[[1]][length(gregexpr("\\.", locationFileName)[[1]])]-1)
  
  write.table(data, file=outputFile, append=FALSE,quote=TRUE,sep=" ", dec = ".",row.names=FALSE,col.names=TRUE)
  
  if(verbose==TRUE){
    print(paste0(outputFile,' sucessfully created!'))
  }
  
}


#' @title Prepare environmental input
#' @description Writes a new environmental file that sambada can work with after having removed too correlated variables. Also calculates population structure from a PCA in SNPRelate and add it at the end of the environmental file
#' @author Solange Duruz, Oliver Selmoni
#' @param envFile char Name of the input environmental file (must be in active directory). Can be .csv or .shp
#' @param outputFile char Name of the output file. Must have a .csv extension.
#' @param maxCorr double A number between 0 and 1 specifying the maximum allowable correlation coefficient between environmental files. If above, one of the variables will be deleted
#' @param idName char Name of the id in the environmental file matching the one of \code{genoFile}
#' @param separator char If envFile is .csv, the separator character. If file created with create_env, separator is ' '
#' @param genoFile char (optional) Name of the input genomic file (must be in active directory). If not null, population variable will be calculated from a PCA relying on the SNPRelate package. Can be .gds, .ped, .bed, .vcf. If different from .gds, a gds file (SNPrelate specific format) will be created
#' @param numPc double If above 1, number of principal components to analyze. If between 0 and 1, automatic detection of number of PC (the program will find the first leap in the proportion of variance where the ratio (difference in variance between PC x and x+1)/(variance of PC x) is greater than NumPc. If 0, PCA and population structure will not be computed: in that case, the \code{genoFile} will only be used to make the sample order in the envFile match the one of the \code{envFile} (necessary for sambada's computation). Set it to null if \code{genoFile} is null 
#' @param mafThresh double A number between 0 and 1 specifying the Major Allele Frequency (MAF) filtering when computing PCA (if null no filtering on MAF will be computed)
#' @param missingnessThresh double A number between 0 and 1 specifying the missing rate filtering when computing PCS(if null no filtering on missing rate will be computed)
#' @param ldThresh double A number between 0 and 1 specifying the linkage disequilibrium (LD) rate filtering before computing the PCA (if null no filtering on LD will be computed)
#' @param numPop integer If not null, clustering based on \code{numPc} first PC will be computed to divide into \code{numPop} populations. If -1 automatic detection of number of cluster (elbow method if \code{clustMethod}='kmeans', maximise branch length if \code{clustMethod}='hclust'). If null, no clustering will be computed: if \code{genoFile} is set, principal component scores will be included as population information in the final file.
#' @param clustMethod char One of 'kmeans' or 'hclust' for K-means and hierarchical clustering respectively. Default 'kmeans'
#' @param interactiveChecks logical If TRUE, plots will show up showing number of populations chosen, and correlation between variables and the user can interactively change the chosen threshold for \code{maxCorr} and \code{numPop} (optional, default value=FALSE)
#' @param includeCol character vector Columns in the environmental file to be considered as variables. If none specified, all numeric variables will be considered as env var except for the id
#' @param excludeCol character vector Columns in the environmental file to exclude in the output (non-variable column). If none specified, all numeric variables will be considered as env var except for the id
#' @param popStrCol character vector Columns in the environmental file describing population structure (ran elsewhere). Those columns won't be excluded when correlated with environmental files
#' @param x character Name of the column corresponding to the x coordinate (or longitude if spherical coordinate). If not null, x column won't be removed even if correlated with other variable. This parameter is also used to display the map of the population structure.
#' @param y character Name of the column corresponding to the y coordinate (or latitude if spherical coordinate). If not null, y column won't be removed even if correlated with other variable. This parameter is also used to display the map of the population structure.
#' @param locationProj integer EPSG code of the projection of x-y coordinate
#' @param verbose boolean If true show information about progress of the process
#' @return None
#' @examples
#' \dontrun{
#' #Calculating PCA-based population structure
#' prepareEnv('myFile-env.csv','myFile-env-export.csv',0.8,'Nom',' ','myFile.gds', 
#'      numPc=0.2, mafThresh=0.05, missingnessThresh=0.1, ldThresh=0.2, numPop=NULL,
#'      x='Longitude', y='Latitude', locationProj=4326, interactiveChecks = TRUE)
#'
#' #Calculating structure membership coefficient based on kmeans clustering
#' prepareEnv('myFile-env.csv','myFile-env-export.csv',0.8,'Nom',' ','myFile.gds', 
#'      numPc=0.2, mafThresh=0.05, missingnessThresh=0.1, ldThresh=0.2, numPop=NULL,
#'      x='Longitude', y='Latitude', locationProj=4326, interactiveChecks = TRUE)
#'
#' #Without calculating population structure.
#' prepareEnv('myFile-env.csv','myFile-env-export.csv',0.8,'Nom',' ', 
#'      x='Longitude',y='Latitude', locationProj=4326, interactiveChecks = TRUE)
#' }
#' @export
prepareEnv=function(envFile, outputFile, maxCorr, idName, separator=' ',genoFile=NULL, numPc=0.5, mafThresh=NULL, missingnessThresh=NULL, ldThresh=NULL, numPop=-1, clustMethod='kmeans', includeCol=NULL, excludeCol=NULL, popStrCol=NULL, x,y,locationProj,interactiveChecks=FALSE, verbose=TRUE){

  ### Check inputs ###
  
  if(is.null(envFile)) stop("envFile argument is required")
  if(typeof(envFile)!='character') stop("envFile argument supposed to be a character string")
  if (!file.exists(envFile)) stop("Input envFile not found.")
  
  if(typeof(outputFile)!='character') stop("outputFile argument supposed to be character")
  extensionO=substr(outputFile,gregexpr("\\.", outputFile)[[1]][length(gregexpr("\\.", outputFile)[[1]])]+1, nchar(outputFile))
  if(extensionO!='csv') stop("outputFile must have a .csv extension")
  
  if(!is.null(maxCorr)){
    if(typeof(maxCorr)!='double') stop("maxCorr argument supposed to be decimal number")
    if(maxCorr>1 | maxCorr<0) stop("maxCorr argument supposed to be between 0 and 1")
  }
  
  if(!is.null(idName)){
    if(typeof(idName)!='character') stop("idName argument supposed to be a character string")
  }
  
  if(!is.null(genoFile)){
    if(typeof(genoFile)!='character') stop("genoFile argument supposed to be a character string")
    if(!file.exists(genoFile)) stop("Input genoFile not found.")
    if(is.null(numPc)) stop("numPc cannot be null if genoFile is not null")
  }
  
  if(!is.null(numPc)){ 
    if(typeof(numPc)!='double')stop("numPc must be a number")
    if(numPc>1 & numPc%%1!=0) stop("If numPc>1, numPc must be an integer")
    if(is.null(genoFile))stop("genoFile must be specified if numPc is not null")
  }
  
  if(!is.null(numPop)){ 
    if(numPop<0 & numPop!=-1) stop("numPop must be positive or equal to -1")
    if(typeof(numPop)!='double') stop("If numPop>1, numPop must be an integer")
    if(numPop>1 & numPop%%1!=0) stop("numPop argument supposed to be an integer")
  }
  
  if(!is.null(mafThresh)){
    if(typeof(mafThresh)!='double') stop("mafThresh argument supposed to be decimal number")
    if((mafThresh>1 | mafThresh<0)) stop("mafThresh argument supposed to be between 0 and 1")
  }
  
  if(!is.null(missingnessThresh)){
    if(typeof(missingnessThresh)!='double') stop("missingnessThresh argument supposed to be decimal number")
    if((missingnessThresh>1 | missingnessThresh<0)) stop("missingnessThresh argument supposed to be between 0 and 1")
  }
  
  if(!is.null(includeCol)){
    if(typeof(includeCol)!='character') stop("includeCol argument supposed to be a character vector")
  }
  if(!is.null(excludeCol)){
    if(typeof(excludeCol)!='character') stop("excludeCol argument supposed to be a character vector")
  }
  
  if(!is.null(x)){
    if(typeof(x)!='character') stop('x supposed to be a string')
    if(is.null(y)) stop('y must be provided if x is provided')
    if(typeof(y)!='character') stop('y supposed to be a string')
    if(is.null(locationProj)) stop('locationProj must be specified when x and y are provided')
    if(typeof(locationProj)!='double')stop('locationProj supposed to be an integer (EPSG code)')
    if(locationProj%%1!=0) stop('locationProj supposed to be an integer (EPSG code)')
  }
  
  if(!is.null(interactiveChecks)){
    if(typeof(interactiveChecks)!='logical') stop('interactiveChecks argument supposed to be logical')
  }
  
  if(!is.null(verbose)){
    if(typeof(verbose)!='logical') stop('verbose argument supposed to be logical')
  }
  
  ### Load required libraries ###
  if(!is.null(genoFile)){
    if (!requireNamespace("SNPRelate", quietly = TRUE)) {
      stop("Package \"SNPRelate\" needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("gdsfmt", quietly = TRUE)) {
      stop("Package \"gdsfmt\" needed for this function to work. Please install it.", call. = FALSE)
    }
  }
  if(interactiveChecks==TRUE & !is.null(x)){
    if (!requireNamespace("sp", quietly = TRUE)) {
      stop("If interactiveChecks is TRUE and x/y is specified, package \"sp\" is needed. Please install it.", call. = FALSE)
    }  
    if (!requireNamespace("rgdal", quietly = TRUE)) {
      stop("If interactiveChecks is TRUE and x/y is specified, package \"rgdal\" is needed. Please install it.", call. = FALSE)
    }  
    if (!requireNamespace("rworldmap", quietly = TRUE)) {
      stop("If interactiveChecks is TRUE and x/y is specified, package \"rworldmap\" is needed. Please install it.", call. = FALSE)
    } 
    if (!requireNamespace("mapplots", quietly = TRUE)) {
      stop("If interactiveChecks is TRUE and x/y is specified, package \"mapplots\" is needed. Please install it.", call. = FALSE)
    }   
  }
  
  
  ### Open Environmental file ###

  if(verbose==TRUE){
    print('Opening envFile')
  }
  
  env=openEnvData(envFile, separator)
  
  #Check that includeCol, excludeCol, popStrCol are in the header of the envFile
  if(!is.null(includeCol)){
    if(sum(includeCol %in% colnames(env))!=length(includeCol)) stop("Not all includeCol present in envFile")
  }
  if(!is.null(excludeCol)){
    if(sum(excludeCol %in% colnames(env))!=length(excludeCol)) stop("Not all excludeCol present in envFile")
  }
  if(!is.null(x)){
    if(sum(x %in% colnames(env))!=1) stop("x column not in envFile")
    if(sum(x %in% colnames(env))!=1) stop("y column not in envFile")
  }
  if(!is.null(popStrCol)){
    if(sum(popStrCol %in% colnames(env))!=length(popStrCol)) stop("Not all popStrCol present in envFile")
  }
  
  ### Check correlation among variables ###
  
  if(verbose==TRUE){
    print('Computing correlation among variables and reducing dataset')
  }
    
  # Find numeric column
  env_colnames=c()
  
  # Only include column present in includeCol
  if(!is.null(includeCol)){
    env2=env[,includeCol]
  }
  else{
    env2=env
  }
  
  for(i in 1:ncol(env2)){
    #If column is numeric and not ID nor in excludeCol, nor popStrCol, nor x/y
    if(is.numeric(env2[,i]) & colnames(env2)[i]!=idName & !(colnames(env2)[i] %in% excludeCol) & !(colnames(env2)[i] %in% popStrCol) & colnames(env2)[i]!=x & colnames(env2)[i]!=y){
      env_colnames=c(env_colnames, colnames(env2)[i])
    }
  }
  env_numeric=env[,env_colnames]
  
  if(interactiveChecks==TRUE){
    #plot: chosen threshold - num variables kept
    kept_env=c()
    corr_thresh=c()
    for(mc in seq(0,1,by=0.01)){
      env_red=redENV(env_numeric, mc)
      kept_env=c(kept_env,length(env_red))
      corr_thresh=c(corr_thresh,mc)
    }
    plot(corr_thresh,kept_env,pch='.',xlab='chosen correlation threshold',ylab='number of variables kept',main='Num of var kept according to correlation threshold')
    lines(corr_thresh,kept_env)
    abline(v=maxCorr,col="red")
    maxCorr2 = readline(prompt="Would you like to choose a different max authorized correlation threshold? (press any letter if no, or enter a new number (between 0-1) indicating the maximum correlation): ")
    if(grepl("[[:digit:]\\.-]",maxCorr2)){
      maxCorr=as.integer(maxCorr2)
    }  
  }
  
  #Keep only variables whose correlation among them is below the maxCorr threshold (ad-hoc function redEnv)
  env_red=redENV(env_numeric, maxCorr) ###
  env_name=names(env_red)
  env_kept=env[,env_name]
  
  ### Calculate propulation structure
  
  if(!is.null(genoFile) & numPc>0){ #If genoFile is not null, calculate population structure

    if(verbose==TRUE){
      print('computing population structure')
    }
   
    #Create GDS file
    gds_file=createGDSfile(genoFile, getwd())
    
    #Open GDS file (and close it on exit)
    gds_obj=SNPRelate::snpgdsOpen(gds_file)
    on.exit(SNPRelate::snpgdsClose(gds_obj))
    
    # Run PCA from SNPRelate package

    if(!is.null(ldThresh)){
      ld_filtered=SNPRelate::snpgdsLDpruning(gds_obj, ld.threshold=ldThresh)
      numvect=min(length(unlist(ld_filtered)),100)
      pca=SNPRelate::snpgdsPCA(gds_obj, snp.id=unlist(ld_filtered),maf=mafThresh, missing.rate=missingnessThresh, eigen.cnt = numvect)
    } else {
      numvect=min(length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.id"))),100)
      pca=SNPRelate::snpgdsPCA(gds_obj,maf=mafThresh, missing.rate=missingnessThresh, eigen.cnt = numvect)
    }
    # Choose best number of PC if numPc<1
    varprop=pca$varprop[1:numvect]
    if(numPc<1){
      #DIFFPC = diff(c(0,varprop))
      #NbPCs: ad-hoc function seaking for the leap in the variance proportion
      numPc=NbPCs2(varprop,numPc)
      if(numPc==length(varprop))stop('Leap in proportion of variance of PCA not found within the first 100 principal component')
    } 
    # Plot VP of first 100 axes and show chosen number of pc
    if(interactiveChecks==TRUE){
      barplot(varprop, main='Variance proportion of first 100 axes', xlab='axis number', ylab='proportion variance explained')
      abline(v=numPc,col="red")
      numPc2 = readline(prompt=paste("Number of PCs chosen: ",numPc,". Would you like to choose a different number of pc? (press any letter if no, or enter a new number): "))
      if(grepl("[[:digit:]\\.-]",numPc2)){
        numPc=as.integer(numPc2)
      }
    }
    
    popvect_tot=cbind(pca$sample.id, pca$eigenvect)
    #ID match between envFile and genoFile
    popvect2=popvect_tot[match(env[,idName],popvect_tot[,1]),]
    pca$eigenvect=popvect2[complete.cases(popvect2),2:ncol(popvect2)]
    class(pca$eigenvect)='numeric'
    #pca$sample.id=popvect2[complete.cases(popvect2),1]
    
    #Retrieve right number of pc
    popvect2=popvect2[,1:(numPc+1)] #First column=ID, keep numPc axis 
    if(numPc==0){
      colnames(popvect2)=c('sampleid')
    } else {
      colnames(popvect2)=c('sampleid',paste0('pop',1:(numPc)))
    }
    
    if(nrow(popvect_tot)!=nrow(pca$eigenvect)){
      stop(paste0('All IDs in envFile and genoFile do not match. Found ',nrow(pca$eigenvect),' match out of ',nrow(popvect_tot),' indiv'))
    }
    

    
    if(!is.null(numPop) & numPc>0){
      changePopNum=TRUE
      while(changePopNum==TRUE){ #The user can interactively change the number of population. Loop begins again 
        
        ### Clustering ###
        
        if(numPop==-1){ #automatic detection of num pop
          
          #hierarchical clustering 
          if(clustMethod=='hclust'){
            
            hiclust=hclust(dist(as.matrix(pca$eigenvect[,1:numPc])))
            #automatic detection of optimal number of cluster (det num clust for each height and keep the most frequent)
            v=seq(0,max(hiclust$height),max(hiclust$height)/100)
            numPop_all = sapply(v, function(height){max(cutree(hiclust, h=height))})
            numPop=as.integer(names(sort(table(numPop_all),decreasing=TRUE)[1])) #most frequent
            
            #Prompt the user for a different number of clusters
            plot(hiclust)
            abline(h=v[numPop_all==numPop][1],col="red")
            numPop2=readline(paste0("Number of suggested cluster ",numPop,". Would you like to change it? (Press a new number or any letter to continue): "))
            if(grepl("[[:digit:]\\.-]",numPop2)){
              numPop=as.integer(numPop2)
            } 
            clust_cut=cutree(hiclust,k=numPop)
            clust=list(cluster=clust_cut)
          } else{ # kmeans
            
            #Kmeans and elbow method
            wss <- sapply(1:20, function(k){kmeans(pca$eigenvect[,1:numPc], k, nstart=50,iter.max = 15 )$tot.withinss})
            numPop=NbPopElbow(diff(wss))
            
            #Prompt the user for a different number of clusters
            plot(1:20, wss)
            abline(v=numPop, col="red")
            numPop2=readline(paste0("Number of suggested cluster ",numPop,". Would you like to change it? (Press a new number or any letter to continue): "))
            if(grepl("[[:digit:]\\.-]",numPop2)){
              numPop=as.integer(numPop2)
            }   
            clust=kmeans(pca$eigenvect[,1:numPc],numPop)        
          }
        } else if (numPop>1){
          if(clustMethod=='hclust'){
            hiclust=hclust(dist(as.matrix(pca$eigenvect[,1:numPc])))
            clust_cut=cutree(hiclust,k=numPop)
            clust=list(cluster=clust_cut)           
          } else {
            clust=kmeans(pca$eigenvect[,1:numPc],numPop)
          }
        }
  
        #   sil=vector(length=19)
        #   scale_mm=function(x){(x-min(x))/(max(x)-min(x))} 
        #   for(i in 2:20){
        #   clust=kmeans(pca$eigenvect[,1:2],i)
        #   ss=silhouette(clust$cluster, dist(as.matrix(pca$eigenval[1]*scale_mm(pca$eigenvect[,1]),pca$eigenval[2]*scale_mm(pca$eigenvect[,2]),pca$eigenval[3]*scale_mm(pca$eigenvect[,3]))))
        #   sil[i-1] <- mean(ss[, 3])
        #   # m = ncol(clust$centers)
        #   # n = length(clust$cluster)
        #   # k = nrow(clust$centers)
        #   # #D = clust$tot.withinss
        #   # D=sum(-clust$size * log(clust$withinss)/2 )
        #   # BIC[i] = -2*D + 2*log(n)*m*k
        # }
        # numclust=which(sil==max(sil))+1
        
        #Biplot
        plot(pca$eigenvect[,2], pca$eigenvect[,1], col=clust$cluster, xlab="eigenvector 2", ylab="eigenvector 1")  
        numPop2 = readline(prompt="Would you like to choose a different number of population? (press any letter if no, -1 to reshow the tree/elbow, or enter a new number): ")
        clustMethod2 = readline(prompt="Would you like to choose a different clustering alogrithm? (press any letter if no, kmeans or hclust to change algorithm): ")
        if(grepl("[[:digit:]\\.-]",numPop2) | clustMethod2=='hclust' | clustMethod2=='kmeans'){
          if(grepl("[[:digit:]\\.-]",numPop2)){
            numPop=as.integer(numPop2)
          }
          if(clustMethod2=='hclust' | clustMethod2=='kmeans'){
            clustMethod=clustMethod2
          }
          next
        }
        
        #Calculate distance to centroid
        for(j in 1:numPop){
          #dist_j=sqrt((pca$eigenvect[,1]-mean(pca$eigenvect[t(clust$cluster==j),1]))^2+(pca$eigenvect[,2]-mean(pca$eigenvect[t(clust$cluster==j),2]))^2)
          dist_j=0
          for(kk in 1:numPc){
            dist_j=dist_j+(pca$eigenvect[,kk]-mean(pca$eigenvect[t(clust$cluster==j),kk]))^2
          }
          dist_j=sqrt(dist_j)
          if(j==1){
            dist=dist_j
          }else{
            dist=cbind(dist, dist_j)
          }
        }
        
        #Pie chart => lock ratio x/y + better define pie size
        if(interactiveChecks==TRUE & !is.null(x)){
            #Transform env to a spatial object, set projection and reproject to mercator
            location=env
            sp::coordinates(location)=c(x,y)
            sp::proj4string(location)=paste0('+init=epsg:',locationProj)
            sp::spTransform(location, '+init=epsg:4326')
            country=data('wrld_simpl', package='maptools', envir=environment())
            plot(country,xlim=c(min(sp::coordinates(location)[,'Longitude']),max(sp::coordinates(location)[,'Longitude'])),ylim=c(min(sp::coordinates(location)[,'Latitude']),max(sp::coordinates(location)[,'Latitude'])))
            rad=abs(max(sp::coordinates(location)[,'Longitude'])-min(sp::coordinates(location)[,'Longitude']))/50
            for(i in 1:nrow(location)){
              mapplots::add.pie(z=1/dist[i,1:(numPop)],x=sp::coordinates(location[i,1])[,'Longitude'], y=sp::coordinates(location[i,1])[,'Latitude'], labels=NA, radius=rad)
            }
            numPop2 = readline(prompt="Would you like to choose a different number of population? (press any letter if no, -1 to reshow the tree/elbow, or enter a new number): ")
            clustMethod2 = readline(prompt="Would you like to choose a different clustering alogrithm? (press any letter if no, kmeans or hclust to change algorithm): ")
            if(grepl("[[:digit:]\\.-]",numPop2) | clustMethod2=='hclust' | clustMethod2=='kmeans'){
              if(grepl("[[:digit:]\\.-]",numPop2)){
                numPop=as.integer(numPop2)
              }
              if(clustMethod2=='hclust' | clustMethod2=='kmeans'){
                clustMethod=clustMethod2
              }
              next
            } else {
              changePopNum=FALSE
            }
        } else {
          changePopNum=FALSE
        }
              
      }
    }
  }
  
  ### Check correlation between kept env variables and population var ###
  
  if((!is.null(genoFile) & numPc>0) | !is.null(popStrCol)){
    
    if(verbose==TRUE){
      print('Checking correlation between kept env variables and population var. If correlation > 70%, the variable will be printed here')
    }    
    if(!is.null(popStrCol)){
      #If population structure already in env file
      numvar=(length(popStrCol))
      popvect2=as.matrix(env[,popStrCol])
      if(numvar==1){ #looses column name
        colnames(popvect2)=popStrCol
      }
    } else if(!is.null(numPop)){
      popvect2=as.matrix(dist[,1:(numPop-1)])
      colnames(popvect2)=c(paste0('pop',1:(numPop-1)))
      numvar=numPop-1
    } else {
      popvect2=as.matrix(popvect2[,2:ncol(popvect2)])
      numvar=numPc
      if(numvar==1){ #looses column name
        colnames(popvect2)='pop1'
      }
    }
    
    for(i in 1:(numvar)){
      for(j in 1:ncol(env_kept)){
        
        corr_value=cor(env_kept[,j],as.numeric(popvect2[,i]), use='complete.obs')*100
        if(abs(corr_value)>70){
          # Inform user if population and env var are more than 70% correlated
          print(paste0('!! Variable ',colnames(env_kept)[j],' is correlated at ',round(corr_value,1),' percent with population ',colnames(popvect2)[i],' !!'))
        }
        
      }
    }
  }
  
  ### bind kept env variables and population and print file ###
  
  if(!is.null(x)){
    if(!is.null(genoFile) & numPc>0){
      total=cbind(env[,c(idName, x, y)],env_kept,popvect2)
      colnames(total)=c(idName,x,y,colnames(env_kept),colnames(popvect2)[1:numvar])
      #ID match between envFile and genoFile
      total=total[match(pca$sample.id,total[,idName]),]
    } else if(!is.null(popStrCol)){
      total=cbind(env[,c(idName,x,y)],env_kept,env[,popStrCol])
      colnames(total)=c(idName,x,y,colnames(env_kept),popStrCol) 
    } else {
      total = cbind(env[,c(idName,x,y)],env_kept)
      colnames(total)=c(idName,x,y,colnames(env_kept))
    }
  } else { #not geographic column
    if(!is.null(genoFile) & numPc>0){
      total=cbind(env[,idName],env_kept,popvect2)
      colnames(total)=c(idName,colnames(env_kept),colnames(popvect2)[1:numvar])
      #ID match between envFile and genoFile
      total=total[match(pca$sample.id,total[,idName]),]
    } else if(!is.null(popStrCol)){
      total=cbind(env[,idName],env_kept,env[,popStrCol])
      colnames(total)=c(idName,colnames(env_kept),popStrCol)   
    } else {
      total = cbind(env[,idName],env_kept)
      colnames(total)=c(idName,colnames(env_kept))
    }    
  }
  envFilename_short=substr(envFile,1,gregexpr("\\.", envFile)[[1]][length(gregexpr("\\.", envFile)[[1]])]-1)
  
  write.table(total, file=outputFile, append=FALSE,quote=FALSE,sep=" ", dec = ".",row.names=FALSE,col.names=TRUE)
  if(verbose==TRUE){
    print(paste0('File ',outputFile,' successfully created'))
  }
}

#Open environmental data
openEnvData = function(envFile, separator){
  if(typeof(envFile)!='character') stop("envFile argument supposed to be decimal number")
  if (!file.exists(envFile)) stop("Input envFile not found.")
  envextension=substr(envFile,gregexpr("\\.", envFile)[[1]][length(gregexpr("\\.", envFile)[[1]])]+1, nchar(envFile))
  if(envextension=='csv'){
    env=read.csv(envFile, header=TRUE, sep=separator)
  } else if(envextension=='shp'){
    env=raster::shapefile(envFile)
  } else {
    stop("Extension not supported! Should be either .csv or .shp")
  } 
  return(env)
}

#Extract raster var from location file
extractVar = function(raster, locations){
  if(raster::projection(locations)!=raster::projection(raster)){ #marche pas forcement car peut avoir init=epsg en plus
    locations2=sp::spTransform(locations, sp::CRS(raster::projection(raster)))
  } 
  #Extract value of wordclim at point location
  extractdata=raster::extract(raster, locations)
  return(extractdata)
}





#To choose the optimal number of populations
NbPCs = function(DIFFPC, cutoff=0.1) {
  
  co=1
  while (abs((DIFFPC[co]-DIFFPC[co+1])/DIFFPC[co])>cutoff) {
    co=co+1
    if(co==length(DIFFPC)){
      break()
    }
  }
  
  print(paste0('Number of PC suggested for describing pop structure: ',co-1))
  return(co-1)
  
}

#To choose the optimal number of populations
NbPCs2 = function(PC, cutoff=0.5) {
  
  co=1
  while (abs((PC[co+1]-PC[co])/PC[co+1])<cutoff) {
    co=co+1
    if(co==length(PC)){
      return(0)
    }
  }
  
  #print(paste0('Number of PC suggested for describing pop structure: ',co-1))
  return(co)
  
}

NbPopElbow = function(DIFFPC, cutoff=0.3) {

  co=length(DIFFPC)-1
  while ((DIFFPC[co+1]/DIFFPC[co])>cutoff) {
    co=co-1
    if(co==1){
      break()
    }
  }
  
  print(paste0('Number of clusters to describe your populations: ',co+1))
  return(co+1)
  
}

#Reduce environmental dataset
redENV = function(ienv, corcutoff=0.7) {
  
  listvar = colnames(ienv)
  
  co=1
  
  olist=list()
  
  while(length(listvar)>0) {
    
    C = cor(ienv)[listvar[co],listvar[-co]]  
    group=names(C[C>corcutoff]) 
    
    olist[listvar[co]]=paste(group, collapse = ', ')
    
    listvar=listvar[-c(1,which(listvar%in%group))]
    
  }
  
  return(olist)
  print(paste0('Number of variables kept: ',length(olist)))
  
}

createGDSfile=function(filename, outputDir){
  filename_short=substr(filename,1,gregexpr("\\.", filename)[[1]][length(gregexpr("\\.", filename)[[1]])]-1)
  extension=substr(filename,gregexpr("\\.", filename)[[1]][length(gregexpr("\\.", filename)[[1]])]+1, nchar(filename))
  gds_file=file.path(outputDir,paste0(filename_short,'.gds'))
  if (extension=='gds') {
    gds_file=filename
    
  } else if (extension=='ped') {
    
    if (!file.exists(paste(filename_short,'.map',sep=''))) stop(".map input file not found. Same name as .ped mandatory")
    SNPRelate::snpgdsPED2GDS(filename,map.fn=paste(filename_short,'.map',sep=''),gds_file, verbose=FALSE)
    
  } else if (extension=='bed') {
    
    if (!file.exists(paste(filename_short,'.bim',sep=''))) stop(".bim input file not found. Same name as .bed mandatory")
    if (!file.exists(paste(filename_short,'.fam',sep=''))) stop(".fam input file not found. Same name as .bed mandatory")
    SNPRelate::snpgdsBED2GDS(bed.fn=filename,fam.fn=paste(filename_short,'.fam',sep=''), bim.fn=paste(filename_short,'.bim',sep=''),gds_file, verbose=FALSE)
  } else if (extension=='vcf') {
    SNPRelate::snpgdsVCF2GDS(vcf.fn=filename,out.fn=gds_file, verbose=FALSE)
  } else {
    stop("File extension not supported")
  }
  return(gds_file)
}

