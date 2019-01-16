#' @title Prepare genomic input
#' @description Writes a new genomic file that sambada can work with after having applied the selected genomic filtering options. The output file has the same name as the input file but with a .csv extension
#' @author Solange Gaillard, Oliver Selmoni
#' @param filename char Name of the input file (must be in active directory). Can be .gds, .ped, .bed, .vcf. If different from .gds, a gds file (SNPrelate specific format) will be created unless no filtering options are chosen
#' @param maf.thresh double A number between 0 and 1 specifying the Major Allele Frequency (MAF) filtering (if null no filtering on MAF will be computed)
#' @param missingness.thresh double A number between 0 and 1 specifying the missing rate filtering (if null no filtering on missing rate will be computed)
#' @param LD.thresh double A number between 0 and 1 specifying the linkage desiquilibrium (LD) rate filtering (if null no filtering on LD will be computed)
#' @param mgf.thresh double A number between 0 and 1 specifying the Major Genotype Frequency (MGF) rate filtering (if null no filtering on MGF will be computed). NB: sambada computations relie on gentoypes
#' @param directory char The directory where binaries of sambada are saved. This parameter is not necessary if directoy path is permanently stored in the PATH environmental variable or if a function invoking samabada executable (prepareGeno or sambada_parallel) has been already run in the R active session.
#' @param interactiveChecks logical If TRUE, plots will show up showing distribution of allele frequency etc... and the user can interactively change the chosen threshold for maf.thresh, missingness.thresh, mgf.thresh (optional, default value=FALSE)
#' @return None
#' @examples
#' prepareGeno('myPlinkFile.ped',maf.tresh=0.05, missingness.thresh=0.05,LD.thresh=0.2,mgf.thresh=0.8,interactiveChecks=TRUE)
#' prepareGeno('myGDSFile.gds',maf.tresh=0.05, missingness.thresh=0.05,LD.thresh=0.2,mgf.thresh=0.8,interactiveChecks=FALSE)
#' @export
prepareGeno=function(filename,maf.thresh=NULL, missingness.thresh=NULL,LD.thresh=NULL,mgf.thresh=NULL, directory=NULL, interactiveChecks=FALSE, verbose=FALSE){

  ### Check inputs ###
  
  if(typeof(filename)!='character') stop("filename argument supposed to be character")
  if (!file.exists(filename)) stop("Input file not found.")
  if(!is.null(maf.thresh)){
    if(typeof(maf.thresh)!='double') stop("maf.thresh argument supposed to be decimal number")
    if(maf.thresh>1 | maf.thresh<0) stop("maf.thresh argument supposed to be between 0 and 1")    
  }
  if(!is.null(missingness.thresh)){
    if(typeof(missingness.thresh)!='double') stop("missingness.thresh argument supposed to be decimal number")
    if(missingness.thresh>1 | missingness.thresh<0) stop("missingness.thresh argument supposed to be between 0 and 1")
  }

  if(!is.null(mgf.thresh)){
    if(typeof(mgf.thresh)!='double') stop("mgf.thresh argument supposed to be decimal number")
    if(mgf.thresh>1 | mgf.thresh<0) stop("mgf.thresh argument supposed to be between 0 and 1")
  }
  if(!is.null(LD.thresh)){
    if(typeof(LD.thresh)!='double') stop("LD.thresh argument supposed to be decimal number")
    if(LD.thresh>1 | LD.thresh<0) stop("LD.thresh argument supposed to be between 0 and 1")
  }
  if(typeof(interactiveChecks)!='logical') stop('interactiveChecks argument supposed to be logical')
  if(typeof(verbose)!='logical') stop('verbose argument supposed to be logical')
  
  # Get file extension
  filename_short=substr(filename,1,gregexpr("\\.", filename)[[1]][length(gregexpr("\\.", filename)[[1]])]-1)
  extension=substr(filename,gregexpr("\\.", filename)[[1]][length(gregexpr("\\.", filename)[[1]])]+1, nchar(filename))
  
  if(!is.null(directory)){
    changePath(directory)
  }

  ### If ped file with no filters => no need to code to GDS ###
  
  if(extension=='ped' & is.null(maf.thresh) & is.null(missingness.thresh) & is.null(mgf.thresh) & is.null(LD.thresh)){
    
    if (!file.exists(paste(filename_short,'.map',sep=''))) stop(".map input file not found. Same name as .ped mandatory")
    numSNP=nrow(read.table(paste(filename_short,'.map',sep=''), colClasses=c("character", rep("NULL",3)),sep="\t"))
    #numIndiv=nrow(read.table(filename, colClasses=c(rep("character",2), rep("NULL",4+numSNP*2)),sep=" "))
    #Calculate numIndiv by counting number of new line in pef file
    f = file(filename, open="rb")
    numIndiv = 0L
    while (length(chunk <- readBin(f, "raw", 20000)) > 0) {
      numIndiv = numIndiv + sum(chunk == as.raw(10L))
    }
    close(f)
    
    #recodePlink(numIndiv,numSNP,filename_short,paste(filename_short,'.csv'))
    system(paste('recode-plink',numIndiv,numSNP,filename_short,paste0(filename_short,'.csv')))
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
  
  if(verbose==TRUE){
    print('Creating and opening GDS file')
    print('===================================')
  }
  #Creates GDS file
  if(interactiveChecks==TRUE & extension!='gds' & file.exists(paste0(filename_short,'.gds'))){
    useGDS = readline(prompt="A .gds file has been found with the same name. Would you like to use it (press y) or rewrite a new file (press any other key): ")
    if(useGDS=='y'){
      gds_file=paste0(filename_short,'.gds')
    } else{
      gds_file=createGDSfile(filename)
    }
  } else{
    gds_file=createGDSfile(filename)
  }
  
  #Open gds file
  gds_obj=SNPRelate::snpgdsOpen(gds_file)
  on.exit(SNPRelate::snpgdsClose(gds_obj))
  
  ### Filtering ###
  
  if(verbose==TRUE){
    print('Filtering using SNPRelate in process')
    print('===================================')  
  }
  # Interactive histograms of filtering
  if(interactiveChecks==TRUE){ #If true, shows histogram of pruning values and asks the user to check chosen thresholds
    snpSummary=SNPRelate::snpgdsSNPRateFreq(gds_obj, with.snp.id = TRUE)
    # Creates MAF histograms
    if(is.null(maf.thresh)==FALSE){
      hist(snpSummary$MinorFreq, breaks=100, xlab='Minor allele Frequency', main='Histogram of minor allele frequency') 
      abline(v=maf.thresh,col="red")
      maf.thresh2 = readline(prompt="Would you like to change your maf threshold? (press n if no, or enter a new threshold): ")
      if(grepl("[[:digit:]\\.-]",maf.thresh2)){
        if(as.numeric(maf.thresh2)>1 | as.numeric(maf.thresh2)<0) stop("maf.thresh argument supposed to be between 0 and 1")
        maf.thresh=as.numeric(maf.thresh2)
      }      
    }
    # Creates Missingness histograms
    if(is.null(missingness.thresh)==FALSE){
      hist(snpSummary$MissingRate, breaks=100, xlab='Missingness',main='Histogram of missingness') 
      abline(v=missingness.thresh,col="red")
      missingness.thresh2 = readline(prompt="Would you like to change your missingness threshold? (press n if no, or enter a new threshold): ")
      if(grepl("[[:digit:]\\.-]",missingness.thresh2)){
        if(as.numeric(missingness.thresh2)>1 | as.numeric(missingness.thresh2)<0) stop("missingness.thresh argument supposed to be between 0 and 1")
        missingness.thresh=as.numeric(missingness.thresh2)
      }
    }
    # Number of pruned SNP by LD
    if(is.null(LD.thresh)==FALSE){
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
        abline(v=LD.thresh*length(seq(0,1,0.05))+1-LD.thresh, col="red")
        LD.thresh2 = readline(prompt="Would you like to change your LD threshold? (press n if no, or enter a new threshold): ")
        if(grepl("[[:digit:]\\.-]",LD.thresh2)){
          if(as.numeric(LD.thresh2)>1 | as.numeric(LD.thresh2)<0) stop("LD.thresh argument supposed to be between 0 and 1")
          LD.thresh=as.numeric(LD.thresh2)
        }
      }
    }
    # Major Genotype Frequency
    if(is.null(mgf.thresh)==FALSE){
      geno=SNPRelate::snpgdsGetGeno(gds_obj)
      mgf_freq=MGF(geno,-1)
      hist(mgf_freq/100, breaks=100, main='Histogram of MGF')
      abline(v=mgf.thresh, col='red')
      mgf.thresh2 = readline(prompt="Would you like to change your MGF threshold? (press n if no, or enter a new threshold): ")
      if(grepl("[[:digit:]\\.-]",mgf.thresh2)){
        if(as.numeric(mgf.thresh2)>1 | as.numeric(mgf.thresh2)<0) stop("mgf.thresh argument supposed to be between 0 and 1")
        mgf.thresh=as.numeric(mgf.thresh2)
      }
    }
  }
  
  #LD, MAF and Missing rate pruning
  if(is.null(maf.thresh)){
    maf.thresh=NaN
  } 
  if(is.null(missingness.thresh)){
    missingness.thresh=NaN
  } 
  if(is.null(LD.thresh)){ #Filtrer à la main or use snpgdsSelectSNP
    #if(!exists('snpSummary')){
    #  snpSummary=SNPRelate::snpgdsSNPRateFreq(gds_obj, with.snp.id=TRUE)
    #}
    #maf_filtered=snpSummary$snp.id[snpSummary$MinorFreq>maf.thresh] #List of snps with maf > than threshold (pass the filter)
    #miss_filtered=snpSummary$snp.id[snpSummary$MissingRate<missingness.thresh] #List of snps with missingness < than threshold (pass the filter)
    #snp_filtered=maf_filtered[maf_filtered %in% miss_filtered] #List of snps that pass both tests
    snp_filtered=SNPRelate::snpgdsSelectSNP(gds_obj, autosome.only=FALSE, maf=maf.thresh, missing.rate=missingness.thresh, verbose=FALSE)
  } else {
    gds_pruned=SNPRelate::snpgdsLDpruning(gds_obj, maf=maf.thresh, missing.rate=missingness.thresh, ld.threshold=LD.thresh, verbose=FALSE)
    snp_filtered=unlist(gds_pruned)
  }
  #Filter out insertion
  list_snptype=SNPRelate::snpgdsSNPList(gds_obj)
  snp_filtered2=list_snptype$snp.id[!(list_snptype$allele %in% c('A/C','C/A','A/G','G/A','A/T','T/A','C/G','G/C','C/T','T/C','G/T','T/G','A/0','C/0','G/0','T/0'))]
  snp_filtered=snp_filtered[!(snp_filtered %in% snp_filtered2)]
  
  
  #Major genotype frequency pruning
  if(!is.null(mgf.thresh)){
    if(interactiveChecks==FALSE){ #If true, geno already calculated
      geno=SNPRelate::snpgdsGetGeno(gds_obj)
    }
    snpMGF=MGF(geno,mgf.thresh)
    snpMGF=snpMGF[snpMGF>0] #End full of zero
    #snp present in general and MGF filters
    list_snp=snpMGF[snpMGF %in% snp_filtered] 
  }
  else{
    list_snp=snp_filtered
  }
  
  if(interactiveChecks==TRUE & file.exists(paste0(filename_short,'.csv'))){
    cont=readline(prompt=paste0(filename_short,".csv already exists and will be replaced. Do you want to continue? (press x to exit, any other key to continue): "))
    if(cont=='x'){
      print('Function ended on user input')
      return(NA)
    }
  } 
  
  #Write ped file 
  SNPRelate::snpgdsGDS2PED(gds_obj, ped.fn=paste(filename_short,'_filtered',sep=''), snp.id=list_snp, verbose=FALSE)

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
  system(paste('recode-plink',numIndiv,numSNP,paste(filename_short,'_filtered',sep=''),paste0(filename_short,'.csv')))
}


#' @title Create env file from raster file(s) and/or glabal database present in the raster r package
#' @description Create env file as an input for SamBada (it is recommended to run prepare_env function before running samBada) raster file(s) and/or glabal database present in the raster r package
#' @author Solange Gaillard
#' @param locationfilename !! char Name of the file containing location of individuals. Must be in the active directory. Supported extension are .csv, .shp. All columns present in this file will also be present in the output file
#' @param x char Name of the x (or longitude if not projected coordinate system) column in the \code{locationfilename}. Required if \code{locationfilename} extension is .csv
#' @param x char Name of the y (or latitude if not projected coordinate system) column in the \code{locationfilename}. Required if \code{locationfilename} extension is .csv
#' @param locationproj integer Coordinate system EPSG code of the \code{locationfilename}. If \code{locationfilename} is already georeferenced, this argument will be skipped. Required if \code{locationfilename} extension is csv.
#' @param rastername char or list Name or list of name of raster files to import. Supported format are the one of raster package. If \code{directory} is TRUE then the path to the directory
#' @param rasterproj integer or list of integer Coordinate system EPSG code of the rasterlayer. If rasterlyer is already georeferenced, this argument will be skipped. If \code{rastername} is a list, can be either a single number if all projections are the same or a list of projection for all files if different. If \code{directory} is TRUE, can only contain one number (all projections must be equal or rasters must be georeferenced)
#' @param directory logical If true, all .tif, .gtiff, .img, .sdat, . present in \code{rastername} will be loaded
#' @param worldclim logical If TRUE worldclim bio, tmin, tmax and prec variables will be downloaded at a resolution of 0.5 minutes of degree (the finest resolution). Rely rgdal and gdalUtils R package to merge the tiles. The downloaded tiles will be stored in the (new) wc0.5 directory of the active directory
#' @param srtm logical If TRUE the SRTM (altitude) variables will be downloaded at a resolution ... Rely rgdal and gdalUtils R package to merge the tiles. The downloaded tiles will be stored in the (new) wc0.5 directory of the active directory
#' @param interactiveChecks logical If TRUE, shows loaded rasters and point locations
#' @param verbose logical If TRUE, indication on process will be shown
#' @return None 
#' @examples
#' createEnv(rastername=c('prec.tif','tmin.sdat'),locationfilename='MyFile.shp',rasterproj=c(4326,21781), worldclim=TRUE,interactiveChecks=TRUE)
#' createEnv(locationfilename='MyFile.csv',x='Longitude',y='Latitude',locationproj=4326, worldclim=TRUE,interactiveChecks=FALSE)
#' @export
createEnv=function(locationfilename, x,y,locationproj, separator=',', rastername, rasterproj,directory=FALSE, worldclim=TRUE, srtm=FALSE, interactiveChecks, verbose=TRUE){
  
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

  #library(dismo)??

  
  
  #locationfilename='ADAPTmap2-env.csv'
  #x='Longitude'
  #y='Latitude'
  #rastername='prec_1.bil'
  #locationproj=4326
  
  ### Check inputs ###
  
  if(!is.null(rastername)){
    if(typeof(rastername)!='character') stop("rastername is supposed to be character or vector of char")
  }
  if(!is.null(rasterproj)){
    if(typeof(rasterproj)!='double' & sum(rasterproj%%1)!=0) stop("rasterproj is supposed to be integer or vector of integers")
    if(length(rasterproj)>1 & length(rastername)!=length(rasterproj)) stop('The size of rasterproj is supposed to be 1 or equals size of rastername')
  }
  if(!is.null(directory)){
    if(typeof(directory)!='logical') stop("directory is supposed to be logical")
  }
  if(!is.null(locationfilename)){
    if(typeof(locationfilename)!='character') stop("locationfilename is supposed to be a character")
  }
  if(!is.null(x)){
    if(typeof(x)!='character') stop("x is supposed to be a character")
  }
  if(!is.null(y)){
    if(typeof(y)!='character') stop("y is supposed to be a character")
  }
  if(!is.null(locationproj)){
    if(typeof(locationproj)!='double' & locationproj%%1!=0) stop("locationproj is supposed to be an integer")
    if(length(locationproj)>1) stop('rasterproj is not supposed to be a vector')
  }
  
  if(!is.null(worldclim)){
    if(typeof(worldclim)!='logical') stop("worldclim is supposed to be logical")
  }
  if(!is.null(srtm)){
    if(typeof(srtm)!='logical') stop("srtm is supposed to be logical")
  }
  if(!is.null(interactiveChecks)){
    if(typeof(interactiveChecks)!='logical') stop("interactiveChecks is supposed to be logical")
  }
  if(is.null(rastername) & worldclim==FALSE & srtm==FALSE) stop('Either provide a rastername or set worldclim or srtm to true')
  if(directory==TRUE & is.null(rastername)) stop('The name of the directory must be indicated in the rastername argument if directory=TRUE')
  if(directory==FALSE & !is.null(rastername)){
    for(i in 1:length(rastername)){
      if(!file.exists(rastername[i])) stop(paste0(rastername[i],' not found'))
    }
  }
  if(!file.exists(locationfilename))stop("locationfilename not found")

  locationextension=substr(locationfilename,gregexpr("\\.", locationfilename)[[1]][length(gregexpr("\\.", locationfilename)[[1]])]+1, nchar(locationfilename))
  if(locationextension=='csv'){
    if(is.null(x))stop("x must be provided if locationfilename is a .CSV")
    if(is.null(y))stop("y must be provided if locationfilename is a .CSV")
    if(is.null(locationproj))stop("locationproj must be provided if locationfilename is a .CSV")
  }
  
  #Save active directory (will be changed)
  active_dir=getwd()
  
  
  ### Open location file and set projection ###
  
  if(verbose==TRUE){
    print('loading location file')
  }
  
  
  locations=openEnvData(locationfilename, separator = separator)
  #Create data to be exported (original columns of locationfile for now)
  data=locations
  #Transform locations into a spatial object.! If shapefile
  if(locationextension=='csv'){
    sp::coordinates(locations) = c(x,y)
  }
  
  #If location file is not georeferenced and locationproj is given, assign projection
  if(is.na(raster::projection(locations)) & !is.null(locationproj)){
    raster::projection(locations) = sp::CRS(paste0('+init=epsg:',locationproj))
  }
  
  ### Rasters provided by the user
  
  if(!is.null(rastername)){
    
    if(verbose==TRUE){
      print('loading rasters provided in rastername')
    }
    
    if(directory==TRUE){
      
      #Find all raster files in the directory
      setwd(rastername)
      files=list.files(pattern = "\\.tif$")
      files=c(files, list.files(pattern = "\\.gtiff$")) 
      files=c(files, list.files(pattern = "\\.img$")) 
      files=c(files, list.files(pattern = "\\.sdat$")) 
      files=c(files, list.files(pattern = "\\.ascii$")) 
      files=c(files, list.files(pattern = "\\.grd$")) 
      files=c(files, list.files(pattern = "\\.bil$")) 
      rastername=files
    
    }
    
    for(i in 1:length(rastername)){
      
      rastername2=rastername[i]
      rastername2_short=substr(rastername2,1,gregexpr("\\.", rastername2)[[1]][length(gregexpr("\\.", rastername2)[[1]])]-1)
      if(!file.exists(rastername2))stop(paste(rastername2,"not found"))
      #Load raster using raster library
      raster=raster::raster(rastername2)
      
      #Set projection if not georeferenced
      if(is.na(raster::projection(raster)) & is.null(rasterproj)) stop(paste0("Rasters must be georeferenced or rasterproj must be provided! Check file ",rastername2))
      if(is.na(raster::projection(raster)) & !is.null(rasterproj)){
        if(length(rasterproj)>1){ #one proj per raster
          rasterproj2=rasterproj[i]
        }
        else{ # one proj for all rasters
          rasterproj2=rasterproj
        }
        raster::projection(raster) = sp::CRS(paste0('+init=epsg:',rasterproj2))
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
      colnames(data)=c(colnames_data, rastername2_short)
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
  setwd(active_dir)
  locationfilename_short=substr(locationfilename,1,gregexpr("\\.", locationfilename)[[1]][length(gregexpr("\\.", locationfilename)[[1]])]-1)
  
  write.table(data, file=paste0(locationfilename_short,'-env.csv'), append=FALSE,quote=TRUE,sep=" ", dec = ".",row.names=FALSE,col.names=TRUE)
  
  if(verbose==TRUE){
    print(paste0(locationfilename_short,'-env.csv',' sucessfully created!'))
  }
  
}


#' @title Prepare environmental input
#' @description Writes a new environmental file that sambada can work with after having removed too correlated variables. Also calculates population structure from a PCA in SNPRelate and add it at the end of the environmental file
#' @author Solange Gaillard, Oliver Selmoni
#' @param envfile char Name of the input environmental file (must be in active directory). Can be .csv or .shp
#' @param maxcorr double A number between 0 and 1 specifying the maximum allowable correlation coefficient between environmental files. If above, one of the variables will be deleted
#' @param idname char Name of the id in the environmental file matching the one of genofile
#' @param separator char If envfile is .csv, the separator character. If file created with create_env, separator is ' '
#' @param genofile char (optional)Name of the input genomic file (must be in active directory). If not null, population variable will be calculated from a PCA relying on the SNPRelate package. Can be .gds, .ped, .bed, .vcf. If different from .gds, a gds file (SNPrelate specific format) will be created
#' @param numpc double If above 1, number of principal components to analyze. If between 0 and 1, automatic detection of number of PC. If 0, PCA and population structure will not be computed: in that case, the \code{genofile} will only be used to make the sample order in the envFile match the one of the \code{enfile} (necessary for sambada's computation). Set it to null if \code{genofile} is null 
#' @param maf.thresh double A number between 0 and 1 specifying the Major Allele Frequency (MAF) filtering when computing PCA (if null no filtering on MAF will be computed)
#' @param missingness.thresh double A number between 0 and 1 specifying the missing rate filtering when computing PCS(if null no filtering on missing rate will be computed)
#' @param LD.thresh double A number between 0 and 1 specifying the linkage desiquilibrium (LD) rate filtering before computing the PCA (if null no filtereing on LD will be computed)
#' @param numpop integer If not null, clustering based on \code{numpc} first PC will be computed to devide into \code{numpop} populations. If -1 automatic detection of number of cluster (elbow method if \code{ClustMethod}='kmeans', maximise branch length if \code{ClustMethod}='hclust'). If null, no clustering will be computed: if \code{genofile} is set, principal component scores will be included as population information in the final file.
#' @param ClustMethod char One of 'kmeans' or 'hclust' for K-means and hierarchical clustering respectively. Default 'kmeans'
#' @param interactiveChecks logical If TRUE, plots will show up showing number of populations chosen, and correlation between variables and the user can interactively change the chosen threshold for maxcorr and numpop (optional, default value=FALSE)
#' @param includeCol character vector Columns in the envrionmental file to be considered as variables. If none specified, all numeric variables will be considered as env var except for the id
#' @param excludeCol character vector Columns in the envrionmental file to exclude in the output (non-variable column). If none specified, all numeric variables will be considered as env var except for the id
#' @param popstrCol character vector Columns in the envrionmental file describing population structure (ran elsewhere). Those columns won't be excluded when correlated with environmental files
#' @param x character Name of the column corresponding to the x coordinate (or longitude if sperical coordinate). If not null, x column won't be removed even if correlated with other variable. This parameter is also used to display the map of the population structure.
#' @param y character Name of the column corresponding to the y coordinate (or latitude if sperical coordinate). If not null, y column won't be removed even if correlated with other variable. This parameter is also used to display the map of the population structure.
#' @param locationProj integer EPSG code of the projection of x-y coordinate
#' @param verbose boolean If true show information about progress of the process
#' @return None
#' @examples
#' prepareEnv('myFile-env.csv',0.8,'Nom',' ',x='Longitude',y='Latitude', locationproj=4326, interactiveChecks = TRUE)
#' @export
prepareEnv=function(envfile, maxcorr, idname, separator=',',genofile=NULL, numpc=NULL, maf.thresh=NULL, missingness.thresh=NULL, LD.thresh=NULL, numpop=-1, ClustMethod='kmeans', includeCol=NULL, excludeCol=NULL, popstrCol=NULL, x,y,locationproj,interactiveChecks=FALSE, verbose=TRUE){
  #setwd("/home/lasigadmin/R/test2/src")
  #envfile='ADAPTmap2-env.csv'
  #separator=';'
  
  ### Check inputs ###
  
  if(is.null(envfile)) stop("envfile argument is required")
  if(typeof(envfile)!='character') stop("envfile argument supposed to be a character string")
  if (!file.exists(envfile)) stop("Input envfile not found.")
  
  if(!is.null(maxcorr)){
    if(typeof(maxcorr)!='double') stop("maxcorr argument supposed to be decimal number")
    if(maxcorr>1 | maxcorr<0) stop("maxcorr argument supposed to be between 0 and 1")
  }
  
  if(!is.null(idname)){
    if(typeof(idname)!='character') stop("idname argument supposed to be a character string")
  }
  
  if(!is.null(genofile)){
    if(typeof(genofile)!='character') stop("genofile argument supposed to be a character string")
    if(!file.exists(genofile)) stop("Input genofile not found.")
    if(is.null(numpc)) stop("numpc cannot be null if genofile is not null")
  }
  
  if(!is.null(numpc)){ 
    if(typeof(numpc)!='double')stop("Numpc must be a number")
    if(numpc>1 & numpc%%1!=0) stop("If numpc>1, numpc must be an integer")
  }
  
  if(!is.null(numpop)){ 
    if(numpop<0 & numpop!=-1) stop("Numpop must be positive or equal to -1")
    if(typeof(numpop)!='double') stop("If numpop>1, numpop must be an integer")
    if(numpop>1 & numpop%%1!=0) stop("numpop argument supposed to be an integer")
  }
  
  if(!is.null(maf.thresh)){
    if(typeof(maf.thresh)!='double') stop("maf.thresh argument supposed to be decimal number")
    if((maf.thresh>1 | maf.thresh<0)) stop("maf.thresh argument supposed to be between 0 and 1")
  }
  
  if(!is.null(missingness.thresh)){
    if(typeof(missingness.thresh)!='double') stop("missingness.thresh argument supposed to be decimal number")
    if((missingness.thresh>1 | missingness.thresh<0)) stop("missingness.thresh argument supposed to be between 0 and 1")
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
    if(is.null(locationproj)) stop('locationproj must be specified when x and y are provided')
    if(typeof(locationproj)!='double')stop('locationproj supposed to be an integer (EPSG code)')
    if(locationproj%%1!=0) stop('locationproj supposed to be an integer (EPSG code)')
  }
  
  if(!is.null(interactiveChecks)){
    if(typeof(interactiveChecks)!='logical') stop('interactiveChecks argument supposed to be logical')
  }
  
  if(!is.null(verbose)){
    if(typeof(verbose)!='logical') stop('verbose argument supposed to be logical')
  }
  
  ### Load required libraries ###
  if(!is.null(genofile)){
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
    print('Opening envfile')
  }
  
  env=openEnvData(envfile, separator)
  
  #Check that includeCol, excludeCol, popstrCol are in the header of the envfile
  if(!is.null(includeCol)){
    if(sum(includeCol %in% colnames(env))!=length(includeCol)) stop("Not all includeCol present in envfile")
  }
  if(!is.null(excludeCol)){
    if(sum(excludeCol %in% colnames(env))!=length(excludeCol)) stop("Not all excludeCol present in envfile")
  }
  if(!is.null(x)){
    if(sum(x %in% colnames(env))!=1) stop("x column not in envfile")
    if(sum(x %in% colnames(env))!=1) stop("y column not in envfile")
  }
  if(!is.null(popstrCol)){
    if(sum(popstrCol %in% colnames(env))!=length(popstrCol)) stop("Not all popstrCol present in envfile")
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
    #If column is numeric and not ID nor in excludeCol, nor popstrCol, nor x/y
    if(is.numeric(env2[,i]) & colnames(env2)[i]!=idname & !(colnames(env2)[i] %in% excludeCol) & !(colnames(env2)[i] %in% popstrCol) & colnames(env2)[i]!=x & colnames(env2)[i]!=y){
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
    abline(v=maxcorr,col="red")
    maxcorr2 = readline(prompt="Would you like to choose a different max authorized correlation threshold? (press any letter if no, or enter a new number (between 0-1) indicating the maximum correlation): ")
    if(grepl("[[:digit:]\\.-]",maxcorr2)){
      maxcorr=as.integer(maxcorr2)
    }  
  }
  
  #Keep only variables whose correlation among them is below the maxcorr threshold (ad-hoc function redEnv)
  env_red=redENV(env_numeric, maxcorr) ###
  env_name=names(env_red)
  env_kept=env[,env_name]
  
  ### Calculate propulation structure
  
  if(!is.null(genofile) & numpc>0){ #If genofile is not null, calculate population structure

    if(verbose==TRUE){
      print('computing population structure')
    }
   
    #Create GDS file
    gds_file=createGDSfile(genofile)
    
    #Open GDS file (and close it on exit)
    gds_obj=SNPRelate::snpgdsOpen(gds_file)
    on.exit(SNPRelate::snpgdsClose(gds_obj))
    
    # Run PCA from SNPRelate package

    #numIndiv=length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "sample.id")))
    #if(numIndiv>100){
    #  numvect=100
    #} else {
    #  numvect=numIndiv
    #}
    if(!is.null(LD.thresh)){
      ld_filtered=SNPRelate::snpgdsLDpruning(gds_obj, ld.threshold=LD.thresh)
      numvect=min(length(unlist(ld_filtered)),100)
      pca=SNPRelate::snpgdsPCA(gds_obj, snp.id=unlist(ld_filtered),maf=maf.thresh, missing.rate=missingness.thresh, eigen.cnt = numvect)
    } else {
      numvect=min(length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.id"))),100)
      pca=SNPRelate::snpgdsPCA(gds_obj,maf=maf.thresh, missing.rate=missingness.thresh, eigen.cnt = numvect)
    }
    
    # Choose best number of PC if numpc<1
    varprop=pca$varprop[1:numvect]
    if(numpc<1){
      DIFFPC = diff(c(0,varprop))
      #NbPCs: ad-hoc function seaking for the leap in the variance proportion
      numpc=NbPCs(DIFFPC,numpc)
      if(numpc==length(varprop))stop('Leap in proportion of variance of PCA not found within the first 100 principal component')
    } 
    
    # Plot VP of first 100 axes and show chosen number of pc
    if(interactiveChecks==TRUE){
      barplot(varprop, main='Variance proportion of first 100 axes', xlab='axis number', ylab='proportion variance explained')
      abline(v=numpc,col="red")
      numpc2 = readline(prompt="Would you like to choose a different number of pc? (press any letter if no, or enter a new number): ")
      if(grepl("[[:digit:]\\.-]",numpc2)){
        numpc=as.integer(numpc2)
      }
    }
    
    popvect_tot=cbind(pca$sample.id, pca$eigenvect)
    #ID match between envfile and genofile
    popvect2=popvect_tot[match(env[,idname],popvect_tot[,1]),]
    pca$eigenvect=popvect2[complete.cases(popvect2),2:ncol(popvect2)]
    class(pca$eigenvect)='numeric'
    #pca$sample.id=popvect2[complete.cases(popvect2),1]
    
    #Retrieve right number of pc
    popvect2=popvect2[,1:(numpc+1)] #First column=ID, keep numpc axis 
    colnames(popvect2)=c('sampleid',paste0('pop',1:(numpc)))
    
    if(nrow(popvect_tot)!=nrow(pca$eigenvect)){
      stop(paste0('All IDs in envfile and genofile do not match. Found ',nrow(popvect2),' out of ',nrow(popvect),' indiv'))
    }
    

    
    if(!is.null(numpop)){
      changePopNum=TRUE
      while(changePopNum==TRUE){ #The user can interactively change the number of population. Loop begins again 
        
        ### Clustering ###
        
        if(numpop==-1){ #automatic detection of num pop
          
          #hierarchical clustering 
          if(ClustMethod=='hclust'){
            
            hiclust=hclust(dist(as.matrix(pca$eigenvect[,1:numpc])))
            #automatic detection of optimal number of cluster (det num clust for each height and keep the most frequent)
            v=seq(0,max(hiclust$height),max(hiclust$height)/100)
            numpop_all = sapply(v, function(height){max(cutree(hiclust, h=height))})
            numpop=as.integer(names(sort(table(numpop_all),decreasing=TRUE)[1])) #most frequent
            
            #Prompt the user for a different number of clusters
            plot(hiclust)
            abline(h=v[numpop_all==numpop][1],col="red")
            numpop2=readline(paste0("Number of suggested cluster ",numpop,". Would you like to change it? (Press a new number or any letter to continue): "))
            if(grepl("[[:digit:]\\.-]",numpop2)){
              numpop=as.integer(numpop2)
            } 
            clust_cut=cutree(hiclust,k=numpop)
            clust=list(cluster=clust_cut)
          } else{ # kmeans
            
            #Kmeans and elbow method
            wss <- sapply(1:20, function(k){kmeans(pca$eigenvect[,1:numpc], k, nstart=50,iter.max = 15 )$tot.withinss})
            numpop=NbPopElbow(diff(wss))
            
            #Prompt the user for a different number of clusters
            plot(1:20, wss)
            abline(v=numpop, col="red")
            numpop2=readline(paste0("Number of suggested cluster ",numpop,". Would you like to change it? (Press a new number or any letter to continue): "))
            if(grepl("[[:digit:]\\.-]",numpop2)){
              numpop=as.integer(numpop2)
            }   
            clust=kmeans(pca$eigenvect[,1:numpc],numpop)        
          }
        } else if (numpop>1){
          if(ClustMethod=='hclust'){
            hiclust=hclust(dist(as.matrix(pca$eigenvect[,1:numpc])))
            clust_cut=cutree(hiclust,k=numpop)
            clust=list(cluster=clust_cut)           
          } else {
            clust=kmeans(pca$eigenvect[,1:numpc],numpop)
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
        numpop2 = readline(prompt="Would you like to choose a different number of population? (press any letter if no, -1 to reshow the tree/elbow, or enter a new number): ")
        ClustMethod2 = readline(prompt="Would you like to choose a different clustering alogrithm? (press any letter if no, kmeans or hclust to change algorithm): ")
        if(grepl("[[:digit:]\\.-]",numpop2) | ClustMethod2=='hclust' | ClustMethod2=='kmeans'){
          if(grepl("[[:digit:]\\.-]",numpop2)){
            numpop=as.integer(numpop2)
          }
          if(ClustMethod2=='hclust' | ClustMethod2=='kmeans'){
            ClustMethod=ClustMethod2
          }
          next
        }
        
        #Calculate distance to centroid
        for(j in 1:numpop){
          dist_j=sqrt((pca$eigenvect[,1]-mean(pca$eigenvect[t(clust$cluster==j),1]))^2+(pca$eigenvect[,2]-mean(pca$eigenvect[t(clust$cluster==j),2]))^2)
          if(j==1){
            dist=dist_j
          }else{
            dist=cbind(dist, dist_j)
          }
        }
        
        #Pie chart => bloquer ratio x/y + mieux définir la taille des pies
        if(interactiveChecks==TRUE & !is.null(x)){
            #Transform env to a spatial object, set projection and reproject to mercator
            location=env
            sp::coordinates(location)=c(x,y)
            sp::proj4string(location)=paste0('+init=epsg:',locationproj)
            sp::spTransform(location, '+init=epsg:4326')
            country=data('wrld_simpl', package='maptools')
            plot(country,xlim=c(min(sp::coordinates(location)[,'Longitude']),max(sp::coordinates(location)[,'Longitude'])),ylim=c(min(sp::coordinates(location)[,'Latitude']),max(sp::coordinates(location)[,'Latitude'])))
            rad=abs(max(sp::coordinates(location)[,'Longitude'])-min(sp::coordinates(location)[,'Longitude']))/50
            for(i in 1:nrow(location)){
              mapplots::add.pie(z=1/dist[i,1:(numpop)],x=sp::coordinates(location[i,1])[,'Longitude'], y=sp::coordinates(location[i,1])[,'Latitude'], labels=NA, radius=rad)
            }
            numpop2 = readline(prompt="Would you like to choose a different number of population? (press any letter if no, -1 to reshow the tree/elbow, or enter a new number): ")
            ClustMethod2 = readline(prompt="Would you like to choose a different clustering alogrithm? (press any letter if no, kmeans or hclust to change algorithm): ")
            if(grepl("[[:digit:]\\.-]",numpop2) | ClustMethod2=='hclust' | ClustMethod2=='kmeans'){
              if(grepl("[[:digit:]\\.-]",numpop2)){
                numpop=as.integer(numpop2)
              }
              if(ClustMethod2=='hclust' | ClustMethod2=='kmeans'){
                ClustMethod=ClustMethod2
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
  
  if((!is.null(genofile) & numpc>0) | !is.null(popstrCol)){
    
    if(verbose==TRUE){
      print('Checking correlation between kept env variables and population var. If correlation > 70%, the variable will be printed here')
    }    
    if(!is.null(popstrCol)){
      #If population structure already in env file
      numvar=(length(popstrCol))
      popvect2=as.matrix(env[,popstrCol])
      if(numvar==1){ #looses column name
        colnames(popvect2)=popstrCol
      }
    } else if(!is.null(numpop)){
      popvect2=as.matrix(dist[,1:(numpop-1)])
      colnames(popvect2)=c(paste0('pop',1:(numpop-1)))
      numvar=numpop-1
    } else {
      popvect2=as.matrix(popvect2[,2:ncol(popvect2)])
      numvar=numpc
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
    if(!is.null(genofile) & numpc>0){
      total=cbind(env[,c(idname, x, y)],env_kept,popvect2)
      colnames(total)=c(idname,x,y,colnames(env_kept),colnames(popvect2)[1:numvar])
      #ID match between envfile and genofile
      total=total[match(pca$sample.id,total[,idname]),]
    } else if(!is.null(popstrCol)){
      total=cbind(env[,c(idname,x,y)],env_kept,env[,popstrCol])
      colnames(total)=c(idname,x,y,colnames(env_kept),popstrCol) 
    } else {
      total = cbind(env[,c(idname,x,y)],env_kept)
      colnames(total)=c(idname,x,y,colnames(env_kept))
    }
  } else { #not geographic column
    if(!is.null(genofile) & numpc>0){
      total=cbind(env[,idname],env_kept,popvect2)
      colnames(total)=c(idname,colnames(env_kept),colnames(popvect2)[1:numvar])
      #ID match between envfile and genofile
      total=total[match(pca$sample.id,total[,idname]),]
    } else if(!is.null(popstrCol)){
      total=cbind(env[,idname],env_kept,env[,popstrCol])
      colnames(total)=c(idname,colnames(env_kept),popstrCol)   
    } else {
      total = cbind(env[,idname],env_kept)
      colnames(total)=c(idname,colnames(env_kept))
    }    
  }
  envfilename_short=substr(envfile,1,gregexpr("\\.", envfile)[[1]][length(gregexpr("\\.", envfile)[[1]])]-1)
  
  write.table(total, file=paste0(envfilename_short,'-export.csv'), append=FALSE,quote=FALSE,sep=" ", dec = ".",row.names=FALSE,col.names=TRUE)
  if(verbose==TRUE){
    print(paste0('File ',envfilename_short,'-export.csv',' successfully created'))
  }
}

#Open environmental data
openEnvData = function(envfile, separator){
  if(typeof(envfile)!='character') stop("envfile argument supposed to be decimal number")
  if (!file.exists(envfile)) stop("Input envfile not found.")
  envextension=substr(envfile,gregexpr("\\.", envfile)[[1]][length(gregexpr("\\.", envfile)[[1]])]+1, nchar(envfile))
  if(envextension=='csv'){
    env=read.csv(envfile, header=TRUE, sep=separator)
  } else if(envextension=='shp'){
    env=raster::shapefile(envfile)
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
  while (((DIFFPC[co]-DIFFPC[co+1])/DIFFPC[co])>cutoff) {
    co=co+1
    if(co==length(DIFFPC)){
      break()
    }
  }
  
  print(paste0('Number of PC suggested for describing pop structure: ',co-1))
  return(co-1)
  
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

createGDSfile=function(filename){
  filename_short=substr(filename,1,gregexpr("\\.", filename)[[1]][length(gregexpr("\\.", filename)[[1]])]-1)
  extension=substr(filename,gregexpr("\\.", filename)[[1]][length(gregexpr("\\.", filename)[[1]])]+1, nchar(filename))
  gds_file=paste0(filename_short,'.gds')
  if (extension=='gds') {
    gds_file=filename
    
  } else if (extension=='ped') {
    
    if (!file.exists(paste(filename_short,'.map',sep=''))) stop(".map input file not found. Same name as .ped mandatory")
    gds_file=paste(filename_short,'.gds',sep='')
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

