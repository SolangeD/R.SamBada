#' @title Prepare output (useful for all postprocessing analysis)
#' @author Solange Duruz, Sylvie Stucki
#' @description Read sambada's output and prepare it by retrieving the snp position and chromosome (useful for plotting manhattan)
#' @param sambadaname char The name of the genofile without extension name given to sambada (or outputfile of sambada without the ending -Out-Dim.csv)
#' @param dimMax integer The maximum number of dimension given in sambada
#' @param gdsFile char Name of the gds file associated with sambada's input file. If null, will try with \code{sambadaname}.gds
#' @param popStr logical Indicates whether sambada was run using the POPSTRVAR parameter (i.e. population structure was taken into account). Default false
#' @param nrows integer Specifies the number of line to read from the input file. Useful if \code{saveType} 'END ALL' was used in \code{sambadaParallel} and that the number of models run is large so that the reading and processing is too slow. The \code{saveType} 'END' parameter ensures that most significant models are located at the top of the file.
#' @param interactiveChecks logical
#' @return a list containing a) \code{$sambadaOutput} a matrix containing the output from sambada with 3 additional column: corresponding snp, chromosome and position of the marker b) \code{$chrSNPNum} The total number of SNPs in each chromosome c) \code{$chrMaxPos} The highest position found in each chromosome
#' @examples
#' # Example with data from the package
#' # First copy needed files into the temporary directory
#' file.copy(system.file("extdata", "uganda-subset-mol-Out-2.csv", package = "R.SamBada"), 
#'      file.path(tempdir(),'uganda-subset-mol-Out-2.csv'), overwrite=TRUE)
#' file.copy(system.file("extdata", "uganda-subset-mol-storey.csv", package = "R.SamBada"), 
#'      file.path(tempdir(),'uganda-subset-mol-storey.csv'), overwrite=TRUE)
#' if(Sys.info()['sysname']=='Windows'){
#'   file.copy(system.file("extdata", "uganda-subset-mol_windows.gds", package = "R.SamBada"),
#'       file.path(tempdir(),'uganda-subset-mol.gds'), overwrite=TRUE) #If you run Windows
#' } else {
#'   file.copy(system.file("extdata", "uganda-subset-mol_unix.gds", package = "R.SamBada"),
#'       file.path(tempdir(),'uganda-subset-mol.gds'), overwrite=TRUE)
#' }
#' ###################
#' # Run prepareOutput
#' ###################
#' prep=prepareOutput(file.path(tempdir(),'uganda-subset-mol'),2,interactiveChecks=FALSE)
#' @export
prepareOutput = function(sambadaname, dimMax, gdsFile=NULL, popStr=FALSE, nrows=NULL, interactiveChecks=TRUE){
  
  ### Checks ###
  #check input type
  if(typeof(sambadaname)!='character') stop("sambadaname argument supposed to be character")
  if(typeof(dimMax)!='double') stop("dimMax argument supposed to be integer")
  if(dimMax%%1!=0) stop("dimMax argument supposed to be integer")
  if(!is.null(nrows)){
    if(typeof(nrows)!='double') stop("nrows argument supposed to be integer")
    if(nrows%%1!=0) stop("nrows argument supposed to be integer")
  }
  if(typeof(popStr)!='logical') stop("popStr argument supposed to be logical")
  if(typeof(interactiveChecks)!='logical') stop("interactiveChecks argument supposed to be logical")
  
  # All output file exists
  # for(i in 1:dimMax){
  if (!file.exists(paste0(sambadaname,'-Out-',dimMax,'.csv'))) stop(paste0('File ',sambadaname,'-Out-',dimMax,'.csv not found. Check input sambadaname and dimMax'))
  #}
  # Histogram output file exists
  if (!file.exists(paste0(sambadaname,'-storey.csv'))) stop(paste0('File ',sambadaname,'-storey.csv does not exists. This is required to calculate the qvalues'))
  # provided gdsFile exists
  if(!is.null(gdsFile)){
    if (!file.exists(gdsFile)) stop (paste0('File ',gdsFile,' specified in the gdsFile not found. Please enter a valid name'))
  }

  # Check dependencies
  
  sambadaShortName=gsub('_recode','',sambadaname)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package \"data.table\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("Package \"SNPRelate\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("Package \"gdsfmt\" needed for this function to work. Please install it.", call. = FALSE)
  }
  
  if(is.null(nrows)){
    nrows=Inf
  }
  output = tryCatch({data.table::fread(paste0(sambadaname,"-Out-",dimMax,".csv"), nrows=nrows, h=T, colClass=c('character', rep('character',dimMax), rep('double', 3), 'integer', rep('double', (9+(dimMax-1)*3))))}, error=function(e){stop(paste0('Cannot open ',sambadaname,'-Out-',dimMax,'.csv. Maybe the format does not correspond to sambadas strandard and you have an unexpected number of columns'))}) 
  # Check that required columns are present
  if(!('WaldScore' %in% colnames(output))) stop("Column WaldScore not present in your outputfile!")
  if(!('Gscore' %in% colnames(output))) stop("Column Gscore not present in your outputfile!")
  if(!('Marker' %in% colnames(output))) stop("Column Marker not present in your outputfile!")
  
  ### calculate p- and q-value: 
  #Read histogram
  storeyTot=read.table(paste0(sambadaname,'-storey.csv'))

  #p-value and start of qvalue
  if(popStr==TRUE){
    pvalueG=pchisq(getElement(output,'GscorePop'), 1, lower.tail=F)
    pvalueW=pchisq(getElement(output,'WaldScorePop'), 1, lower.tail=F)
    m=sum(storeyTot[nrow(storeyTot),2:ncol(storeyTot)]) #Number of models
    #calculate lambda for G and Wald score
    pi_lambdaG=cumsum(t(as.vector(storeyTot[(nrow(storeyTot)-1),2:ncol(storeyTot)])))/(m*(1-t(as.vector(storeyTot[1,2:ncol(storeyTot)]))))
    pi_lambdaW=cumsum(t(as.vector(storeyTot[nrow(storeyTot),2:ncol(storeyTot)])))/(m*(1-t(as.vector(storeyTot[1,2:ncol(storeyTot)]))))
  } else {
    pvalueG=pchisq(getElement(output,'Gscore'), 1, lower.tail=F)
    pvalueW=pchisq(getElement(output,'WaldScore'), 1, lower.tail=F)  
    m=sum(storeyTot[(3+(dimMax-1)*4),2:ncol(storeyTot)]) #Number of models
    #calculate lambda for G and Wald score
    pi_lambdaG=cumsum(t(as.vector(storeyTot[(3+(dimMax-1)*4),2:ncol(storeyTot)])))/(m*(1-t(as.vector(storeyTot[1,2:ncol(storeyTot)]))))
    pi_lambdaW=cumsum(t(as.vector(storeyTot[(5+(dimMax-1)*4),2:ncol(storeyTot)])))/(m*(1-t(as.vector(storeyTot[1,2:ncol(storeyTot)]))))
  }
  #Bonferroni correction
  pvalueG_Bon=pvalueG*m
  pvalueW_Bon=pvalueW*m

  #Qvalue
  
  #Qvalue based on Gscore
  #splineG=splinefun(t(storeyTot[1,2:ncol(storeyTot)]), pi_lambdaG) #in qvalue package, use spline.smooth (and predict)
  #pi0G=splineG(1)
  splineG <- stats::smooth.spline(t(storeyTot[1,2:ncol(storeyTot)]), pi_lambdaG, df = 3)
  #Estimate pi0 for Gscore
  pi0G <- getElement(stats::predict(splineG, x = t(storeyTot[1,2:ncol(storeyTot)])),'y')
  pi0G <- min(pi0G[1], 1)
  if(interactiveChecks==TRUE){
    #plot histo
    plot(t(storeyTot[1,2:ncol(storeyTot)]), pi_lambdaG, xlab='pValue_G', ylim=c(min(pi0G-0.05, min(pi_lambdaG)),1.01))
    graphics::abline(h=pi0G, col='red')
    cont=readline('Do you want to continue? (press x to exit, any other letter to continue): ')
    if(cont=='x'){
      return(NA)
    }
  }
  
  #Qvalue based on Walsdscore
  #splineW=splinefun(storeyTot[1,2:ncol(storeyTot)], pi_lambdaW)
  #pi0W=splineW(1)
  splineW <- stats::smooth.spline(t(storeyTot[1,2:ncol(storeyTot)]), pi_lambdaW, df = 3)
  #Estimate pi0 for Waldscore
  pi0W <- getElement(stats::predict(splineW, x = t(storeyTot[1,2:ncol(storeyTot)])),'y')
  pi0W <- min(pi0W[1], 1)
  if(interactiveChecks==TRUE){
    ####plot histo + estimated pi0
    plot(t(storeyTot[1,2:ncol(storeyTot)]), pi_lambdaW, xlab='pValue_W', ylim=c(min(pi0W-0.05, min(pi_lambdaW)),1.01))
    graphics::abline(h=pi0W, col='red')
    cont=readline('Do you want to continue? (press x to exit, any other letter to continue): ')
    if(cont=='x'){
      return(NA)
    }
    
  }
  
  #Calculate qvalue
  i=length(pvalueG):1
  qvalueG=pi0G * pmin(1, cummin(pvalueG[order(pvalueG, decreasing=TRUE)] * m /i ))[order(order(pvalueG, decreasing=TRUE))] #from qvalue package
  qvalueW=pi0W * pmin(1, cummin(pvalueW[order(pvalueW, decreasing=TRUE)] * m /i ))[order(order(pvalueW, decreasing=TRUE))]
  
  #Add calculated value to output
  output=cbind(output, "pvalueG"=pvalueG, "pvalueW"=pvalueW, "pvalueG_Bon"=pvalueG_Bon, "pvalueW_Bon"=pvalueW_Bon,"qvalueG"=qvalueG, "qvalueW"=qvalueW)
  
  ### get snp chromosome and position from gds file
  # open gds file
  if(is.null(gdsFile)){
    gdsFile=paste0(sambadaname,'.gds')
  }
  if(!file.exists(gdsFile)){
    stop("A .gds file from package SNPRelate must exist for this function to work. Please provide it in the argument gdsFile if different from sambadaname. Use function prepare_geno from this package or one of the function of SNPRelate")
  }
  gds_obj=SNPRelate::snpgdsOpen(gdsFile)
  on.exit(SNPRelate::snpgdsClose(gds_obj))
  
  # extract chromosome and position
  all_chr=gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.chromosome"))
  all_pos=gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.position"))
  #Take rs.id if exists, id otherwise
  all_snp=tryCatch({gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.rs.id"))}, error=function(e){gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.id"))})
  map2=cbind.data.frame(all_chr, all_snp, as.integer(all_pos))
  colnames(map2)=c('chr','snp','pos')
  
  ### Merge both information (p/qvalue and snp position)

  #Merge acts like join in SQL
  output2=merge(cbind.data.frame('snp'=substr(getElement(output,'Marker'), start = 1, stop = nchar(getElement(output,'Marker')) - 3), output), map2[,c('chr','snp','pos')], by='snp', sort=FALSE, all.x=TRUE)
  
  # Number of SNP by chromosome
  numSNP=table(all_chr)
  # Reorder chromosome: number first in numerical order, letter last in alphabetical order
  numSNP=suppressWarnings(numSNP[c(sort(as.integer(names(numSNP))), names(numSNP)[is.na(as.integer(names(numSNP)))])])
  
  
  # Max position by chromosome
  maxPos=aggregate(all_pos, by=list(all_chr), FUN='max')
  maxPos=suppressWarnings(maxPos[match(c(sort(as.numeric(getElement(maxPos,'Group.1'))), getElement(maxPos,'Group.1')[is.na(as.numeric(getElement(maxPos,'Group.1')))]), getElement(maxPos,'Group.1')),])
  
  #Final output
  finalOutput=list("sambadaOutput"=output2, "chrSNPNum"=numSNP, "chrMaxPos"=maxPos)
  return(finalOutput)
  
}



#' @title Interactive plotting of results
#' @description Plots the manhattan plot for a given environmental variable. The plot is interactive and a map of the distribution of the marker can be retrieved as well as nearby genes listed in Ensembl.
#' @author Solange Duruz
#' @param preparedOutput char The prepared output list from prepare_output function
#' @param varEnv char The name of the environmental variable one wish to study (as in the header of \code{envFile})
#' @param envFile char The file containing the input environmental variable of sambada. 
#' @param species char The abbreviated latin name of the species without capitals nor punctuation (e.g. btaurus, chircus,...). Can be set to null if species not present in ensembl database
#' @param pass integer Number of BP around a SNP in which to look for an annotation in Ensembl. Set to null if species is null
#' @param x char The name of the column corresponding to the x-coordinate in the envFile. Can be set to null if unknown, in this case the maps will not be available
#' @param y char The name of the column corresponding to the y-coordinate in the env file. Can be set to null if x is null.
#' @param valueName char Name of the p- or q-value one wish to plot the manhattan on. This can be either pvalueG, pvalueW, qvalueG, qvalueW for G- or Waldscore respectively.
#' @param chromo char/integer Name or vector of name of the chromosome to investigate. If all is chosen (default), all numerical chromosome will be mapped. If your sambada output is large (typically if you are working with more than 50K genomic file), you should probably map a subset of your dataset (e.g. \code{chromo}=1)
#' @param gdsFile char The GDS file created in the preprocessing of sambada. If null, will try with envFile(without -env.csv or -env-export.csv) and .gds
#' @param IDCol char The name of the column in \code{envFile} corresponding to the ID of the individual. If provided, hover on the output map will give the id of the animal
#' @param popStrCol char The name or vector of name of column(s) in \code{envFile} describing population structure. If provided, additional layers on the map will be available representing population structure.
#' @return None 
#' @examples
#' \dontrun{
#' # Example with data from the package
#' # First copy needed files into the temporary directory
#' file.copy(system.file("extdata", "uganda-subset-mol-Out-2.csv", package = "R.SamBada"), 
#'      file.path(tempdir(),'uganda-subset-mol-Out-2.csv'), overwrite=TRUE)
#' file.copy(system.file("extdata", "uganda-subset-mol-storey.csv", package = "R.SamBada"), 
#'      file.path(tempdir(),'uganda-subset-mol-storey.csv'), overwrite=TRUE)
#' file.copy(system.file("extdata", "uganda-subset-env-export.csv", package = "R.SamBada"), 
#'      file.path(tempdir(),'uganda-subset-env-export.csv'), overwrite=TRUE)
#' if(Sys.info()['sysname']=='Windows'){
#'   file.copy(system.file("extdata", "uganda-subset-mol_windows.gds", package = "R.SamBada"),
#'       file.path(tempdir(),'uganda-subset-mol.gds'), overwrite=TRUE) #If you run Windows
#' } else {
#'   file.copy(system.file("extdata", "uganda-subset-mol_unix.gds", package = "R.SamBada"),
#'       file.path(tempdir(),'uganda-subset-mol.gds'), overwrite=TRUE)
#' }
#' # Run prepareOutput
#' prep=prepareOutput(file.path(tempdir(),'uganda-subset-mol'),2,interactiveChecks=FALSE)
#' ###################
#' # Run plotResultInteractive
#' ###################
#' plotResultInteractive(prep,'bio1','uganda-subset-env-export.csv',species='btaurus',
#'      pass=25000,x='longitude',y='latitude', gdsFile='uganda-subset-mol.gds',
#'      IDCol='short_name',popStrCol='pop1')
#' }
#' @export
plotResultInteractive = function(preparedOutput, varEnv, envFile,species=NULL, pass=NULL,x=NULL,y=NULL,  valueName='pvalueG',chromo='all',gdsFile=NULL, IDCol=NULL, popStrCol=NULL){

  ### Checks
  #preparedOutput
  if(!('sambadaOutput' %in% names(preparedOutput))) stop('preparedOutput should have a component named samabadaOutput. Use the result of the function prepare_output')
  if(!('chrSNPNum' %in% names(preparedOutput))) stop('preparedOutput should have a component named chrSNPNum. Use the result of the function prepare_output')  
  if(!('chrMaxPos' %in% names(preparedOutput))) stop('preparedOutput should have a component named chrMaxPos. Use the result of the function prepare_output')  
  #species
  if(!is.null(species)){
    if(typeof(species)!='character') stop('species argument should be a character')
    #pass
    if(typeof(pass)!='double') stop("pass argument supposed to be a positive integer")
    if(pass%%1>0) stop("pass argument supposed to be a positive integer")
    if(pass<0) stop("pass argument supposed to be a positive integer")
  }
  #envFile and corresponding column
  if(!file.exists(envFile)) stop("envFile not found!")
  envDataTest = read.csv(envFile, header=TRUE, sep=" ", nrows=1)
  if(!is.null(x)){
    if(is.null(y))stop('If x argument is provided, y argument should be provided as well')
    if(!(x %in% names(envDataTest)))stop('x-argument should be part of the header of the envFile')
    if(!(y %in% names(envDataTest)))stop('y-argument should be part of the header of the envFile')
  }
  if(!(varEnv %in% names(envDataTest)))stop('varEnv-argument should be part of the header of the envFile')
  if(!(IDCol %in% names(envDataTest)))stop('IDCol-argument should be part of the header of the envFile')
  if(sum(popStrCol %in% names(envDataTest))!=length(popStrCol))stop('All elements of popStrCol-argument should be part of the header of the envFile')
  #valueName
  if(!(valueName %in% names(preparedOutput$sambadaOutput)))stop('valueName should be a component of preparedOutput$sambadaOutput. Use the result of the function prepare_output as preparedOutput')

  # Test if required libraries are installed
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("Package \"SNPRelate\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package \"httr\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package \"shiny\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package \"plotly\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package \"biomaRt\" needed for this function to work. Please install it.", call. = FALSE)
  }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.", call. = FALSE)
  }
  
  #Read envFile and output
  envData = read.csv(envFile, header=TRUE, sep=" ")
  sambadaOutput = preparedOutput$sambadaOutput
  chrSNPNum = preparedOutput$chrSNPNum
  chrMaxPos = preparedOutput$chrMaxPos
    
  #Connection to ensembl database
  ensemblOutput = ensembl_connection(species, TRUE)
  snp = ensemblOutput$snp
  ensembl = ensemblOutput$ensembl
  
  #Prepare Manhattan 
  #Only model involving chosen varenv
  subset=sambadaOutput[which(sambadaOutput[,'Env_1']==varEnv),]
  if(nrow(subset)==0) stop('No records found in sambada output corresponding to the chosen environmental variable (varEnv argument)')
  #Only models in specified chromosome
  if(length(chromo)>1){
    subset=subset[with(subset, which(chr %in% chromo)),]
  } else if(chromo != 'all' ){
    subset=subset[with(subset, which(chr %in% chromo)),]
  }
  if(nrow(subset)==0) stop('No record found in sambada output corresponding to the chosen chromosome (chromo argument)')
  #..valueName=NULL #so that R check doesn't complain
  subset$pval=-log10(getElement(subset,valueName))
  #Calculate final position (if several chromosome, shift the position, otherwise take initial position)
  prevPos=data.frame("chr"=chrMaxPos$Group.1, "maxPos"=cumsum(as.numeric(chrMaxPos$x))-chrMaxPos$x)
  prevPos=data.frame(prevPos, "chrPos"=rownames(prevPos))
  subset=merge(subset,prevPos,by='chr', sort=FALSE, all.x=TRUE)
  subset$color=colors()[as.integer(as.character(subset$chrPos))%%2+6] 
  if(length(chromo)>1 | chromo=='all'){
    subset$xcoord=subset$maxPos+subset$pos
  } else {
    subset$xcoord=subset$pos
  }
  
  #Define and open GDS file
  if(is.null(gdsFile)){
    gdsFile = paste0(gsub('-env-export.csv','',envFile),'.gds')
    if(!file.exists(gdsFile)){
      gdsFile = paste0(gsub('-env.csv','',envFile),'.gds')
      if(!file.exists(gdsFile)){
        stop("A gds file is needed for this function to work. Specify the name in the input of the function if it is already created, or create it with the package SNPRelate or prepareGeno from this package")
      }
    }
  }
  gds_obj=SNPRelate::snpgdsOpen(gdsFile)

  #Define layout of webpage
  ui <- shiny::shinyUI(shiny::fluidPage(
    shiny::fluidRow(
      shiny::column(
        width = 6,
        shiny::h4("Manhattan"),
        plotly::plotlyOutput("manhattan", height = 350)
      ),
      shiny::column(
        width = 6,
        shiny::h4("Map"),
        plotly::plotlyOutput("map", height = 350)
      )
    ),
    shiny::fluidRow(
      shiny::column(
        width = 8,
        shiny::h4("Information on SNP"),
        shiny::verbatimTextOutput("event2"),
        shiny::h4("Ensembl genes within determined window"),
        shiny::verbatimTextOutput("event"),
        shiny::h4("Variant consequence"),
        shiny::verbatimTextOutput("event4"),
        shiny::h4("Same SNP, different env var"),
        shiny::verbatimTextOutput("event3")
      ),
      shiny::column(
        width = 2,
        shiny::h4("Boxplot"),
        shiny::plotOutput("boxplot", height = 350, width=200)
      )
      
    )

  ))
  
  server <- function(input, output, session) {
    #Prepare manhattan
    p <- ggplot2::ggplot(data=subset, ggplot2::aes_string(x='xcoord', y='pval', colour='color', label='snp', text='pos'), showlegend=FALSE) 
    p <- p + ggplot2::geom_point(size=1)
    p <- p + ggplot2::theme(legend.position="none")
    p <- p + ggplot2::scale_color_manual(values=c("#8B8378", "#CDC0B0"))
    p <- p + ggplot2::scale_y_continuous(name ="-log10(pval)")
    if(length(chromo)>1){
      p <- p + ggplot2::scale_x_continuous(name ="Chromosome", breaks=prevPos$maxPos[chromo], labels=as.character(prevPos$chr[chromo]))
    } else if (chromo == 'all'){
      p <- p + ggplot2::scale_x_continuous(name ="Chromosome", breaks=prevPos$maxPos, labels=as.character(prevPos$chr))
    } else {
      p <- p + ggplot2::scale_x_continuous(name ="Position")
    }

    # Manhattan
    output$manhattan <- plotly::renderPlotly({
      plotly::ggplotly(p = p,tooltip = c("y", "label","x","text"))
    })
    
    #When clicked: query Ensembl
    if(!is.null(species)){
      output$event <- shiny::renderPrint({
        f <- plotly::event_data("plotly_click")
        if (is.null(f)) {
          "Select a point!" 
        }else {
          if(length(chromo)>1 | chromo=='all'){
            selectSNP=subset[which(subset$pos+subset$maxPos==f$x),c('pos','chr')]
          } else {
            selectSNP=subset[which(subset$pos==f$x),c('pos','chr')]
          }
          selectBP=selectSNP[1]$pos
          selectCHR=selectSNP[1]$chr
          #Query ensembl database to get nearby genes
          #c=tryCatch({biomaRt::getBM(attributes=c('chromosome_name','ensembl_gene_id','wikigene_name','start_position','end_position','description'), filters=c("chromosome_name","start","end"), values=list(selectCHR,selectBP-pass,selectBP+pass), mart=ensembl)}, error=function(e){"no gene found!"})
          #c
          biomaRt::getBM(attributes=c('chromosome_name','ensembl_gene_id','wikigene_name','start_position','end_position','description'), filters=c("chromosome_name","start","end"), values=list(selectCHR,selectBP-pass,selectBP+pass), mart=ensembl)
          #Google link => doesn't work because open a local webpage
          #     test=biomaRt::getBM(attributes=c('description'), filters=c("chromosome_name","start","end"), values=list(selectCHR,selectBP-pass,selectBP+pass), mart=ensembl)
          #     test=substr(test,1,regexpr(pattern ='\\[',test)[1]-1)
          #     test=gsub(" ", "+", test)
          #     #renderUI({tagList("URL link:", url)})
          #     #renderUI({tags$a(href = 'www.google.com/search?q=carbohydrate+sulfotransferase','www.google.com/search?q=carbohydrate+sulfotransferase')})
          #     tags$div(class="header", checked=NA,
          #              tags$a(href=paste0('www.google.com/search?q=',test), test))
        }
          
      })
    }
      #When clicked: query Ensembl VCE to get the variant consequence
      if(!is.null(species)){
        output$event4 <- shiny::renderPrint({
          f <- plotly::event_data("plotly_click")
          if (is.null(f)) {
            "Select a point!" 
          }else {
            if(length(chromo)>1 | chromo=='all'){
              selectSNP=subset[which(subset$pos+subset$maxPos==f$x),c('pos','chr')]
            } else {
              selectSNP=subset[which(subset$pos==f$x),c('pos','chr')]
            }
            selectSNP=selectSNP[1]
            server <- "https://rest.ensembl.org"
            #Get reference allele from gds
            snp_id=which(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.chromosome"))==selectSNP$chr  & gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.position"))==selectSNP$pos)
            alleles=gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.allele"), start=c(snp_id), count=c(1))
            alternative_allele=substr(alleles,3,3)
            base_allele=substr(alleles,1,1)
            ext <- paste0("/vep/",species,"/region/",selectSNP$chr,":",selectSNP$pos,"/")
            r <- httr::GET(paste(server, ext,alternative_allele, sep = ""), httr::content_type("application/json"))
            #Try both allele as reference allele
            r<-tryCatch({httr::stop_for_status(r)} , error=function(e){return(httr::GET(paste(server, ext, base_allele,sep = ""), httr::content_type("application/json")))}, finally={r})
            httr::stop_for_status(r)  
            cont=httr::content(r)
            cont[[1]]$most_severe_consequence
          }
        })
      }
    
    #When clicked: give position of SNP
    output$event2<- shiny::renderPrint({
      f <- plotly::event_data("plotly_click")
      if (is.null(f)) {
        "Select a point!" 
      }else {
        if(length(chromo)>1 | chromo=='all'){
          selectSNP=subset[which(subset$pos+subset$maxPos==f$x),c('chr','pos','snp')]
        } else {
          selectSNP=subset[which(subset$pos==f$x),c('chr','pos','snp')]
        }
        selectSNP=selectSNP[1]
        infoOutput=data.frame('snp'=selectSNP$snp, 'chr'=selectSNP$chr, 'BP'=selectSNP$pos)
        infoOutput
      }
    })
    
    #When clicked: get all marker in output of same SNP
    output$event3<- shiny::renderPrint({
      f <- plotly::event_data("plotly_click")
      if (is.null(f)) {
        "Select a point!" 
      }else {
        if(length(chromo)>1 | chromo=='all'){
          selectSNP=subset[which(subset$pos+subset$maxPos==f$x),'snp']
        } else {
          selectSNP=subset[which(subset$pos==f$x),'snp']
        }
        selectSNP=selectSNP[1]
        otherVar=sambadaOutput[sambadaOutput$snp==selectSNP$snp,]
        otherVar=data.frame('Marker'=otherVar$Marker, 'Var'=otherVar$Env_1, 'p/q-value'=otherVar[[valueName]])
        otherVar
      }
    })
    
    #Map
    if(!is.null(x)){
      output$map <- plotly::renderPlotly({
        g <- plotly::event_data("plotly_click")
        if(is.null(g)) {
          plotly::plot_ly(type='scatter', mode='markers')
        } else{ 
          x=envData[,x]
          y=envData[,y]
          varenv=envData[,varEnv]
          popCol=envData[,popStrCol]
          ID=envData[,IDCol]
          # Get marker and snp
          if(length(chromo)>1 | chromo=='all'){
            selectSNP=subset[which(subset$pos+subset$maxPos==g$x),c('chr','pos','Marker')]
          } else {
            selectSNP=subset[which(subset$pos==g$x),c('chr','pos','Marker')]
          }
          selectSNP=selectSNP[1]
          #Retrieve genotype
          snp_id=which(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.chromosome"))==selectSNP$chr  & gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.position"))==selectSNP$pos)
          
          # #Check SNP in LD in the area, not yet implemented
          # snp1 = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "genotype"), start=c(1,snp_id), count=c(-1,1))
          # #loop on nearby snps, break when ld too small
          # snp2 = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "genotype"), start=c(1,snp_id+1), count=c(-1,1))
          # SNPRelate::snpgdsLDpair(snp1, snp2, method = "corr")
          
          geno=SNPRelate::snpgdsGetGeno(gds_obj, snp.id=snp_id)
          pres=genoToMarker(gds_obj, selectSNP$Marker)
          xy=data.frame(x,y,varenv,ID, geno, pres, popCol)
         
          pal1 <- c("antiquewhite2", "black")
          
          graph <- plotly::plot_ly()
          #Plot marker distribution
          graph <- plotly::add_trace(graph,data=xy, x=x,y=y, type='scatter',mode='markers', color=pres, marker=list(showscale=FALSE), name='marker', colors=pal1, text=ID,hoverinfo = c("color"))
          graph <- plotly::hide_colorbar(graph)  
          #Plot environmental variable
          graph <- plotly::add_markers(graph,data=xy, inherit=FALSE, mode='markers', x=~x,y=~y, name='varenv', marker=list(color=~varenv,colorscale = 'YlOrRd', showscale=TRUE, colorbar=list(len=0.3, title='varenv',y=1)))
          #Plot population structure
          if(!is.null(popStrCol)){
            if(length(popStrCol)==1){
              graph <- plotly::add_markers(graph,data=xy, inherit=FALSE, mode='markers',x=~x,y=~y, name='pop1', marker=list(color=~popCol,colorscale = 'Greens', showscale=TRUE, colorbar=list(len=0.3, title='pops', y=0.7)))
            } else {  
              minpop=min(xy[,popStrCol])
              maxpop=max(xy[,popStrCol])
              graph <- plotly::add_markers(graph,data=xy, inherit=FALSE, mode='markers',x=~x,y=~y, name='pop1', marker=list(color=~get(popStrCol[1]),colorscale = 'Greens', cmin=minpop, cmax=maxpop, showscale=TRUE, colorbar=list(len=0.2, title='pops', y=0.75)))
              for(p in 2:length(popStrCol)){
                graph <- plotly::add_markers(graph,data=xy, inherit=FALSE, mode='markers',x=~x,y=~y, name=paste0('pop',p), marker=list(color=~get(popStrCol[p]),colorscale = 'Greens',cmin=minpop, cmax=maxpop,  showscale=FALSE))
              }
            }
          }
          graph
        }
        
      })
    } 
    
    #When clicked shows boxplot
    output$boxplot <- shiny::renderPlot({
      g <- plotly::event_data("plotly_click")
      if(is.null(g)) {
        #plotly::plot_ly(type='scatter', mode='markers')
      } else{
        varenv=envData[,varEnv]
        popCol=envData[,popStrCol]
        # Get marker and snp
        if(length(chromo)>1 | chromo=='all'){
          selectSNP=subset[which(subset$pos+subset$maxPos==g$x),c('chr','pos','Marker')]
        } else {
          selectSNP=subset[which(subset$pos==g$x),c('chr','pos','Marker')]
        }
        selectSNP=selectSNP[1]
        #Retrieve genotype
        snp_id=which(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.chromosome"))==selectSNP$chr  & gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.position"))==selectSNP$pos)
        pres=genoToMarker(gds_obj, selectSNP$Marker)
        xy=data.frame(varenv, pres, popCol)
        boxplot(varenv~pres, data=xy)
      }
    })
    
    #When close browser, stop the app to stop the function and close the gds file
    session$onSessionEnded(function() {
      shiny::stopApp()
      SNPRelate::snpgdsClose(gds_obj)
    })
  }
  
  shiny::shinyApp(ui, server)

}

#' @title Manhattan plot 
#' @description Plot the manhattan plot for a given environmental data
#' @author Solange Duruz
#' @param preparedOutput char The prepared output list from prepare_output function
#' @param varEnv char The name of the environmental variable one wish to study. Can be a vector of char if you want to plot several \code{varEnv} at a row. If \code{saveType} is NULL, the program prompts to continue. If \code{saveType} is png or pdf, several files are saved
#' @param valueName char Name of the p- or q-value one wish to plot the manhattan on. This can be either pvalueG, pvalueW, qvalueG, qvalueW for G- or Waldscore respectively.
#' @param chromo char/integer Name or vector of name of the chromosome to investigate. If all is chosen (default), all numerical chromosome will be mapped. If your sambada output is large (typically if you are working with more than 50K genomic file), you should probably map a subset of your dataset (e.g. \code{chromo}=1)
#' @param saveType char One of NULL, 'png' or 'pdf'. If NULL is set, the plot will be shown in the R plotting window. Otherwise, it will be saved in the specified format in your working directory with the name 'manhattan-' followed by varEnv.
#' @param threshold double A digit number indicating a value to draw a threshold line
#' @param highlight char Name of the genotype to highlight in red on plot (should be SNPName_Genotype e.g. 'ARS-BFGL-NGS-106879_AA')
#' @return The last plot object (if several \code{varEnv} are specified, only the last one is returned)
#' @examples
#' # Example with data from the package
#' # First copy needed files into the temporary directory
#' file.copy(system.file("extdata", "uganda-subset-mol-Out-2.csv", package = "R.SamBada"), 
#'     file.path(tempdir(),'uganda-subset-mol-Out-2.csv'), overwrite=TRUE)
#' file.copy(system.file("extdata", "uganda-subset-mol-storey.csv", package = "R.SamBada"), 
#'     file.path(tempdir(),'uganda-subset-mol-storey.csv'), overwrite=TRUE)
#' if(Sys.info()['sysname']=='Windows'){
#'   file.copy(system.file("extdata", "uganda-subset-mol_windows.gds", package = "R.SamBada"),
#'       file.path(tempdir(),'uganda-subset-mol.gds'), overwrite=TRUE) #If you run Windows
#' } else {
#'   file.copy(system.file("extdata", "uganda-subset-mol_unix.gds", package = "R.SamBada"),
#'       file.path(tempdir(),'uganda-subset-mol.gds'), overwrite=TRUE) #If you run Unix
#' }
#' # Run prepareOutput
#' prep=prepareOutput(file.path(tempdir(),'uganda-subset-mol'),2,interactiveChecks=FALSE)
#' ###################
#' # Run plotManhattan
#' ###################
#' plotManhattan(prep, c('bio1'),chromo='all',valueName='pvalueG')
#' \donttest{
#' # Example with several environmental variables
#' plotManhattan(prep,c('bio1','bio2'),'pvalueG')
#' }
#' @export
plotManhattan=function(preparedOutput, varEnv, valueName, chromo='all',saveType=NULL, threshold=NULL, highlight=NULL){
  ### Checks
  #preparedOutput
  if(!('sambadaOutput' %in% names(preparedOutput))) stop('preparedOutput should have a component named samabadaOutput. Use the result of the function prepare_output')
  if(!('chrSNPNum' %in% names(preparedOutput))) stop('preparedOutput should have a component named chrSNPNum. Use the result of the function prepare_output')  
  if(!('chrMaxPos' %in% names(preparedOutput))) stop('preparedOutput should have a component named chrMaxPos. Use the result of the function prepare_output')  

  #valueName
  if(!(valueName %in% names(preparedOutput$sambadaOutput)))stop('valueName should be a component of preparedOutput$sambadaOutput. Use the result of the function prepare_output as preparedOutput')

  # Test if required libraries are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.", call. = FALSE)
  }
  
  sambadaOutput = preparedOutput$sambadaOutput
  chrSNPNum = preparedOutput$chrSNPNum
  chrMaxPos = preparedOutput$chrMaxPos
  
  for(i in 1:length(varEnv)){
    #Only chosen varEnv
    subset=sambadaOutput[which(sambadaOutput[,'Env_1']==varEnv[i]),]
    if(nrow(subset)==0) stop('No records found in sambada output corresponding to the chosen environmental variable (varEnv[i] argument)')
    #Only chosen chromosome
    if(length(chromo)>1){
      subset=subset[with(subset, which(chr %in% chromo)),]
    } else if(chromo != 'all' ){
      subset=subset[with(subset, which(chr %in% chromo)),]
    }
    if(nrow(subset)==0) stop('No record found in sambada output corresponding to the chosen chromosome (chromo argument)')
    #..valueName=NULL
    subset$pval=-log10(getElement(subset,valueName))

    #Calculate final position (if several chromosome, shift the position, otherwise take initial position)
    prevPos=data.frame("chr"=chrMaxPos$Group.1, "maxPos"=cumsum(as.numeric(chrMaxPos$x))-chrMaxPos$x)
    prevPos=data.frame(prevPos, "chrPos"=rownames(prevPos))
    subset=merge(subset,prevPos,by='chr', sort=FALSE, all.x=TRUE)
    subset$color=as.character(as.integer(as.character(subset$chrPos))%%2)
    if(!is.null(highlight)){
      subset[subset$Marker %in% highlight,]$color='3'
    }
    if(length(chromo)>1 | chromo=='all'){
      subset$xcoord=subset$maxPos+subset$pos
    } else {
      subset$xcoord=subset$pos
    }
    
    #Prepare manhattan
    p <- ggplot2::ggplot(data=subset, showlegend=FALSE)
    p <- p + ggplot2::geom_point(ggplot2::aes_string(x='xcoord', y='pval', colour='color'),size=1)
    if(!is.null(highlight)){
      p <- p + ggplot2::geom_point(data = subset(subset, subset$color == '3'),ggplot2::aes_string(x='xcoord', y='pval', colour='color'),size=1)
    }
    p <- p + ggplot2::theme(legend.position="none")
    if(!is.null(highlight) & (length(chromo)>1 | chromo=='all')){
      #If SNP(s) must be highlighted + at least two chromosome, define 3 colors (2 greys, 1 red)
      p <- p + ggplot2::scale_color_manual(values=c("#8B8378", "#CDC0B0", '#FF0000'))
    } else if (!is.null(highlight) & length(chromo)==1){
      #If no SNP must be highlighted + at least two chromosome, define 2 colors (2 greys)
      p <- p + ggplot2::scale_color_manual(values=c("#8B8378", '#FF0000'))
    } else {
      #If SNP(s) must be highlighted + one chromosome, define 2 colors (1 grey, 1 red)
      p <- p + ggplot2::scale_color_manual(values=c("#8B8378", "#CDC0B0"))
    }
    p <- p + ggplot2::scale_y_continuous(name =paste0("-log10(",valueName,")"))
    p <- p + ggplot2::ggtitle(varEnv[i])
    
    #x-label: either chromosome number if several chromosomes, or bp
    if(length(chromo)>1){
      p <- p + ggplot2::scale_x_continuous(name ="Chromosome", breaks=prevPos$maxPos[chromo], labels=as.character(prevPos$chr[chromo]))
    } else if (chromo == 'all'){
      p <- p + ggplot2::scale_x_continuous(name ="Chromosome", breaks=prevPos$maxPos, labels=as.character(prevPos$chr))
    } else {
      p <- p + ggplot2::scale_x_continuous(name ="Position")
    }
    if(!is.null(threshold)){
      p <- p + ggplot2::geom_hline(yintercept=-log10(threshold), colour='red')
    }

    if(!is.null(saveType)){
      if(saveType=='png'){
        png(paste0('manhattan-',varEnv[i],'.png'))
      }
      if(saveType=='pdf'){
        pdf(paste0('manhattan-',varEnv[i],'.pdf'))
      }
    }
    # Manhattan
    plot(p)
    
    if(!is.null(saveType)){
      dev.off()
    } else {
      if(i<length(varEnv)){
        cont=readline('Would you like to continue? Press x to exit, any other letter to continue: ')
        if(cont=='x'){
          return(NA)
        }
      }
    }
  }
  #Only returns the last plot
  return(p)
}


#' @title Plotting of maps 
#' @description Plots several kinds of maps (environmental variable distribution, population structure, marker absence or presence, autocorrelation of marker). Unlike \code{\link{plotResultInteractive}}, the resulting maps are non-interactive. The function can handle several marker/variables at once and create separate output files.
#' @author Solange Duruz
#' @param envFile char The file containing the input environmental variable of sambada. 
#' @param x char The name of the column corresponding to the x-coordinate in the \code{envFile}. Can be set to null if unknown, in this case the maps will not be available
#' @param y char The name of the column corresponding to the y-coordinate in the env file. Can be set to null if x is null.
#' @param gdsFile char The GDS file created in the preprocessing of sambada. If null, will try with \code{envFile}(without -env.csv) and .gds
#' @param popStrCol char The name or vector of name of column(s) in \code{envFile} describing population structure. If provided, additional layers on the map will be available representing population structure.
#' @param locationProj integer EPSG code of the geographical projection in the \code{envFile}
#' @param markerName name of the marker to be plotter if \code{mapType} is 'marker' or 'AS'. \code{markerName} can be found in preparedOutput$sambadaOutput[,''] where preparedOutput would be the result of the function \code{prepareOutput}
#' @param mapType char A string or vector of string containing one or several of 'marker' (presence/absence of marker), 'env' (environmental variable distribution), 'popStr' (population variable on continuous scale), 'popPieChart' (belonging to a population in pie charts), 'AS' (autocorrelation of the marker). Note that the background of all maps, if found, will be the raster of the environmental variable. Thus the 'env' \code{mapType} is preferred when no raster is provided. For the 'AS' type, it is calculated on the fly for the markers provided and not the one possibly calculated by sambada.
#' @param varEnvName char Name of the environmental variable. If a raster of the variable is located in your working directory, you can provide \code{varEnvName} even for \code{mapType} such as 'marker' or 'AS'. The function will scan the folder of your working directory for raster with the same name as \code{varEnvName} (and commonly used extension for raster) and put it as background.
#' @param SAMethod char If \code{mapType} contains 'AS', then you must specify the method for setting the weights of neighbours. Can be one of 'knn' (k-nearest neighbours) or 'distance' 
#' @param SAThreshold char If \code{mapType} contains 'AS' and \code{SAMethod} is 'knn' then the number of neighbours. If \code{SAMethod} is 'distance' then the distance in map-unit (unless you use a spherical projection (latitude/longitude), in which case you should use km)
#' @param saveType char One of NULL, 'png' or 'pdf'. If NULL is set, the maps will be shown in the R plotting window. Otherwise, it will be saved in the specified format in your working directory.
#' @param rasterName char If a raster file with the environmental variable distribution exists with a different name than \code{varEnvName}, provide it here (including extension)
#' @param simultaneous boolean If TRUE and \code{mapType} contains several kinds of maps, all maps corresponding to the same marker will be plotted on the same window. The resulting maps can be very small.
#' @return None 
#' @examples
#' # Define right GDS file according to your OS
#' if(Sys.info()['sysname']=='Windows'){
#'   gdsFile=system.file("extdata", "uganda-subset-mol_windows.gds", package = "R.SamBada")
#' } else {
#'   gdsFile=system.file("extdata", "uganda-subset-mol_unix.gds", package = "R.SamBada")
#' }
#' #############
#' # Run plotMap
#' #############
#' plotMap(envFile=system.file("extdata", "uganda-subset-env-export.csv", package = "R.SamBada"), 
#'      x='longitude', y='latitude', locationProj=4326,  popStrCol='pop1', gdsFile=gdsFile, 
#'      markerName='Hapmap28985-BTA-73836_GG', mapType='marker', varEnvName='bio1', simultaneous=FALSE)
#' \donttest{
#' # Maps of marker and population structure (two subplot)
#' plotMap(envFile=system.file("extdata", "uganda-subset-env-export.csv", package = "R.SamBada"),
#'      'longitude','latitude', locationProj=4326,  popStrCol='pop1', 
#'      gdsFile=gdsFile, markerName='Hapmap28985-BTA-73836_GG', 
#'      mapType=c('marker', 'popStr'), varEnvName='bio1', simultaneous=TRUE)
#' }
#' @export
plotMap = function(envFile, x, y, locationProj,  popStrCol, gdsFile, markerName, mapType, varEnvName, SAMethod=NULL, SAThreshold=NULL, saveType=NULL, rasterName=NULL, simultaneous=FALSE){

  # Test if required libraries are installed
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("Package \"SNPRelate\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("Package \"gdsfmt\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package \"sp\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("packcircles", quietly = TRUE)) {
    stop("Package \"packcircles\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Package \"raster\" needed for this function to work. Please install it.", call. = FALSE)
  }  
  if (!requireNamespace("mapplots", quietly = TRUE)) {
    stop("Package \"mapplots\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if('AS' %in% mapType){
    if (!requireNamespace("spdep", quietly = TRUE)) {
	  stop("Package \"spdep\" needed for this function to work. Please install it.", call. = FALSE)
    }  
  }
  ### Check inputs ###
  if(typeof(envFile)!='character') stop('envFile argument should be a character')
  if(!file.exists(envFile)){
    stop("envFile not found")
  }
  if(typeof(x)!='character' | typeof(y)!='character') stop('x and y arguments should be a character')
  firstLine=read.csv(envFile, header=TRUE, sep=" ", nrows=1)
  if(!(x %in% names(firstLine)) | !(y %in% names(firstLine))) stop('x and y should be in the header line of envFile')
  if(typeof(locationProj)!='double') stop('locationProj should be an integer (EPSG code)')
  t=tryCatch(sp::CRS(paste0("+init=epsg:",locationProj)),error = function(e) {stop('Invalid EPSG code in locationProj')})
  if(sum(popStrCol %in% names(firstLine))!=length(popStrCol)) stop('popStrCol should all be included in the header line of envFile')
  if(is.null(gdsFile)){
    gdsFile = paste0(gsub('-env.csv','',envFile),'.gds')
  }
  if(!file.exists(gdsFile) & ('AS' %in% mapType | 'marker' %in% mapType)){
    stop("A gds file is needed for this function to work if mapType contains 'marker' or 'AS'. Specify the name in the input of the function if it is already created, or create it with the package SNPRelate or prepareGeno from this package")
  }
  if(typeof(varEnvName)!='character') stop('markerName argument should be a character')
  #Test if all SNPs in GDS
  gds_obj=SNPRelate::snpgdsOpen(gdsFile)
  on.exit(SNPRelate::snpgdsClose(gds_obj))
  if(!is.null(markerName)){
    for(m in 1:length(markerName)){
      if(length(which(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.rs.id"))==substr(markerName[m], 1, nchar(markerName[m])-3)))==0) stop(paste0(markerName[m]," not found in gds"))
    }
  }
  if(sum(mapType %in% c('marker','env','AS','popStr', 'popPieChart'))!=length(mapType)) stop("mapType should be one, or several of, 'marker','env','AS','popStr', 'popPieChart'")
  if(('popPieChart' %in% mapType) & (length(popStrCol)<2)) stop('popStrCol should have a length >= 2 if mapType is popPieChart')
  if(typeof(varEnvName)!='character') stop('varEnvName argument should be a character')
  if(!is.null(SAMethod)) if(!(SAMethod %in% c('knn','distance'))) stop('SAMethod should be one of knn or distance')
  if(!is.null(SAThreshold)) if(typeof(SAThreshold)!='double') stop ('SAThreshold should be numeric')
  if(('AS' %in% mapType | 'marker' %in% mapType) & is.null(markerName)) stop('markerName cannot be null if mapType is "AS" or "marker"')
  if('AS' %in% mapType & is.null(SAMethod)) stop('SAMethod must be not null if mapType include AS')
  if('AS' %in% mapType & is.null(SAThreshold)) stop('SAThreshold must be not null if mapType include AS')
  if('AS' %in% mapType) if(SAMethod=='knn' & SAThreshold%%1!=0) stop ("SAThreshold should be an integer if SAMethod is 'knn'")
  if(!is.null(saveType)) if(!(saveType %in% c('png','pdf'))) stop("saveType should be one of NULL,'png','pdf'")
  if(!is.null(rasterName)) if(typeof(rasterName)!='character') stop('rasterName argument should be a character')
  if(!is.null(rasterName)) if(!file.exists(rasterName)) stop('rasterName file not found')
  
  #Open Env file
  envData=read.csv(envFile, header=TRUE, sep=" ")
  envData2=envData
  sp::coordinates(envData) = c(x,y)
  sp::proj4string(envData)=sp::CRS(paste0("+init=epsg:",locationProj))
  
  #Define layout
  numSimPlot=0 #number of simultaneous plots
  if('marker' %in% mapType) numSimPlot=numSimPlot+1
  if('AS' %in% mapType) numSimPlot=numSimPlot+1
  if('env' %in% mapType) numSimPlot=numSimPlot+1
  if('popStr' %in% mapType) numSimPlot=numSimPlot+length(popStrCol)
  if('popPieChart' %in% mapType) numSimPlot=numSimPlot+1
  #Define the order in which plots should appear
  allMapType=c('marker','AS','env','popStr','popPieChart')
  mapTypeb=allMapType[allMapType %in% mapType]
  if('popStr' %in% mapTypeb){
    mapTypeb=c(mapTypeb, rep('popStr',length(popStrCol)-1))
  }
  mapType2=rep(mapTypeb, max(1,length(markerName)))
  
  if(numSimPlot>1 & simultaneous==TRUE){
    numrow=2*numSimPlot
    matrixLayout=rep(c(1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,3),numSimPlot)+rep(seq(from=0, by=3, length=numSimPlot),each=20)
  } else {
    numrow=2
    matrixLayout=c(1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,3)
  }
  tryCatch({dev.off()}, error=function(e){})
  layout(matrix(matrixLayout, nrow = numrow, ncol = 10, byrow = TRUE))
  par(mar=c(2,2,2,2))
  
  #plot.new()
  #layout(matrix(c(1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,3), nrow = 2, ncol = 10, byrow = TRUE))
  #par(mar=c(2,2,2,2))
  
  #Jitter the points to avoid overlapp 
  size=0.1*max((max(envData@coords[,x])-min(envData@coords[,x]))/dev.size(units='cm')[1],(max(envData@coords[,y])-min(envData@coords[,y]))/dev.size(units='cm')[2])
  data_df=data.frame(size=rep(size,nrow(envData)),x=envData@coords[,x],y=envData@coords[,y])
  scattered_point=packcircles::circleRepelLayout(data_df,NULL,NULL,xysizecols=c('x','y','size'),sizetype='radius',maxiter = 50)
  sp::coordinates(scattered_point$layout)=c('x','y')

  #Draw plots
  numMark=0
  numEnv=0
  numPopCol=0
  for(t in 1:length(mapType2)){
    if(mapType2[t]=='marker' | (mapType2[t]=='AS' & mapType2[max(1,t-1)]!='marker')){
      numMark=numMark+1
    } else if(mapType2[t]=='env'){
      numEnv=numEnv+1
    } else if(mapType2[t]=='popStr'){
      if(numPopCol==length(popStrCol)){
        numPopCol=0
      }
      numPopCol=numPopCol+1
    } 
    numEnv=max(numEnv, numMark)
    
    if(simultaneous==TRUE & !is.null(saveType)){ #####!!!!!!!!!!!!!!!!!!!!!!!!!############
      if('marker' %in% mapType2 | 'AS' %in% mapType2){
        mapName=paste0(markerName[numMark],'_map.',saveType)
      } else if ('env' %in% mapType2){
        mapName=paste0(varEnvName[numEnv],'_map.',saveType)
      } else {
        mapName=paste0('popStr_map.',saveType)
      }
      
    } else if (simultaneous==FALSE & !is.null(saveType)){
      if('marker' %in% mapType2[t]){
        mapName=paste0(markerName[numMark],'_map.',saveType)
      } else if('AS' %in% mapType2[t]){
        mapName=paste0(markerName[numMark],'_AS_map.',saveType)
      } else if ('env' %in% mapType2[t]){
        mapName=paste0(varEnvName[numEnv],'_map.',saveType)
      } else {
        mapName=paste0('popStr_map.',saveType)
      }      
    }
    #for(m in 1:loopNum){
      #Get Marker info
      if('marker' %in%  mapType | 'AS' %in% mapType ){
        #Open gds File
        pres=genoToMarker(gds_obj, markerName[numMark])
      }
      #Try to find corresponding raster
      allowedExtension=c('bil','tif')
      if(exists('rasterName')) {
        rm(rasterName)
      }
      if(length(varEnvName)>1){
        varEnvName2=varEnvName[numEnv]
      } else {
        varEnvName2=varEnvName
      }
      
      if(regexpr('bio',varEnvName2)>0){
        if(file.exists(paste0('wc0.5/',varEnvName2,'.tif'))){
          rasterName=paste0('wc0.5/',varEnvName2,'.tif')
        }
      } else if (regexpr('bio',varEnvName2)>0) {
        if(file.exists(paste0('srtm/',varEnvName2,'.tif'))){
          rasterName=paste0('srtm/',varEnvName2,'.tif')
        }
      } else {
        for(aE in 1:length(allowedExtension))
          if(file.exists(paste0(varEnvName2,'.',allowedExtension(aE)))){
            rasterName=paste0(varEnvName2,'.',allowedExtension(aE))
            break
          }
      }
      #Open raster
      if(exists('rasterName')){
        raster=raster::raster(rasterName)
        #Get real raster data
        raster_df=as.data.frame(raster::sampleRegular(raster, size=1e5, asRaster=FALSE), xy=TRUE)
      }
      
      par(mar=c(2,2,2,2), xpd=FALSE)
      #Draw background
      if(exists('rasterName')){
        #If raster found, put it as background
        #Put coordinates of scattered point or envData???
        raster::image(raster, asp=1, maxpixels=10000000000,  col=terrain.colors(100),xlim = c(min(envData@coords[,x]), max(envData@coords[,x])), ylim = c(min(envData@coords[,y]), max(envData@coords[,y])))
        
      }else {
        #If raster not found, put countries as background
        country=data('wrld_simpl', package='maptools', envir=environment())
        raster::plot(country,xlim=c(min(sp::coordinates(envData)[,x]),max(sp::coordinates(envData)[,x])),ylim=c(min(sp::coordinates(envData)[,y]),max(sp::coordinates(envData)[,y])))
      }
      
      #Draw lines between original location and scattered one (if not scattered, the lines will be masked by the points)
      for(i in 1:nrow(envData)){
        lines(c(envData@coords[i,x],scattered_point$layout@coords[i,'x']),c(envData@coords[i,y],scattered_point$layout@coords[i,'y']), col='antiquewhite3')
      }
      
      #Draw points
      if(mapType2[t]=='marker'){
        raster::plot(scattered_point$layout,col=colors()[pres*19+5],pch=16,add=TRUE)
      } else if(mapType2[t]=='env'){
        #Define color palette
        new.pal=colorRampPalette(c("yellow", "orange","red"))( 100 )
        raster::plot(scattered_point$layout, col=new.pal[round((envData@data[,varEnvName]-min(envData@data[,varEnvName]))/(max(envData@data[,varEnvName])-min(envData@data[,varEnvName]))*100)],pch=16,add=TRUE)
      } else if(mapType2[t]=='popStr'){
          new.pal=colorRampPalette(c("white","black"))( 100 )
          raster::plot(scattered_point$layout, col=new.pal[round((envData@data[,popStrCol[numPopCol]]-min(envData@data[,popStrCol[numPopCol]]))/(max(envData@data[,popStrCol[numPopCol]])-min(envData@data[,popStrCol[numPopCol]]))*100)],pch=16,add=TRUE)
      } else if(mapType2[t]=='popPieChart'){
        for(i in 1:nrow(scattered_point$layout)){ 
          mapplots::add.pie(z=as.double(abs(1/envData@data[i,popStrCol])),x=scattered_point$layout@coords[i,'x'], scattered_point$layout@coords[i,'y'], labels=NA, col=terrain.colors(3), radius=size*2)
        }
      } else if(mapType2[t]=='AS'){
        ### autocorrelation
        #Calculate autocorrelation
        #knearneigh, dnearneigh
        if(SAMethod=='knn'){
          knn=spdep::knearneigh(envData,5)
          nb=spdep::knn2nb(knn)
        } else if (SAMethod=='distance'){
          nb=spdep::dnearneigh(envData,0,20)
        }
        nblist=spdep::nb2listw(nb, zero.policy=TRUE) 
        #iglobal=spdep::moran.test(as.vector(pres),nblist, na.action=na.omit)
        ilocal=spdep::localmoran(as.vector(pres+1),nblist, zero.policy=TRUE,na.action=na.exclude)
        #ilocalPval=spdep::moran.mc(as.vector(pres),nblist, nSim, zero.policy=TRUE,na.action=na.exclude)
        #Define color for map
        color=vector('character',nrow(ilocal))
        color[ilocal[,'Ii']<(-0.5)]='blue4'
        color[ilocal[,'Ii']<(-0.1) & ilocal[,'Ii']>=(-0.5)]='lightblue3'
        color[ilocal[,'Ii']<(0.1) & ilocal[,'Ii']>=(-0.1)]='wheat'
        color[ilocal[,'Ii']<(0.5) & ilocal[,'Ii']>=(0.1)]='salmon'
        color[ilocal[,'Ii']>=0.5]='red3'
        color[is.na(ilocal[,'Ii'])]='white'
        #Define pch of points in map according to significance
        pointch=vector('integer',nrow(ilocal))
        pointch[ilocal[,'Pr(z > 0)']<=(0.05)]=19
        pointch[ilocal[,'Pr(z > 0)']>(0.05)]=21
        pointch[is.na(ilocal[,'Pr(z > 0)'])]=19
        # pointch[ilocalPval$p.value<=(0.05)]=19
        # pointch[ilocalPval$p.value>(0.05)]=21
        # pointch[is.na(ilocalPval$p.value)]=19
        #Draw points
        raster::plot(scattered_point$layout,col=color,pch=pointch,cex=1.2,bg='grey',add=TRUE)
      }
      
      #Draw legend
      if(exists('rasterName')){
        #Raster legend
        par(mar=c(2,1,3,2), xpd=NA)
        #raster.pal=colorRampPalette(c("yellow", "orange","red"))( 100 )
        raster.pal=terrain.colors( 100 )
        
        image(1, 1:100, t(seq_along(1:100)), col=raster.pal, axes=FALSE , xlab="", ylab="")
        axis(4, at=(pretty(raster_df[,varEnvName2])[2:(length(pretty(raster_df[,varEnvName2]))-1)]-min(raster_df[,varEnvName2], na.rm=TRUE))/(max(raster_df[,varEnvName2], na.rm=TRUE)-min(raster_df[,varEnvName2], na.rm=TRUE))*100, labels=pretty(raster_df[,varEnvName2])[2:(length(pretty(raster_df[,varEnvName2]))-1)])
        text(1,107,varEnvName2)
      } else {
        plot.new()
        par(mar=c(2,1,3,2), xpd=NA)
      }

      if(mapType2[t]=='popPieChart'){
        # Point legend
        #par(mar=c(2,1,3,2), xpd=NA)
          points(rep(0,length(popStrCol)),seq(from=-20, by=-10, length.out=length(popStrCol)),pch=19, col=terrain.colors(length(popStrCol)))
          text(rep(0.1,length(popStrCol)),seq(from=-20, by=-10, length.out=length(popStrCol)),popStrCol, pos=4)
          text(1,-10,'Population')          
      } else if(mapType[t]=='env'){
        # Point legend
        par(mar=c(2,1,3,2), xpd=NA)
        image(1, 1:100, t(seq_along(1:100)), col=new.pal, axes=FALSE, xlab="", ylab="")
        axis(4, at=(pretty(envData@data[,varEnvName2])[2:(length(pretty(envData@data[,varEnvName2]))-1)]-min(envData@data[,varEnvName2], na.rm=TRUE))/(max(envData@data[,varEnvName2], na.rm=TRUE)-min(envData@data[,varEnvName2], na.rm=TRUE))*100, labels=pretty(envData@data[,varEnvName2])[2:(length(pretty(envData@data[,varEnvName2]))-1)])
        text(1,107,varEnvName2)
      } else if(mapType2[t]=='popStr'){
        par(mar=c(2,1,3,2), xpd=NA)
        pop.pal=colorRampPalette(c("white", "black"))( 100 )
        image(1, 1:100, t(seq_along(1:100)), col=pop.pal, axes=FALSE , xlab="", ylab="")
        axis(4, at=(pretty(envData@data[,popStrCol[numPopCol]])[2:(length(pretty(envData@data[,popStrCol[numPopCol]]))-1)]-min(envData@data[,popStrCol[numPopCol]], na.rm=TRUE))/(max(envData@data[,popStrCol[numPopCol]], na.rm=TRUE)-min(envData@data[,popStrCol[numPopCol]], na.rm=TRUE))*100, labels=pretty(envData@data[,popStrCol[numPopCol]])[2:(length(pretty(envData@data[,popStrCol[numPopCol]]))-1)])
        text(1,107,'Population')
      } else if(mapType2[t]=='AS'){
        # Point legend
        #par(mar=c(2,1,3,2), xpd=NA)
        plot.new()
        points(rep(0,10),c(seq(from=1, by=-0.1, length.out=5),seq(from=0.4, by=-0.1, length.out=5)),pch=c(19,19,19,19,19,21,21,21,21,21), col=c('blue4','lightblue3','wheat','salmon','red3','blue4','lightblue3','wheat','salmon','red3'), bg='gray')
        text(rep(0.1,10),c(seq(from=1, by=-0.1, length.out=5),seq(from=0.4, by=-0.1, length.out=5)),c('< -0.5','-0.5 - -0.1','-0.1 - 0.1','0.1 - 0.5','> 0.5'), pos=4)
        text(0,1.1,'I (signif)', pos=4)
        text(0, 0.5,'I (not signif)', pos=4)
      } else if(mapType2[t]=='marker'){
        plot.new()
        points(c(0,0),c(1,0.8),pch=19, col=c(colors()[5],colors()[24]))
        text(c(0.1,0.1),c(1,0.8),c('Absent','Present'),pos=4)
      }
      if((simultaneous==FALSE & is.null(saveType) & t<length(mapType2)) | (t%%numSimPlot==0 & t<length(mapType2))){
        continue = readline(prompt="Would you like to continue? (press x to exit, any other letter to continue): ")
        if(continue=='x'){
          return(NA)
        } 
      } else if (simultaneous==FALSE & !is.null(saveType)){
        grDevices::dev.copy(get(saveType), mapName)  
        dev.off()
      }
    #}
    if(simultaneous==TRUE & !is.null(saveType)){
      grDevices::dev.copy(get(saveType), mapName)  
      dev.off()        
    }
  }
  
}



ensembl_connection = function(species, interactiveChecks){
  
  ensembl_dataset = biomaRt::useMart('ensembl')
  ensembl_dataset=biomaRt::listDatasets(ensembl_dataset)
  dataset_found=ensembl_dataset[grepl(species, ensembl_dataset[,'dataset']),'dataset']
  if(length(dataset_found)==0) stop('Species not found in ensembl database. Either change the name of the species (latin naming e.g. btaurus for cattle) or set it to NULL')
  
  snp_dataset = biomaRt::useMart('ENSEMBL_MART_SNP')
  snp_dataset=biomaRt::listDatasets(snp_dataset)
  snp_found=snp_dataset[grepl(species, snp_dataset[,'dataset']),'dataset']
  
  if(interactiveChecks==TRUE){
    n=readline(paste0('Dataset found :', dataset_found,'. Would you like to continue? Press x to exit, any other key to continue '))
    if(n=='x'){
      print('Function ended on user input')
      return(NA)
    }
    n=readline(paste0('SNP Data found :', snp_found,'. Would you like to continue? Press x to exit, any other key to continue '))
    if(n=='x'){
      print('Function ended on user input')
      return(NA)
    }  
  }
  ensembl = biomaRt::useMart("ensembl",dataset=dataset_found[[1]]) #connection to ensembl database 
  snp <- biomaRt::useMart("ENSEMBL_MART_SNP", dataset=snp_found[[1]]) 
  
  ensemblOutput=list("ensembl"=ensembl, "snp"=snp)
  return(ensemblOutput)
}

genoToMarker = function(gds_obj, selectMarker){
  
  snp_name=substr(selectMarker, 1, nchar(selectMarker)-3)
  snp_id=which(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.rs.id"))==snp_name)
  geno=SNPRelate::snpgdsGetGeno(gds_obj, snp.id=snp_id)
  allele_comb=substr(selectMarker, nchar(selectMarker)-1, nchar(selectMarker))
  if(substr(allele_comb,1,1)!=substr(allele_comb,2,2)){ #Heterozygote
    geno[geno == 2] <- 0 
  } else { #Homozygote
    alleles=gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.allele"), start=c(snp_id), count=c(1)) # Ex 'A/G'
    ref_allele=substr(alleles,1,1) #First allele is reference allele
    if(substr(allele_comb,1,1)==ref_allele){ #Homozygote, of ref allele 0-1 => 0, 2 => 1
      geno[geno == 1] <- 0 
      geno[geno == 2] <- 1 
    } else { #Homozygote, of other allele, 1-2 => 0 , 0 => 1
      geno[geno == 1] <- 2
      geno[geno == 0] <- 1
      geno[geno == 2] <- 0
    }
  }
  return(geno)
}



