#' @title Prepare output (usefull for all postprocessing analysis)
#' @description Read sambada's output and prepare it by retrieving the snp position and chromosome (usefull for plotting manhattan)
#' @param sambadaname char The name of the genofile without extension name given to sambada (or outputfile of sambada without the ending -Out-Dim.csv)
#' @param dimMax integer The maximum number of dimension given in sambada
#' @param gdsFile char Name of the gds file associated with sambada's input file. If null, will try with \code{sambadaname}.gds
#' @param popStr logical Indicates whether sambada was run using the POPSTRVAR parameter (i.e. population structure was taken into account). Default false
#' @param nrows integer Specifies the number of line to read from the input file. Useful if saveType END ALL was used in sambada and that the number of models run is large so that the reading and processing is too slow. The saveType END parameter ensures that most significant models are located at the top of the file.
#' @param interactiveChecks logical
#' @return a list containing a) $sambadaOutput a matrix containing the output from sambada with 3 additional column: corresponding snp, chromosome and position of the marker b) chrSNPNum The total number of SNPs in each chromosme c) $chrMaxPos The highest position found in each chromosome
#' @examples
#' prepare_output('myFile',1)
#' @export
prepareOutput = function(sambadaname, dimMax, gdsFile=NULL, popStr=FALSE, nrows=NULL, interactiveChecks=TRUE){
  #setwd('/home/lasigadmin/R/test2/test_data')
  #
  #SNPRelate::snpgdsBED2GDS(bed.fn='ADAPTmap2.bed',fam.fn='ADAPTmap2.fam',bim.fn='ADAPTmap2.bim',out.gdsfn='ADAPTmap2.gds',compress.annotation = "")
  #setwd("D:/sambada/test2/test_data")
  #sambadaname='ADAPTmap2_recode'
  #gdsFile=NULL
  #dimMax=1
  
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
  for(i in 1:dimMax){
    if (!file.exists(paste0(sambadaname,'-Out-',i,'.csv'))) stop(paste0('File ',sambadaname,'-Out-',i,'.csv not found. Check input sambadaname and dimMax'))
  }
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
  # Check that it looks correct!
  
  # blabla
  
  m=sum(storeyTot[3,2:ncol(storeyTot)]) #Number of models
  #p-value
  if(popStr==TRUE){
    pvalueG=pchisq(output$GscorePop, 1, lower.tail=F)
    pvalueW=pchisq(output$WaldScorePop, 1, lower.tail=F)    
  } else {
    pvalueG=pchisq(output$Gscore, 1, lower.tail=F)
    pvalueW=pchisq(output$WaldScore, 1, lower.tail=F)    
  }

  #Qvalue
  #estimate pi0
  #pi0 function
  pi_lambdaG=cumsum(t(as.vector(storeyTot[3,2:ncol(storeyTot)])))/(m*(1-t(as.vector(storeyTot[1,2:ncol(storeyTot)]))))
  splineG=splinefun(t(storeyTot[1,2:ncol(storeyTot)]), pi_lambdaG) #in qvalue package, use spline.smooth (and predict)
  pi0G=splineG(1)
  if(interactiveChecks==TRUE){
    ####plot histo + estimated pi0
    
  }
  pi_lambdaW=cumsum(t(as.vector(storeyTot[5,2:ncol(storeyTot)])))/(m*(1-t(as.vector(storeyTot[1,2:ncol(storeyTot)]))))
  splineW=splinefun(storeyTot[1,2:ncol(storeyTot)], pi_lambdaW)
  pi0W=splineW(1)
  if(interactiveChecks==TRUE){
    ####plot histo + estimated pi0
    
  }
  #Calculate qvalue
  i=length(pvalueG):1
  qvalueG=pi0G * pmin(1, cummin(pvalueG[length(pvalueG):1] * m /i ))[length(pvalueG):1] #from qvalue package
  qvalueW=pi0W * pmin(1, cummin(pvalueW[order(pvalueW, decreasing=TRUE)] * m /i ))[order(pvalueW, decreasing=TRUE)]
  #Add it to the output vector
  output=cbind(output, "pvalueG"=pvalueG, "pvalueW"=pvalueW, "qvalueG"=qvalueG, "qvalueW"=qvalueW)
  
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
  
  # Merge both information
  # add q.value
  #Merge acts like join in SQL
  output2=merge(cbind.data.frame('snp'=substr(output[, 1]$Marker, start = 1, stop = nchar(c(output[,1])$Marker) - 3), output), map2[,c('chr','snp','pos')], by='snp', sort=FALSE, all.x=TRUE)
  
  # Number of SNP by chromosome
  numSNP=table(all_chr)
  # Reorder chromosome: number first in numerical order, letter last in alphabetical order
  numSNP=suppressWarnings(numSNP[c(sort(as.integer(names(numSNP))), names(numSNP)[is.na(as.integer(names(numSNP)))])])
  
  
  # Max position by chromosome
  maxPos=aggregate(all_pos, by=list(all_chr), FUN='max')
  maxPos=suppressWarnings(maxPos[match(c(sort(as.numeric(maxPos$Group.1)), maxPos$Group.1[is.na(as.numeric(maxPos$Group.1))]), maxPos$Group.1),])
  
  #Final output
  finalOutput=list("sambadaOutput"=output2, "chrSNPNum"=numSNP, "chrMaxPos"=maxPos)
  return(finalOutput)
  
}



#' @title Interactive plotting of results
#' @description Plots the manhattan plot for a given environmental variable. The plot is interactive and a map of the distribution of the marker can be retrieved as well as nearby genes listed in Ensembl.
#' @author Solange Gaillard
#' @param preparedOutput char The prepared output list from prepare_output function
#' @param varEnv char The name of the environmental variable one wish to study (as in the header of \code{envFIle})
#' @param envFile char The file containing the input environmental variable of sambada. 
#' @param species char The abbreviated latin name of the species without capitals nor punctuation (e.g. btaurus, chircus,...). Can be set to null if species not present in ensembl database
#' @param pass integer Number of BP around a SNP in which to look for an annotation in Ensembl. Set to null if species is null
#' @param x char The name of the column corresponding to the x-coordinate in the envFile. Can be set to null if unknown, in this case the maps will not be available
#' @param y char The name of the column corresponding to the y-coordinate in the env file. Can be set to null if x is null.
#' @param valueName char Name of the p- or q-value one wish to plot the manhattan on. This can be either pvalueG, pvalueW, qvalueG, qvalueW for G- or Waldscore respectively.
#' @param chromo char/integer Name or vector of name of the chromosome to investigate. If all is chosen (default), all numerical chromosome will be mapped. If your sambada output is large (typically if you are working with more than 50K genomic file), you should probably map a subset of your dataset (e.g. chr=1)
#' @param gdsFile char The GDS file created in the preprocessing of sambada. If null, will try with envFile(without -env.csv) and .gds
#' @param IDCol char The name of the column in envFile corresponding to the ID of the individuall. If provided, hover on the output map will give the id of the animal
#' @param popStrCol char The name or vector of name of column(s) in envFile describing population structure. If provided, additional layers on the map will be available reprenting population structure.
#' @return None 
#' @examples
#' plotResultInteractive('myFile','chircus',1,'Longitude','Latitude','bio1',c('1','2'))
#' @export
plotResultInteractive = function(preparedOutput, varEnv, envFile,species=NULL, pass=NULL,x=NULL,y=NULL,  valueName='pvalueG',chromo='all',gdsFile=NULL, IDCol=NULL, popstrCol=NULL){

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
  if(sum(popstrCol %in% names(envDataTest))!=length(popstrCol))stop('All elements of popstrCol-argument should be part of the header of the envFile')
  #valueName
  if(!(valueName %in% names(preparedOutput$sambadaOutput)))stop('valueName should be a component of preparedOutput$sambadaOutput. Use the result of the function prepare_output as preparedOutput')
  #chromo
  # if('all' %in% chromo){
  # }  


  # Test if required libraries are installed
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("Package \"SNPRelate\" needed for this function to work. Please install it.", call. = FALSE)
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
  
  #Define and open GDS file
  if(is.null(gdsFile)){
    gdsFile = paste0(gsub('-env.csv','',envFile),'.gds')
    if(!file.exists(gdsFile)){
      stop("A gds file is needed for this function to work. Specify the name in the input of the function if it is already created, or create it with the package SNPRelate or prepare_geno from this package")
    }
  }
  gds_obj=SNPRelate::snpgdsOpen(gdsFile)
  #on.exit(SNPRelate::snpgdsClose(gds_obj))
  
  #Read envFile and output
  envData = read.csv(envFile, header=TRUE, sep=" ")
  sambadaOutput = preparedOutput$sambadaOutput
  chrSNPNum = preparedOutput$chrSNPNum
  chrMaxPos = preparedOutput$chrMaxPos
    
  #Connection to ensembl database

  ensemblOutput = ensembl_connection(species, TRUE)
  snp = ensemblOutput$snp
  ensembl = ensemblOutput$ensembl
  
  #Prepare Manhattan (to be changed)
  subset=sambadaOutput[which(sambadaOutput[,'Env_1']==varEnv),]
  if(nrow(subset)==0) stop('No records found in sambada output corresponding to the chosen environmental variable (varEnv argument)')
  if(length(chromo)>1){
    subset=subset[with(subset, which(chr %in% chromo)),]
  } else if(chromo != 'all' ){
    subset=subset[with(subset, which(chr %in% chromo)),]
  }
  if(nrow(subset)==0) stop('No record found in sambada output corresponding to the chosen chromosome (chromo argument)')
  subset$pval=-log10(subset[,get(valueName)])
  
  prevPos=data.frame("chr"=chrMaxPos$Group.1, "maxPos"=cumsum(as.numeric(chrMaxPos$x))-chrMaxPos$x)
  prevPos=data.frame(prevPos, "chrPos"=rownames(prevPos))
  subset=merge(subset,prevPos,by='chr', sort=FALSE, all.x=TRUE)
  subset$color=colors()[as.integer(as.character(subset$chrPos))%%2+6] 

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
        #shiny::h4("Link to google search of description of found genes"),
        #shiny::htmlOutput("event4"),
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
    rhg_cols=c("#CCCC99","#999966")
    p <- ggplot2::ggplot(data=subset, ggplot2::aes(x=maxPos+pos, y=pval, colour=color, label=snp, text=pos), showlegend=FALSE) 
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
    # Problem when serveral chromo and chromosomes not side by side


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
          selectBP=subset[which(subset$pos+subset$maxPos==f$x),'pos']
          selectBP=selectBP[1]
          selectCHR=subset[which(subset$pos+subset$maxPos==f$x),'chr']
          selectCHR=selectCHR[1]
          #Query ensembl database to get nearby genes
          #c=tryCatch({biomaRt::getBM(attributes=c('chromosome_name','ensembl_gene_id','wikigene_name','start_position','end_position','description'), filters=c("chromosome_name","start","end"), values=list(selectCHR,selectBP-pass,selectBP+pass), mart=ensembl)}, error=function(e){"no gene found!"})
          #c
          biomaRt::getBM(attributes=c('chromosome_name','ensembl_gene_id','wikigene_name','start_position','end_position','description'), filters=c("chromosome_name","start","end"), values=list(selectCHR,selectBP-pass,selectBP+pass), mart=ensembl)
        }
        
      })
      # output$event4 <- shiny::renderPrint({
      #   f <- plotly::event_data("plotly_click")
      #   if (is.null(f)) {
      #     "Select a point!" 
      #   }else {
      #     selectBP=subset[which(subset$pos+subset$maxPos==f$x),'pos']
      #     selectBP=selectBP[1]
      #     selectCHR=subset[which(subset$pos+subset$maxPos==f$x),'chr']
      #     selectCHR=selectCHR[1]
      #     #Query ensembl database to get nearby genes
      #     #c=tryCatch({biomaRt::getBM(attributes=c('chromosome_name','ensembl_gene_id','wikigene_name','start_position','end_position','description'), filters=c("chromosome_name","start","end"), values=list(selectCHR,selectBP-pass,selectBP+pass), mart=ensembl)}, error=function(e){"no gene found!"})
      #     #c
      #     #biomaRt::getBM(attributes=c('chromosome_name','ensembl_gene_id','wikigene_name','start_position','end_position','description'), filters=c("chromosome_name","start","end"), values=list(selectCHR,selectBP-pass,selectBP+pass), mart=ensembl)
      #     test=biomaRt::getBM(attributes=c('description'), filters=c("chromosome_name","start","end"), values=list(selectCHR,selectBP-pass,selectBP+pass), mart=ensembl)
      #     test=substr(test,1,regexpr(pattern ='\\[',test)[1]-1)
      #     test=gsub(" ", "+", test)
      #     #renderUI({tagList("URL link:", url)})
      #     #renderUI({tags$a(href = 'www.google.com/search?q=carbohydrate+sulfotransferase','www.google.com/search?q=carbohydrate+sulfotransferase')})
      #     tags$div(class="header", checked=NA,
      #              tags$a(href=paste0('www.google.com/search?q=',test), test))
      #   }
      #   
      # })
    }
    
    #When clicked: give position of SNP
    output$event2<- shiny::renderPrint({
      f <- plotly::event_data("plotly_click")
      if (is.null(f)) {
        "Select a point!" 
      }else {
        selectBP=subset[which(subset$pos+subset$maxPos==f$x),'pos']
        selectBP=selectBP[1]
        selectCHR=subset[which(subset$pos+subset$maxPos==f$x),'chr']
        selectCHR=selectCHR[1]
        selectSNP=subset[which(subset$pos+subset$maxPos==f$x),'snp']
        selectSNP=selectSNP[1]
        infoOutput=data.frame('snp'=selectSNP, 'chr'=selectCHR, 'BP'=selectBP)
        infoOutput
      }
    })
    
    #When clicked: get all marker in output of same SNP
    output$event3<- shiny::renderPrint({
      f <- plotly::event_data("plotly_click")
      if (is.null(f)) {
        "Select a point!" 
      }else {
        selectSNP=subset[which(subset$pos+subset$maxPos==f$x),'snp']
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
          popCol=envData[,popstrCol]
          ID=envData[,IDCol]
          # Get marker and snp
          selectSNP=subset[which(subset$pos+subset$maxPos==g$x),'snp']
          selectSNP=selectSNP[1]
          selectMarker=subset[which(subset$pos+subset$maxPos==g$x),'Marker']
          selectMarker=selectMarker[1]$Marker
          #Retrieve genotype
          snp_id=which(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.rs.id"))==selectSNP$snp)
          
          # #Check SNP in LD in the area, too be implemented
          # snp1 = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "genotype"), start=c(1,snp_id), count=c(-1,1))
          # #loop on nearby snps, break when ld too small
          # snp2 = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "genotype"), start=c(1,snp_id+1), count=c(-1,1))
          # SNPRelate::snpgdsLDpair(snp1, snp2, method = "corr")
          
          geno=SNPRelate::snpgdsGetGeno(gds_obj, snp.id=snp_id)
          pres=genoToMarker(gds_obj, selectMarker)
          xy=data.frame(x,y,varenv,ID, geno, pres, popCol)
         
          graph <- plotly::plot_ly()
          graph <- plotly::add_trace(graph,data=xy, x=x,y=y, type='scatter',mode='markers', color=pres, marker=list(showscale=FALSE), name='marker', colors=pal1, text=ID,hoverinfo = c("color"))
          graph <- plotly::hide_colorbar(graph)  
          graph <- plotly::add_markers(graph,data=xy, inherit=FALSE, mode='markers', x=~x,y=~y, name='varenv', marker=list(color=~varenv,colorscale = 'YlOrRd', showscale=TRUE, colorbar=list(len=0.3, title='varenv',y=1)))
          if(length(popstrCol)==1){
            graph <- plotly::add_markers(graph,data=xy, inherit=FALSE, mode='markers',x=~x,y=~y, name='pop1', marker=list(color=~popCol,colorscale = 'Greens', showscale=TRUE, colorbar=list(len=0.3, title='pops', y=0.7)))
          } else {  
            minpop=min(xy[,popstrCol])
            maxpop=max(xy[,popstrCol])
            graph <- plotly::add_markers(graph,data=xy, inherit=FALSE, mode='markers',x=~x,y=~y, name='pop1', marker=list(color=~get(popstrCol[1]),colorscale = 'Greens', cmin=minpop, cmax=maxpop, showscale=TRUE, colorbar=list(len=0.2, title='pops', y=0.75)))
            for(p in 2:length(popstrCol)){
              graph <- plotly::add_markers(graph,data=xy, inherit=FALSE, mode='markers',x=~x,y=~y, name=paste0('pop',p), marker=list(color=~get(popstrCol[p]),colorscale = 'Greens',cmin=minpop, cmax=maxpop,  showscale=FALSE))
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
          popCol=envData[,popstrCol]
          # Get marker and snp
          selectSNP=subset[which(subset$pos+subset$maxPos==g$x),'snp']
          selectSNP=selectSNP[1]
          selectMarker=subset[which(subset$pos+subset$maxPos==g$x),'Marker']
          selectMarker=selectMarker[1]$Marker
          #Retrieve genotype
          snp_id=which(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj, "snp.rs.id"))==selectSNP$snp)
          pres=genoToMarker(gds_obj, selectMarker)
          xy=data.frame(varenv, pres, popCol)
          boxplot(varenv~pres, data=xy)
        }
      })
      
      session$onSessionEnded(function() {
        stopApp()
        SNPRelate::snpgdsClose(gds_obj)
      })
  }
  
  shiny::shinyApp(ui, server)
  #shiny::stopApp()
}


#' @title Plotting of maps 
#' @description Plots several kinds of maps (environmental variable distribution, population structure, marker absence or presence, autocorrelation of marker). Unlike \code{plotResultInteractive}, the resulting maps are non-interactive. The function can handle several marker/variables at once and create separate outputfiles.
#' @author Solange Gaillard
#' @param envFile char The file containing the input environmental variable of sambada. 
#' @param x char The name of the column corresponding to the x-coordinate in the envFile. Can be set to null if unknown, in this case the maps will not be available
#' @param y char The name of the column corresponding to the y-coordinate in the env file. Can be set to null if x is null.
#' @param gdsFile char The GDS file created in the preprocessing of sambada. If null, will try with envFile(without -env.csv) and .gds
#' @param popStrCol char The name or vector of name of column(s) in envFile describing population structure. If provided, additional layers on the map will be available reprenting population structure.
#' @param locationProj integer EPSG code of the geographical projection in the envFile
#' @param markerName name of the marker to be plotter if \code{mapType} is 'marker' or 'AS'. \code{markerName} can be found in preparedOutput$sambadaOutput[,''] where preparedOutput would be the result of the function \code{prepareOutput}
#' @param mapType char A string or vector of string containing one or several of 'marker' (presence/absence of marker), 'env' (envrionmental variable distribution), 'popStr' (appartnance to a population in pie charts), 'AS' (autocorrelation of the marker). Note that the background of all maps, if found, will be the raster of the environmental variable. Thus the 'env' \code{mapType} is preferred when no raster is provided. For the 'AS' type, it is calculated on the fly for the markers provided and not the one possibly calculated by sambada.
#' @param varEnvName char Name of the environmental variable. If a raster of the variable is located in your working directory, you can provide \code{varEnvName} even for \code{mapType} such as 'marker' or 'AS'. The function will scan the folder of your working directory for raster with the same name as \code{varEnvName} (and commonly used extension for raster) and put it as background.
#' @param SAMethod char If \code{mapType} contains 'AS', then you must specify the method for setting the weights of neighbours. Can be one of 'knn' (k-nearest neighbours) or 'distance' 
#' @param SAMThreshold char If \code{mapType} contains 'AS' and \code{SAMethod} id 'knn' then the number of neighbours. If \code{SAThreshold} is 'distance' then the distance in map-unit (unless you use a spherical projection (latitude/longitude), in which case you should use km)
#' @param saveType char One of NULL, 'png' or 'pdf'. To be implemented... If NULL is set, the maps will be shown in the R plotting window. Otherwise, it will be saved in the specified format in your working directory.
#' @param rasterName char If a raster file with the environmental variable distribution exists with a different name than \code{varEnvName}, provide it here (including extension)
#' @param simultaneous boolean If TRUE and \code{mapType} contains several kinds of maps, all maps corresponding to the same marker will be plotted on the same window. The resulting maps can be very small.
#' @return None 
#' @examples
#' plotMap('EnvFile.csv','longitude','latitude', locationProj=4326,  popStrCol='pop1', gdsFile='GDSFile.gds', markerName='ARS-BFGL-NGS-106879_AA', mapType=c('marker'), varEnvName='bio1')
#' @export
plotMap = function(envFile, x, y, locationProj,  popStrCol, gdsFile, markerName, mapType, varEnvName, SAMethod=NULL, SAThreshold=NULL, saveType=NULL, rasterName=NULL, simultaneous=FALSE){

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
    stop("A gds file is needed for this function to work if mapType contains 'marker' or 'AS'. Specify the name in the input of the function if it is already created, or create it with the package SNPRelate or prepare_geno from this package")
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
  if(sum(mapType %in% c('marker','env','AS','popStr'))!=length(mapType)) stop("mapType should be one, or several of, 'marker','env','AS','popStr'")
  if(typeof(varEnvName)!='character') stop('varEnvName argument should be a character')
  if(!is.null(SAMethod)) if(!(SAMethod %in% c('knn','distance'))) stop('SAMethod should be one of knn or distance')
  if(!is.null(SAThreshold)) if(typeof(SAThreshold)!='double') stop ('SAThreshold should be numeric')
  if('AS' %in% mapType & is.null(SAMethod)) stop('SAMethod must be not null if mapType include AS')
  if('AS' %in% mapType & is.null(SAThreshold)) stop('SAThreshold must be not null if mapType include AS')
  if('AS' %in% mapType) if(SAMethod=='knn' & SAThreshold%%1!=0) stop ("SAThreshold should be an integer if SAMethod is 'knn'")
  if(!is.null(saveType)) if(saveType %in% c('png','pdf')) stop("saveType should be one of NULL,'png','pdf'")
  if(!is.null(rasterName)) if(typeof(rasterName)!='character') stop('rasterName argument should be a character')
  if(!is.null(rasterName)) if(!file.exists(rasterName)) stop('rasterName file not found')
  
  #Open Env file
  envData=read.csv(envFile, header=TRUE, sep=" ")
  envData2=envData
  sp::coordinates(envData) = c(x,y)
  sp::proj4string(envData)=sp::CRS(paste0("+init=epsg:",locationProj))
  
  #Define layout
  if(length(mapType)>1 & simultaneous==TRUE){
    numrow=2*length(mapType)
    matrixLayout=rep(c(1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,3),length(mapType))+rep(seq(from=0, by=3, length=length(mapType)),each=20)
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
  # Est-ce que dev.size va marcher???
  size=0.1*max((max(envData@coords[,x])-min(envData@coords[,x]))/dev.size(units='cm')[1],(max(envData@coords[,y])-min(envData@coords[,y]))/dev.size(units='cm')[2])
  data_df=data.frame(size=rep(size,nrow(envData)),x=envData@coords[,x],y=envData@coords[,y])
  scattered_point=packcircles::circleRepelLayout(data_df,NULL,NULL,xysizecols=c('x','y','size'),sizetype='radius',maxiter = 50)
  sp::coordinates(scattered_point$layout)=c('x','y')

  
  #Draw plots
  for(m in 1:length(markerName)){
    #Get Marker info
    if('marker' %in%  mapType | 'AS' %in% mapType ){
      #Open gds File
      pres=genoToMarker(gds_obj, markerName)
      
    }
    #Try to find corresponding raster
    if(m>1)
    allowedExtension=c('bil','tif')
    if(exists('rasterName')) {
      rm(rasterName)
    }
    if(regexpr('bio',varEnvName[m])>0){
      if(file.exists(paste0('wc0.5/',varEnvName[m],'.tif'))){
        rasterName=paste0('wc0.5/',varEnvName[m],'.tif')
      }
    } else if (regexpr('bio',varEnvName[m])>0) {
      if(file.exists(paste0('srtm/',varEnvName[m],'.tif'))){
        rasterName=paste0('srtm/',varEnvName[m],'.tif')
      }
    } else {
      for(aE in 1:length(allowedExtension))
        if(file.exists(paste0(varEnvName[m],'.',allowedExtension(aE)))){
          rasterName=paste0(varEnvName[m],'.',allowedExtension(aE))
          break
        }
    }
    #Open raster
    if(exists('rasterName')){
      raster=raster::raster(rasterName)
      #Get real raster data
      raster_df=as.data.frame(raster::sampleRegular(raster, size=1e5, asRaster=FALSE), xy=TRUE)
    }
    
    for(t in 1:length(mapType)){
      par(mar=c(2,2,2,2), xpd=FALSE)
      #if(t==1){
      #Draw background
      if(exists('rasterName')){
        #If raster found, put it as background
        #Attention mettre les coordonnées de scattered point ou envData???
        raster::image(raster, asp=1, maxpixels=10000000000,  col=terrain.colors(100),xlim = c(min(envData@coords[,x]), max(envData@coords[,x])), ylim = c(min(envData@coords[,y]), max(envData@coords[,y])))
        
      }else {
        #If raster not found, put countries as background
        country=data('wrld_simpl', package='maptools')
        raster::plot(wrld_simpl,xlim=c(min(sp::coordinates(envData)[,x]),max(sp::coordinates(envData)[,y])),ylim=c(min(sp::coordinates(envData)[,x]),max(sp::coordinates(envData)[,y])))
      }
      
      #Draw lines between original location and scattered one (if not scattered, the lines will be masked by the points)
      for(i in 1:nrow(envData)){
        lines(c(envData@coords[i,x],scattered_point$layout@coords[i,'x']),c(envData@coords[i,y],scattered_point$layout@coords[i,'y']), col='antiquewhite3')
      }
      
      #Draw points
      if(mapType[t]=='marker'){
        raster::plot(scattered_point$layout,col=colors()[pres*19+5],pch=16,add=TRUE)
      } else if(mapType[t]=='env'){
        #Define color palette
        new.pal=colorRampPalette(c("yellow", "orange","red"))( 100 )
        raster::plot(scattered_point$layout, col=new.pal[round((envData@data[,varEnvName]-min(envData@data[,varEnvName]))/(max(envData@data[,varEnvName])-min(envData@data[,varEnvName]))*100)],pch=20,add=TRUE)
      } else if(mapType[t]=='popStr'){
        if(length(popStrCol)>1){
          for(i in 1:nrow(scattered_point$layout)){ 
            mapplots::add.pie(z=as.double(abs(1/envData@data[i,popStrCol])),x=scattered_point$layout@coords[i,'x'], scattered_point$layout@coords[i,'y'], labels=NA, col=terrain.colors(3), radius=size*2)
          }
        } else {
          new.pal=colorRampPalette(c("white","black"))( 100 )
          raster::plot(scattered_point$layout, col=new.pal[round((envData@data[,popStrCol]-min(envData@data[,popStrCol]))/(max(envData@data[,popStrCol])-min(envData@data[,popStrCol]))*100)],pch=16,add=TRUE)
        }
      
      } else if(mapType[t]=='AS'){
        ### autocorrelation
        #Calculate autocorrelation
        #knearneigh, dnearneigh
        if(method=='knn'){
          knn=spdep::knearneigh(envData,5)
          nb=spdep::knn2nb(knn)
        } else if (method=='distance'){
          nb=spdep::dnearneigh(envData,0,20)
        }
        nblist=spdep::nb2listw(nb, zero.policy=TRUE) 
        #iglobal=spdep::moran.test(as.vector(pres),nblist, na.action=na.omit)
        ilocal=spdep::localmoran(as.vector(pres),nblist, zero.policy=TRUE,na.action=na.exclude)
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
        axis(4, at=(pretty(raster_df[,varEnvName])[2:(length(pretty(raster_df[,varEnvName]))-1)]-min(raster_df[,varEnvName], na.rm=TRUE))/(max(raster_df[,varEnvName], na.rm=TRUE)-min(raster_df[,varEnvName], na.rm=TRUE))*100, label=pretty(raster_df[,varEnvName])[2:(length(pretty(raster_df[,varEnvName]))-1)])
        text(1,107,'Raster')
      } else {
        if(mapType[t]=='env'){
          # Point legend
          par(mar=c(2,1,3,2), xpd=NA)
          image(1, 1:100, t(seq_along(1:100)), col=new.pal, axes=FALSE, xlab="", ylab="")
          axis(4, at=(pretty(envData@data[,varEnvName])[2:(length(pretty(envData@data[,varEnvName]))-1)]-min(envData@data[,varEnvName], na.rm=TRUE))/(max(envData@data[,varEnvName], na.rm=TRUE)-min(envData@data[,varEnvName], na.rm=TRUE))*100, label=pretty(envData@data[,varEnvName])[2:(length(pretty(envData@data[,varEnvName]))-1)])
          text(1,107,'Points')
          plot.new()
          next
        }
      }
      if(mapType[t]=='popStr'){
        # Point legend
        #par(mar=c(2,1,3,2), xpd=NA)
        if(length(popStrCol)>1){
          points(rep(0,length(popStrCol)),seq(from=-20, by=-10, length.out=length(popStrCol)),pch=19, col=terrain.colors(length(popStrCol)))
          text(rep(0.1,length(popStrCol)),seq(from=-20, by=-10, length.out=length(popStrCol)),popStrCol, pos=4)
          text(1,-10,'Population')          
        } else {
          par(mar=c(2,1,3,2), xpd=NA)
          pop.pal=colorRampPalette(c("white", "black"))( 100 )
          image(1, 1:100, t(seq_along(1:100)), col=pop.pal, axes=FALSE , xlab="", ylab="")
          axis(4, at=(pretty(envData@data[,popStrCol])[2:(length(pretty(envData@data[,popStrCol]))-1)]-min(envData@data[,popStrCol], na.rm=TRUE))/(max(envData@data[,popStrCol], na.rm=TRUE)-min(envData@data[,popStrCol], na.rm=TRUE))*100, label=pretty(envData@data[,popStrCol])[2:(length(pretty(envData@data[,popStrCol]))-1)])
          text(1,107,'Poulation')
        }

      } else if(mapType[t]=='AS'){
        # Point legend
        #par(mar=c(2,1,3,2), xpd=NA)
        plot.new()
        points(rep(0,10),c(seq(from=1, by=-0.1, length.out=5),seq(from=0.4, by=-0.1, length.out=5)),pch=c(19,19,19,19,19,21,21,21,21,21), col=c('blue4','lightblue3','wheat','salmon','red3','blue4','lightblue3','wheat','salmon','red3'), bg='gray')
        text(rep(0.1,10),c(seq(from=1, by=-0.1, length.out=5),seq(from=0.4, by=-0.1, length.out=5)),c('< -0.5','-0.5 - -0.1','-0.1 - 0.1','0.1 - 0.5','> 0.5'), pos=4)
        text(0,1.1,'I (signif)', pos=4)
        text(0, 0.5,'I (not signif)', pos=4)
      } else if(mapType[t]=='marker'){
        plot.new()
        points(c(0,0),c(1,0.8),pch=19, col=c(colors()[5],colors()[24]))
        text(c(0.1,0.1),c(1,0.8),c('Absent','Present'),pos=4)
      }
    #}
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



