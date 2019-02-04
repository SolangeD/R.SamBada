#' @title Download samBada
#' @description Downloads from GitHub the version of samBada that corresponds to your OS. Unzips the folder and adds the path to the binary folder to the environmental path variable. This operation is only valid for the current R session. You must run change_path for every new R session. Alternatively, you can manually edit your "PATH" environmental variable permanently on your OS so that it entails the path to the binaries folder of sambada (this procedure different for every OS).
#' @author Solange Gaillard
#' @param directory character The directory where sambada should be downloaded. If null, downloads in a (new) folder named sambada in the active directory.
#' @examples
#' \dontrun{
#' downloadSamada('D:/Sambada')
#' }
#' @export
downloadSambada = function(directory = NULL){
  if(is.null(directory)){
    directory=paste0(getwd(),'/sambada')
  }
  
	#Get OS name
	env=Sys.info()['sysname']
	#Create tempfile
	temp <- tempfile()
	#Download compressed folder from GitHub according to right OS
	if(env=='Windows'){
		utils::download.file("https://github.com/Sylvie/sambada/releases/download/v0.7.1/sambada-0.7.1-windows.zip",temp)
	  #Unzip the file
	  utils::unzip(temp, exdir = directory)
	  changePath(paste0(directory,'/binaries'))
	}
	else if(env=='Linux'){
		utils::download.file("https://github.com/Sylvie/sambada/releases/download/v0.7.1/sambada-0.7.1-ubuntu.tar.gz",temp)
	  #Untar the file => Needs twice?
	  utils::untar(temp, exdir = directory)
	  changePath(paste0(directory,'/sambada-0.7.1-ubuntu/binaries'))
	}
	else if(env=='MacOS'){
		utils::download.file("https://github.com/Sylvie/sambada/releases/download/v0.7.1/sambada-0.7.1-osx.tar.gz",temp)
	  #Untar the file
	  utils::untar(temp, exdir = directory)	
	  changePath(paste0(directory,'/sambada-0.7.1-osx/binaries'))
	}
	else{
		stop('Unknown operating system. Please download sambada manually from GitHub')
	}
	
	
	#Delete temp file
	unlink(temp)

	#Set current folder in PATH so that sambada can be called without specifying its path for the current session
	
}

#' @title Adds folder to the PATH environmental variable
#' @description Adds the \code{directory} path to the environmental PATH variable. This operation is only valid for the current R session. You must run change_path for every new R session. Alternatively, you can permantently edit your "PATH" environmental variable on your OS so that it entails the path to the binaries folder of sambada.
#' @author Solange Gaillard
#' @param directory character The path to samBada binaries folder 
changePath=function(directory){
  #check if directory already in path
  if(grepl(gsub('\\','/',directory, fixed=TRUE),gsub('\\','/',Sys.getenv('PATH'), fixed=TRUE))==FALSE){
    #Add directory to path
    env=Sys.info()['sysname']
    if(env=='Windows'){
      Sys.setenv(
        PATH = paste(
          Sys.getenv("PATH"), 
          directory, 
          sep = ";"
        )
      )
    } else {
      Sys.setenv(
        PATH = paste(
          Sys.getenv("PATH"), 
          directory, 
          sep = ":"
        )
      )      
    }
  
  }
}
