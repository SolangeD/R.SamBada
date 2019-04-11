#' @title Download samBada
#' @description Downloads from GitHub the version of samBada that corresponds to your OS. Unzips the folder and adds the path to the binary folder to the environmental path variable. This operation is only valid for the current R session. You must run change_path for every new R session. Alternatively, you can manually edit your "PATH" environmental variable permanently on your OS so that it entails the path to the binaries folder of sambada (this procedure different for every OS).
#' @author Solange Duruz
#' @param directory character The directory where sambada should be downloaded. If null, downloads in a (new) folder named sambada in the active directory.
#' @examples
#' # Downloads SamBada to temporary folder (tempdir)
#' downloadSambada(tempdir())
#' @export
downloadSambada = function(directory = NULL){
  if(is.null(directory)){
    directory=paste0(getwd(),'/sambada')
  }
  continue=readline("This funcion will download Sambada in the specified or the current directory. Press any letter to continue or x to exit the process ")
  if(continue=='x'){
    return(NA)
  }
  
	#Get OS name
	env=Sys.info()['sysname']
	#Create tempfile
	temp <- tempfile()
	#Download compressed folder from GitHub according to right OS
	if(env=='Windows'){
		utils::download.file("https://github.com/Sylvie/sambada/releases/download/v0.8.0/sambada-0.8.0-windows.zip",temp)
	  #Unzip the file
	  utils::unzip(temp, exdir = directory)
	  changePath(paste0(directory,'/binaries'))
	}
	else if(env=='Linux'){
		utils::download.file("https://github.com/Sylvie/sambada/releases/download/v0.8.0/sambada-0.8.0-ubuntu.tar.gz",temp)
	  #Untar the file => Needs twice?
	  utils::untar(temp, exdir = directory)
	  changePath(paste0(directory,'/sambada-0.8.0-ubuntu/binaries'))
	}
	else if(env=='MacOS'){
		utils::download.file("https://github.com/Sylvie/sambada/releases/download/v0.8.0/sambada-0.8.0-osx.tar.gz",temp)
	  #Untar the file
	  utils::untar(temp, exdir = directory)	
	  changePath(paste0(directory,'/sambada-0.8.0-osx/binaries'))
	}
	else{
		stop('Unknown operating system. Please download sambada manually from GitHub')
	}
	
	
	#Delete temp file
	unlink(temp)

	#Set current folder in PATH so that sambada can be called without specifying its path for the current session
	
}

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
