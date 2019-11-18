# R.SamBada 0.1.2

## Bug fixes
- prepareGeno: suppressed warnings "arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only"
- createEnv: checks that the csv file has x and y column inside + a third column
- prepareEnv: suppress variables correlated above a threshold IN ABSOLUTE VALUE
- sambadaParallel: - bug fix when keepAllFiles=TRUE (files previously saved in temporary directory)
	- bug fix when wordDelim=,
	- bug fix when cores=1 to save storey file
	- removes SamBada error log that is created when SamBada's installation is tested
	- improved parameter checks
	- suppressed warnings "arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only"
- prepareOutput: bug fix when reading the file when popStr=FALSE
	- bug fix when opening storey file when wordDelim was different from space in sambadaParallel
- plotResultInteractive: added axis-label to the boxplot
- plotMap: - bug fix when rasterName was not null 
	- bug fix when environmental variable was not a bioclimatic variable

## Update in the documentation
- better explanation on plots shown with interactiveChecks=TRUE  or if file given in relative path in prepareOutput 
- better indication on which species can be handled and of the details of plots shown in plotResultInteractive
- indication on how to install sambada in sambadaParallel and prepareGeno
- indication of which variable is deleted when two variables are too highly correlated

# R.SamBada 0.1.1

## Bug fixes
- downloadSambada: bug fix for MacOS
- prepareGeno, sambadaParallel: error message when sambada is not installed
- plotResultInteractive: bug fix when click on Manhattan plot
## Update of the vignette and documentation
- Simplification of the vignette. Some details have been moved from the vignette to the documentation. Reference to files used in the vignette is simplified.
## Updates in the test dataset
- Update of the file 'uganda-subset-mol-Out-2.csv'
