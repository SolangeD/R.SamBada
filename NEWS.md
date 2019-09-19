# R.SamBada 0.1.2

## Bug fixes
- prepareGeno: suppressed warnings "arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only"
- sambadaParallel: - bug fix when keepAllFiles=TRUE 
	- bug fix when wordDelim=,
	- bug fix when cores=1 to save storey file
	- removes SamBada error log that is created when SamBada's installation is tested
	- improved parameter checks
	- suppressed warnings "arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only"
- prepareOutput: bug fix when reading the file when popStr=FALSE
	- bug fix when opening storey file when wordDelim was different from space in sambadaParallel

## Update in the documentation
- better indication on which species can be handled in plotResultInteractive
- indication on how to install sambada in sambadaParallel and prepareGeno

# R.SamBada 0.1.1

## Bug fixes
- downloadSambada: bug fix for MacOS
- prepareGeno, sambadaParallel: error message when sambada is not installed
- plotResultInteractive: bug fix when click on Manhattan plot
## Update of the vignette and documentation
- Simplification of the vignette. Some details have been moved from the vignette to the documentation. Reference to files used in the vignette is simplified.
## Updates in the test dataset
- Update of the file 'uganda-subset-mol-Out-2.csv'
