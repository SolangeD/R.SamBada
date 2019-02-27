RSambada R-package

* This R-package is designed to 
	-Install sambada, 
	-Preprocess data to comply to sambada standard, 
	-Run sambada on multiple cores
	-Post-process the results and show interactive plots
	
* Sambada https://github.com/Sylvie/sambada is a landscape genomic software designed to compute correlations between genomic and environmental variables to detect signatures of selection

* To install this R-package
install.packages('devtools')
require('devtools')
install_github('SolangeD/RSambada', build_vignettes=TRUE)

* to start with 
Read RSambada-manual.pdf
Get documentation from sambada https://github.com/Sylvie/sambada 
Once installed, read the vignette with
	vignette('RSambada')