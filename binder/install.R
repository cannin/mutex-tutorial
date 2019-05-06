system("R CMD javareconf -e")

install.packages("devtools")
install.packages("stringr")

library(devtools)

install_cran("rmarkdown")
install_cran("rJava")

source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c("clusterProfiler", "GSEABase", "enrichplot"))

install_bitbucket("mil2041/netboxr")
install_github("BioPAX/paxtoolsr")

#source("https://gist.githubusercontent.com/cannin/6b8c68e7db19c4902459/raw/installPackages.R")
#installPackages("r-requirements.txt")
