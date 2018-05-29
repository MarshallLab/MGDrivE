# Mosquito Gene Drive Explorer

**MGDrivE** is a model designed to be a reliable testbed where various gene drive interventions for mosquito-borne diseases control. It is being developed to accommodate the use of various mosquito-specific gene drive systems within a population dynamics framework that allows migration of individuals between patches in landscape.

For more information take a look at our <a href="https://marshalllab.github.io/MGDrivE/">website</a>; where you'll find a thorough description, documentation and examples.


[![Demo](https://marshalllab.github.io/MGDrivE/images/crispr.jpg)](https://www.youtube.com/watch?time_continue=3&v=sZXuUtToszw)
_Click the image to watch a video description._

<hr>


## Installing

**MGDrivE** depends upon the following third-party packages:

* [RCPP](https://cran.r-project.org/web/packages/Rcpp/index.html)
* [RCPPArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [R6](https://cran.r-project.org/web/packages/R6/index.html)
* [Rdpack](https://cran.r-project.org/web/packages/Rdpack/index.html)
* [roxygen2](https://cran.r-project.org/web/packages/roxygen2/index.html)

which can be installed with the following command (in the _R_ console):

```R
install.packages(c("Rcpp", "RcppArmadillo","data.table","R6","Rdpack","roxygen2"))
```

Additionally **Windows** users need to install [RTools](https://cran.r-project.org/bin/windows/Rtools/).

### Github

To install the package directly from this repository make sure you have [devtools](https://cran.r-project.org/web/packages/devtools/index.html) installed:

```R
install.packages("devtools")
```

Having installed "devtools", run the following commands in the _R_ terminal:

```R
devtools::install_github("MarshallLab/MGDrivE",subdir="MGDrivE")
```

And, finally, run the following command to test if the package installed correctly:

```R
library(MGDrivE)
```

### Download

Alternatively, __MGDrivE__ can be installed by downloading the source code and following these instructions:

1. Download and unzip the current _tar_ or _zip_ from our repository.
2. Open the file _./MGDrivE/MGDrivE.Rproj_ on [RStudio](https://www.rstudio.com/).
3. Click on the _Build>Document_ toolbar.
4. Click on the _Build>Install and Restart_ toolbar.

<hr>

## Authors

* Lead: <a href="https://chipdelmal.github.io/">Héctor M. Sánchez C.</a>,<br>
* Core Development: <a href="https://slwu89.github.io/">Sean L. Wu</a>,Jared Bennett<br>
* Spatial Analysis: Biyonka Liang, Sarafina Smith, Sabrina Wong<br>
* Movement Kernels: Partow Imani<br>
* PI: <a href="http://www.marshalllab.com/">John M. Marshall</a>

## Projects

This model is being supported and developed as part of the following projects:

* [UCI Malaria Initiative](https://news.uci.edu/7517/05/08/uci-establishes-malaria-initiative-to-fight-deadly-disease-in-africa/)
* [DARPA: Safe Genes](https://www.darpa.mil/program/safe-genes)

## Research Centers

<img src="https://marshalllab.github.io/MGDrivE/images/berkeley.jpg" height="40px" align="middle"><img src="https://marshalllab.github.io/MGDrivE/images/UCI.png" height="40px" align="middle"><img src="https://marshalllab.github.io/MGDrivE/images/UCD.jpg" height="40px" align="middle"><img src="https://marshalllab.github.io/MGDrivE/images/UCSD.png" height="40px" align="middle"><img src="https://marshalllab.github.io/MGDrivE/images/UCLA.png" height="40px" align="middle"><img src="https://marshalllab.github.io/MGDrivE/images/JPL.png" height="40px" align="middle"><img src="https://marshalllab.github.io/MGDrivE/images/IGI.png" height="40px" align="middle"><img src="https://marshalllab.github.io/MGDrivE/images/DARPA.jpg" height="40px" align="middle">
