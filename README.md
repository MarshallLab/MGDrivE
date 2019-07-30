# [Mosquito Gene Drive Explorer](https://marshalllab.github.io/MGDrivE/)


**MGDrivE** is a model designed to be a reliable testbed where various gene drive interventions for mosquito-borne diseases control. It is being developed to accommodate the use of various mosquito-specific gene drive systems within a population dynamics framework that allows migration of individuals between patches in landscape.

For more information take a look at our <a href="https://marshalllab.github.io/MGDrivE/">website</a>; where you'll find a thorough description, documentation and examples.


[![Demo](https://marshalllab.github.io/MGDrivE/images/crispr.jpg)](https://www.youtube.com/watch?time_continue=3&v=sZXuUtToszw)
_Click the image to watch a video description._

<hr>


## Installation Instructions

**MGDrivE** depends upon the following third-party packages: [RCPP](https://cran.r-project.org/web/packages/Rcpp/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [R6](https://cran.r-project.org/web/packages/R6/index.html), and [Rdpack](https://cran.r-project.org/web/packages/Rdpack/index.html). The documentation and vignettes depend on [roxygen2](https://cran.r-project.org/web/packages/roxygen2/index.html), [knitr](https://cran.r-project.org/web/packages/knitr/index.html), and [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html); which can be installed with the following command (in the _R_ console):

```R
install.packages(c("Rcpp","data.table","R6","Rdpack","roxygen2","knitr","rmarkdown"))
```

Once these packages are installed, the package can be installed with:

```R
install.packages("MGDrivE")
```

and loaded:

```R
library(MGDrivE)
```

<hr>

## License

This software is freely available under the [GPL3 License](https://www.gnu.org/licenses/gpl-3.0.en.html).


## Projects

This model is being supported and developed as part of the following research projects: [UCI Malaria Initiative](https://news.uci.edu/7517/05/08/uci-establishes-malaria-initiative-to-fight-deadly-disease-in-africa/), [DARPA: Safe Genes](https://www.darpa.mil/program/safe-genes), [IGI](https://innovativegenomics.org/)

<hr>

##  [MoNeT_MGDrivE](https://pypi.org/project/MoNeT-MGDrivE/)

To facilitate the data analysis of the results produced by [MGDrivE](https://marshalllab.github.io/MGDrivE/), we provide a python package [MoNeT_MGDrivE](https://pypi.org/project/MoNeT-MGDrivE/) installable through pip.

<img src="https://marshalllab.github.io/MGDrivE/images/Homing_01Cb.png" align="middle">

<hr>

##  Funders and Collaborators

<img src="https://marshalllab.github.io/MGDrivE/images/berkeley.jpg" height="40px" align="middle"><img src="https://marshalllab.github.io/MGDrivE/images/UCD.jpg" height="40px" align="middle"> <img src="https://marshalllab.github.io/MGDrivE/images/UCI.png" height="40px" align="middle"> <img src="https://marshalllab.github.io/MGDrivE/images/UCLA.png" height="40px" align="middle"> <img src="https://marshalllab.github.io/MGDrivE/images/UCR.jpg" height="40px" align="middle"> <br><br> <img src="https://marshalllab.github.io/MGDrivE/images/UCSD.png" height="40px" align="middle"> <img src="https://marshalllab.github.io/MGDrivE/images/JPL.png" height="40px" align="middle"> <img src="https://marshalllab.github.io/MGDrivE/images/IGI.png" height="40px" align="middle"> <img src="https://marshalllab.github.io/MGDrivE/images/DARPA.jpg" height="40px" align="middle"> <img src="https://marshalllab.github.io/MGDrivE/images/nvidia.jpg" height="40px" align="middle"> <img src="https://marshalllab.github.io/MGDrivE/images/UCIMI.png" height="40px" align="middle">
