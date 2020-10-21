# Comparing Correction Methods for Misclassification Bias: Dashboard

## Description
This GitHub page contains additional information for the paper *Comparing Correction Methods for Misclassification Bias*, written for the [BNAIC/BENELEARN Conference 2020](https://bnaic.liacs.leidenuniv.nl/). At this page, you can find the paper (including the appendix with mathematical derivations and proofs) and a interactive dashboard. This dashboard makes it possible to apply the theory presented in the paper in an accessible way.

## Installation
In order to use the dashboard, you need to install the following software packages:
  * [R](https://cran.r-project.org/bin/windows/base/) (version 4.0.2 or higher)
  * [RStudio](https://rstudio.com/products/rstudio/download/)

Within R, you will need the following packages:
  * [abind](https://www.rdocumentation.org/packages/abind)
  * [ggpubr](https://www.rdocumentation.org/packages/ggpubr)
  * [plotly](https://www.rdocumentation.org/packages/plotly)
  * [plyr](https://www.rdocumentation.org/packages/plyr)
  * [RColorBrewer](https://www.rdocumentation.org/packages/RColorBrewer)
  * [shiny](https://www.rdocumentation.org/packages/shiny)
  * [shinydashboard](https://www.rdocumentation.org/packages/shinydashboard)
  * [tidyverse](https://www.rdocumentation.org/packages/tidyverse)
  
After installing the software, the app can be opened by pressing **Run App** at the top right of the R Script or running the whole script.
  
## Usage 
The dashboard contains two pages: *Descriptive One Point* and *Descriptive Curve*. At the page *Descriptive One Point*, you can find the RMSE and a boxplot of the data points of all estimators that are introduced in the paper. At the page *Descriptive Curve*, you can find curve of the RMSE and a overview which estimators performs the best given the input.

### Descriptive One Point
Page that shows the RMSE and a boxplot of estimates of the five estimators. The new data can be applied by pressing **Update Boxplot**.

#### Input
  * Sensitivity 
  * Specificity
  * Size of test set
  * Size of unlabeled data set
  * Proportion of objects in class 1 (alpha)
  * Amount of runs in simulation
 
#### Output
  * Boxes with RMSE of the estimators
  * Boxplot that shows the distribution of the estimates
  
### Descriptive Curve
Page that shows an interactive curve of the RMSE and a heat map which selected method has the lowest RMSE. Updating the input values to the plots can be done by pressing **Update Plots**.

#### Input
  * Range of Sensitivity 
  * Range of Specificity
  * Size of test set
  * Step size of graphs
  * Proportion of objects in class 1 (alpha)
  * Selection of estimators

#### Output
  * Interactive plot of the RMSE
  * Heat map of method with lowest RMSE
