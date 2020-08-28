## app.R ##

#################################### Libraries ##########################################
library(abind)
library(shiny)
library(shinydashboard)
library(plotly)
library(tidyverse)

#################################### Functions ##########################################

prop.CI <- function(p, a, n){
  z <- qnorm(a/2)
  lb <- p + (z*sqrt(p*(1-p)/n) + 1/(2*n))
  ub <- p - (z*sqrt(p*(1-p)/n) + 1/(2*n))
  return(c(lb,ub))
}

estimate.cali <- function(n00,n10,n01,n11,a.hat.star){
  p00 <- n00 / (n00 + n01)
  p11 <- n11 / (n10 + n11)
  n <- sum(n00,n10,n01,n11)
  val.dat <- matrix(c(n00,n10,n01,n11), nrow = 2)

  alphas.hat <- c(1 - a.hat.star ,a.hat.star)
  p.est.cali <- validation.to.calibration(val.dat)
  est <- p.est.cali %*% alphas.hat
  return(est[2])
}

rmse.thr.prob <- function(n, p00, p11, a1){
  val1 <- (1-a1)*p00*(1-p00)*(1+a1/(n*(1-a1)))
  val2 <- a1*p11*(1-p11)*(1 + (1-a1)/(n*a1))
  val3 <- n*(p00+p11-1)^2

  tot.var <- (val1 + val2)/val3

  bias.prob <- function(n,p00,p11,a){
    p1 <- (a*(p00+p11-1) - (p00-1))/(n*(p00+p11-1)^3) * (p00*(1-p00)/(1-a) + p11*(1-p11)/a)
    p2 <- (p00-1)*p11/(n*(p00+p11-1)^3) * ((1-p11)/a + p00/(1-a))
    return(p1+p2)
  }
  tot.bias <- bias.prob(n,p00,p11,a1)

  return(sqrt(tot.var + tot.bias^2))
}

rmse.thr.naive <- function(N, p00, p11, a1){
  #vari <- (p11*(1-p11)*a1 + p00*(1-p00)*(1-a1)) / N
  bias <- (p11 - 1)*a1 + (1 - p00)*(1 - a1)
  #tot.mse <- bias^2 + vari
  tot.mse <- bias^2
  return(sqrt(tot.mse))
}

rmse.thr.vali <- function(n, alpha){
  vari <- alpha * (1 - alpha) / n
  return(sqrt(vari))
}

rmse.thr.cali <- function(n, p00, p11, a1){
  a.hat.star <- p11*a1 + (1-p00)*(1-a1)
  q00 <- ((1-a1)*p00) / ((1-a1)*p00 + a1*(1-p11))
  q11 <- (a1*p11) / (a1*p11 + (1-a1)*(1-p00))

  val1 <- (a.hat.star/n + (1 - a.hat.star)/n^2) * (q11 - q11^2)
  val2 <- (a.hat.star/n^2 + (1 - a.hat.star)/n) * (q00 - q00^2)

  return(sqrt(val1 + val2))
}

rmse.thr.esbi <- function(n, p00, p11, a1){
  e.a.hat.star <- p11*a1 + (1-p00)*(1-a1)
  bias.esbi <- a1 - (e.a.hat.star * (3-p11-p00) - (1-p00))

  v.p00 <- p00*(1-p00) / (n * (1 - a1)) * (1 + a1/(n*(1-a1)))
  v.p11 <- p11*(1-p11) / (n * a1) * (1 + (1- a1) /(n*a1))

  var.esbi <- v.p00 * (e.a.hat.star - 1)^2 + v.p11 * e.a.hat.star^2

  mse.esbi <- bias.esbi^2 + var.esbi
  return(sqrt(mse.esbi))
}

predictions.2classes <- function(runs, n, N, p00, p11, alpha1){
  p.real <- matrix(c(p00, 1 - p11, 1 - p00, p11), nrow = 2)
  alphas <- c(1 - alpha1,alpha1)

  # Get "real data set"
  pop.dat <- populationdata.nclasses(runs, p.real, alphas, N)
  # Take sample to validate
  val.dat <- sample.from.populationdata(pop.dat, n)
  # Get non-corrected estimates of alpha1
  alphas.hat <- apply(pop.dat, 3, function(x) {sum(x[,2]) / sum(x)})
  alphas.v <- apply(val.dat, 3, function(x) {sum(x[2,]) / n})
  # Compute Inversed Contigency Matrix and calibration matrix
  p.est.prob <- validation.to.probability(val.dat)
  p.est.cali <- validation.to.calibration(val.dat)


  # Probabilities as a vector
  p00.hat <- p.est.prob[1,1,]
  p11.hat <- p.est.prob[2,2,]

  c01.hat <- p.est.cali[2,1,]
  c11.hat <- p.est.cali[2,2,]
  #Compute estimations per method and calculate RMSE
  est.prob <- alphas.hat/(p00.hat+p11.hat-1) + (p00.hat-1)/(p00.hat+p11.hat-1)
  est.cali <- c01.hat * (1-alphas.hat) + c11.hat * alphas.hat
  est.esbi <- alphas.hat - (p11.hat - 1)*alphas.hat - (1 - p00.hat)*(1-alphas.hat)


  dfr <- data.frame("C.Validation.Data" = alphas.v,
                    "B.Contingency.Matrix" = est.prob,
                    "A.Calibration.Matrix"= est.cali,
                    "E.Naive.Estimation" = alphas.hat,
                    "D.Bias.Correction" = est.esbi)
  dfr <- pivot_longer(dfr,
                      cols = c("C.Validation.Data", "B.Contingency.Matrix",
                               "A.Calibration.Matrix", "D.Bias.Correction",
                               "E.Naive.Estimation"),
                      names_to = "Method",
                      values_to = "Value")
  return(dfr)
}


expected.populationdata.nclasses <- function(p.mat, a.mat, N){
  if(nrow(p.mat) != length(a.mat)){stop("Length alphas unequal to dimension p-matrix.")}
  if(nrow(p.mat) != ncol(p.mat)){stop("Inversed Contigency Matrix not correctly defined.")}
  # Obtain dimensions
  dims <- nrow(p.mat)
  # Obtain probability per cell
  val.dat <- p.mat * rep(a.mat, times = dims) * N

  # Obtain all integers
  val.dat.integers <- floor(val.dat)
  # And their decimals
  val.dat.decimals <- val.dat %% 1

  # Check how many numbers are left to fill in
  rest.numbers <- N - sum(val.dat.integers)

  # And assign them to the ones with the highest decimals
  highest.indices <- order(val.dat.decimals, decreasing = T)[seq_len(rest.numbers)]
  val.dat.integers[highest.indices] =  val.dat.integers[highest.indices] + 1

  ## Nice output
  rownames(val.dat.integers) <- paste0("True:", 0:(dims-1))
  colnames(val.dat.integers) <- paste0("Obs:", 0:(dims-1))

  ## Output
  return(val.dat.integers)
}

populationdata.nclasses <- function(runs, p.mat, a.mat, N){
  pop.dat <- expected.populationdata.nclasses(p.mat,a.mat,N)
  n <- nrow(p.mat)

  Sampler <- function(n){
    samps <- lapply(1:n, function(x) {
      sample(1:n, size = sum(pop.dat[x,]), prob = p.mat[x,], replace = T)
    })
    out <- lapply(samps, tabulate, nbins = n) %>% unlist() %>% matrix(nrow = n, byrow = T)
  }

  out <- replicate(runs, Sampler(n))


  return(out)
}


sample.from.populationdata <- function(val.dat, n){
  #sample from a dataset (without replacement)
  Sampler <- function(val.dat, n){
    vec <- as.vector(val.dat)
    len <- length(vec)
    samples <- sample(rep(1:len, times = vec), size = n, replace = F)
    counts <- sapply(1:len, function(x) {sum(samples == x)})
    counts <- matrix(counts, nrow = sqrt(len))
    return(counts)
  }

  #repeat a certain amount of times for simulation purposes
  arr <- apply(val.dat, 3, function(x) {Sampler(x, n)})
  arr <- array(arr, dim=dim(val.dat))
  #create nice output
  dim3 <- paste0("Iteration ", 1:dim(val.dat)[3])
  return(arr)
}

validation.to.probability <- function(val.dat){
  if(length(dim(val.dat)) == 2){
    pmat <- val.dat / rowSums(val.dat)
    return(pmat)
  }
  else{
    #create probability matrix
    pmat <- apply(val.dat,3,function(x) x/rowSums(x))
    #set dimensions right for output
    dim <- sqrt(nrow(pmat))
    runs <- ncol(pmat)
    pmat <- array(pmat, dim = c(dim, dim, runs), dimnames = dimnames(val.dat))
    #return matrix
    return(pmat)
  }
}

validation.to.calibration <- function(val.dat){
  if(length(dim(val.dat)) == 2){
    pc.mat <- t(t(val.dat) / colSums(val.dat))
    return(pc.mat)
  }
  else{
    #create probability matrix
    pc.mat <- apply(val.dat,3,function(x) t(t(x)/colSums(x)))
    #set dimensions right for output
    dim <- sqrt(nrow(pc.mat))
    runs <- ncol(pc.mat)
    pc.mat <- array(pc.mat, dim = c(dim, dim, runs), dimnames = dimnames(val.dat))
    #return matrices
    return(pc.mat)
  }
}

data.rmseplot <- function(p00_left, p00_right, p11_left, p11_right, n, N, alpha, methods, steps){
  p00 <- seq(p00_left, p00_right, length.out = steps)
  p11 <- seq(p11_left, p11_right, length.out = steps)
  p.grid <- expand.grid(p00, p11)

  data.prob <- data.cali <- data.vali <- data.naiv <- data.esbi <- data.esbi2 <- NULL

  if("probability" %in% methods){
    data.prob <- mapply(rmse.thr.prob, n, p.grid[,1], p.grid[,2], alpha)
    data.prob <- matrix(data.prob, nrow = length(p00), byrow = T)
  }
  if("calibration" %in% methods){
    data.cali <- mapply(rmse.thr.cali, n, p.grid[,1], p.grid[,2], alpha)
    data.cali <- matrix(data.cali, nrow = length(p00), byrow = T)
  }
  if("validation" %in% methods){
    data.vali <- matrix(sqrt(alpha*(1-alpha)/n), nrow = length(p00), ncol = length(p11))
  }
  if("naive" %in% methods){
    data.naiv <- mapply(rmse.thr.naive, N, p.grid[,1], p.grid[,2], alpha)
    data.naiv <- matrix(data.naiv, nrow = length(p00), byrow = T)
  }
  if("estbias" %in% methods){
    data.esbi <- mapply(rmse.thr.esbi, n, p.grid[,1], p.grid[,2], alpha)
    data.esbi <- matrix(data.esbi, nrow = length(p00), byrow = T)
  }
  return(list(data.naiv = data.naiv,
              data.esbi = data.esbi,
              data.vali = data.vali,
              data.prob = data.prob,
              data.cali = data.cali,
              n = n,
              N = N,
              p00 = c(p00_left, p00_right),
              p11 = c(p11_left, p11_right),
              alpha = alpha,
              methods = methods,
              steps = steps))
}

dash.rmseplot <- function(p00_left, p00_right, p11_left, p11_right, n, N, alpha, methods, steps){

  p00 <- seq(p00_left, p00_right, length.out = steps)
  p11 <- seq(p11_left, p11_right, length.out = steps)
  p.grid <- expand.grid(p00, p11)

  dat <- data.rmseplot(p00_left, p00_right, p11_left, p11_right, n, N, alpha, methods, steps)

  color1 <- rep(0, length(p00) * length(p11))
  dim(color1) <- dim(dat$data.prob)
  color2 <- color1 + 1/4
  color3 <- color1 + 2/4
  color4 <- color1 + 3/4
  color5 <- color1 + 1

  # create plot
  p <- plot_ly(x = ~p00, y = ~p11, showscale = F)
  if ("probability" %in% methods){
    p <- p %>% add_surface(z = ~dat$data.prob, surfacecolor = color1,
                           cauto = F, cmax = 1, cmin = 0,
                           showscale = F, name = "Inversed Contigency Matrix")
  }
  if("validation" %in% methods){
    p <- p %>% add_surface(z = ~dat$data.vali, surfacecolor = color2,
                           cauto = F, cmax = 1, cmin = 0,
                           showscale = F, name = "Validation Data")
  }
  if("calibration" %in% methods){
    p <- p %>% add_surface(z = ~dat$data.cali, surfacecolor = color3,
                           cauto = F, cmax = 1, cmin = 0,
                           showscale = F, name = "Calibration Matrix")
  }
  if("naive" %in% methods){
    p <- p %>% add_surface(z ~ dat$data.naiv, surfacecolor = color4,
                           cauto = F, cmax = 1, cmin = 0,
                           showscale = F, name = "Naive Estimation")
  }
  if("estbias" %in% methods){
    p <- p %>% add_surface(z ~ dat$data.esbi, surfacecolor = color5,
                           cauto = F, cmax = 1, cmin = 0,
                           showscale = F, name = "Estimated Bias: Biased")
  }
  p <- p %>% layout(
    title = paste("RMSE with alpha = ", alpha, ", n = ", n ,sep = ""),
    scene = list(
      xaxis = list(title = "p00", showgrid = FALSE),
      yaxis = list(title = "p11", showgrid = FALSE),
      zaxis = list(title = "RMSE", showgrid = FALSE)
    ))

  return(p)
}

################################ Dashboard #############################################
header <- dashboardHeader(title = "Bias Correction Dashboard")
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Descriptive One Point", tabName = "descriptives1p", icon = icon("dashboard")),
    menuItem("Descriptive Curve", tabName = "descriptives", icon = icon("chart-line")),
    menuItem("Prediction", tabName = "prediction", icon = icon("lightbulb"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "prediction",
        fluidRow(
          box(numericInput(inputId = "n00",
                           label = "Number of objects correctly classfied in class 0:",
                           value = 10),
              numericInput(inputId = "n01",
                           label = "Number of objects incorrectly classfied in class 0:",
                           value = 10),
              width = 4),
          box(numericInput(inputId = "n11",
                           label = "Number of objects correctly classfied in class 1:",
                           value = 10),
              numericInput(inputId = "n10",
                           label = "Number of objects incorrectly classfied in class 1:",
                           value = 10),
              width = 4),
          box(numericInput(inputId = "alpha.hat.star",
                           label = "Naive estimation of percentage class 1 (%)",
                           value = 50),
              width = 4)),

        fluidRow(
          infoBoxOutput("sensitivityBox"),
          infoBoxOutput("specificityBox")
        ),

        fluidRow(
          infoBoxOutput("naiveEstimate"),
          infoBoxOutput("observedEstimate"),
          infoBoxOutput("correctedEstimate")
        )
    ),
    tabItem(tabName = "descriptives1p",
        fluidRow(
          box(numericInput(inputId = "p00",
                           label = "Probability of objects in class 0 correctly classfied (p00):",
                           value = 0.80),
              numericInput(inputId = "p11",
                           label = "Probability of objects in class 1 correctly classfied (p11):",
                           value = 0.90),
              width = 4),
          box(numericInput(inputId = "n",
                           label = "Sample size of the validation set (n):",
                           value = 300),
              numericInput(inputId = "N",
                           label = "Population size (N):",
                           value = 300000),
              width = 4),
          box(numericInput(inputId = "alpha",
                           label = "True proportion of objects in class 1 (alpha)",
                           value = 0.85),
              numericInput(inputId = "runs",
                           label = "Amount of runs in simulation",
                           value = 10),
              width = 4)),
          fluidRow(
            valueBoxOutput("naiveMSE", width = 6),
            valueBoxOutput("estbiasMSE", width = 6)),
          fluidRow(
            valueBoxOutput("validationMSE", width = 4),
            valueBoxOutput("probabilityMSE", width = 4),
            valueBoxOutput("calibrationMSE", width = 4)),
          fluidRow(
            plotOutput("boxplot")
          )
    ),
    tabItem(tabName = "descriptives",
        fluidRow(
          box(sliderInput(inputId = "p00Range",
                          label = "Range for probabilities of objects in class 0 correctly classfied (p00)",
                          min = 0, max = 1, value = c(0.6,1)),
              sliderInput(inputId = "p11Range",
                          label = "Range for probabilities of objects in class 1 correctly classfied (p11)",
                          min = 0, max = 1, value = c(0.6,1)),
              width = 4,
              height = "20em"),
          box(numericInput(inputId = "n2",
                           label = "Sample size of the validation set (n):",
                           value = 300),
              sliderInput(inputId = "steps",
                          label = "Steps in the graph. Higher number -> more accuracy, but longer computation time:",
                          min = 0, max = 1000, value = 500),
              numericInput(inputId = "alpha2",
                           label = "True proportion of objects in class 1 (alpha)",
                           value = 0.85),
              width = 4,
              height = "20em"),
          box(checkboxGroupInput(inputId = "methods",
                             label = "Which methods are shown in the plot?",
                             choiceNames = list("Naive Method",
                                                "Estimated Bias Method",
                                                "Validation Data",
                                                "Inversed Contigency Matrix",
                                                "Calibration Matrix"),
                             choiceValues = list("naive",
                                                 "estbias",
                                                 "validation",
                                                 "probability",
                                                 "calibration"),
                             selected = "calibration"),
              actionButton("submit", "Update Plots"),
              width = 4,
              height = "20em"),
          ),
        fluidRow(
          splitLayout(cellWidths = c("50%", "50%"),
                      plotlyOutput(outputId = "d3plot"),
                      plotOutput(outputId = "mini"))
        )
    )
  )
)
ui <- dashboardPage(
  skin = "green",
  header,
  sidebar,
  body
)


server <- function(input, output) {
  output$sensitivityBox <- renderInfoBox({
    infoBox("Accuracy in predicting Class 0 (95% CI)",
            paste0(round(100 * input$n00 / (input$n00 + input$n01),3),
                   "% (",
                   100 * round(prop.CI(input$n00 / (input$n00 + input$n01), 0.05, (input$n00 + input$n01))[1],5),
                   "% - ",
                   100 * round(prop.CI(input$n00 / (input$n00 + input$n01), 0.05, (input$n00 + input$n01))[2],5),
                   "%)"),
            icon = icon("dice", lib = "font-awesome"),
            color = "blue",
            width = 3,
            fill = T)
  })
  output$specificityBox <- renderInfoBox({
    #Hoe ver afgerond?
    infoBox("Accuracy in predicting Class 1 (95% CI)",
            paste0(round(100 * input$n11 / (input$n10 + input$n11),3),
                   "% (",
                   100 * round(prop.CI(input$n11 / (input$n10 + input$n11), 0.05, (input$n10 + input$n11))[1],5),
                   "% - ",
                   100 * round(prop.CI(input$n11 / (input$n10 + input$n11), 0.05, (input$n10 + input$n11))[2],5),
                   "%)"),
            icon = icon("dice", lib = "font-awesome"),
            color = "blue",
            width = 3,
            fill = T)
  })
  output$naiveEstimate <- renderInfoBox({
    infoBox("Naive estimate of proportion class 1",
            paste0(100 * round(sum(input$n01, input$n11)/sum(input$n01, input$n11, input$n00, input$n10), 5),
                   "%"),
            color = "yellow")
  })
  output$observedEstimate <- renderInfoBox({
    infoBox("Observed estimate of proportion class 1 in validation set",
            paste0(100 * round(sum(input$n10, input$n11)/sum(input$n01, input$n11, input$n00, input$n10), 5),
                   "%"),
            color = "green")
  })
  output$correctedEstimate <- renderInfoBox({
    infoBox("Corrected estimate of proportion class 1",
            paste0(100 * round(estimate.cali(input$n00, input$n10, input$n01,
                                             input$n11, input$alpha.hat.star ),5),
                   "%"),
            color = "red")
  })
  output$naiveMSE <- renderValueBox({
    valueBox(subtitle = "RMSE with Naive Estimation",
            value = rmse.thr.naive(input$N, input$p00 , input$p11 , input$alpha) %>% signif(digits = 4),
            color = "yellow",
            icon = icon("laptop-code"))
  })
  output$estbiasMSE <- renderValueBox({
    valueBox(subtitle = "RMSE with Estimated Bias Correction: Biased",
             value = rmse.thr.esbi(input$n, input$p00 , input$p11 , input$alpha) %>% signif(digits = 4),
             color = "red",
             icon = icon("laptop-code"))
  })
  output$validationMSE <- renderValueBox({
    valueBox(subtitle = "RMSE with Validation Data",
            value = rmse.thr.vali(input$n, input$alpha ) %>% signif(digits = 4),
            color = "green",
            icon = icon("laptop-code"))
  })
  output$probabilityMSE <- renderValueBox({
    valueBox(subtitle = "RMSE with Inversed Contigency Matrix",
            value = rmse.thr.prob(input$n, input$p00 , input$p11 , input$alpha ) %>% signif(digits = 4),
            color = "light-blue",
            icon = icon("laptop-code"))
  })
  output$calibrationMSE <- renderValueBox({
    valueBox(subtitle = "RMSE with Calibration Matrix",
            value = rmse.thr.cali(input$n, input$p00 , input$p11 , input$alpha ) %>% signif(digits = 4),
            color = "purple",
            icon = icon("laptop-code"))
  })

  output$boxplot <- renderPlot({
    d <- predictions.2classes(input$runs, input$n, input$N, input$p00, input$p11, input$alpha)
    p <- ggplot(d) +
         geom_boxplot(aes(x = Value , y = Method, fill = Method)) +
         scale_fill_manual(values = c("purple", "lightblue", "green", "red", "yellow")) +
         theme(axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank()) +
         xlab('Estimated alpha') +
      geom_vline(aes(xintercept = input$alpha), lty = 2, lwd = 2, color = "black")
    p
  })
  dat <- eventReactive(input$submit, {
    data.rmseplot(input$p00Range[1], input$p00Range[2], input$p11Range[1], input$p11Range[2],
                  input$n2, 300000, input$alpha2, input$methods, input$steps)})
  output$d3plot <- renderPlotly({
    p00 <- seq(dat()$p00[1], dat()$p00[2], length.out = dat()$steps)
    p11 <- seq(dat()$p11[1], dat()$p11[2], length.out = dat()$steps)
    p.grid <- expand.grid(p00, p11)

    color1 <- rep(0, length(p00) * length(p11))
    dim(color1) <- c(length(p00), length(p11))
    color2 <- color1 + 1/4
    color3 <- color1 + 2/4
    color4 <- color1 + 3/4
    color5 <- color1 + 1
    # create plot
    p <- plot_ly(x = ~p00, y = ~p11, showscale = F)
    if (!is.null(dat()$data.prob)){
      p <- p %>% add_surface(z = ~dat()$data.prob, surfacecolor = color1,
                             cauto = F, cmax = 1, cmin = 0,
                             showscale = F, name = "Inversed Contigency Matrix")
    }
    if(!is.null(dat()$data.vali)){
      p <- p %>% add_surface(z = ~dat()$data.vali, surfacecolor = color2,
                             cauto = F, cmax = 1, cmin = 0,
                             showscale = F, name = "Validation Data")
    }
    if(!is.null(dat()$data.cali)){
      p <- p %>% add_surface(z = ~dat()$data.cali, surfacecolor = color3,
                             cauto = F, cmax = 1, cmin = 0,
                             showscale = F, name = "Calibration Matrix")
    }
    if(!is.null(dat()$data.naiv)){
      p <- p %>% add_surface(z ~ dat()$data.naiv, surfacecolor = color4,
                             cauto = F, cmax = 1, cmin = 0,
                             showscale = F, name = "Naive Estimation")
    }
    if (!is.null(dat()$data.esbi2)){
      p <- p %>% add_surface(z = ~dat()$data.esbi2, surfacecolor = color6,
                             cauto = F, cmax = 1, cmin = 0,
                             showscale = F, name = "Estimated Bias")
    }
    p <- p %>% layout(
      title = paste("RMSE with alpha = ", dat()$alpha, ", n = ", dat()$n ,sep = ""),
      scene = list(
        xaxis = list(title = "p00", showgrid = F),
        yaxis = list(title = "p11", showgrid = F),
        zaxis = list(title = "RMSE", showgrid = F)
      ))

    p
  })
  output$mini <- renderPlot({
    meth <- c("naive", "estbias", "validation", "contingency", "calibration")
    fac <- meth[which(meth %in% dat()$methods)]
    vals <- abind(dat()[1:5], along = 3)
    min <- apply(vals, c(1,2), which.min)
    p00.seq <- seq(dat()$p00[1],dat()$p00[2], length.out = dat()$steps)
    p11.seq <- seq(dat()$p11[1],dat()$p11[2], length.out = dat()$steps)
    dfr <- data.frame(p00 = rep(p00.seq, each = length(p11.seq)),
                      p11 = rep(p11.seq, times = length(p00.seq)),
                      minimum.rmse = factor(min, labels = fac, levels = 1:length(fac)))
    ggplot(dfr, aes(x = p00, y = p11, fill = minimum.rmse)) +
      geom_raster() +
      ggtitle(paste("Method with lowest RMSE for alpha =",
                    dat()$alpha, "and n =", dat()$n, sep = " "))
  })
}

####################################### RUN #######################################
shinyApp(ui, server)
