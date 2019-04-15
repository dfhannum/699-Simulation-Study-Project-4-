#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Power Calculation"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         numericInput("n",
                     "Sample Size:",
                     min = 20,
                     max = 300,
                     value = 100),
         numericInput("p",
                     "Disease Rate",
                     min = 0.05,
                     max = 0.45, 
                     value = 0.40),
         checkboxInput('c',
                       "Censoring of Non-Compliers",
                       value = T),
         checkboxInput('i',
                       "Interim Analysis", 
                       value = T),
         submitButton()
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         textOutput("power")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$power <- renderText({
           library(ggplot2)
           library(ReIns)
           library(survival)
           library(ggfortify)
           library(svMisc)
           library(survminer)
           library(GGally)
           library(truncnorm)
      # generate bins based on input$bins from ui.R
           lambda <- 0.0203821
           
           #Gives a number of random observations (x) from the truncated exponential distr
           rdist <- function(x = 1){
                   rtexp(x, lambda, endpoint = 100)
           }
           
           #Random number from the truncated normal distribution
           rdistn <- function(x = 1){
                   rtruncnorm(x, a = 0, b = 100, mean = 50, sd = 25)
           }
           
           #Dataframe with no censoring -> Intent to treat dataframe
           dataf <- function(n,p){
                   #n is the number of samples being tested
                   #p is the disease rate
                   
                   #Setting up the null dataframe
                   df <- as.data.frame(matrix(ncol = 2, nrow = n, rep(c(200,0), each = n)))
                   colnames(df) <- c('Time','Event')
                   
                   #Using a random binomial to determine if an event has occured
                   df$Event <- rbinom(n,1,p)
                   #Using the truncated exponential distr. to get the time for events
                   df[df$Event == 1,]$Time <- rdist(sum(df$Event))
                   
                   return (df)
           }
           
           #A dataframe using censoring to model the non-compliance of people to the trtmnt
           datafc <- function(n,p){
                   #n is the number of samples being tested
                   #p is the disease rate
                   df <- as.data.frame(matrix(ncol = 4, nrow = n,
                                              rep(c(200,0), each = n)))
                   colnames(df) <- c('Time1','Event1', 'Time2','Event2')
                   
                   #Event one is getting the GVHD
                   df$Event1 <- rbinom(n,1,p)
                   df[df$Event1 == 1,]$Time1 <- rdist(sum(df$Event1))
                   
                   #Event two is getting censored for non-compliance
                   df$Event2 <- rbinom(n,1,.8)
                   df[df$Event2 == 0,]$Time2 <- rdistn(n-sum(df$Event2))
                   
                   #Set the actual time for the minimum time between time 1 and 2
                   df$Time <- NA
                   df$Time <- ifelse(df$Time1 < df$Time2, df$Time <- df$Time1, 
                                     df$Time <- df$Time2)
                   df$Event <- ifelse(df$Time ==df$Time1, df$Event1, df$Event2)
                   
                   return (df)
           }
           #Seeing if the survival curve contains 0.50 in its CI at time = 100
           surv_sig <- function (df, confint = 0.95){
                   #Confint variable lets me vary the alpha to get the desired alpha
                   #of 0.05
                   
                   #Fitting a survival curve
                   fit <- survfit(Surv(Time, Event) ~1, data = df, conf.int = confint)
                   #Seeing if the 95% CI passes historic threshold = 0.50
                   l <- summary(fit, times = 100)$lower
                   return(l > 0.50)
           }
           #Power calculation with censoring
           sample_power_c <- function(p = c(.40,.35,.3), samples = 10, l = 40 , u = 100,
                                      by = 10, confint = 0.94, c = TRUE){
                   # p is the disease rates to be tested
                   # samples refers to the number of simulations
                   # l is the lower bound of the sample size
                   # u is the upper bound of the sample size
                   # by is the interval between sample sizes to be tested
                   # confint is the confidence interval applied to the logrank test
                   
                   
                   n <- seq(l, u, by = by)
                   
                   data <- as.data.frame(matrix(ncol = 3, nrow = length(n)*length(p),0))
                   
                   colnames(data) <- c('Sample Size', 'Rate', 'Power')
                   data$`Sample Size` <- rep(n,length(p))
                   data$Rate <- as.factor(rep(p, each = length(n)))
                   
                   for (k in 1:length(p)){
                           for (i in 1:length(n)){
                                   
                                   vect <- rep(NA,samples)
                                   for (j in 1:samples){
                                           if (c == T){
                                                   df <- datafc(n[i],p[k])
                                           } else {
                                                   df <- dataf(n[i], p[k])
                                           }
                                           vect[j] <- surv_sig(df, confint)
                                   }
                                   data[length(n)*(k-1) + i, 3] <- sum(vect)/samples
                           }
                   }
                   return (data)
           }
           
           #Sample power with an interim analysis
           sample_power_int <- function(p = .5, samples = 10, l = 40 , u = 100,
                                        by = 10, confint1 = 0.85, confint2 = 0.85, c = TRUE){
                   # two different variables to control the alpha for the interim and
                   # for the final analysis
                   n <- seq(l, u, by = by)
                   
                   data <- as.data.frame(matrix(ncol = 5, nrow = length(n)*length(p),0))
                   
                   colnames(data) <- c('Sample Size', 'Rate', 'Power','Rint',
                                       'Rfin')
                   
                   data$`Sample Size` <- rep(n,length(p))
                   data$Rate <- as.factor(rep(p, each = length(n)))
                   
                   for (k in 1:length(p)){
                           for (i in 1:length(n)){
                                   
                                   vect <- rep(NA,samples)
                                   
                                   #Creating empty vectors to track where trials fial
                                   Rfin <- 0
                                   Rint <- 0
                                   
                                   #Interim setting n = 2
                                   for (j in 1:samples){
                                           #Testing the data halfway through enrollment
                                           if (c == T){
                                                   df <- datafc(n[i]/2,p[k])
                                           } else {
                                                   df <- dataf(n[i]/2, p[k])
                                           }
                                           
                                           sig <- surv_sig(df, confint1)
                                           
                                           #Seeing if the interim fails
                                           if (sig == 0){
                                                   vect[j] <- sig
                                                   Rint <- Rint + 1
                                           }
                                           # If the interim does not fail then do the rest
                                           else{
                                                   if (c == T){
                                                           df <- rbind(df,
                                                                       datafc(n[i]/2,p[k]))
                                                   } else {
                                                           df <- rbind(df,
                                                                       dataf(n[i]/2, p[k]))
                                                   }
                                                   sig <- surv_sig(df,confint2)
                                                   vect[j] <- sig
                                           }
                                   }
                                   data[length(n)*(k-1) + i, 3] <- sum(vect)/samples
                                   data[length(n)*(k-1) + i, 4] <- Rint
                                   data[length(n)*(k-1) + i, 5] <- samples - Rint - sum(vect)
                           }
                   }
                   return (data)
           }
           power_only <- function (n, p = .5, samples = 1000, c = F, i = F){
                   #Decreasing by the 0.15 that will dropout
                   n <- n * .85
                   
                   if (i == T){
                           if (c == T){
                                   sp <- sample_power_int(p = p, samples = samples, 
                                                          l = n, u = n, 
                                                          by = 1, c = T)
                           } else{
                                   sp <- sample_power_int(p = p, samples = samples, 
                                                          l = n, u = n, 
                                                          by = 1, c = F)
                           }
                   } else{
                           if (c == T){
                                   sp <- sample_power_c(p = p, samples = samples, 
                                                        l = n, u = n, 
                                                        by = 1, c = T)
                           } else{
                                   sp <- sample_power_c(p = p, samples = samples, 
                                                        l = n, u = n, 
                                                        by = 1, c = F)
                           }
                   }
                   return(sp$Power)
           }
      # draw the histogram with the specified number of bins
      paste0("Power = ", round(power_only(n = input$n, p =input$ p, 
                                          c = input$c, i = input$i), 2))
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

