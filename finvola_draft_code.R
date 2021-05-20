#Volatility Modeling and NOPE
#Jonas Schmitten, Moritz Köhler, Giovanni Magagnin, Noah Angara
#May 2021

#Packages
library(fGarch)
library(tidyquant)
library(fBasics)
library(vars)
library(tseries)
library(dplyr)


# Data Prep & Cleaning ----------------------------------------------------

#set wd according to who is executing the code
Paths = c("/Users/jonasschmitten/Desktop/FS 2021/Economics in Practice/Clean Data", 
          "/Users/noahangara/Documents/Master's/8th Semester/Financial Volatility")
names(Paths) = c("jonasschmitten", "noahangara")
setwd(Paths[Sys.info()[7]])

### Loading Option and SPX Data
data_option = read.csv("options_data.csv")
sp500 = tq_get('^GSPC', get = "stock.prices", from ='2019-01-01', to = '2020-12-31')

### Not interested in High or Low Price
data_spx = sp500 %>%
  subset(select = -c(high, low))
rm(sp500)

data_option = data_option %>%
  subset(select = -c(secid, index_flag, issue_type, issuer, exercise_style))

### Dataframe with SPX Returns
spx_log_returns = data.frame("Returns" = diff(log(data_spx$close)))
  
# Exploring Data ----------------------------------------------------------

### Calculating basic statistics (from fBasics)
basic_stats = basicStats(spx_log_returns)

### Plotting histogram for the distribution
hist(spx_log_returns)
# investigating the kurtosis for the tails
basic_stats
# already from here we can say that it is not normal, the kurtosis is too high but let's do the test

### Jarque-Bera test for normality
# first, turn log returns into percentages
percentage_returns=log(spx_log_returns+1)*100 
# run the test
normalTest(percentage_returns, method = "jb")
# reject normality

### Testing the hypothesis of a white noise
# professor does this with monthly return, I am doing with daily
acf(spx_log_returns,lag.max = 50,main = "S&P500")
pacf(spx_log_returns, lag.max = 50, main = "S&P500")
# we can reject the hypothesis of white noise 

### Testing for stationarity
# qualitatively assessing stationarity with the plot of returns
plot(spx_log_returns)
# Augmented Dickey Fuller Test
adf.test(spx_log_returns)
# it seems stationary, but let's run other tests to be sure
# Phillip Perron Unit Root Test
pp.test(spx_log_returns)
# ok it is stationary

### Plotting returns, acf, acf for squared returns, pacf
plot(spx_log_returns, main = "S&P500 log returns")
acf(spx_log_returns, main = "S&P500 log returns")
acf(spx_log_returns^2, main = "S&P500 squared log returns")
pacf(spx_log_returns^2, main = "S&P500 squared log returns")

### Ljung-Box test for autocorrelation
Box.test(spx_log_returns,lag=30,type="Ljung")
Box.test(spx_log_returns^2,lag=30,type="Ljung")
# randomness hypothesis rejected for both the cases

### Idea behind the ARCH / GARCH: we are accounting for the volatility in addition to the mean
# We have to estimate conditional variance and conditional mean jointly
# We can use ARMA for the conditional mean (since there's some information flow)


# Fitting Volatility Models -----------------------------------------------------


### Fitting the model
# ARCH(3) - trying this randomly
arch3 = garchFit(spx_log_returns ~ garch(3,0), data = spx_log_returns, trace=F)
summary(arch3)
# all the alphas seem statistically significant
# trying now different ARCH

# ARCH(5)
arch5 = garchFit(spx_log_returns ~ garch(5,0), data = spx_log_returns, trace=F)
summary(arch5)
# alpha 5 is not significant, probably arch4 is the way to go

# ARCH(4)
arch4 = garchFit(spx_log_returns ~ garch(4,0), data = spx_log_returns, trace=F)
summary(arch4)
# looks ok
# predicting for next 90 days
predict(arch4, 90)

# trying now arch4 with student t for the error distribution
arch4_student = garchFit(spx_log_returns ~ garch(4,0), data = spx_log_returns, 
                         trace=F, cond.dist = "std")
summary(arch4_student)
# predicting
predict(arch4_student, 90)
# plotting random stuff according to the prof
acf(arch4_student@residuals/arch4_student@sigma.t,main="ARCH(4) residuals")
acf(arch4_student@residuals^2/arch4_student@sigma.t^2,main="ARCH(4) squared residuals")
acf(abs(arch4_student@residuals/arch4_student@sigma.t),main="ARCH(4) absolute residuals")


################################################################################
### TRYING WITH GARCH

### Testing for the GARCH effects in the daily log returns
LM=function(x,h)
{
  n=length(x)
  x2=x^2-mean(x^2)
  dat=matrix(,n-h,h+1)
  for (i in 1:(h+1))
  {
    dat[,i]=x2[(h+2-i):(n-i+1)]
  }
  a=lm(dat[,1]~dat[,2:(h+1)])
  r2=summary(a)$r.squared
  print(r2 * n)
  print(1-pchisq(r2*n,h))
}
# calling the function, 4 lags
LM(spx_log_returns,4)
# we can reject the hypothesis, there are GARCH effects in the returns

### Fitting the model
garch1 = garchFit(spx_log_returns ~ garch(1,1), data = spx_log_returns, trace=F)
summary(garch1)
# for the second moment m2 we should do alpha1 + beta1
# second moment
(garch1@fit$matcoef[3,1]+garch1@fit$matcoef[4,1])
# fourth moment
(garch1@fit$matcoef[3,1]+garch1@fit$matcoef[4,1])^2+
  (mean(garch1@residuals^4/garch1@sigma.t^4)-1)*garch1@fit$matcoef[3,1]^2
# predicting
predict(garch1, 90)
# for the plots, check the profs code





