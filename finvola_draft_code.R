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
library(tidyr)
library(tibble)
library(rugarch)
library(MSGARCH)


# Data Prep & Cleaning ----------------------------------------------------

#set wd according to who is executing the code
Paths = c("/Users/jonasschmitten/Downloads/NOPE and Volatility/Data", 
          "/Users/noahangara/Documents/Master's/8th Semester/Financial Volatility",
          "C:/Users/magag/Documents/UNIVERSITA/_St.GALLEN - Master in Quant Economics and Finance/Subjects/2^ Semester/Financial volatility/project")
names(Paths) = c("jonasschmitten", "noahangara", "magag")
setwd(Paths[Sys.info()[7]])

### Loading Option and SPX Data
#data_option = read.csv("options_data.csv", nrows = 1000000)
sp500 = tq_get('^GSPC', get = "stock.prices", from = as.Date('2019-01-01')-732, to = '2020-12-31')

### Not interested in High or Low Price
data_spx = sp500 %>%
  subset(select = -c(high, low))


### Dataframe with SPX Returns
spx_log_returns = data.frame("Returns" = diff(log(data_spx$close)))
spx_log_returns = spx_log_returns %>%
  mutate(data.frame(sp500[2:nrow(sp500), "date"])) %>%
  xts(spx_log_returns, order.by = spx_log_returns$date)
spx_log_returns = spx_log_returns[,-2]
rm(sp500)
  
#remove not needed columns 
data_option = data_option %>%
  subset(select = -c(secid, index_flag, issue_type, issuer, exercise_style, optionid, contract_size))

#convert date and expiration date to date value
data_option[c('date', 'exdate')] =  lapply(data_option[c('date', 'exdate')], function(x) as.Date(as.character(x), "%Y%m%d"))

#convert strike price, best bid, bid ask to dollars instead of cents
data_option[c('strike_price', 'best_bid', 'best_offer')] = data_option[c('strike_price', 'best_bid', 'best_offer')]/1000 

#order data
data_option = data_option %>% 
  arrange(date, exdate, cp_flag, strike_price)

#Check if there is 0 delta
unique(data_option$delta == 0)

#Remove rows with NA delta 
data_option = data_option %>% 
  drop_na(delta)

#Get deltas to 100 and -100 (Saw this on NOPECord)
data_option$delta = data_option$delta*100

#NOPE
data_option$NO = group_by(data_option, date, cp_flag)$volume * data_option$delta

NOPE =  group_by(data_option, date, cp_flag) %>% 
  summarise(NOPE = sum(NO))

NOPE = as.data.frame(rowsum(NOPE$NOPE, as.integer(gl(nrow(NOPE), 2, nrow(NOPE)))))
NOPE = add_column(NOPE, date = unique(data_option$date), .before = 1)

NOPE$volume = sp500[1:nrow(NOPE),which(colnames(sp500)=="volume")]

NOPE$NOPE = (NOPE$V1/NOPE$volume)*100


#Gamma exposure GEX
data_option$GEX = ifelse(data_option$cp_flag == "C",
                         data_option$gamma*data_option$open_interest*100,data_option$gamma*data_option$open_interest*-100)




#Black-Scholes Option Pricing
#Has to be made dynamic

# inputs
K_call = 3700
K_put = 3400
m = tau/252

# formula (21.3) & (21.4) from lecture notes to calculate option prices
black_scholes <- function(S, K, y, m, sig, call = TRUE) {
  d1 <- (log(S/K) + (y + (sig^2)/2) * m) / (sig * sqrt(m))
  d2 <- d1 - sig * sqrt(m)
  if (call == TRUE) { # CALL OPTION
    C <- S * pnorm(d1) - exp(-y*m) * K * pnorm(d2)  
  } else { # PUT OPTION
    C <- exp(-y*m) * K * pnorm(-d2) - S * pnorm(-d1)
  }
  return(C)
}


#define risk-free as 3-month T-Bill
y = 0.0001

#Apply BS to each row 
#not working yet
#stock price has to be exchanged
#price cannot be negative 
data_option$BS = ifelse(data_option$cp_flag=='C', black_scholes(3000, data_option$strike_price, y, as.numeric(unlist(data_option['exdate'] - data_option['date']))/252, data_option$impl_volatility, call =T), 
                        black_scholes(3000, data_option$strike_price, y, as.numeric(unlist(data_option['exdate'] - data_option['date']))/252, data_option$impl_volatility, call =F))







# GREEKS MANUALLY
d1_call <- (log(s0/K_call) + (risk_free_rate + (expected_vola^2)/2) * m) / (expected_vola * sqrt(m))
d2_call <- d1_call - expected_vola * sqrt(m)

# DELTA
delta_call <- pnorm(d1_call)
# GAMMA
gamma_call <- dnorm(d1_call) / (s0 * expected_vola * sqrt(m))
# THETA
theta_call <- (-(risk_free_rate) * K_call * exp(-risk_free_rate * m) * pnorm(d2_call) - (1/(2*sqrt(m))) * s0 * dnorm(d1_call) * expected_vola) / 365
# VEGA
vega_call <- s0 * dnorm(d1_call) * sqrt(m) / 100
# RHO
rho_call <- m * K_call * exp(-risk_free_rate * m) * pnorm(d2_call) / 100

# GREEKS WITH PACKAGE
delta_call
gamma_call
theta_call
vega_call
rho_call

library(derivmkts)
greeks(bscall(s = s0, K_call, expected_vola, risk_free_rate, m, 0))















# Exploring Data ----------------------------------------------------------

### Calculating basic statistics (from fBasics)
basic_stats = basicStats(spx_log_returns)

### Plotting histogram for the distribution
hist(spx_log_returns$Returns,breaks = 200)
# investigating the kurtosis for the tails
basic_stats
# already from here we can say that it is not normal, the kurtosis is too high but let's do the test

### Jarque-Bera test for normality
# first, turn log returns into percentages
#percentage_returns=log(spx_log_returns+1)*100 
# run the test
normalTest(log(spx_log_returns+1)*100 , method = "jb")
# reject normality

### Testing the hypothesis of a white noise
# professor does this with monthly return, I am doing with daily
acf(spx_log_returns,lag.max = 50,main = "S&P500")
pacf(spx_log_returns, lag.max = 50, main = "S&P500")
# we can reject the hypothesis of white noise 

### Testing for stationarity
# qualitatively assessing stationarity with the plot of returns
plot(y=spx_log_returns$Returns, x = sp500$date[2:nrow(sp500)], ty)
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



# rolling forecast for GARCH(1,1) with forecast length of 504 days
spec_garch = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                         mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                         distribution.model = "norm")

mod_garch = ugarchroll(spec_garch, data = spx_log_returns, n.ahead = 1, 
                 n.start = 504,  refit.every = 1, window.size= 504, refit.window=c('moving'), 
                 solver = "hybrid", fit.control = list(), keep.coef = TRUE)

# rolling forecast for GJR-GARCH with forecast length of 504 days
spec_gjrgarch = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)), 
                           mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                           distribution.model = "std")

mod_gjrgarch = ugarchroll(spec_gjrgarch, data = spx_log_returns, n.ahead = 1, 
                       n.start = 504,  refit.every = 1, window.size= 504, refit.window=c('moving'), 
                       solver = "hybrid", fit.control = list(), keep.coef = TRUE)

# rolling forecast for Markov-Switching GARCH with forecast length of 504 days
n.ots <- 504
n.its <- 502
k.update <- 100

y.ots <- matrix(NA, nrow = n.ots, ncol = 1) #pre-allocate memory
model.fit <- vector(mode = "list", length = length(models)) #pre-allocate memory
MS_Vola <- matrix(NA, nrow = n.ots, ncol = 1) ##pre-allocate memory
ms2.garch.s <- CreateSpec(variance.spec = list(model = "sGARCH"),
                          distribution.spec = list(distribution = "std"),
                          switch.spec = list(K = 2))
models <- list(ms2.garch.s)

#loop to create rolling forecast
for (i in 1:n.ots) {
  cat("Backtest - Iteration: ", i, "\n")
  y.its <- as.numeric(spx_log_returns[i:(n.its + i - 1)])
  y.ots[i] <- spx_log_returns[n.its + i]
  for (j in 1:length(models)) {
    if (k.update == 1 || i %% k.update == 1) {
      cat("Model", j, "is reestimated\n")
      model.fit[[j]] <- FitML(spec = models[[j]], data = y.its,
                                ctr = list(do.se = FALSE))
    }
    MS_Vola[i] = predict(model.fit[[j]]$spec, par = model.fit[[j]]$par,
                           newdata = y.its)
    }
  }





