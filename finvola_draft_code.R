#Volatility Modeling and NOPE
#Jonas Schmitten, Moritz Köhler, Giovanni Magagnin, Noah Angara
#May 2021

# Packages ----------------------------------------------------------------
library(fGarch)
library(tidyquant)
library(fBasics)
library(strucchange)
library(vars)
library(tseries)
library(dplyr)
library(tidyr)
library(tibble)
library(rugarch)
library(MSGARCH)
library(data.table)
library(stargazer)
library(hrbrthemes)
library(here)

# Data Prep & Cleaning ----------------------------------------------------

#set wd according to who is executing the code
Paths = c("/Users/jonasschmitten/Downloads/NOPE and Volatility/Data", 
          "/Users/noahangara/Documents/Master's/8th Semester/Financial Volatility",
          "C:/Users/magag/Documents/UNIVERSITA/_St.GALLEN - Master in Quant Economics and Finance/Subjects/2^ Semester/Financial volatility/project",
          "/Users/MK/Desktop/HSG/master/sem_3/fin_vola")
names(Paths) = c("jonasschmitten", "noahangara", "magag", "MK")
setwd(Paths[Sys.info()[7]])

### Loading Option and SPX Data
data_option = read.csv("options_data.csv", nrows = 100000)
sp500 = tq_get('^GSPC', get = "stock.prices", from = as.Date('2019-01-01')-750, to = '2020-12-31')

### Not interested in High or Low Price
data_spx = sp500 %>%
  subset(select = -c(high, low))


### Dataframe with SPX Returns
spx_log_returns = data.frame("Returns" = diff(log(data_spx$close)))
spx_log_returns = spx_log_returns %>%
  mutate(data.frame(sp500[2:nrow(sp500), "date"]))
spx_log_returns = xts(spx_log_returns, order.by = spx_log_returns$date)
spx_log_returns = spx_log_returns[,-2]
spx_log_returns$Returns <- apply(spx_log_returns$Returns,2,as.numeric)

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

#Get deltas ang gammas to 100 and -100 (Saw this on NOPECord)
data_option$delta = data_option$delta*100
data_option$gamma = data_option$gamma*100

#NOPE
data_option$NO = group_by(data_option, date, cp_flag)$volume * data_option$delta

NOPE =  group_by(data_option, date, cp_flag) %>% 
  summarise(NOPE = sum(NO))

NOPE = as.data.frame(rowsum(NOPE$NOPE, as.integer(gl(nrow(NOPE), 2, nrow(NOPE)))))
NOPE = add_column(NOPE, date = unique(data_option$date), .before = 1)

NOPE$volume = sp500[1:nrow(NOPE),which(colnames(sp500)=="volume")]

NOPE$NOPE = (NOPE$V1/NOPE$volume)*100

#Gamma exposure GEX
get_gex = function(yourdata){
  df = yourdata
  dates = unique(df$date)
  gammas = c()
  gammas_xday = c()
  for (i in dates) {
    sub_df = df %>% 
      filter(date == i)
    gamma_i = ifelse(sub_df$cp_flag == "C",
                     sub_df$gamma*sub_df$open_interest*100,
                     sub_df$gamma*sub_df$open_interest*-100)
    gammas = c(gammas, gamma_i)
    gammas_xday = c(gammas_xday, sum(gamma_i))
  }
  df$GEX = gammas
  df_gammas = data.frame(dates, gammas_xday)
  names(df_gammas) = c("date", "GEX")
  return(list(df, df_gammas))
}


data_option = get_gex(data_option)[[1]]
gex_xday = get_gex(data_option)[[2]]

#add price of SPX to options (first change date format)
data_option = data_option %>%
  mutate(time_to_exp = as.numeric(data_option$exdate-data_option$date)) %>%
  mutate(spx_price = lapply(as.numeric(rownames(data_option)), 
                            function(x) pull(data_spx[which(data_option$date[x] == data_spx$date), "close"])))
# Pre-Tests SPX Data ------------------------------------------------------

### Calculating basic statistics (from fBasics)
basic_stats = basicStats(as.numeric(spx_log_returns))
basic_stats = basic_stats %>%
  slice(-c(5:6, 9:12))
stargazer(t(basic_stats[1:5,]))
stargazer(t(basic_stats[6:10,]))

### Plotting histogram for the distribution
hist(as.numeric(spx_log_returns$Returns),breaks = 200, xlab = "", main = "")
abline(v = mean(as.numeric(spx_log_returns)), col="red", lwd=2)

plot(y = as.numeric(spx_log_returns$Returns), type = 'l', x = sp500[2:nrow(sp500), 'date'], ylab = "SP500 Daily Return", xlab = "")

### Normality 
normalTest(log(as.numeric(spx_log_returns)+1)*100 , method = "jb")
#QQ-plots
qqnorm(as.numeric(spx_log_returns))
qqline(as.numeric(spx_log_returns), distribution = qnorm, colour = 'r')

### Correlogram
acf(as.numeric(spx_log_returns),lag.max = 50,main = "S&P500 Log Returns")
pacf(as.numeric(spx_log_returns), lag.max = 50, main = "S&P500 Log Returns")
#Squared returns 
acf(as.numeric(spx_log_returns)^2,lag.max = 50,main = "S&P500 Log Returns Squared")
pacf(as.numeric(spx_log_returns)^2, lag.max = 50, main = "S&P500 Log Returns Squared")

### Ljung-Box test for autocorrelation
Box.test(as.numeric(spx_log_returns),lag=30,type="Ljung")
Box.test(as.numeric(spx_log_returns)^2,lag=30,type="Ljung")

### Stationarity
#Augmented-Dickey-Fuller
adf.test(as.numeric(spx_log_returns))

### Asymmetries and the Leverage Effect

### Testing for structural breaks
d = cbind(spx_log_returns, stats::lag(spx_log_returns))
fs = Fstats(as.numeric(d$Returns) ~ as.numeric(d$Returns.1))
plot(fs)
lines(breakpoints(fs))

# Volatility Forecasts ----------------------------------------------------

# rolling forecast for GARCH(1,1) with forecast length of 504 days
spec_garch = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                        distribution.model = "norm")

mod_garch = ugarchroll(spec_garch, data = spx_log_returns, n.ahead = 1, 
                       n.start = 504,  refit.every = 1, window.size= 504, refit.window=c('moving'), 
                       solver = "hybrid", fit.control = list(), keep.coef = TRUE)

fit_garch = ugarchfit(data = spx_log_returns, spec = spec_garch)
garch_news = newsimpact(z = NULL, fit_garch)
#BIC and AIC
-2*mean(mod_garch@model[["loglik"]])+log(504)*4
-2*mean(mod_garch@model[["loglik"]])+2*length(mod_garch@model[["coef"]][[1]])

# rolling forecast for GJR-GARCH with forecast length of 504 days
spec_gjrgarch = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)), 
                           mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                           distribution.model = "std")

mod_gjrgarch = ugarchroll(spec_gjrgarch, data = spx_log_returns, n.ahead = 1, 
                          n.start = 504,  refit.every = 1, window.size= 504, refit.window=c('moving'), 
                          solver = "hybrid", fit.control = list(), keep.coef = TRUE)

fit_gjrgarch = ugarchfit(data = spx_log_returns, spec = spec_gjrgarch)
gjrgarch_news = newsimpact(z = NULL, fit_gjrgarch)

-2*mean(mod_gjrgarch@model[["loglik"]])+log(504)*4
-2*mean(mod_gjrgarch@model[["loglik"]])+2*length(mod_gjrgarch@model[["coef"]][[1]])

#news impact curve
plot(gjrgarch_news$zx, gjrgarch_news$zy, ylab=gjrgarch_news$yexpr, 
     xlab=gjrgarch_news$xexpr, type="l", main = "News Impact Curve" ,ylim=c(0,0.01), lwd = 2)
lines(garch_news$zx, garch_news$zy, xlab=garch_news$xexpr, ylab=garch_news$yexpr
     , type="l", main = "News Impact Curve" ,ylim=c(0,0.01), lwd = 2, col = "red")
legend(0.03, 0.008, legend=c("GARCH", "GJR-GARCH"),
       col=c("red", "black"), lty=1, cex=0.8)

# rolling forecast for Markov-Switching GARCH with forecast length of 504 days
n.ots <- 504 #window size for forecast
n.its <- nrow(spx_log_returns)-n.ots #number of forecasts we can produce given by sample size - window size
k.update <- 1

y.ots <- matrix(NA, nrow = n.ots, ncol = 1) #pre-allocate memory
model.fit <- vector(mode = "list", length = length(models)) #pre-allocate memory
MS_Vola <- matrix(NA, nrow = n.ots, ncol = 1) ##pre-allocate memory
MS_LogLik <- matrix(NA, nrow = n.ots, ncol = 1) ##pre-allocate memory
ms2.garch.s <- CreateSpec(variance.spec = list(model = c("sGARCH","gjrGARCH")),
                          distribution.spec = list(distribution = c("norm", "std")),
                          switch.spec = list(do.mix = FALSE))
models <- list(ms2.garch.s)


#loop to create rolling forecast
for (i in 1:n.ots) {
  #indicate which i-step ahead forecast is produced
  cat("Backtest - Iteration: ", i, "\n")
  y.its <- as.numeric(spx_log_returns[i:(n.its + i - 1)])
  y.ots[i] <- spx_log_returns[n.its + i]
  for (j in 1:length(models)) {
    if (k.update == 1 || i %% k.update == 1) {
      #indicate which model is re-estimated
      cat("Model", j, "is reestimated\n")
      #estimate MS-GARCH on data
      model.fit[[j]] <- FitML(spec = models[[j]], data = y.its,
                              ctr = list(do.se = FALSE))
    }
    #add conditional vola forecast to list with MS-GARCH model spec
    MS_Vola[i] = predict(model.fit[[j]]$spec, par = model.fit[[j]]$par,
                         newdata = y.its)
    MS_LogLik[i] = model.fit[[j]][["loglik"]]
  }
}

#BIC and AIC
-2*mean(MS_LogLik) + 2*length(model.fit[1][[1]][["par"]])
-2*mean(MS_LogLik) + log(n.ots)*length(model.fit[1][[1]][["par"]])


# Black-Scholes Option Pricing --------------------------------------------
# plot delta and gamma vs strike price
plot_df = data_option[which(data_option$cp_flag == "C" & data_option$date == "2019-01-02" &
                           data_option$exdate == "2020-01-17"), c("gamma", "delta", "strike_price")]
plot_df$gamma = plot_df$gamma*1000
plot_df = melt(plot_df, "strike_price")
greeks_plot = ggplot(data = plot_df)+
  geom_smooth(aes(y = value, x = strike_price, color = variable), se = F) + 
  geom_vline(xintercept = as.numeric(data_spx[data_spx$date == "2019-01-02", "close"]), linetype = "dashed", color = "red") + 
  scale_colour_manual(values = c(gamma = "#69b3a2" , delta = rgb(0.2, 0.6, 0.9, 1)))+
  scale_y_continuous(
    
    # Features of the first axis
    name = "Option Delta",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~./1000, name="Option Gamma")
  ) + labs(x = "Strike Price", color = "") + theme_ipsum() + theme(legend.position="bottom", 
                                                                   axis.title.x = element_text(hjust=0.5), axis.title.y.left = element_text(hjust=0.5, color = rgb(0.2, 0.6, 0.9, 1)),
                                                                   axis.title.y.right = element_text(hjust=0.5, color = "#69b3a2"), legend.box.margin=margin(-15,0, 0, 0),
                                                                   axis.text.y.right = element_text(color = "#69b3a2"), axis.text.y.left = element_text(color = rgb(0.2, 0.6, 0.9, 1)))
ggsave(greeks_plot, file="greeks_plot.png", width = 14, height = 10, units = "cm")


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
y = 0.0026 #LIBOR as of 21/05/2021
data_option$risk_free = y*data_option$time_to_exp/365

#Apply BS to each row using implied volatility
data_option = data_option %>%
  mutate(BS = as.numeric(lapply(as.numeric(rownames(data_option)), 
                     function(x) black_scholes(
                       S = as.numeric(data_option$spx_price[x]),
                       K = data_option$strike_price[x],
                       y = data_option$risk_free[x],
                       m = data_option$time_to_exp[x]/365,
                       sig = data_option$impl_volatility[x],
                       call = ifelse(data_option$cp_flag[x] == "C", T, F)
                     )))/1000)

#investigate differences between mid of bid/offer and our computed BS price
data_option$acc_price = (data_option$best_offer-data_option$best_bid)/2 + data_option$best_bid
data_option$diff = as.numeric(data_option$BS) - data_option$acc_price
plot(y = data_option$diff, x = data_option$time_to_exp)

#create df of vola forecasts to add to options df
garch_vola = data.frame("garch_vola" = tail(mod_garch@forecast[["density"]][["Sigma"]], 504),
                        "date" = tail(mod_garch@model[["index"]], 504)) #504 because this is number of days from start date of options df to end
gjr_garch_vola = data.frame("gjr_garch_vola" = tail(mod_gjrgarch@forecast[["density"]][["Sigma"]], 504),
                            "date" = tail(mod_gjrgarch@model[["index"]], 504))
ms_garch_vola = data.frame("ms_garch_vola" = unlist(MS_Vola), #MS Garch doesn't appear to save dates so can't retrieve it directly
                           "date" = tail(mod_gjrgarch@model[["index"]], 504))

#joining vola forecasts into option df
data_option = data_option %>%
  mutate(garch_vola = as.numeric(lapply(as.numeric(rownames(data_option)), 
                             function(x) garch_vola[which(data_option$date[x] == garch_vola$date), "garch_vola"]))) %>%
  mutate(gjr_garch_vola = as.numeric(lapply(as.numeric(rownames(data_option)), 
                          function(x) gjr_garch_vola[which(data_option$date[x] == gjr_garch_vola$date), "gjr_garch_vola"]))) %>%
  mutate(ms_garch_vola = as.numeric(lapply(as.numeric(rownames(data_option)), 
                                 function(x) ms_garch_vola[which(data_option$date[x] == ms_garch_vola$date), "ms_garch_vola"])))

#BS Price with GARCH, GJR-GARCH and MS-GARCH Volatility
vola_models = c("garch", "gjr_garch", "ms_garch")
for (i in vola_models) {
data_option = data_option %>%
  mutate("BS_{i}" := as.numeric(lapply(as.numeric(rownames(data_option)), 
                                function(x) black_scholes(
                                  S = as.numeric(data_option$spx_price[x]),
                                  K = data_option$strike_price[x],
                                  y = data_option$risk_free[x],
                                  m = data_option$time_to_exp[x]/365,
                                  sig = data_option[x, paste0(i, "_vola")]*sqrt(data_option$time_to_exp[x]),#scale vola over maturity of option
                                  call = ifelse(data_option$cp_flag[x] == "C", T, F)
                                )))/1000)
}

#Calculate Option Greeks with forecasted B&S Prices
delta_fct <- function(S, K, y, m, sig, call = T){
  d_1 <- (log(S/K) + (y + (sig^2)/2) * m) / (sig * sqrt(m))
  if (call == T){#Call Option
    option_delta <- pnorm(d_1)
  }
  else {#Put Option
    option_delta <- (pnorm(d_1)-1)
  }
  return(option_delta)
}

#function to compute option gamma (doesn't matter if put or call)
gamma_fct <- function(S, K, y, m, sig){
  d_1 <- (log(S/K) + (y + (sig^2)/2) * m) / (sig * sqrt(m))
  option_gamma <- dnorm(d_1)/(S*sig*sqrt(m))#density not cdf
}

#join option delta into option df
for (i in vola_models) {
  data_option = data_option %>%
    mutate("delta_{i}" := as.numeric(lapply(as.numeric(rownames(data_option)), 
                                         function(x) delta_fct(
                                           S = as.numeric(data_option$spx_price[x]),
                                           K = data_option$strike_price[x],
                                           y = data_option$risk_free[x],
                                           m = data_option$time_to_exp[x]/365,
                                           sig = data_option[x, paste0(i, "_vola")]*sqrt(data_option$time_to_exp[x]),#scale vola over maturity of option
                                           call = ifelse(data_option$cp_flag[x] == "C", T, F)
                                         )))*100)
}

#join option gamma into option df
for (i in vola_models) {
  data_option = data_option %>%
    mutate("gamma_{i}" := as.numeric(lapply(as.numeric(rownames(data_option)), 
                                            function(x) gamma_fct(
                                              S = as.numeric(data_option$spx_price[x]),
                                              K = data_option$strike_price[x],
                                              y = data_option$risk_free[x],
                                              m = data_option$time_to_exp[x]/365,
                                              sig = data_option[x, paste0(i, "_vola")]*sqrt(data_option$time_to_exp[x]))
                                            ))/10)
}

#join option delta with implied vola into option df
data_option = data_option %>%
    mutate("delta_implied_vola" = as.numeric(lapply(as.numeric(rownames(data_option)), 
                                            function(x) delta_fct(
                                              S = as.numeric(data_option$spx_price[x]),
                                              K = data_option$strike_price[x],
                                              y = data_option$risk_free[x],
                                              m = data_option$time_to_exp[x]/365,
                                              sig = data_option[x, "impl_volatility"],#scale vola over maturity of option
                                              call = ifelse(data_option$cp_flag[x] == "C", T, F)
                                            )))*100)


#join option gamma with implied vola into option df
data_option = data_option %>%
    mutate("gamma_implied_vola" = as.numeric(lapply(as.numeric(rownames(data_option)), 
                                            function(x) gamma_fct(
                                              S = as.numeric(data_option$spx_price[x]),
                                              K = data_option$strike_price[x],
                                              y = data_option$risk_free[x],
                                              m = data_option$time_to_exp[x]/365,
                                              sig = data_option[x, "impl_volatility"])))/10)

#for visual analysis
plot(x = data_option[which(data_option$exdate == "2019-01-04" & data_option$cp_flag == "C"), "strike_price"], 
     y = data_option[which(data_option$exdate == "2019-01-04" & data_option$cp_flag == "C"), "delta_ms_garch"], col = "red")

plot(x = data_option[which(data_option$exdate == "2019-01-04" & data_option$cp_flag == "C"), "strike_price"], 
     y = data_option[which(data_option$exdate == "2019-01-04" & data_option$cp_flag == "C"), "delta_implied_vola"])

