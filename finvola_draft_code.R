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
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(tibble)
library(rugarch)
library(MSGARCH)
library(data.table)
library(stargazer)
library(hrbrthemes)
library(here)

# ----------------------------------------------------------------- #
### Data Prep & Cleaning
# ----------------------------------------------------------------- #

#set wd according to who is executing the code
Paths = c("/Users/jonasschmitten/Downloads/NOPE and Volatility/Data",
          "/Users/noahangara/Documents/Master's/8th Semester/Financial Volatility",
          "C:/Users/magag/Documents/UNIVERSITA/_St.GALLEN - Master in Quant Economics and Finance/Subjects/2^ Semester/Financial volatility/project",
          "/Users/MK/Desktop/HSG/master/sem_3/fin_vola")
names(Paths) = c("jonasschmitten", "noahangara", "magag", "MK")
setwd(Paths[Sys.info()[7]])


### Loading Option and SPX Data
# ----------------------------------------------------
data_option = read_csv("options_data.csv")

## CLEANING OPTIONS
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
names(data_option)[which(names(data_option)=="delta")] <- "delta_given"
names(data_option)[which(names(data_option)=="gamma")] <- "gamma_given"

### LOADING SP500
# ----------------------------------------------------
sp500 = tq_get('^GSPC', get = "stock.prices", from = as.Date('2019-01-01')-750, to = '2021-01-05')
sp500 <- sp500 %>% subset(select = -c(high, low))
# add log-returns (close-to-close)
sp500$log_returns <- c(NA,diff(log(sp500$close)))


# SPX returns as xts-object
spx_log_returns = data.frame("Returns" = diff(log(sp500$close)))
spx_log_returns = spx_log_returns %>%
  mutate(data.frame(sp500[2:nrow(sp500), "date"]))
spx_log_returns = xts(spx_log_returns, order.by = spx_log_returns$date)
spx_log_returns = spx_log_returns[,-2]
spx_log_returns$Returns <- apply(spx_log_returns$Returns,2,as.numeric)


### add price of SPX to options (first change date format)
# ----------------------------------------------------

# Start the clock!
start_time <- proc.time()

# TIME: 40 sec
data_option <- data_option %>%
  mutate(time_to_exp = as.numeric(data_option$exdate-data_option$date))

# underlying price pver date --> data_option
data_option$spx_price <- NA

for (i in which(sp500$date == "2019-01-02"):length(sp500$date)) {
  dat <- sp500$date[i]
  data_option$spx_price[which(data_option$date == dat)] <- sp500$close[i]
}

# Stop the clock
proc.time() - start_time



# ----------------------------------------------------------------- #
# Pre-Tests SPX Data
# ----------------------------------------------------------------- #

### Calculating basic statistics (from fBasics)
basic_stats = basicStats(as.numeric(spx_log_returns$Returns))
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


# ----------------------------------------------------------------- #
# Volatility Forecasts
# ----------------------------------------------------------------- #

### GARCH FORMULA
# ----------------------------------------------------

# create data frame with the error-type
models_error <- data.frame(error = c("BIC","AIC"))

custom_garch <- function(yourdata, error_df, vola_model = "sGARCH", model_order = c(1,1), arma_order = c(0,0), wind_size = 505) {
  
  # find model distribution
  mod_dist <- ifelse(vola_model == "sGARCH", "norm", "std")
  
  # model_sepecifics
  spec_model <- ugarchspec(variance.model = list(model = vola_model, garchOrder = model_order),
                           mean.model = list(armaOrder = arma_order, include.mean = FALSE),
                           distribution.model = mod_dist)
  
  # train the model
  mod <- ugarchroll(spec_model, data = yourdata, n.ahead = 1,
                    n.start = wind_size,  refit.every = 1, window.size= wind_size, refit.window=c('moving'),
                    solver = "hybrid", fit.control = list(), keep.coef = TRUE)
  
  # fit the model and pull news impact
  fit_model <- ugarchfit(data = yourdata, spec = spec_model)
  mod_news <- newsimpact(z = NULL, fit_model)
  
  # add error-vector for the model
  error_df[vola_model] <- c(-2 * mean(mod@model[["loglik"]]) + log(wind_size) * 4,
                            -2*mean(mod@model[["loglik"]])+2*length(mod@model[["coef"]][[1]]))
  
  # assign global variables for the model, news and the new error data frame
  assign(paste("mod", vola_model, sep="_"), mod,  envir = .GlobalEnv)
  assign(paste("news", vola_model, sep="_"), mod_news,  envir = .GlobalEnv)
  assign("models_error", error_df,  envir = .GlobalEnv)
}

### GARCH (1,1) - rolling forecast with length 504 days
# ----------------------------------------------------

# assigns 'mod_sGARCH' and 'news_sGARCH'
custom_garch(spx_log_returns, models_error, vola_model = "sGARCH")


### GJR-GARCH - rolling forecast with length 504 days
# ----------------------------------------------------

# assigns 'mod_gjrGARCH' and 'news_gjrGARCH'
custom_garch(spx_log_returns, models_error, vola_model = "gjrGARCH")


### IMPACT CURVE - GARCH vs gjr-GARCH -----------------

plot(news_gjrGARCH$zx, news_gjrGARCH$zy, ylab=news_gjrGARCH$yexpr,
     xlab=news_gjrGARCH$xexpr, type="l", main = "News Impact Curve" ,ylim=c(0,0.01), lwd = 2)
lines(news_sGARCH$zx, news_sGARCH$zy, xlab=news_sGARCH$xexpr, ylab=news_sGARCH$yexpr
      , type="l", main = "News Impact Curve" ,ylim=c(0,0.01), lwd = 2, col = "red")
legend(0.03, 0.008, legend=c("GARCH", "GJR-GARCH"),
       col=c("red", "black"), lty=1, cex=0.8)


### MARKOV-SWITHCHING
# ----------------------------------------------------

# rolling forecast for Markov-Switching GARCH with forecast length of 504 days
n_ots <- 505 #window size for forecast
n_its <- nrow(spx_log_returns)-n_ots #number of forecasts we can produce given by sample size - window size

y_ots <- c() #pre-allocate memory
MS_Vola <- c() ##pre-allocate memory
MS_LogLik <- c() ##pre-allocate memory

ms2_garch_s <- CreateSpec(variance.spec = list(model = c("gjrGARCH","gjrGARCH")),
                          distribution.spec = list(distribution = c("std", "std")),
                          switch.spec = list(do.mix = FALSE))

msGARCH <- list(ms2_garch_s)[[1]]



# loop to create rolling forecast
for (i in 1:n_ots) {
  # indicate which i-step ahead forecast is produced 
  cat("Backtest - Iteration:", i, "/ 505 \n")
  y_its <- as.numeric(spx_log_returns[i:(n_its + i - 1)])
  y_ots <- c(y_ots, spx_log_returns[n_its + i])
  
  # fit the model
  model_fit <- FitML(spec = msGARCH, data = y_its, ctr = list(do.se = FALSE))
  
  #add conditional vola forecast to list with MS-GARCH model spec
  MS_Vola = c(MS_Vola, predict(model_fit$spec, par = model_fit$par, newdata = y_its))
  
  MS_LogLik = c(MS_LogLik, model_fit$loglik)
}


# Markov Switching Errors
models_error$markov <- c(-2 * mean(MS_LogLik) + 2 * length(model_fit$par),
                         -2 * mean(MS_LogLik) + log(n_ots) * length(model_fit$par))

# see error summary of all models
print(models_error)

### VOLA FORECAST DATA FRAME
# ----------------------------------------------------

vola_forecasts = data.frame("date" = tail(mod_sGARCH@model[["index"]], 505), #504 because this is number of days from start date of options df to end
                            "garch_vola" = tail(mod_sGARCH@forecast[["density"]][["Sigma"]], 505),
                            "gjr_garch_vola" = tail(mod_gjrGARCH@forecast[["density"]][["Sigma"]], 505),
                            "ms_garch_vola" = as.numeric(unlist(MS_Vola)))

# ----------------------------------------------------------------- #
# Black-Scholes Option Pricing
# ----------------------------------------------------------------- #

# plot delta and gamma vs strike price
plot_df = data_option[which(data_option$cp_flag == "C" & data_option$date == "2019-01-02" &
                              data_option$exdate == "2020-01-17"), c("gamma_given", "delta_given", "strike_price")]
plot_df$gamma_given = plot_df$gamma_given*1000
colnames(plot_df) = c("gamma", "delta", "strike")
plot_df = melt(plot_df, "strike")
greeks_plot = ggplot(data = plot_df)+
  geom_smooth(aes(y = value, x = strike, color = variable), se = F) +
  geom_vline(xintercept = as.numeric(sp500[sp500$date == "2019-01-02", "close"]), linetype = "dashed", color = "red") +
  scale_colour_manual(values = c(gamma = "red" , delta = "black"))+
  scale_y_continuous(
    
    # Features of the first axis
    name = "Option Delta",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~./1000, name="Option Gamma")
  ) + labs(x = "Strike Price", color = "") + theme_ipsum() + theme(legend.position= c(0.8, 0.8),
                                                                   axis.title.x = element_text(hjust=0.5), axis.title.y.left = element_text(hjust=0.5),
                                                                   axis.title.y.right = element_text(hjust=0.5), legend.box.margin=margin(-15,0, 0, 0),
                                                                   axis.text.y.right = element_text(), axis.text.y.left = element_text())

ggsave(greeks_plot, file="greeks_plot.png", width = 14, height = 10, units = "cm")


greeks_plot



# ----------------------------------------------------------------- #
# BLACK & SCHOLES
# ----------------------------------------------------------------- #

### BS FROMULA
# ----------------------------------------------------
black_scholes <- function(S, K, y, m, sig, call = TRUE) {
  d1 <- (log(S/K) + (y + (sig^2)/2) * m) / (sig * sqrt(m))
  d2 <- d1 - sig * sqrt(m)
  
  # option price dependding on CALL/PUT
  C <- ifelse(call ==TRUE, 
              S * pnorm(d1) - exp(-y*m) * K * pnorm(d2), 
              exp(-y*m) * K * pnorm(-d2) - S * pnorm(-d1))
  
  return(C)
}


#define risk-free as 3-month T-Bill
y <- 0.0026 #LIBOR as of 21/05/2021
data_option$risk_free <- y*data_option$time_to_exp/365


### DELTA FORMULA
# ----------------------------------------------------

delta_fct <- function(S, K, y, m, sig, call = T){
  d1 <- (log(S/K) + (y + (sig^2)/2) * m) / (sig * sqrt(m))
  
  # DELTA depending on CALL/PUT
  option_delta <- ifelse(call == T, pnorm(d1), pnorm(d1)-1)
  
  return(option_delta)
}


### GAMMA FORMULA
# ----------------------------------------------------

gamma_fct <- function(S, K, y, m, sig){
  d1 <- (log(S/K) + (y + (sig^2)/2) * m) / (sig * sqrt(m))
  option_gamma <- dnorm(d1)/(S*sig*sqrt(m))#density not cdf
}


### BS-PRICES with IMPLIED VOLATILITY --> data_option
# ----------------------------------------------------
# TIME: 2 sec
data_option$BS <- black_scholes(S = as.numeric(data_option$spx_price), 
                                K = data_option$strike_price, 
                                y = data_option$risk_free, 
                                m = data_option$time_to_exp/365,
                                sig = data_option$impl_volatility, 
                                call = ifelse(data_option$cp_flag == "C", T, F)) / 1000


# #investigate differences between mid of bid/offer and our computed BS price
# data_option$acc_price = (data_option$best_offer-data_option$best_bid)/2 + data_option$best_bid
# data_option$diff = as.numeric(data_option$BS) - data_option$acc_price
# plot(y = data_option$diff, x = data_option$time_to_exp)



### VOLATILITY INTO data_option
# ----------------------------------------------------

# TIME : 1.5 min
data_option$vola_sGARCH <- NA
data_option$vola_gjrGARCH <- NA
data_option$vola_msGARCH <- NA

for (i in 1:length(vola_forecasts$date)) {
  dat <- vola_forecasts$date[i]
  data_option$vola_sGARCH[which(data_option$date == dat)] <- vola_forecasts$garch_vola[i]
  data_option$vola_gjrGARCH[which(data_option$date == dat)] <- vola_forecasts$gjr_garch_vola[i]
  data_option$vola_msGARCH[which(data_option$date == dat)] <- vola_forecasts$ms_garch_vola[i]
}

### NEW OPTION PRICES --> data_option
# ----------------------------------------------------

# TIME: 6 sec
# BS Price with GARCH, GJR-GARCH and MS-GARCH Volatility (MISSING!!!!!!!!!!!)
vola_models = c("sGARCH", "gjrGARCH", "msGARCH")

# define row numbers with right date
row_nums <- which(data_option$date != "2019-01-02")

for (i in vola_models) {
  
  # create empty column
  data_option[paste0("BS_", i)] <- NA
  
  # fill the empty column only for the right dates
  data_option[[paste0("BS_", i)]] <- black_scholes(S = as.numeric(data_option$spx_price),
                                                   K = data_option$strike_price, 
                                                   y = data_option$risk_free, 
                                                   m = data_option$time_to_exp/365,
                                                   sig = data_option[[paste0("vola_", i)]] * sqrt(data_option$time_to_exp), # scale volatility over maturity of option
                                                   call = ifelse(data_option$cp_flag == "C", T, F)) / 1000
  }


### DELTA (models) --> data_option
# ----------------------------------------------------

for (i in vola_models) {
  
  # create empty column
  data_option[paste0("delta_", i)] <- NA
  
  # fill the empty column only for the right dates
  data_option[[paste0("delta_", i)]] <- delta_fct(S = as.numeric(data_option$spx_price), 
                                                            K = data_option$strike_price, 
                                                            y = data_option$risk_free, 
                                                            m = data_option$time_to_exp/365,
                                                            sig = data_option[[paste0("vola_", i)]] * sqrt(data_option$time_to_exp), #scale vola over maturity of option
                                                            call = ifelse(data_option$cp_flag == "C", T, F)) * 100
  
 
}


### DELTA (implied volatility) --> data_option
# ----------------------------------------------------

#join option delta with implied vola into option df
data_option$delta_impl_vola <- delta_fct(S = as.numeric(data_option$spx_price), 
                                                   K = data_option$strike_price, 
                                                   y = data_option$risk_free, 
                                                   m = data_option$time_to_exp/365,
                                                   sig = data_option$impl_volatility,
                                                   call = ifelse(data_option$cp_flag == "C", T, F)) * 100





### GAMMA (models) --> data_option
# ----------------------------------------------------

# TIME: 4 sec
# gammas with GARCH, GJR-GARCH and MS-GARCH Volatility
for (i in vola_models) {
  # create empty column
  data_option[paste0("gamma_", i)] <- NA
  
  # fill the empty column only for the right dates
  data_option[[paste0("gamma_", i)]] <- gamma_fct(S = as.numeric(data_option$spx_price), 
                                                            K = data_option$strike_price, 
                                                            y = data_option$risk_free, 
                                                            m = data_option$time_to_exp/365,
                                                            sig = data_option[[paste0("vola_", i)]] * sqrt(data_option$time_to_exp)) * 100
}

### GAMMA (implied volatility) --> data_option
# ----------------------------------------------------

data_option$gamma_impl_vola <- gamma_fct(S = as.numeric(data_option$spx_price), 
                                                   K = data_option$strike_price, 
                                                   y = data_option$risk_free, 
                                                   m = data_option$time_to_exp/365,
                                                   sig = data_option$impl_volatility) * 100


### NOPE & GEX
# ----------------------------------------------------

## NOPE & GEX
get_nope_gex <- function(opt_data, underl_data) {
  
  # initialize df_NOPE
  df_NOPE <- underl_data[c("date", "volume")][which(underl_data$date %in% unique(opt_data$date)),]
  
  ### NOPE
  # find delta columns and isolate model
  delta_cols <- names(opt_data)[which(str_detect(names(data_option), "delta"))]
  models <- str_split(delta_cols, "_")
  
  # loop through the different delta columns
  for (col in delta_cols) {
    mod <-  models[[which(delta_cols == col)]][2]
    
    # required hedging position
    opt_data[[paste0("NO_", mod)]] <- opt_data$volume * opt_data[[col]]
    
    # sum up the hedging positions per day
    sub_df <- group_by(opt_data, date) %>%
      summarise(nope_vol = sum(!!as.name(paste0("NO_", mod))))
    
    # calculate the NOPE per day and store it as column in 'df_NOPE'
    df_NOPE[paste0("NOPE_", mod)] <- 100 * (sub_df$nope_vol/df_NOPE$volume)
  }
  
  ### GEX
  gamma_cols <- names(opt_data)[which(str_detect(names(data_option), "gamma"))]
  models <- str_split(gamma_cols, "_")
  for (col in gamma_cols) {
    mod <- models[[which(gamma_cols == col)]][2]
    
    opt_data[[paste0("GEX_", mod)]] <- opt_data$open_interest * opt_data[[col]]
    opt_data[[paste0("GEX_", mod)]][which(opt_data$cp_flag == "P")] <- -1 *
                opt_data[[paste0("GEX_", mod)]][which(opt_data$cp_flag == "P")]
    
    sub_df <- group_by(opt_data, date) %>%
      summarise(!!as.name(paste0("GEX_", mod)) := sum(!!as.name(paste0("GEX_", mod))))
    
    df_NOPE[[paste0("GEX_", mod)]] <- sub_df[[paste0("GEX_", mod)]]
  }
  
  df_NOPE$close_close <- underl_data$log_returns[which(underl_data$date %in% df_NOPE$date)+1]
  
  # returns the options data frame with new NOPE column
  # and a data frame with the date and corresponding NOPE
  return(list(opt_data, df_NOPE))
}

data_option <- get_nope_gex(data_option, sp500)[[1]]
NOPE_GEX_Rt <- get_nope_gex(data_option, sp500)[[2]]

### CORRELATION & REGRESSION
# ----------------------------------------------------
corr_nope <- data.frame(cor(NOPE_GEX_Rt[3:13], use="complete.obs"))[11,1:5]
names(corr_nope) <- c("SqueezeMetrics", "GARCH", "gjrGARCJ", "msGARCH", "ImpliedVolatility")
rownames(corr_nope) <- c("1-Day Close-Close")

corr_gex <- data.frame(cor(NOPE_GEX_Rt[3:13], use="complete.obs"))[11,6:10]
names(corr_gex) <- c("SqueezeMetrics", "GARCH", "gjrGARCJ", "msGARCH", "ImpliedVolatility")
rownames(corr_gex) <- c("1-Day Close-Close")

# stargazer
stargazer(corr_nope,  summary = FALSE)
stargazer(corr_gex,  summary = FALSE)

# regressions
# separate
summary(lm(close_close ~ GEX_given + GEX_impl + GEX_msGARCH + GEX_sGARCH + GEX_gjrGARCH, data=NOPE_GEX_Rt))
summary(lm(close_close ~ NOPE_impl+ NOPE_given + NOPE_sGARCH + NOPE_gjrGARCH + NOPE_msGARCH, data=NOPE_GEX_Rt))

# mixed
summary(lm(close_close ~ NOPE_impl + GEX_impl, data=NOPE_GEX_Rt))
summary(lm(close_close ~ NOPE_sGARCH + GEX_sGARCH, data=NOPE_GEX_Rt))
summary(lm(close_close ~ NOPE_msGARCH + GEX_msGARCH, data=NOPE_GEX_Rt))
summary(lm(close_close ~ NOPE_gjrGARCH + GEX_gjrGARCH, data=NOPE_GEX_Rt))

# ----------------------------------------------------------------- #
# PLOTS
# ----------------------------------------------------------------- #

### NOPE & GEX
# ----------------------------------------------------
# ggplot function
custom_scatter <- function(yourdata, model){
  ggplot(data = yourdata) + geom_point(aes(x = !!as.name(model), y=close_close), size = 0.5) +
    labs(x = unlist(str_split(model, "_"))[1], y = "1-Day Close-Close Log-Return", color = "") + 
    theme_bw() + theme(panel.grid = element_blank()) + 
    geom_vline(xintercept = 0, color = "red", size = 0.1)
}

# NOPE
custom_scatter(NOPE_GEX_Rt, "NOPE_given")
custom_scatter(NOPE_GEX_Rt, "NOPE_sGARCH")
custom_scatter(NOPE_GEX_Rt, "NOPE_gjrGARCH")
ggsave(custom_scatter(NOPE_GEX_Rt, "NOPE_msGARCH"), file="NOPE_msGARCH.png", width = 14, height = 10, units = "cm")
ggsave(custom_scatter(NOPE_GEX_Rt, "NOPE_impl"), file="NOPE_implied.png", width = 14, height = 10, units = "cm")

# GEX
custom_scatter(NOPE_GEX_Rt, "GEX_given")
custom_scatter(NOPE_GEX_Rt, "GEX_sGARCH")
custom_scatter(NOPE_GEX_Rt, "GEX_gjrGARCH")
ggsave(custom_scatter(NOPE_GEX_Rt, "GEX_msGARCH"), file="GEX_msGARCH.png", width = 14, height = 10, units = "cm")
ggsave(custom_scatter(NOPE_GEX_Rt, "GEX_impl"), file="GEX_implied.png", width = 14, height = 10, units = "cm")


### GREEKS
# ----------------------------------------------------

summary(data_option$delta_impl_vola)
#plotting delta with implied volatility
plot(x = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
     y = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_impl_vola"], col = "red", type = "l", lwd = 2, ylab = "Delta", xlab = "")
lines(x = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_impl_vola"], col = "blue", type = "l", lwd = 2)
lines(x = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_impl_vola"], col = "black", type = "l", lwd = 2)
abline(v = sp500[sp500$date == "2019-01-02", "close"], col = "green")
legend(3200, 80, legend=c("2019-01-18", "2019-03-29", "2020-01-17", "SPX Close"),
       col=c("red", "blue", "black", "green"), lty=1, cex=0.8, lwd = 2)

#plotting delta with MS-GARCH volatility
plot(x = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
     y = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_msGARCH"], col = "red", type = "l", lwd = 2, ylab = "Delta", xlab = "")
lines(x = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_msGARCH"], col = "blue", type = "l", lwd = 2)
lines(x = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_msGARCH"], col = "black", type = "l", lwd = 2)
abline(v = data_spx[data_spx$date == "2019-01-02", "close"], col = "green")
legend(3200, 80, legend=c("2019-01-18", "2019-03-29", "2020-01-17", "SPX Close"),
       col=c("red", "blue", "black", "green"), lty=1, cex=0.8, lwd = 2)

#plotting gamma with implied volatility
plot(x = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
     y = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "gamma_implied_vola"], col = "red", type = "l", lwd = 2, ylab = "Gamma", xlab = "")
lines(x = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "gamma_implied_vola"], col = "blue", type = "l", lwd = 2)
lines(x = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "gamma_implied_vola"], col = "black", type = "l", lwd = 2)
abline(v = data_spx[data_spx$date == "2019-01-02", "close"], col = "green")
legend(3400, 0.0002, legend=c("2019-01-18", "2019-03-29", "2020-01-17", "SPX Close"),
       col=c("red", "blue", "black", "green"), lty=1, cex=0.8, lwd = 2)

#plotting gamma with MS-GARCH volatility
plot(x = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
     y = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "gamma_ms_garch"], col = "red", type = "l", lwd = 2, ylab = "Gamma", xlab = "")
lines(x = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "gamma_ms_garch"], col = "blue", type = "l", lwd = 2)
lines(x = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "gamma_ms_garch"], col = "black", type = "l", lwd = 2)
abline(v = data_spx[data_spx$date == "2019-01-02", "close"], col = "green")
legend(3400, 0.001, legend=c("2019-01-18", "2019-03-29", "2020-01-17", "SPX Close"),
       col=c("red", "blue", "black", "green"), lty=1, cex=0.8, lwd = 2)

plot(x = data_option[which(data_option$exdate == "2019-01-14" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
     y = data_option[which(data_option$exdate == "2019-01-14" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "impl_volatility"], col = "red", type = "l", lwd = 2, ylab = "Implied Volatility", xlab = "")
lines(x = data_option[which(data_option$exdate == "2019-02-08" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2019-02-08" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "impl_volatility"], col = "blue", type = "l", lwd = 2)
lines(x = data_option[which(data_option$exdate == "2019-06-21" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2019-06-21" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "impl_volatility"], col = "black", type = "l", lwd = 2)
legend("topleft", legend=c("2019-01-14", "2019-02-08", "2019-06-21"),
       col=c("red", "blue", "black"), lty=1, cex=0.8, lwd = 2)

#plot MS forecasts vs Realized
plot(y = unlist(MS_Vola), x = tail(index(spx_log_returns), 504), type = "l", col = "black", ylab = "Volatility", xlab = "", lwd = 2)
plot(y = unlist(MS_Vola)-tail(unlist(list(rollapplyr(spx_log_returns, 2, sd, na.rm = T))), 504), x = tail(index(spx_log_returns), 504),
     ylab = "Est. Difference to Realized Volatility", xlab = "", type ="l", col = "black", lwd = 2)
