# Volatility Modeling and NOPE
# Jonas Schmitten, Moritz Köhler, Giovanni Magagnin, Noah Angara
# May 2021

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

# set working directory according to who is executing the code
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
# remove not needed columns
data_option = data_option %>%
  subset(select = -c(secid, index_flag, issue_type, issuer, exercise_style, optionid, contract_size))

# convert date and expiration date to date value
data_option[c('date', 'exdate')] =  lapply(data_option[c('date', 'exdate')], function(x) as.Date(as.character(x), "%Y%m%d"))

# convert strike price, best bid, bid ask to dollars instead of cents
data_option[c('strike_price', 'best_bid', 'best_offer')] = data_option[c('strike_price', 'best_bid', 'best_offer')]/1000

# order data
data_option = data_option %>%
  arrange(date, exdate, cp_flag, strike_price)

# check if there is 0 delta
unique(data_option$delta == 0)

# remove rows with NA delta
data_option = data_option %>%
  drop_na(delta)

# get deltas and gammas to 100 (Saw this on NOPECord)
data_option$delta = data_option$delta*100
data_option$gamma = data_option$gamma*100
# rename for later (NOPE-GEX-calculations)
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

# TIME: 40 sec
data_option <- data_option %>%
  mutate(time_to_exp = as.numeric(data_option$exdate-data_option$date))

# underlying price pver date --> data_option
data_option$spx_price <- NA

for (i in which(sp500$date == "2019-01-02"):length(sp500$date)) {
  dat <- sp500$date[i]
  data_option$spx_price[which(data_option$date == dat)] <- sp500$close[i]
}


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
models_error <- data.frame(error = c("AIC", "BIC"))

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
  
  # add error-vector to the data frame for the model
  error_df[vola_model] <- c(-2 * mean(mod@model[["loglik"]])+ 2 * length(mod@model[["coef"]][[1]]),
                            -2 * mean(mod@model[["loglik"]]) + log(wind_size) * 4)
  
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


### MARKOV-SWITHCHING
# ----------------------------------------------------

# TIME: 17 min
# rolling forecast for Markov-Switching GARCH with forecast length of 504 days

# window size for forecast
n_ots <- 505
# number of forecasts we can produce given by (sample_size - window_size)
n_its <- nrow(spx_log_returns)-n_ots 

# pre-allocate memory
y_ots <- c()
MS_Vola <- c()
MS_LogLik <- c()

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

# stargazer errors
rownames(models_error) <- models_error$error
models_error <- models_error[-1]
models_error
stargazer(models_error, summary = FALSE, align = TRUE, no.space = TRUE, column.labels = c("GARCH", "GJR GARCH", "MS GARCH"))



### VOLA FORECAST DATA FRAME
# ----------------------------------------------------

vola_forecasts = data.frame("date" = head(tail(mod_sGARCH@model[["index"]], 506), 505), #504 because this is number of days from start date of options df to end
                            "garch_vola" = head(tail(mod_sGARCH@forecast[["density"]][["Sigma"]], 506), 505),
                            "gjr_garch_vola" = head(tail(mod_gjrGARCH@forecast[["density"]][["Sigma"]], 506), 505),
                            "ms_garch_vola" = as.numeric(unlist(MS_Vola)))


### VOLATILITY INTO data_option
# ----------------------------------------------------

# TIME : 1.5 min
# initiate new volatility columns in main data frame
data_option$vola_sGARCH <- NA
data_option$vola_gjrGARCH <- NA
data_option$vola_msGARCH <- NA

# loop through the unique dates to overwrite NAs
for (i in 1:length(vola_forecasts$date)) {
  dat <- vola_forecasts$date[i]
  # assign multiple rows at once
  data_option$vola_sGARCH[which(data_option$date == dat)] <- vola_forecasts$garch_vola[i]
  data_option$vola_gjrGARCH[which(data_option$date == dat)] <- vola_forecasts$gjr_garch_vola[i]
  data_option$vola_msGARCH[which(data_option$date == dat)] <- vola_forecasts$ms_garch_vola[i]
}

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


### NEW OPTION PRICES --> data_option
# ----------------------------------------------------

# BS Price with GARCH, GJR-GARCH and MS-GARCH Volatility
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

## function that calculates the NOPEs and the GEX for the different models 
## and creates a summary table for easy plotting
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

# ----------------------------------------------------------------- #
# PLOTS
# ----------------------------------------------------------------- #

### NEWS IMPACT CURVE - GARCH vs GJR-GARCH 
# ----------------------------------------------------
news_plot = ggplot() +
  geom_line(aes(x = news_gjrGARCH$zx, y = news_gjrGARCH$zy, color = "black"), show.legend = T) + 
  geom_line(aes(x = news_sGARCH$zx, y = news_sGARCH$zy, color = "red"), show.legend = T) + 
  labs(x = news_gjrGARCH$xexpr, y = news_gjrGARCH$yexpr, color = "") + 
  scale_colour_manual(labels = c("GJR GARCH", "GARCH"), values = c("red", "black"))+
  theme_bw() + theme(panel.grid = element_blank(), legend.position= c(0.8, 0.8))

ggsave(news_plot, file = "newsimpact.png", width = 14, height = 10, units = "cm")



### GREEKS --> delta & gamma (from Option Metrics) vs. strike price
# ----------------------------------------------------
# reshape data
plot_df = data_option[which(data_option$cp_flag == "C" & data_option$date == "2019-01-02" &
                              data_option$exdate == "2020-01-17"), c("gamma_given", "delta_given", "strike_price")]
plot_df$gamma_given = plot_df$gamma_given*1000
colnames(plot_df) = c("gamma", "delta", "strike")
plot_df = melt(plot_df, "strike")

# plotting
greeks_plot = ggplot(data = plot_df) +
  geom_smooth(aes(y = value, x = strike, color = variable), se = F) +
  geom_vline(xintercept = as.numeric(sp500[sp500$date == "2019-01-02", "close"]), linetype = "dashed", color = "red") +
  scale_colour_manual(values = c(gamma = "red" , delta = "black")) +
  scale_y_continuous(
    # Features of the first axis
    name = "Option Delta",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~./1000, name="Option Gamma")) + 
  labs(x = "Strike Price", color = "") + 
  theme_ipsum() + theme(legend.position= c(0.8, 0.8),
                        axis.title.x = element_text(hjust=0.5), axis.title.y.left = element_text(hjust=0.5),
                        axis.title.y.right = element_text(hjust=0.5), legend.box.margin=margin(-15,0, 0, 0),
                        axis.text.y.right = element_text(), axis.text.y.left = element_text())

# save the plot
ggsave(greeks_plot, file="greeks_plot.png", width = 14, height = 10, units = "cm")


### NOPE & GEX
# ----------------------------------------------------
# ggplot function for scatter plots
custom_scatter <- function(yourdata, model, clr = "black"){
  ggplot(data = yourdata) + geom_point(aes(x = !!as.name(model), y=close_close), size = 0.5, color = clr) +
    labs(x = unlist(str_split(model, "_"))[1], y = "1-Day Close-Close Log-Return", color = "") + 
    theme_bw() + theme(panel.grid = element_blank()) + 
    geom_vline(xintercept = 0, color = "red", size = 0.1)
}

# NOPE
custom_scatter(NOPE_GEX_Rt, "NOPE_given")
custom_scatter(NOPE_GEX_Rt, "NOPE_sGARCH")
custom_scatter(NOPE_GEX_Rt, "NOPE_gjrGARCH")
ggsave(custom_scatter(NOPE_GEX_Rt, "NOPE_msGARCH", clr = "blue"), file="NOPE_msGARCH.png", width = 14, height = 10, units = "cm")
ggsave(custom_scatter(NOPE_GEX_Rt, "NOPE_impl"), file="NOPE_implied.png", width = 14, height = 10, units = "cm")

# GEX
custom_scatter(NOPE_GEX_Rt, "GEX_given")
custom_scatter(NOPE_GEX_Rt, "GEX_sGARCH")
custom_scatter(NOPE_GEX_Rt, "GEX_gjrGARCH")
ggsave(custom_scatter(NOPE_GEX_Rt, "GEX_msGARCH", clr = "blue"), file="GEX_msGARCH.png", width = 14, height = 10, units = "cm")
ggsave(custom_scatter(NOPE_GEX_Rt, "GEX_impl"), file="GEX_implied.png", width = 14, height = 10, units = "cm")


# HISTOGRAM
# pull relevant columns vom NOPE_GEX_Rt
histo_df <- NOPE_GEX_Rt %>% 
  select(close_close, NOPE_impl, NOPE_msGARCH) %>%
  arrange(NOPE_impl, close_close)

# function that divides the returns up by range of the NOPEs
get_interval <- function(yourdata, col, low_lim = -0.3, up_lim = 0.2, step_size=0.1){
  # initiate list with returns above and blow the limits
  interval <- list(histo_df$close_close[histo_df[col] < low_lim], histo_df$close_close[histo_df[col] >= up_lim])
  sequence <- c(seq(low_lim + step_size, up_lim, step_size))
  
  # loop through the range windows and storre respective returns in 'interval' as sublist
  for (i in 1:length(sequence)) {
    interval[[i+2]] <- c(histo_df$close_close[(sequence[i]-step_size) <= histo_df[col] & histo_df[col] < sequence[i]])
  }
  
  # create data frame with a range indicator, share of positive returns 
  # & number of observations per bin
  bin_df <- data.frame("bin" = c(-0.35, 0.25, seq(low_lim+0.05, up_lim-0.05, 0.1)), 
                            "positives" = unlist(lapply(interval, function(x) sum(x > 0) / length(x))), 
                            "n"= unlist(lapply(interval, function(x) length(x))))
  
  # sort data frame by 'bin'-value
  bin_df <- bin_df %>%
    arrange(bin)
  
  # add column with custom names (for the plot)
  bin_df$bin_names <- c("<-.3", "[-.3, -.2)", "[-.2, -.1)", "[-.1, 0)", "[0, .1)", "[.1, .2)", ">.2")
  
  return(bin_df)
}

bin_df_impl <- get_interval(histo_df, "NOPE_impl")
bin_df_msGARCH <- get_interval(histo_df, "NOPE_msGARCH")

# plot bins
ggsave(ggplot(data=bin_df_impl, aes(x = bin_names, y = positives)) + 
         geom_col(width = 0.5, fill="black", alpha = 0.7) +
         coord_cartesian(ylim=c(0.45, 0.65)) +
         labs(x = "NOPE", y = "Share Positive 1-Day Close-Close") + 
         theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = .714) + 
         scale_x_discrete(limits=bin_df_msGARCH$bin_names)
       , file = "hist_impl.png", width = 14, height = 10, units = "cm")
  


ggsave(ggplot(data=bin_df_msGARCH, aes(x = bin_names, y = positives)) + 
         geom_col(width = 0.5, fill = "steelblue", alpha = 1) + 
         coord_cartesian(ylim=c(0.45, 0.65)) +
         labs(x = "NOPE", y = "Share Positive 1-Day Close-Close") + 
         theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = .714) + 
         scale_x_discrete(limits=bin_df_msGARCH$bin_names)
       , file = "hist_msGARCH.png", width = 14, height = 10, units = "cm")



### GREEKS
# ----------------------------------------------------

# CUSTOM PLOT FUNCTION --> DELTA + GAMMA
custom_deltagamma_plot <- function(yourdata, greeks){
  
  # reshape the data
  yourdata = melt(yourdata, c("exdate", "strike_price"))
  # plotting 
  ggplot(data = yourdata) + 
    geom_line(aes(x = strike_price, y = value, color = factor(exdate))) + 
    geom_vline(aes(xintercept = pull(sp500[sp500$date == "2019-01-02", "close"]), 
                   color = "green"), show.legend = F) +
    labs(x = "Strike Price", y = greeks , color = "") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position= c(0.8, 0.8)) + 
    scale_colour_manual(labels = c("2019-01-18", "2019-03-29", "2020-01-17", "SPX Close"), 
                        values = c("red", "black", "blue", "green"))

}

# DELTA (implied volatility)
ggsave(custom_deltagamma_plot(data_option[which(data_option$exdate == c("2019-01-18", "2019-03-29", "2020-01-17")  & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), 
                                          c("strike_price", "delta_impl_vola", "exdate")], "delta")
                                        , file = "implied_delta.png", width = 14, height = 10, units = "cm")

# DELTA (MS-GARCH volatility)
ggsave(custom_deltagamma_plot(data_option[which(data_option$exdate == c("2019-01-18", "2019-03-29", "2020-01-17")  & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), 
                                          c("strike_price", "delta_msGARCH", "exdate")], "delta")
                                        , file = "msgarch_delta.png", width = 14, height = 10, units = "cm")

# GAMMA (implied volatility)
ggsave(custom_deltagamma_plot(data_option[which(data_option$exdate == c("2019-01-18", "2019-03-29", "2020-01-17")  & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), 
                                          c("strike_price", "gamma_impl_vola", "exdate")], "gamma")
                                        , file = "implied_gamma.png", width = 14, height = 10, units = "cm")

# GAMMA (MS-GARCH volatility)
ggsave(custom_deltagamma_plot(data_option[which(data_option$exdate == c("2019-01-18", "2019-03-29", "2020-01-17")  & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), 
                                          c("strike_price", "gamma_msGARCH", "exdate")], "gamma")
                                        , file = "gamma_msgarch.png", width = 14, height = 10, units = "cm")

### VOLA SMILE
# ----------------------------------------------------
# function for the vola smile plot
vola_smile_plot <- function(yourdata){
  yourdata = melt(yourdata, c("exdate", "strike_price"))
  ggplot(data = yourdata) + 
    geom_line(aes(x = strike_price, y = value, color = factor(exdate))) +
    labs(x = "Strike Price", y = "Implied Volatility" , color = "") + theme_bw() + theme(panel.grid = element_blank(), legend.position= c(0.8, 0.8)) + 
    scale_colour_manual(labels = c("2019-01-14", "2019-02-08", "2020-06-21"), values = c("red", "black", "blue"))
}

ggsave(vola_smile_plot(data_option[which(data_option$exdate == c("2019-01-14", "2019-02-08", "2019-06-21")  & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), 
                                   c("strike_price", "impl_volatility", "exdate")])
                                , file = "vola_smile.png", width = 14, height = 10, units = "cm")


# plot MS forecasts vs Realized
ms_vola_plot = ggplot() + geom_line(aes(y = unlist(MS_Vola), x = tail(index(spx_log_returns), 505))) + 
  labs(x = "", y = "Volatility") + theme_bw() + theme(panel.grid = element_blank())
ggsave(ms_vola_plot, file = "ms_vola.png", width = 14, height = 10, units = "cm")


### CORRELATION & REGRESSION
# ----------------------------------------------------
## CORRELATION of NOPEs & GEX with close-close-returns
corr_nope_gex <- data.frame(NOPE = cor(NOPE_GEX_Rt[3:13])[11, 1:5], GEX = cor(NOPE_GEX_Rt[3:13])[11, 6:10])
rownames(corr_nope_gex) <- c("Option Metrics", "GARCH", "GJR GARCH", "MS GARCH", "Implied Volatility")

# stargazer correlation
stargazer(corr_nope_gex, align = TRUE, no.space = TRUE, summary=FALSE)

## REGRESSIONS

# NOPEs & GEX
summary(lm(close_close ~ GEX_given + GEX_impl + GEX_msGARCH + GEX_sGARCH + GEX_gjrGARCH, data=NOPE_GEX_Rt))
summary(lm(close_close ~ NOPE_impl+ NOPE_given + NOPE_sGARCH + NOPE_gjrGARCH + NOPE_msGARCH, data=NOPE_GEX_Rt))

# mixed (NOPEs & GEX)
summary(lm(close_close ~ NOPE_impl + GEX_impl, data=NOPE_GEX_Rt))
summary(lm(close_close ~ NOPE_sGARCH + GEX_sGARCH, data=NOPE_GEX_Rt))
summary(lm(close_close ~ NOPE_msGARCH + GEX_msGARCH, data=NOPE_GEX_Rt))
summary(lm(close_close ~ NOPE_gjrGARCH + GEX_gjrGARCH, data=NOPE_GEX_Rt))
stargazer(lm(close_close ~ NOPE_impl + GEX_impl, data=NOPE_GEX_Rt),
          lm(close_close ~ NOPE_msGARCH + GEX_msGARCH, data=NOPE_GEX_Rt), 
          no.space = TRUE, align = TRUE, covariate.labels = 
            c("NOPE Implied Volatility","GEX Implied Volatility", "NOPE MS GARCH","GEX MS GARCH", "Constant"))

# Mincer Zarnowitz Regressions
stargazer(lm(tail(as.numeric(spx_log_returns)^2, 505) ~  vola_forecasts$ms_garch_vola),
          lm(tail(as.numeric(spx_log_returns)^2, 505) ~  vola_forecasts$gjr_garch_vola),
          lm(tail(as.numeric(spx_log_returns)^2, 505) ~  vola_forecasts$garch_vola), align = T, no.space = T)
