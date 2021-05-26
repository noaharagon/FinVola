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


### LOADING SP500
# ----------------------------------------------------
sp500 = tq_get('^GSPC', get = "stock.prices", from = as.Date('2019-01-01')-750, to = '2021-01-01')
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

# TIME: 7 min!
data_option <- data_option %>%
  mutate(time_to_exp = as.numeric(data_option$exdate-data_option$date)) %>%
  mutate(spx_price = lapply(as.numeric(rownames(data_option)),
                            function(x) sp500$close[which(data_option$date[x] == sp500$date)]))

# Stop the clock
proc.time() - start_time


### NOPE & GEX
# ----------------------------------------------------

## NOPE
get_nope <- function(opt_data, underl_data) {
  # required hedging position per contract
  opt_data$NO <- opt_data$volume * opt_data$delta
  
  # caclulate the NOPE per day
  df_NOPE <- group_by(opt_data, date) %>%
    summarise(nope_vol = sum(NO))
  df_NOPE$volume <- underl_data$volume[which(underl_data$date %in% df_NOPE$date)]
  df_NOPE$NOPE <- 100 * (df_NOPE$nope_vol/df_NOPE$volume)
  
  # returns the options data frame with new NOPE column
  # and a data frame with the date and corresponding NOPE
  return(list(opt_data, df_NOPE[, -which(names(df_NOPE) %in% c("nope_vol", "volume"))]))
}

## GEX
get_gex <- function(opt_data){
  # pull the unique dates from the option data
  dates <- unique(opt_data$date)
  gammas <- c()
  gex_xday <- c()
  for (i in dates) {
    sub_df <- opt_data %>%
      filter(date == i)
    gamma_i <- ifelse(sub_df$cp_flag == "C",
                      sub_df$gamma*sub_df$open_interest*100,
                      sub_df$gamma*sub_df$open_interest*-100)
    gammas <- c(gammas, gamma_i)
    gex_xday <- c(gex_xday, sum(gamma_i))
  }
  
  opt_data$GEX <- gammas
  
  # returns the options data frame with new GEX column
  # and a vector with the calculated gammas per day
  return(list(opt_data, gex_xday))
}


## NOPE + GEX + CLOSE-CLOSE-RETURNS
NOPE_GEX_Rt <- function(opt_data, underl_data){
  
  # NOPE
  opt_df <- get_nope(opt_data, underl_data)[[1]]
  nope_gex_df <- get_nope(opt_data, underl_data)[[2]]
  
  # GEX
  opt_df <- get_gex(opt_df)[[1]]
  nope_gex_df$GEX <- get_gex(opt_df)[[2]]
  
  # CLOSE-CLOSE
  nope_gex_df$close_close <- underl_data$log_returns[which(underl_data$date %in% nope_gex_df$date)+1]
  
  # return the new option data and NOPE_GEX
  return(list(opt_df, nope_gex_df))
}

# Start the clock!
start_time <- proc.time()

# TIME: 2.5 min
output_list <- NOPE_GEX_Rt(data_option, sp500)
data_option <- output_list[[1]]
sum_nope_gex <- output_list[[2]]

# Stop the clock
proc.time() - start_time


# plot NOPE -> close-close
plot(sum_nope_gex$NOPE, sum_nope_gex$close_close)

# plot GEX -> close-close
plot(sum_nope_gex$GEX, sum_nope_gex$close_close)


# ----------------------------------------------------------------- #
# GRAVEYARD NOPE - START
# ----------------------------------------------------------------- #

### Not interested in High or Low Price
# data_spx = sp500 %>%
#   subset(select = -c(high, low))
#

# # OLD NOPE
# data_option$NO = group_by(data_option, date, cp_flag)$volume * data_option$delta
#
# NOPE =  group_by(data_option, date, cp_flag) %>%
#   summarise(NOPE = sum(NO))
#
# NOPE = as.data.frame(rowsum(NOPE$NOPE, as.integer(gl(nrow(NOPE), 2, nrow(NOPE)))))
# NOPE = add_column(NOPE, date = unique(data_option$date), .before = 1)
#
# NOPE$volume = sp500[which(sp500$date %in% dates),which(colnames(sp500)=="volume")]
#
# NOPE$NOPE = (NOPE$V1/NOPE$volume)*100

# ----------------------------------------------------------------- #
# GRAVEYARD END
# ----------------------------------------------------------------- #




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

custom_garch <- function(yourdata, error_df, vola_model = "sGARCH", model_order = c(1,1), arma_order = c(0,0), wind_size = 504) {
  
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
n_ots <- 504 #window size for forecast
n_its <- nrow(spx_log_returns)-n_ots #number of forecasts we can produce given by sample size - window size
k_update <- 1

y_ots <- matrix(NA, nrow = n_ots, ncol = 1) #pre-allocate memory
MS_Vola <- matrix(NA, nrow = n_ots, ncol = 1) ##pre-allocate memory
MS_LogLik <- matrix(NA, nrow = n_ots, ncol = 1) ##pre-allocate memory
ms2_garch_s <- CreateSpec(variance.spec = list(model = c("sGARCH","gjrGARCH")),
                          distribution.spec = list(distribution = c("norm", "std")),
                          switch.spec = list(do.mix = FALSE))

models <- list(ms2_garch_s)
model_fit <- vector(mode = "list", length = length(models)) #pre-allocate memory


# loop to create rolling forecast
for (i in 1:n_ots) {
  #indicate which i-step ahead forecast is produced
  cat("Backtest - Iteration: ", i, "\n")
  y.its <- as.numeric(spx_log_returns[i:(n_its + i - 1)])
  y_ots[i] <- spx_log_returns[n_its + i]
  for (j in 1:length(models)) {
    if (k_update == 1 || i %% k_update == 1) {
      #indicate which model is re-estimated
      cat("Model", j, "is reestimated\n")
      #estimate MS-GARCH on data
      model_fit[[j]] <- FitML(spec = models[[j]], data = y.its,
                              ctr = list(do.se = FALSE))
    }
    #add conditional vola forecast to list with MS-GARCH model spec
    MS_Vola[i] = predict(model_fit[[j]]$spec, par = model_fit[[j]]$par,
                         newdata = y.its)
    MS_LogLik[i] = model_fit[[j]][["loglik"]]
  }
}

# Markov Switching Errors

models_error$markov <- c(-2*mean(MS_LogLik) + 2*length(model_fit[1][[1]][["par"]]),
                         -2*mean(MS_LogLik) + log(n.ots)*length(model_fit[1][[1]][["par"]]))

# see error summary of all models
print(models_error)




# create df of vola forecasts to add to options df
vola_forecasts = data.frame("date" = tail(mod_garch@model[["index"]], 504), #504 because this is number of days from start date of options df to end
                            "garch_vola" = tail(mod_garch@forecast[["density"]][["Sigma"]], 504),
                            "gjr_garch_vola" = tail(mod_gjrgarch@forecast[["density"]][["Sigma"]], 504),
                            "ms_garch_vola" = unlist(MS_Vola)
)

# gjr_garch_vola = data.frame("gjr_garch_vola" = tail(mod_gjrgarch@forecast[["density"]][["Sigma"]], 504),
#                             "date" = tail(mod_gjrgarch@model[["index"]], 504))
# ms_garch_vola = data.frame("ms_garch_vola" = unlist(MS_Vola), #MS Garch doesn't appear to save dates so can't retrieve it directly
#                            "date" = tail(mod_gjrgarch@model[["index"]], 504))




# ----------------------------------------------------------------- #
# GRAVEYARD MODELS - START
# ----------------------------------------------------------------- #
##### GARCH
# spec_garch = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
#                         mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
#                         distribution.model = "norm")
#
# mod_garch = ugarchroll(spec_garch, data = spx_log_returns, n.ahead = 1,
#                        n.start = 504,  refit.every = 1, window.size= 504, refit.window=c('moving'),
#                        solver = "hybrid", fit.control = list(), keep.coef = TRUE)
#
# fit_garch = ugarchfit(data = spx_log_returns, spec = spec_garch)
# garch_news = newsimpact(z = NULL, fit_garch)
# -2*mean(mod_garch@model[["loglik"]])+log(504)*4
# -2*mean(mod_garch@model[["loglik"]])+2*length(mod_garch@model[["coef"]][[1]])



##### gjrGARCH
# spec_gjrgarch = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
#                            mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
#                            distribution.model = "std")
#
# mod_gjrgarch = ugarchroll(spec_gjrgarch, data = spx_log_returns, n.ahead = 1,
#                           n.start = 504,  refit.every = 1, window.size= 504, refit.window=c('moving'),
#                           solver = "hybrid", fit.control = list(), keep.coef = TRUE)
#
# fit_gjrgarch = ugarchfit(data = spx_log_returns, spec = spec_gjrgarch)
# gjrgarch_news = newsimpact(z = NULL, fit_gjrgarch)

# -2*mean(mod_gjrgarch@model[["loglik"]])+log(504)*4
# -2*mean(mod_gjrgarch@model[["loglik"]])+2*length(mod_gjrgarch@model[["coef"]][[1]])


# ----------------------------------------------------------------- #
# GRAVEYARD - END
# ----------------------------------------------------------------- #





# ----------------------------------------------------------------- #
# Black-Scholes Option Pricing
# ----------------------------------------------------------------- #

# plot delta and gamma vs strike price
plot_df = data_option[which(data_option$cp_flag == "C" & data_option$date == "2019-01-02" &
                              data_option$exdate == "2020-01-17"), c("gamma", "delta", "strike_price")]
plot_df$gamma = plot_df$gamma*1000
plot_df = melt(plot_df, "strike_price")
greeks_plot = ggplot(data = plot_df)+
  geom_smooth(aes(y = value, x = strike_price, color = variable), se = F) +
  geom_vline(xintercept = as.numeric(sp500[sp500$date == "2019-01-02", "close"]), linetype = "dashed", color = "red") +
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


greeks_plot



# ----------------------------------------------------------------- #
# BLACK & SCHOLES
# ----------------------------------------------------------------- #

### BS FROMULA
# ----------------------------------------------------
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
y <- 0.0026 #LIBOR as of 21/05/2021
data_option$risk_free <- y*data_option$time_to_exp/365


# Start the clock
start_time <- proc.time()

# TIME:
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

# Stop the clock
proc.time() - start_time

#investigate differences between mid of bid/offer and our computed BS price
data_option$acc_price = (data_option$best_offer-data_option$best_bid)/2 + data_option$best_bid
data_option$diff = as.numeric(data_option$BS) - data_option$acc_price
plot(y = data_option$diff, x = data_option$time_to_exp)



#joining vola forecasts into option df
data_option = data_option %>%
  mutate(garch_vola = as.numeric(lapply(as.numeric(rownames(data_option)),
                                        function(x) garch_vola[which(data_option$date[x] == garch_vola$date), "garch_vola"]))) %>%
  mutate(gjr_garch_vola = as.numeric(lapply(as.numeric(rownames(data_option)),
                                            function(x) gjr_garch_vola[which(data_option$date[x] == gjr_garch_vola$date), "gjr_garch_vola"]))) %>%
  mutate(ms_garch_vola = as.numeric(lapply(as.numeric(rownames(data_option)),
                                           function(x) ms_garch_vola[which(data_option$date[x] == ms_garch_vola$date), "ms_garch_vola"])))




# Start the clock
start_time <- proc.time()


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
                                           sig = data_option[x, paste0(i, "_vola")]*sqrt(data_option$time_to_exp[x]), #scale vola over maturity of option
                                           call = ifelse(data_option$cp_flag[x] == "C", T, F)
                                         )))/1000)
}

# Stop the clock
proc.time() - start_time

### DELTA
# ----------------------------------------------------

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


### GAMMA
# ----------------------------------------------------

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

#plotting delta with implied volatility
plot(x = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
     y = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_implied_vola"], col = "red", type = "l", lwd = 2, ylab = "Delta", xlab = "")
lines(x = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_implied_vola"], col = "blue", type = "l", lwd = 2)
lines(x = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_implied_vola"], col = "black", type = "l", lwd = 2)
abline(v = data_spx[data_spx$date == "2019-01-02", "close"], col = "green")
legend(3200, 80, legend=c("2019-01-18", "2019-03-29", "2020-01-17", "SPX Close"),
       col=c("red", "blue", "black", "green"), lty=1, cex=0.8, lwd = 2)

#plotting delta with MS-GARCH volatility
plot(x = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
     y = data_option[which(data_option$exdate == "2019-01-18" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_ms_garch"], col = "red", type = "l", lwd = 2, ylab = "Delta", xlab = "")
lines(x = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2019-03-29" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_ms_garch"], col = "blue", type = "l", lwd = 2)
lines(x = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "strike_price"],
      y = data_option[which(data_option$exdate == "2020-01-17" & data_option$cp_flag == "C" & data_option$date == "2019-01-02"), "delta_ms_garch"], col = "black", type = "l", lwd = 2)
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
