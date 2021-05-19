library(tidyquant)
library(fBasics)


### Loading data
data_option = read.csv("options_data.csv")
sp500 = getYahooData('^GSPC',start=20190102,end=20201231,freq='daily')

### Removing useless columns
data_spx = subset(sp500, select = -c(High, Low))
data_option = subset(data_option, select = -c(secid, index_flag, issue_type, issuer, exercise_style))

### Calculating spx log-returns, calculating net returns first
spx_net_returns = diff(data_spx$Close)/lag(data_spx$Close)
spx_log_returns = log(1+spx_net_returns)
# removing the first row of NA
spx_log_returns = spx_log_returns[-1,]

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
# professor does this with montly return, I am doing with daily
acf(spx_log_returns,lag.max=50,main="S&P500")
# we can reject the hypothesis of white noise 






