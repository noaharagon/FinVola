################################
#Jonas Schmitten, Moritz Köhler,
#Giovanni Magagnin, Noah Angara
#May 2021
###############################

#Import libraries 
import os 
import pandas as pd
import yfinance as yf

#--------------------------------------------------------------------------------------------------------
#Set working directory
if os.getcwd().split('/')[2] == 'jonasschmitten':
    os.chdir('/Users/jonasschmitten/Downloads/NOPE and Volatility/Data')  
elif os.getcwd().split('/')[2] == 'YOURNAME1':
    os.chdir('YOURPATH1')
elif os.getcwd().split('/')[2] == 'YOURNAME2':
    os.chdir('YOURPATH2')
elif os.getcwd().split('/')[2] == 'YOURNAME3':
    os.chdir('YOURPATH3')
#--------------------------------------------------------------------------------------------------------

#Read in data and check variables
data = pd.read_csv('options_data.csv', nrows=1000000)

#not entirely sure what 'issue_type' means, had only A. 
data.drop(['secid', 'index_flag', 'issue_type', 'issuer', 'exercise_style'], axis=1 , inplace=True)

#Convert date column to date-type
data[['date', 'exdate']] = data[['date', 'exdate']].apply(pd.to_datetime, utc= False, format='%Y%m%d')

#strike price, best bid, best offer seem to be in cents?
data[['strike_price', 'best_bid', 'best_offer']]  = data[['strike_price', 'best_bid', 'best_offer']]/ 1000

#Get SPX data from Yahoo Finance
SPX = yf.download("^GSPC", start="2019-01-01", end="2020-12-31")

