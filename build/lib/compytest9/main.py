"""
Discomfort Severity of Different Comfort Definitions
Python version used: 3.7
Created by: Shide Salimi, Harvard University, February 2021
ssalimi@gsd.harvard.edu,
shide.salimi@gmail.com,
Version 0.0.0
"""
import glob, os
import csv
import pandas as pd
from compytest9.discomfortmethod import *
from compytest9.thermaldefinitions import *
import numpy as np
from datetime import *
import copy
import warnings

# Find the working directory
dir_path = os.path.abspath('')
# Change directory to the location of input file
os.chdir(dir_path)


for file in glob.glob("*.csv"):
    split_name = file.split("_")
    thermal_model = split_name[1].split(".")[0]

               
# Input dataframe 
df = pd.read_csv(file, delimiter = ',')

if thermal_model == 'PMV':
    # PMV Model  
    
    # Changing the headings of the df 
    old_headings = df.columns
    new_headings = ['#', 'Date', 'Time', 'SET', 'OpTem', 'dbTem', 'mrTem', 'v', 'rh', 'met', 'clo']
    df.rename(columns=dict(zip(old_headings, new_headings)), inplace=True)
    
    # Changing the format of Date and Time columns
    df['Date'] = df['Date'].apply(lambda row: datetime.strptime(row, '%d-%b').strftime('%m/%d'))
    df['Time'] = df['Time'].apply(lambda row: datetime.strptime(row, '%H:%M:%S').strftime('%H:%M:%S'))
    df['Time'] = df['Time'].apply(lambda row: row.replace("00:00:00", "24:00:00"))
      
    discomfort_list_pmv = []
    for i in range(df.shape[0]):        
        comfortzone = comfort_zone(df, i)                     
        discomfort = pmv_model(df, i, comfortzone)
        discomfort_list_pmv.append(discomfort)
        
    # Save results as a csv file
    df_discomfort_pmv = pd.DataFrame(discomfort_list_pmv, columns = ['Discomfort Severity']) 
    df_discomfort_pmv.to_csv(dir_path + '\\discomfort_severities_' + thermal_model + '.csv', index=False)
else:
    # Adaptive Model
    
    # Changing the headings of the df 
    old_headings = df.columns
    new_headings = ['#', 'Date', 'Time', 'SET', 'OpTem', 'dbTem', 'mrTem', 'v', 'rh', 'met', 'clo', 'PrevTem']
    df.rename(columns=dict(zip(old_headings, new_headings)), inplace=True)
    
    # Changing the format of Date and Time columns
    df['Date'] = df['Date'].apply(lambda row: datetime.strptime(row, '%d-%b').strftime('%m/%d'))
    df['Time'] = df['Time'].apply(lambda row: datetime.strptime(row, '%H:%M:%S').strftime('%H:%M:%S'))
    df['Time'] = df['Time'].apply(lambda row: row.replace("00:00:00", "24:00:00"))
    
    discomfort_list_adaptive = []
    for i in range(df.shape[0]): 
        t_rm = df['PrevTem'][i]
        discomfort = adaptive_model(df, i, t_rm) 
        discomfort_list_adaptive.append(discomfort)
        
    # Save results as a csv file
    df_discomfort_adaptive = pd.DataFrame(discomfort_list_adaptive, columns = ['Discomfort Severity']) 
    df_discomfort_adaptive.to_csv(dir_path + '\\discomfort_severities_' + thermal_model + '.csv', index=False)

