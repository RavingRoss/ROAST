import pandas as pd
import numpy as np
# Setting paths to all txt files
ks = "measurement files/Measurement/All-KS.txt"
first_three = "measurement files/Measurement/12_14-12_16-KS.txt"
dir0 = "measurement files/TimeSeries/TimeSeries-2459198.8027535-121420UT.txt"
dir1 = "measurement files/TimeSeries/TimeSeries-2459199.6807775-121520UT.txt"
dir2 = "measurement files/TimeSeries/TimeSeries-2459200.7093692-121620UT.txt"
dir3 = "measurement files/TimeSeries/TimeSeries-2459223.7600876-010821UT.txt"
dir4 = "measurement files/TimeSeries/TimeSeries-2459546.8129136-112721UT.txt"
dir5 = "measurement files/TimeSeries/TimeSeries-2459582.736962-010222UT.txt"
dir6 = "measurement files/TimeSeries/TimeSeries-2459583.9094831-010322UT.txt"
dir7 = "measurement files/TimeSeries/TimeSeries-2459584.718099-010422UT.txt"
# Taking the KS data from the All-BJD.txt file
ks = np.loadtxt(ks,skiprows=1)
# Taking the first 3 nights
first_three = np.loadtxt(first_three,skiprows=1)
# Making directories for these files
p0 = np.loadtxt(dir0,skiprows=1)
p1 = np.loadtxt(dir1,skiprows=1)
p2 = np.loadtxt(dir2,skiprows=1)
p3 = np.loadtxt(dir3,skiprows=1)
p4 = np.loadtxt(dir4,skiprows=1)
p5 = np.loadtxt(dir5,skiprows=1)
p6 = np.loadtxt(dir6,skiprows=1)
p7 = np.loadtxt(dir7,skiprows=1)
# Converting into a DataFrame
ks = pd.DataFrame(ks)
first_three = pd.DataFrame(first_three)
df0 = pd.DataFrame(p0)
df1 = pd.DataFrame(p1)
df2 = pd.DataFrame(p2)
df3 = pd.DataFrame(p3)
df4 = pd.DataFrame(p4)
df5 = pd.DataFrame(p5)
df6 = pd.DataFrame(p6)
df7 = pd.DataFrame(p7)

# Concatenating all files
'''frames = [df0,df1,df2,df3,df4,df5,df6,df7]
concat = pd.concat(frames)
concat.columns = ['Kilo_Secs','Relative_Flux','Relative_Flux_Error']
ks.columns = ['Kilo_Secs']
concat = concat.replace([concat['Kilo_Secs']],[ks['Kilo_Secs']])
concat.to_csv('measurement files/concat-TimeSeries-new.csv')'''
# Making Time Series for first three nights
frames = [df0,df1,df2]
concat_first_three = pd.concat(frames)
concat_first_three.columns = ['Kilo_Secs','Relative_Flux','Relative_Flux_Error']
first_three.columns = ['Kilo_Secs']
concat_first_three = concat_first_three.replace([concat_first_three['Kilo_Secs']],[first_three['Kilo_Secs']])
concat_first_three.to_csv('measurement files/concat 12_14 to 12_16-TimeSeries.csv')