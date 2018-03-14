import pandas as pd
import matplotlib.pyplot as plt


temp = pd.read_csv('/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/Info/temperature/temperature-26-06--03-07.csv', delimiter=';', parse_dates = [[0,1]])
temp['moving_avg'] = pd.rolling_mean(temp['UT_lower'], window=5000)
temp['moving_avg'] = temp.moving_avg.shift(-2500)
print temp
plt.figure()
temp.plot(x='Date_Time', y=['UT_lower', 'moving_avg'])
plt.savefig('temp_vs_moving_avg.pdf')
