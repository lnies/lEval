import os

for path, subdirs, files in os.walk('/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/Z1-2'):
    for filename in files:
        if 'positions.csv' in filename:
            print path+'/'+filename
