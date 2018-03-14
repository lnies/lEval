from scipy.stats import poisson
import matplotlib.pyplot as plt
import numpy as np

mu = 0.2643
real_ions_z1 = []

for i in range(100000``):
    measured_event = list(poisson.rvs(0.2643, size=1))[0]
    mcp_efficience = list(poisson.rvs(1/0.3, size=1))[0]

    if measured_event == 1:
        real_ions_z1.append(measured_event*mcp_efficience)


y, x = np.histogram(real_ions_z1, bins=range(-1,max(real_ions_z1)+1))

counts = 0
for i in list(y):
    counts += i

plt.bar([i-1 for i in list(x)[1:]], [float(i)/counts for i in list(y)])
plt.savefig('/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/Z1-1/simulation.pdf')
plt.show()
