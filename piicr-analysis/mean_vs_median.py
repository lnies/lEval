import numpy as np
import matplotlib.pyplot as plt

random = [1,1,1,1,99]
random_error = [5., 5.]
print int(np.mean(random)), int(np.median(random))
while int(np.mean(random)) != int(np.median(random)):
    print 'here'
    random = list(np.random.randint(1000, size=(1, 10000))[0])
    random_error = [i for i in random if i < 600 or i > 950]

plt.figure(1, figsize=(5,5))
plt.hist(random, bins=500)
plt.axvline(x=np.mean(random), linewidth=3, color='r')
plt.axvline(x=np.median(random), linewidth=3, color='y')

plt.savefig('mean_vs_median_random.pdf')

plt.figure(2, figsize=(5,5))
plt.hist(random_error, bins=500)
plt.axvline(x=np.mean(random_error), linewidth=3, color='r')
plt.axvline(x=np.median(random_error), linewidth=3, color='y')

plt.savefig('mean_vs_median_random_cut.pdf')

print np.mean(random)
print np.median(random)
print np.mean(random_error)
print np.median(random_error)
