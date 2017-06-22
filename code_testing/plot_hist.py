import numpy as np
import matplotlib.pyplot as plt

#filename = str(input("Please enter file name: "))
#filename += ".csv"
filename="prob_test.csv"
print filename
file_arr = np.loadtxt(filename, dtype=float, skiprows=1)
plt.hist(file_arr)
plt.gcf()
plt.savefig("hist_test.png", bins=1000)
