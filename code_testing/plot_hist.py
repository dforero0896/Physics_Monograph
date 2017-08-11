import numpy as np
import matplotlib.pyplot as plt

#filename = str(input("Please enter file name: "))
#filename += ".csv"
filename="prob_test.csv"
'''
file_arr = np.loadtxt(filename, dtype=float, skiprows=1)
plt.hist(file_arr)
plt.gcf()
plt.savefig("hist_test.png", bins=1000)
'''

bresenham_test = np.loadtxt('bresenham_test.csv', delimiter=',', dtype=float)
plt.imshow(bresenham_test)
plt.gcf()
plt.savefig("bresenham_test.png", interpolation=False)
