import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy import genfromtxt

# my_data = genfromtxt('frames/frame_mat100.csv', delimiter=' ')
# plt.imshow(my_data, cmap='hot', interpolation='nearest')

tb = genfromtxt('center.csv', delimiter=',')
print(tb)

plt.plot(tb)
plt.show()