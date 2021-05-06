import random
import numpy as np
import fourier
import time
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline

data_size = [4,5,6,7,8,9,10,11,12]
fftTime = []
dftTime = []
error = []



for i in data_size:
    signal = []
    signal = np.random.uniform(low=0.1, high=1, size=(2**(i),))
    dft = []
    start = time.time()
    dft = fourier.DFT(signal)
    end = time.time()
    dftTime.append((end-start) * 10**6)
    fft = []
    start = time.time()
    fft = fourier.FFT(signal)
    end = time.time()
    fftTime.append((end-start) * 10**6)
    err = []
    for i in range(len(fft)):
        err.append((abs(dft[i] - fft[i]))**2)
    error.append(np.average(err))    

# 300 represents number of points to make between T.min and T.max
xnew = np.linspace(min(data_size), max(data_size), 300)  
xnew2 = np.linspace(min(error), max(error), 300)

spl = make_interp_spline(data_size, dftTime, k = 5)
power_smoothD = spl(xnew)
spl = make_interp_spline(data_size, fftTime, k = 3)
power_smoothT = spl(xnew)
spl = make_interp_spline(data_size, error, k = 3)
p = spl(xnew2)
# plt.plot(xnew,power_smooth)
# plt.show()
fig,(ax1,ax2) = plt.subplots(2)
x1 = ax1.plot(xnew, power_smoothT)
x2 = ax1.plot(xnew, power_smoothD)
ax1.set(ylabel = 'Time Complexity',xlabel = 'List Length')
ax1.legend((x2) ,"Tf", loc='upper left', shadow=True)# ax1.se.xlabel('List Length')
ax2.plot(xnew2,p,label = "Mean square error")
plt.grid()
plt.show()