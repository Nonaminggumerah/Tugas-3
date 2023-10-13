# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 18:34:38 2023

@author: ADMIN
"""

print ("FFT 1D Sederhana")
print ("Nadya Putri Arisni")
print ("NRP = 5009211079")

#panggil library cmath dan scipy,fft serta matplotlib
#library cmath untuk melakukan operasi e^
#library scipy.fft untuk melakukan fft


from cmath import exp, pi
from scipy.fft import fft
import matplotlib.pyplot as plt

def oprasiFFT(x):
    N = len(x)
    ans = []
    if N <= 1:
        return x
    else:
        for k in range(N):
            temp = []
            for i in range (N):
                temp.append(x[i] * exp(-2j * pi * k * i / N))
            ans.append(sum(temp))
        return ans
    
# Fungsi fft di atas tanpa library
 
x = [2, 3, 4, 5, 6]
    
myans = oprasiFFT(x)
ans = fft(x)

print(myans)
print(ans)



