import sys
sys.path.append('../../software/models/')
from dftModel import dftAnal, dftSynth
from scipy.signal import get_window
from scipy.fftpack import fft
import numpy as np
from fractions import gcd
import matplotlib.pyplot as plt

def minimizeEnergySpreadDFT(x, fs, f1, f2):
    """
    Inputs:
        x (numpy array) = input signal 
        fs (float) = sampling frequency in Hz
        f1 (float) = frequency of the first sinusoid component in Hz
        f2 (float) = frequency of the second sinusoid component in Hz
    Output:
        The positive half of the DFT spectrum (in dB) of the M sample segment of x. 
                           mX is (M/2)+1 samples long (M is to be computed)
    """
    n1 = fs/f1
    n2 = fs/f2
    M = n1*n2 / gcd(n1,n2)
    
    X = fft(x, M)
    XdB = 20 * np.log10(np.absolute(X))
    mX = XdB[:(M/2)+1]
    return mX

def optimalZeropad(x, fs, f):
    """
    Inputs:
        x (numpy array) = input signal of length M
        fs (float) = sampling frequency in Hz
        f (float) = frequency of the sinusoid in Hz
    Output:
        The positive half of the DFT spectrum of the N point DFT after zero-padding 
        x appropriately (zero-padding length to be computed). mX is (N/2)+1 samples long
    """
    M = len(x)
    m = fs/f
    N = (np.floor(M/m)+1) * m
    print N
    X = np.absolute(fft(x, N))
    low_Xs = X < np.exp(-6)
    X[low_Xs] = 0
    mX = 20.0 * np.log10(X)
    return mX

'''
Example
import numpy as np
import matplotlib.pyplot as plt
x = genSine(1, 250.0, 0, 10000.0, 210.0/10000)
mx = optimalZeropad(x,10000,250)
n= np.arange(len(mx))
plt.plot(n, mx)
plt.show()
'''
"""

def suppressFreqDFTmodel(x, fs, N):
    """
    Inputs:
        x (numpy array) = input signal of length M (odd)
        fs (float) = sampling frequency (Hz)
        N (positive integer) = FFT size
    Outputs:
        The function should return a tuple (y, yfilt)
        y = Output of the dftSynth() without filtering (M samples long)
        yfilt = Output of the dftSynth() with filtering (M samples long)
    The first few lines of the code have been written for you, do not modify it. 
    """
    M = len(x)
    w = get_window('hamming', M)
    outputScaleFactor = sum(w)
    
    ## Your code here
    # compute the dft
    mX, pX = dftAnal(x, w, N)
    # compute the inverse dft
    y = dftSynth(mX, pX, w.size)*sum(w)
    # DFT Filtering
    # mX = np.absolute(20 * log10(mX))
    cutoff_index = np.floor(70*N/fs)
    mX[:cutoff_index] = -120
    yfilt = dftSynth(mX, pX, w.size)*sum(w)
    return (y, yfilt)

'''
Example
import numpy as np
import matplotlib.pyplot as plt
from A3utils import genSine
from A3Part4 import suppressFreqDFTmodel
x = genSine(1, 250.0, 0, 10000.0, 210.0/10000) +  genSine(1, 40.0, 0, 10000.0, 210.0/10000)
y, yfilt = suppressFreqDFTmodel(x, 10000, 256)
n = np.arange(len(yfilt))
plt.plot(n, yfilt)
plt.show()
'''
