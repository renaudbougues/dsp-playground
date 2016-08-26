import numpy as np

def genSine(A, f, phi, fs, t):
    """
    Inputs:
        A (float) =  amplitude of the sinusoid
        f (float) = frequency of the sinusoid in Hz
        phi (float) = initial phase of the sinusoid in radians
        fs (float) = sampling frequency of the sinusoid in Hz
        t (float) =  duration of the sinusoid (is second)
    Output:
        The generated sinusoid (use np.cos())
    """
    ts = np.arange(0, t, 1.0/fs) #timescale
    x = A * np.cos(2 * np.pi * f * ts + phi)
    return x

def complexSine(k, N):
    """
    Inputs:
        k (integer) = frequency index of the complex sinusoid of the DFT (k < N-1)
        N (integer) = length of complex sinusoid in samples
    Output:
        The generated complex sinusoid (length N)
    """
    n = np.array(np.arange(N))
    cSine = np.exp(-1j * 2 * np.pi * k * n / N)
    return cSine
    
def DFT(x):
    """
    Input:
        x (numpy array) = input sequence of length N
    Output:
        The N point DFT of the input sequence x
    """
    X = np.array([])
    N = len(x)
    for k in np.arange(N):	
      X = np.append(X, np.sum(x * complexSine(k, N)))
    return X

def IDFT(X):
    """
    Input:
        X (numpy array) = frequency spectrum (length N)
    Output:
        The N point IDFT of the frequency spectrum X
    """
    x = np.array([])
    N = len(X)
    k = np.arange(N)
    
    for n in np.arange(N):
      # Complex positive exponential
      icSine = np.exp(1j * 2 * np.pi * k * n / N)
      x = np.append(x, np.sum(X * icSine) / N)
    return x
    
def magSpect(x):
    """
    Input:
        x (numpy array) = input sequence of length N
    Output:
        The magnitude spectrum of the input sequence x (length N)
    """
    X = DFT(x)
    magX = np.absolute(X)
    return magX
