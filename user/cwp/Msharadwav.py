#!/usr/bin/env python
'''
Create the time domain range compressed waveform for SHARAD
Created by: Dylan Hickson, Colorado School of Mines
Created on: Mar 9, 2022
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag
window = par.bool('window')  # spectral Hann window flag
delay  = par.float('delay')  # latency delay in transmitted signal in microseconds

# Local functions
def hann(freq, BW, df, BW_l, BW_h):
    N = BW/df
    n = np.arange(0,N,1)
    w = np.sin(np.pi*n/N)**2

    out = np.zeros(len(freq))
    cnt=0
    for i in range(len(freq)):
        if freq[i] >= BW_l and freq[i] <= BW_h:
            if cnt==len(w):
                out[i] = w[cnt-1]
            else:
                out[i] = w[cnt]
            cnt+=1
    return out

# SHARAD chirp parameters
pulse_w = 85e-6             # pulse_width (s)
BW = 10e6                   # chirp bandwidth (Hz)
fc = 20e6                   # chirp centre frequency (Hz)
f0 = fc + (BW / 2)          # beginning frequency of chirp (Hz), use (+) for downswept, (-) for upswept
f1 = fc - (BW / 2)          # ending frequency of chirp (Hz), use (-) for downswept, (+) for upswept

dt = 0.0375e-6                # time sampling (non-baseband)
N = 2268                        # number of samples in transmitted chirp
n = np.arange(0,N,1)

# Create real SHARAD chirp
real_chirp_t = np.sin(2*np.pi*(n*dt)*(f0 - 0.5*(BW/pulse_w)*(n*dt)))

# Pad to 3600
N = 3600
n = np.arange(0,N,1)
time = n*dt
real_chirp_t = np.append(real_chirp_t, np.zeros(3600-2268))

# FFT to frequency
freq = np.fft.fftfreq(N,d=dt)
df = freq[1]-freq[0]
real_chirp_f = np.fft.fft(real_chirp_t)

# Window
if window:
    pos_window = hann(freq, BW, df, 1.67e6, 11.67e6)
    neg_window = hann(freq, BW, df, -11.67e6, -1.67e6)
    window_f = pos_window + neg_window
    real_chirp_f_windowed = real_chirp_f*window_f

# Range compression
if window:
    h_real = np.conj(real_chirp_f_windowed)
    rc_chirp_f = real_chirp_f_windowed*h_real
else:
    h_real = np.conj(real_chirp_f)
    rc_chirp_f = real_chirp_f*h_real

# Shift to latency delay
rc_chirp_f *= np.exp(-1j*2*np.pi*freq*(delay*1e-6))
rc_chirp_t = np.fft.ifft(rc_chirp_f).real
rc_chirp_t /= np.max(rc_chirp_t)

# ------------------------------------------------------------
Fou = rsf.Output()            # output file
Fou.put("n1",N)
Fou.put("o1",0)
Fou.put('d1',dt*1e6)

Fou.write(rc_chirp_t)

Fou.close()
