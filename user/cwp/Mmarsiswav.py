#!/usr/bin/env python
'''
Create the time domain range compressed waveform for MARSIS
Created by: Dylan Hickson, Colorado School of Mines
Created on: Feb 16, 2022
Edited May 4, 2022:
    - add option to zero-pad
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

verb    = par.bool('verb',False) # verbosity flag
window  = par.bool('window')     # spectral Hann window flag
delay   = par.float('delay')     # latency delay in transmitted signal in microseconds
cmplxBB = par.bool('cmplxBB')    # complex baseband chirp
zpad    = par.bool('zpad')       # zero pad to closest 2^n

def hann(freq, BW, df, BW_l, BW_h):
    freq = np.fft.fftshift(freq)
    N = int(BW/df)
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
    return np.fft.ifftshift(out)

def marsis_rc_wav_cmplx(delay, window=False, zpad=True):
    pulse_w = 250e-6             # pulse_width (s)
    BW = 1e6                     # chirp bandwidth (Hz)

    # 490 complex samples recorded in receive window
    Ns = 490
    fs = 1.4e6
    dt = 1/fs
    t = np.arange(0,Ns,1)*dt

    # Create complex MARSIS chirp
    cmplx_chirp_t = np.exp(np.pi * 1j * BW / pulse_w * (t - pulse_w/2)**2)
    cmplx_chirp_t[t > pulse_w] = 0 + 1j*0

    if zpad:
        # Pad to 512
        Np = 512
        cmplx_chirp_t = np.append(cmplx_chirp_t, np.zeros(Np-Ns, dtype=complex))
        Ns = Np

    # FFT to frequency
    freq = np.fft.fftfreq(Ns,d=dt)
    df = freq[1]-freq[0]
    cmplx_chirp_f = np.fft.fft(cmplx_chirp_t)

    # Window 1 MHz BW around DC
    if window:
        window_f = hann(freq, BW, df, -0.5e6, 0.5e6)
        cmplx_chirp_f_windowed = cmplx_chirp_f*window_f

    # Range compression
    if window:
        h_cmplx = np.conj(cmplx_chirp_f_windowed)
        rc_chirp_f = cmplx_chirp_f_windowed*h_cmplx
    else:
        h_cmplx = np.conj(cmplx_chirp_f)
        rc_chirp_f = cmplx_chirp_f*h_cmplx

    # Shift to latency delay and normalize
    rc_chirp_f *= np.exp(-1j*2*np.pi*freq*delay*1e-6)
    rc_chirp_t = np.fft.ifft(rc_chirp_f).real
    rc_chirp_t /= np.max(rc_chirp_t)

    return rc_chirp_t, rc_chirp_f

def marsis_rc_wav_real(delay, window=False, zpad=True):
    pulse_w = 250e-6             # pulse_width (s)
    BW = 1e6                   # chirp bandwidth (Hz)

    # 980 real samples recorded in receive window
    Ns = 980
    fs = 2.8e6
    dt = 1/fs
    t = np.arange(0,Ns,1)*dt
    # down-converted carrier frequency is 0.7 MHz
    f_carrier = 0.7e6
    fl = f_carrier - (BW / 2)          # lowest frequency of chirp (Hz)
    fh = f_carrier + (BW / 2)          # highest frequency of chirp (Hz)

    # Create real MARSIS chirp
    real_chirp_t = np.cos(2*np.pi*f_carrier*t + np.pi * BW / pulse_w * (t - pulse_w/2)**2)
    real_chirp_t[t > pulse_w] = 0

    if zpad:
        # Pad to 1024
        Np = 1024
        real_chirp_t = np.append(real_chirp_t, np.zeros(Np-Ns))
        Ns = Np

    # FFT to frequency
    freq = np.fft.fftfreq(Ns,d=dt)
    df = freq[1]-freq[0]
    real_chirp_f = np.fft.fft(real_chirp_t)

    # Window
    if window:
        pos_window = hann(freq, BW, df, fl, fh)
        neg_window = hann(freq, BW, df, -fh, -fl)
        window_f = pos_window + neg_window
        real_chirp_f_windowed = real_chirp_f*window_f

    # Range compression
    if window:
        h_real = np.conj(real_chirp_f_windowed)
        rc_chirp_f = real_chirp_f_windowed*h_real
    else:
        h_real = np.conj(real_chirp_f)
        rc_chirp_f = real_chirp_f*h_real

    # Shift to latency delay and normalize
    rc_chirp_f *= np.exp(-1j*2*np.pi*freq*delay*1e-6)
    rc_chirp_t = np.fft.ifft(rc_chirp_f).real
    rc_chirp_t /= np.max(rc_chirp_t)

    return rc_chirp_t, rc_chirp_f

if cmplxBB:
    marsisOut_t, marsisOut_f = marsis_rc_wav_cmplx(delay, window=window, zpad=zpad)
    dt = 1/1.4e6
    if zpad:
        Ns = 512
    else:
        Ns = 490
else:
    marsisOut_t, marsisOut_f = marsis_rc_wav_real(delay, window=window, zpad=zpad)
    dt = 1/2.8e6
    if zpad:
        Ns = 1024
    else:
        Ns = 980

# Write to RDF
wavTout = rsf.Output()
wavTout.put('n1',Ns)
wavTout.put('d1',dt*1e6)
wavTout.put('l1','us')
wavTout.put('o1',0)

wavTout.write(marsisOut_t)

wavTout.close()
