'''
Created on Mar 31, 2015

@author: mlai
'''
from __future__ import division  #force floating point division
import numpy as np
from scipy import signal

def createCubedGaussianWhiteNoiseConvolvedWithRickerWavelet(
                                                num_signal,
                                                num_sample,
                                                num_point_ricker_wavelet,
                                                width_parameter_ricker_wavelet):
    gaussian_white_noise = np.random.randn(num_signal,num_sample)
    cubed_noise_by_signal_sample = np.power(gaussian_white_noise, 3)
    cubed_noise_by_signal_sample[:, num_sample * 0 // 10: num_sample*1//10] = 0
    #cubed_noise_by_signal_sample *= np.random.poisson(100./num_sample, (num_signal, num_sample))
    #cubed_noise_by_signal_sample = np.random.poisson(1, (num_signal, num_sample))
    ricker_wavelet = signal.ricker(num_point_ricker_wavelet, 
                         width_parameter_ricker_wavelet)
    output_sample_length = num_sample 
    output_by_signal_sample = np.zeros((num_signal,output_sample_length))
    for idx_signal in range(num_signal):
        output_by_signal_sample[idx_signal,:] = np.convolve(
                                        cubed_noise_by_signal_sample[idx_signal,:], 
                                        ricker_wavelet,
                                        'same') 
        #output_by_signal_sample[idx_signal,17] = 1e-15
        
    return output_by_signal_sample

if __name__ == '__main__':
    pass