import numpy as np
from numpy.fft import fft
from tqdm import tqdm

from utils.Functions import gaussian_kernel

"""
Calculation fourier estimator with non spot volatility by using fft
Input:
    process : 
    x : 
    n : 
    h : 
    L : 
    T : 
"""
def fft_cal_fourier_est(process, x, n, h, L, T):
    
    LR = np.zeros(len(x))
    fourier_coef = np.zeros(len(x))
    fourier_estimator = np.zeros(len(x))

    ret = np.diff(process) # return 
    for i in range(len(x)):
        
        R = np.sqrt( (1 / h) * gaussian_kernel((process[0:n-1] - x[i]) / h) ) * ret
        
        fft_v = fft(R)
        # np.savetxt("fft_v.csv", fft_v, delimiter=",")
        fft_coef_0 = fft_v[0]
        fft_front = fft_v[1:L+1]
        fft_back = np.flipud(np.conj(fft_front))
        fft_def = np.concatenate([fft_back, fft_v[0:L+1]])
        
        fft_coef = np.sum(fft_def * np.flipud(fft_def))
        fourier_coef[i] = fft_coef * (1 / (2*L + 1))

        LR[i] += (1/h) * np.sum(gaussian_kernel((process[0:n] - x[i]) / h)) * (T / n)

        fourier_estimator[i] = fourier_coef[i] / LR[i]
    
    return fourier_estimator, LR

# direct culculation
def cal_fourier_est(obs_process, timestamp, x, T, n, L, h):
    
    LR = np.zeros(len(x))
    fourier_coef = np.zeros(len(x)) # fourier coefficient
    fourier_coef_sub_1 = np.zeros((L, len(x)))
    fourier_coef_sub_2 = np.zeros((L, len(x)))
    fourier_estimator = np.zeros(len(x)) # our estimator
    
    ret = np.diff(obs_process) # return 

    for i in tqdm(range(len(x))):
        s1 = np.sum(np.sqrt((1/h) * gaussian_kernel((obs_process[1:n] - x[i]) / h)) * ret)
        s2 = np.sum(np.sqrt((1/h) * gaussian_kernel((obs_process[1:n] - x[i]) / h)) * ret)
        fourier_coef[i] += s1 * s2 # l = 0

        LR[i] += (1/h) * np.sum(gaussian_kernel((obs_process[1:n] - x[i]) / h)) * (T / n)
        
        # 非対角成分 つまり，l neq 0
        for l in range(1, L):
            exp_variable = -2j * np.pi * l * (timestamp[1:n]) / n
            fourier_coef_sub_1 = np.sum(np.exp(exp_variable) * np.sqrt((1/h) * gaussian_kernel((obs_process[1:n] - x[i]) / h)) * ret)
            fourier_coef_sub_2 = np.conj(fourier_coef_sub_1)

            fourier_coef[i] += 2 * fourier_coef_sub_1 * fourier_coef_sub_2 # 対称性を用いているから，２倍
        # print(fourier_coef[i])
        
        fourier_estimator[i] += (1 / (2*L + 1)) * fourier_coef[i] / LR[i]
    
    return fourier_estimator, LR