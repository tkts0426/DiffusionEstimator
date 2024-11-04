import numpy as np
from numpy.fft import fft

from src.Functions import gaussian_kernel

def cal_FE_day_est(process, x, timestamp, ND, T, n, h, n_obs_intraday):
    
    FE_day_est = np.zeros(len(x))
    vol = np.zeros(ND)
    L = int(np.floor(n_obs_intraday / 2))
    M = 10

    for k in range(ND):
        vol[k] = cal_FE_fft_spot_vol(
            process[k*n_obs_intraday : (k+1) * n_obs_intraday], 
            timestamp[k * n_obs_intraday], 
            T / ND, 
            L,  
            M
        )
    daily_process = process[0 : n : n_obs_intraday]
    for i in range(len(x)):
        num = gaussian_kernel((daily_process[0:ND] - x[i]) / h) * vol
        denorm = gaussian_kernel((daily_process[0:ND] - x[i]) / h)

        FE_day_est[i] = np.sum(num) / np.sum(denorm)
    
    return FE_day_est

def cal_FE_whole_est(process, x, timestamp, ND, T, n, h, n_obs_intraday, N):

    FE_whole_est = np.zeros(len(x))
    L = int(np.floor(N / 2))
    M = 100
    vol = cal_FE_fft_spot_vol(process, timestamp[0 : n-1 : n_obs_intraday], T, L, M)
    # vol = cal_FE_fft_spot_vol(process, timestamp[0 : n-1], T, L, M)

    daily_process = process[0 : n : n_obs_intraday]
    for i in range(len(x)):
        num = gaussian_kernel((daily_process[0:ND] - x[i]) / h) * vol
        denorm = gaussian_kernel((daily_process[0:ND] - x[i]) / h)
        # num = gaussian_kernel((process[0:n-1] - x[i]) / h) * vol
        # denorm = gaussian_kernel((process[0:n-1] - x[i]) / h)

        FE_whole_est[i] = np.sum(num) / np.sum(denorm)

    return FE_whole_est
    
def cal_FE_fft_spot_vol(process, tau, T, L, M):
    
    nv = np.max(tau.size)
    const = 2*np.pi / T
    ret = np.diff(process)
    spot = np.zeros(nv, dtype=complex)
    c_0 = np.sum(ret)
    c_s = np.zeros(2 * M + 1, dtype=complex)

    fft_front = fft(ret)[1 : L + M + 1]
    fft_back = np.flipud(np.conj(fft_front))
    fourier_coef = np.concatenate([fft_back, np.array([c_0]), fft_front])

    c_p = fourier_coef / T
    fact = T / (2*L + 1)
    n_shift = L + M
    
    for k in np.arange(-M, M+1):
        for l in np.arange(-L, L+1):
            c_s[k+M] += fact * (c_p[l + n_shift] * c_p[k - l + n_shift])

    if nv == 1:
        tau = np.array([tau])
        
    for t in range(nv):
        spot[t] = 0
        for k in np.arange(-M, M+1):
            spot[t] = spot[t] + (1 - (abs(k) / M)) * c_s[k+M] * np.exp( 1j * tau[t] * const * k)

    return spot.real

def cal_FE_direct_spot_vol(process, timestamp, tau, delta_t, L, M):
    
    nv = np.max(tau.size)
    const = 2 * np.pi / delta_t
    c_pp = np.zeros(L + M)
    c_p = np.zeros(2 * L + (2 * M + 1))
    ret = np.diff(process)
    spot = np.zeros(nv)
    c_0 = np.sum(ret)
    c_s = np.zeros(2 * M + 1)

    for k in range(L + M):
        c_pp[k] = np.sum(np.exp( -1j * const * (k+1) * timestamp[0:nv]) * ret)

    for j in range(L + M):
        c_p[j] = np.conj(c_pp[L + M - j - 1]) / delta_t
    
    c_p[L + M] = c_0 / delta_t
    
    for j in range(L + M):
        c_p[L + M + 1 + j] = c_pp[j] / delta_t

    fact = delta_t / (2 * L + 1)
    n_shift = L + M

    for k in np.arange(-M, M+1):
        c_s[k + M] = 0
        for l in np.arange(-L, L+1):
            c_s[k + M] += fact * (c_p[l + n_shift - 1] * c_p[k-l + n_shift - 1])

    for t in range(nv):
        spot[t] = 0
        for k in np.arange(-M, M+1):
            spot[t] += (1 - (abs(k) / M)) * c_s[k+M] * np.exp(1j * timestamp[t] * const * k)

    return spot.real