import numpy as np

from utils.Functions import gaussian_kernel

def cal_RV_est(process, x, ND, n, T, n_obs_intraday, h):
    RV_est = np.zeros(len(x))
    vol = np.zeros(ND)
    
    for k in range(ND):
        # vol[k] = np.sum(np.diff(process[k:k+100])**2) / n_obs_intraday
        vol[k] = np.sum(np.diff(process[k*n_obs_intraday + 1 : (k+1) * n_obs_intraday + 1])**2)
    
    for i in range(len(x)):
        num = gaussian_kernel((process[0:n-1:n_obs_intraday] - x[i]) / h) * vol
        denorm = gaussian_kernel((process[0:n-1:n_obs_intraday] - x[i]) / h)

        RV_est[i] = (np.sum(num)*ND / (T*np.sum(denorm)))
    
    # print(RV_est)
    return RV_est