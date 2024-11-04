import numpy as np

from src.Functions import gaussian_kernel
 
def cal_FZ_est(process, x, ND, n, T, n_obs_intraday, h):
    
    N = ND * n_obs_intraday
    daily_process = process[0:n:n_obs_intraday]
    FZ_est = np.zeros(len(x))
    vol = np.diff(daily_process)**2 / (T / ND)
    # vol = np.diff(process[0:n])**2 / (T / N)

    for i in range(len(x)):
        # print((obs_processes[0][1:N:n_obs_intraday].shape))
        num = gaussian_kernel((daily_process[0:ND] - x[i]) / h) * vol
        denorm = gaussian_kernel((daily_process[0:ND] - x[i]) / h)
        # num = gaussian_kernel((process[0:n-1] - x[i]) / h) * vol
        # denorm = gaussian_kernel((process[0:n-1] - x[i]) / h)

        FZ_est[i] = np.sum(num) / np.sum(denorm)

    return FZ_est