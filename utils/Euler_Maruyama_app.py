import numpy as np
import pandas as pd

def _Cham_sigma(eta, gamma, pre_value):
    return eta * pre_value**gamma

def _Cham_mu(beta, alpha, pre_value):
    return beta * (alpha - pre_value)

def Cham(T, N, r_0, alpha, beta, gamma, eta):
    
    r_process = [r_0]
    timestamp = [0]
    dt = T/N
    for i in range(1, N+1):
        # random_term = random.normalvariate(0, np.sqrt(dt)) # normal distribution which mean is 0 and variance is dt
        random_term = np.random.randn(1)[0] * np.sqrt(dt)
        pre_r = r_process[i-1]

        r = pre_r + _Cham_sigma(eta, gamma, pre_r) * random_term + _Cham_mu(beta, alpha, pre_r) * dt
        if pd.isna(r) or pd.isna(_Cham_mu(beta, alpha, pre_r)) or pd.isna(_Cham_sigma(eta, gamma, pre_r)):
            print("nan is included in the process")
            print(f"  mu = {_Cham_mu(beta, alpha, pre_r)}")
            print(f"  sigma = {_Cham_sigma(eta, gamma, pre_r)}")
            print(r)
            break
        r_process.append(r)
        timestamp.append(timestamp[i-1] + dt)
    
    r_process = np.array(r_process)
    time_stamp = np.array(timestamp)

    non_noise_sigma = (eta * r_process**gamma)**2
    
    return r_process, time_stamp, non_noise_sigma