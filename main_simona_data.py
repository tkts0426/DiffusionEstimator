# Library
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from numpy.fft import fft # fast fourier transformation
import warnings
warnings.simplefilter('ignore')

# My code:
from utils.Calculation_Integral import cal_analytical_int, cal_denorm_int
from utils.FZ_Estimator import cal_FZ_est
from utils.RV_Estimator import cal_RV_est
from utils.FE_spot import cal_FE_day_est, cal_FE_whole_est
from utils.FE_nonSpot import fft_cal_fourier_est
from utils.Plot_Functions import hist_plot, cum_dist_plot, approximation_plot

# fix random seed 
# np.random.seed(326)

def main():

    print("********************************************************************************")
    print("*             Diffusion estimation on N estimation points over ND days         *")
    print("*                by means of Realized Volatility type methods                  *")
    print("*                           and Fourier estimator                              *")
    print("*                                                                              *")
    print("*                                                                              *")

    n_traj = 1 # the number of paths
    ND = 250 # the number of days
    T = ND  # the total period
    n_obs_intraday = 100 # the number of trading times per one day
    N = n_obs_intraday * ND # the number of trading times
    D = 1 # serves for possible scattered sampling
    # x = np.arange(0.04, 0.1, 0.0001) # candidates of inital value of price process
    x = np.arange(0.04, 0.51, 0.01) # for simona data

    print("*... calculation                                                               *")
    print("*                                                                              *")
    """
    Estimation of the diffusion coefficient using all intraday prices
    """
    L = int((N/D)/2) # D = 1 => L/n = 1/2
    # L = int(np.round(math.pow(N, 1/3)))
    # Parameter of Chan: ---------------------------------------------
    alpha = 0.079 
    beta=0.093
    gamma=1.474
    eta=0.794
    r_0=0.065
    # -----------------------------------------------------------------
    
    # Generating process
    r_processes = np.zeros((n_traj, N))
    timestamps = np.zeros((n_traj, N))
    obs_processes = np.zeros((n_traj, N))
    non_noise_sigmas = np.zeros((n_traj, N))
    fourier_estimators = np.zeros((n_traj, len(x)))
    
    analytical_ints = np.zeros((n_traj, len(x)))
    denorm_ints = np.zeros((n_traj, len(x)))
    var_of_var = np.zeros((n_traj, len(x)))
    
    FZ_estimators = np.zeros((n_traj, len(x)))
    RV_estimators = np.zeros((n_traj, len(x)))
    FE_day_estimators = np.zeros((n_traj, len(x)))
    FE_whole_estimators = np.zeros((n_traj, len(x)))
    
    # import simona data
    obs_processes = np.genfromtxt("data/data1.csv").reshape(1, N+1)
    timestamp = np.arange(0, N+1, 1/n_obs_intraday)
    total_obs_prices = np.genfromtxt("data/data1.csv")

    n = len(obs_processes[0])

    # plot histgram of total obserbed process
    sorted_total_obs_prices = np.sort(total_obs_prices)
    hist_plot(total_obs_prices)
    cum_dist_plot(sorted_total_obs_prices, "obserbed prices", label="simona_")
    
    """
    Calculate the fourier estimator
    """
    h_1 = 3
    S = np.std(obs_processes, axis=1)
    h_2 = h_1 * S * N**(-1/5) # h = n^{-1/5}

    for k in range(n_traj):
        h = h_2[k]
        fourier_estimator, denorm_int = fft_cal_fourier_est(obs_processes[k], x, n, h, L, T) # fft
        analytical_int = cal_analytical_int(obs_processes[k], x, eta, gamma, n, h, T)
        FZ_estimator = cal_FZ_est(obs_processes[k], x, ND, n,  T, n_obs_intraday, h)
        RV_estimator = cal_RV_est(obs_processes[k], x, ND, n , T, n_obs_intraday, h)
        FE_day_estimator = cal_FE_day_est(obs_processes[k], x, timestamp, ND, T, n, h_2[k], n_obs_intraday)
        FE_whole_estimator = cal_FE_whole_est(obs_processes[k], x, timestamp, ND, T, n, h_2[k], n_obs_intraday, N)
        
        fourier_estimators[k] = fourier_estimator
        denorm_ints[k] = denorm_int
        analytical_ints[k] = analytical_int
        var_of_var[k] = fourier_estimator**2

        FZ_estimators[k] = FZ_estimator
        RV_estimators[k] = RV_estimator
        FE_day_estimators[k] = FE_day_estimator
        FE_whole_estimators[k] = FE_whole_estimator
        
        
    if n_traj == 1:
        fourier_est = fourier_estimators[0]
        var = var_of_var[0]

        FZ_est = FZ_estimators[0]
        RV_est = RV_estimators[0]
        FE_day_est = FE_day_estimators[0]
        FE_whole_est = FE_whole_estimators[0]
    else:
        fourier_est = np.mean(fourier_estimators, axis=0)
        var = (np.sum(var_of_var, axis=0) - n_traj*(fourier_est**2)) / (n_traj - 1) # Specimen Dispersion

        FZ_est = np.mean(FZ_estimators, axis=0)
        RV_est = np.mean(RV_estimators, axis=0)
        FE_day_est = np.mean(FE_day_estimators, axis=0)
        FE_whole_est = np.mean(FE_whole_estimators, axis=0)  

    fourier_mse_error = np.mean(np.abs(fourier_est - (eta*x**gamma)**2)**2)
    FZ_mse_error = np.mean(np.abs(FZ_est - (eta*x**gamma)**2)**2)
    RV_mse_error = np.mean(np.abs(RV_est - (eta*x**gamma)**2)**2)
    FE_day_mse_error = np.mean(np.abs(FE_day_est - (eta*x**gamma)**2)**2)
    FE_whole_mse_error = np.mean(np.abs(FE_whole_est - (eta*x**gamma)**2)**2)
    
    print("*                                                                              *")
    print("* The Error:                                                                   *")
    print("*")
    print(f"* - The Mean Squared Error of Fourier estimator : {fourier_mse_error}")
    print("*")
    print(f"* - The Mean Squared Error of FZ estimator : {FZ_mse_error}")
    print("*")
    print(f"* - The Mean Squared Error of RV estimator : {RV_mse_error}")
    print("*")
    print(f"* - The Mean Squared Error of FE[day by day] : {FE_day_mse_error}")
    print("*")
    print(f"* - The Mean Squared Error of FE[all in all]: {FE_whole_mse_error}")
    print("********************************************************************************")

    
    # approximation_plot(x, [fourier_est, FZ_est, RV_est], ["fourier", "Florens-Zmirou", "Realied volatility"], (eta * x**gamma)**2)
    approximation_plot(
        x, 
        [fourier_est, FZ_est, RV_est,  FE_day_est, FE_whole_est], 
        ["fourier", "Florens-Zmirou", "Realied volatility", "FE day by day", "FE all in all"], 
        (eta * x**gamma)**2, label="simona_"
    )
    # approximation_plot(x, [fourier_est, RV_est], ["fourier", "Realized vol"], (eta * x**gamma)**2)
    # approximation_plot(x, [fourier_est], ["fourier"], (eta * x**gamma)**2)
    # approximation_plot(x, [FE_whole_est], ["FE_whole"], (eta * x**gamma)**2)
    # approximation_plot(x, fourier_est, (eta * x**gamma)**2, "fourier estimator[direct culc]") # direct culculation
    
    
if __name__ == "__main__":
    main()

