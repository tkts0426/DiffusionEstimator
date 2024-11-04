import matplotlib.pyplot as plt

def hist_plot(process_list, labels=[]):
    
    fig, ax = plt.subplots()
    
    for i, process in enumerate(process_list):
        ax.hist(process, bins=200, label=labels[i], alpha=0.5)
        ax.set_title("observed data")
    
    plt.legend()
    fig.savefig(f"figure/Distribution_of_obs_prices.png")

def cum_dist_plot(sorted_process, x_label=None, label=""):
    
    fig, ax = plt.subplots()
    
    n_samples = len(sorted_process)
    prob = [k / n_samples for k in range(n_samples)]
    ax.plot(sorted_process, prob)
    ax.axhline(0.5, ls = "-.", color = "gray")
    ax.axhline(0.95, ls = "-.", color = "gray")

    ax.set_title("Cumulative Distribution")
    ax.set_xlabel(x_label)
    ax.set_ylabel('probability')

    fig.savefig(f"figure/{label}Cumlative_distribution_of_obs_prices.png")

def approximation_plot(x, est_list, est_names, real, label=""):
    
    fig, ax = plt.subplots()
    
    for est, est_name in zip(est_list, est_names):
        ax.plot(x, est, label=est_name, linestyle="-") # plots of estimator
    
    ax.plot(x, real, color = "blue", label="sigma^2", linestyle="--")# plots of real function

    # ax.plot(x, fourier_est + 1.96 * np.sqrt(var / n_traj), color = "gray", linestyle="-", label="upper CI") # upper confidential interval
    # ax.plot(x, fourier_est - 1.96 * np.sqrt(var / n_traj), color = "gray", linestyle="-", label="lower CI") # lower confidential interval
    ax.set_title(f"Compare between estimaters and sigma^2 function")
    plt.legend()
    
    fig.savefig(f"figure/{label}plots of estimators and sigma^2 function.png")