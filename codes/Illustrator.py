### Conceptual figure making

import numpy as np
import matplotlib.pyplot as plt

 
#np.random.seed(42)

fig_size=(7, 5)
 
mean_values = np.arange(22.5, 360, 45)
num_points=len(mean_values)

num_curves = 2 # Number of neruons per stimulus 
# Initialize an array to store Gaussian curves
gaussian_curves = np.zeros((num_points, num_points * num_curves))
 

fig, ax = plt.subplots(1, figsize=fig_size)
 
colors=  ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf']

for idx, mean in enumerate(mean_values):
    for curve_num in range(num_curves): 
        std_dev = np.random.uniform(40, 60) 
        peak_amplitude = np.random.uniform(1, 10) 
        x = mean_values 
        y = peak_amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2)) 
        gaussian_curves[:, idx * num_curves + curve_num] = y 
        ax.scatter(x, y, color=colors[idx], marker='o', s=25) 
        ax.plot(x, y, color=colors[idx], alpha=0.8,lw=3)
 
 
ax.set_xticks([]) 
for mean in mean_values:
    ax.axvline(x=mean, color='gray', linestyle='--', alpha=0.5) 
 
ax.set_ylim(-0.1,10) 
 
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['bottom', 'left']].set_linewidth(3)
ax.tick_params(axis='both', which='major', left=False, right=False, labelleft=False)

fig.tight_layout(pad=2)  
plt.show()
fig.savefig("~/Desktop/Homo_neurons.tiff",dpi=300)
 

given_stimulus=3
row_to_plot = gaussian_curves[given_stimulus]   
averaged_row = np.mean(row_to_plot.reshape(-1, 2), axis=1)
reshaped_values = row_to_plot.reshape(-1, 2)
flattened_values = reshaped_values.flatten()

x_positions = np.arange(0, 8)


fig, ax = plt.subplots(1, figsize=fig_size)
# Plot the scatter plot
for k in range(len(x_positions)):
	ax.scatter(np.repeat(x_positions[k], 2), flattened_values[2*k:2*k+2], color=colors[k], marker='o')
ax.plot(averaged_row, color='k', linestyle='-', linewidth=3,marker='o',alpha=0.3)

ax.set_xticks([])
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['bottom', 'left']].set_linewidth(3)
ax.tick_params(axis='both', which='major', left=False, right=False, labelleft=False)
ax.set_ylim(-0.1,10)

fig.tight_layout(pad=2)  
plt.show()
fig.savefig("~/Desktop/Homo_PopTuning.tiff",dpi=300)


fig, ax = plt.subplots(1, figsize=fig_size)
for given_stimulus in range(8):
 
    row_to_plot = gaussian_curves[given_stimulus]   
    averaged_row = np.mean(row_to_plot.reshape(-1, 2), axis=1)
    reshaped_values = row_to_plot.reshape(-1, 2)
    flattened_values = reshaped_values.flatten()
    
    x_positions = np.arange(0, 8) 
    # Plot the scatter plot
    #for k in range(len(x_positions)):
    #	ax.scatter(np.repeat(x_positions[k], 2), flattened_values[2*k:2*k+2], color=colors[k], marker='o')
    ax.plot(averaged_row, color='k', linestyle='-', linewidth=3,marker='o',alpha=0.3)
    
ax.set_xticks([])
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['bottom', 'left']].set_linewidth(3)
ax.tick_params(axis='both', which='major', left=False, right=False, labelleft=False)
ax.set_ylim(-0.1,10)
    
fig.tight_layout(pad=2)  
plt.show()