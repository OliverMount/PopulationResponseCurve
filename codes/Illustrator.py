### Illustration of population response curve
import numpy as np
import matplotlib.pyplot as plt 
 
np.random.seed(2)  # for reproducibility (vary this and play around)

###################################################
# PLOT 1. Neurons tuning curves
###################################################

fig_size=(7, 5) 
mean_values = np.arange(22.5, 360, 45)  # Stimulus values (here mean motion direction)
num_points=len(mean_values)


colors=  ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf']
neurons_per_stimulus = 2 # Number of neruons per stimulus 
# Initialize an array to store individual neuron tuning curves
tuning_curves = np.zeros((num_points, num_points * neurons_per_stimulus))
 
# Plotting tuning curves of inidividual neurons 
fig, ax = plt.subplots(1, figsize=fig_size) 

for idx, mean in enumerate(mean_values):
    for curve_num in range(neurons_per_stimulus): 
        std_dev = np.random.uniform(40, 60) 
        peak_amplitude = np.random.uniform(1, 10) 
        x = mean_values 
        y = peak_amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2)) 
        tuning_curves[:, idx * neurons_per_stimulus + curve_num] = y 
        ax.scatter(x, y, color=colors[idx], marker='o', s=25) 
        ax.plot(x, y, color=colors[idx], alpha=0.8,lw=3)
 
 
ax.set_xticks(x,range(num_points))
ax.set_xlabel("Stimulus index",fontsize=20)
for mean in mean_values:
    ax.axvline(x=mean, color='gray', linestyle='--', alpha=0.5)  
ax.set_ylim(-0.1,10) 
 
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['bottom', 'left']].set_linewidth(3)
ax.tick_params(axis='y', which='major', left=False, right=False, labelleft=False)
ax.tick_params(axis='x', which='major',labelsize=20,pad=8,length=10,width=2,direction='out')

fig.tight_layout(pad=2)  
plt.show()
fig.savefig("/home/olive/Desktop/NeuronsTuningCurves.png",dpi=300) 


###################################################
# PLOT 2. Construction of population tuning curve for a give stimulus (here it is 3)
###################################################

given_stimulus=3
row_to_plot = tuning_curves[given_stimulus]  
reshaped_values = row_to_plot.reshape(-1, neurons_per_stimulus) 
flattened_values = reshaped_values.flatten()
averaged_row = np.mean(reshaped_values, axis=1)

x_positions = np.arange(0, num_points)


fig, ax = plt.subplots(1, figsize=fig_size)
# Plot the scatter plot
for k in range(len(x_positions)):
	ax.scatter(np.repeat(x_positions[k], neurons_per_stimulus), flattened_values[neurons_per_stimulus*k:neurons_per_stimulus*k+neurons_per_stimulus], color=colors[k], marker='o')
ax.plot(x_positions,averaged_row, color='k', linestyle='-', linewidth=3,marker='o',alpha=0.3)

ax.set_xlabel("Stimulus index",fontsize=20)
ax.set_xticks(x_positions)
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['bottom', 'left']].set_linewidth(3)
ax.tick_params(axis='y', which='major', left=False, right=False, labelleft=False)
ax.tick_params(axis='x', which='major',labelsize=20,pad=8,length=10,width=2,direction='out')
ax.set_ylim(-0.1,10)


fig.tight_layout(pad=2)  
plt.show()
fig.savefig("/home/olive/Desktop/Stimulus3_PopRespCur.png",dpi=300) 


###################################################
# PLOT 3. An array of population response curves for all stimulus
################################################### 

fig, ax = plt.subplots(1, figsize=fig_size)
for given_stimulus in range(num_points):
 
    row_to_plot = tuning_curves[given_stimulus]   
    averaged_row = np.mean(row_to_plot.reshape(-1, neurons_per_stimulus), axis=1)
    reshaped_values = row_to_plot.reshape(-1, neurons_per_stimulus)
    flattened_values = reshaped_values.flatten()
    
    x_positions = np.arange(0, num_points)  
    ax.plot(averaged_row, color=colors[given_stimulus], linestyle='-', linewidth=3,marker='o')
    
ax.set_xticks(x_positions)
ax.set_xlabel("Stimulus index",fontsize=20)
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['bottom', 'left']].set_linewidth(3)
ax.tick_params(axis='y', which='major', left=False, right=False, labelleft=False)
ax.tick_params(axis='x', which='major',labelsize=20,pad=8,length=10,width=2,direction='out')
ax.set_ylim(-0.1,10)
    
fig.tight_layout(pad=2)  
plt.show()
fig.savefig("/home/olive/Desktop/All_PopRespCurs.png",dpi=300) 
