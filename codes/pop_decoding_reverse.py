# Population decoding  of calcium imaging data

# How to run?
# 1. Make sure the individual data files of the mouse are stored.
# 2. P-values stored in the respective directory
# 3. Run the PrefDirection.R to create the csv files that end with _prefer.csv for task and _prefer_passive.csv
# 4. Now run the pop_decoding.py (upto montaging the slopes plot)
# 5. Run in R slopesSummary.R to create final.csv
# 6. Now run the rest of the pop_decoding.py 

import sys
import os
 
import numpy as np 
import pandas as pd
import multiprocessing as mp  # For parallel processing
 
import re	
import time	


from scipy.ndimage import gaussian_filter1d
from scipy import interpolate
from scipy.io import loadmat,savemat   
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d 
from scipy.stats import zscore 

import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt  
matplotlib.rcParams['axes.linewidth'] = 2

# Load local modules 

os.chdir('/media/olive/Research/oliver/IEMdecodingForCalciumData/scripts')
 
sys.path.append(os.getcwd())
from utils import *	 # it imports population response curve making methods

from mne.stats import permutation_cluster_test,permutation_cluster_1samp_test
import mne 

# Data path and parameter settings  
data_path =  '/media/olive/Research/oliver/data_down/' 
#Pvalues and Pref. Direction plots
pval_pref_path='/media/olive/Research/oliver/prefDir'

paradigm = 'task'   
data_path_task =  os.path.join(data_path,paradigm)  

paradigm = 'passive' 
data_path_passive =   os.path.join(data_path,paradigm)

decoding_res_path = '/media/olive/Research/oliver/decoding/pop_decoding_reverse/' 
decoding_res_data_path=os.path.join(decoding_res_path,'tuning_curves')
decoding_res_fig_path=os.path.join(decoding_res_path ,'plots')  
decoding_res_slopes_path=os.path.join(decoding_res_path ,'slopes') 

task_save_path=os.path.join(decoding_res_data_path,'task')
passive_save_path=os.path.join(decoding_res_data_path,'passive')

create_dir(decoding_res_path) 
create_dir(os.path.join(decoding_res_data_path))
create_dir(task_save_path)
create_dir(passive_save_path)
create_dir(decoding_res_fig_path) 
create_dir(decoding_res_slopes_path) 

os.chdir(decoding_res_path)
decoding_res_path=os.getcwd() 
decoding_res_data_path=os.path.join(decoding_res_path,'tuning_curves')
decoding_res_fig_path=os.path.join(decoding_res_path ,'plots')  
decoding_res_slopes_path=os.path.join(decoding_res_path ,'slopes') 

task_save_path=os.path.join(decoding_res_data_path,'task')
passive_save_path=os.path.join(decoding_res_data_path,'passive')
 
# Adjustable parameter settings for decoding  
fs=20	 # sampling frequency 
ts=1/fs   # sampling time

#Data parameter settings
analyse_time=[-1,5] # sec 
time_values=np.arange(analyse_time[0],analyse_time[1],ts)  
nt=len(time_values)

angles=np.arange(22.5,360,45) # Stimulus angles
ns=len(angles)   # no. of stimulus values 
# Parameters for slope computation
tt=np.arange(-1,5,ts)
trun=2
sig=1
wrap_around=5
slope_angles=[-180,-135,-90,-45,0]
  
Gaussian_smoothening_needed=False
time_resolved_tuning_desired=False
iteration_included=False
zscoring_rawdata_desired=True

baseline_correction_needed=True

if baseline_correction_needed:
	base_idx=np.arange(17,20)

nt=120   # 120 time-points for 20 Hz

# parameters for cluster computation
nperm,tail,cluster_alpha=5000,1,0.05

# parameters for plotting figures
xmin=-0.5
xmax=4

ymin=-0.05
ymax=0.3 


# Significant time points placer
# for task
#first_sig_task=0.0950
#second_sig_task=0.0925
#diff_sig_task=0.0900
# for passive
#first_sig_passive=0.0850
#second_sig_passive=0.0825
#diff_sig_passive=0.0800
 
  

ROIs_hetero=['V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135']

for k in ROIs_hetero:
	create_dir(os.path.join(decoding_res_data_path,'task',k))
	create_dir(os.path.join(decoding_res_data_path,'passive',k)) 
 
# Percentage of cells to use for decoding (tuned + % of untuned cells)
#percent_data=[0,10,20,40,60,100]	 
#percent_data=[0,10,25,50,100]
#percent_data=[0,10,50,100]
pvals_threshold=[0.9,0.5,0.05,0]
#########################################################################################
############################  ANALYSIS FOR TASK DATA #####################################
#########################################################################################
 
# Custom sorting key
def extract_number(animal):
	return int(animal[5:-4]) # Assuming the format is 'animalX'

  
paradigm='task' 
# mouse name order
#V1_task=['1R.mat','2L.mat','newV1.mat','Ylbd.mat','Ylbe.mat','Ylbf.mat','Ylbg.mat','2.mat',	'3.mat','4.mat']
#PPC_task=['3L.mat','3R.mat','YLbc.mat','YLcj.mat','YLck.mat','YLcl.mat','YLcm.mat']



print_status('Will begin decoding with the following setting ')
print_status('Gaussian_smoothening_needed: ' + str(Gaussian_smoothening_needed )) 
print_status('zscoring_rawdata_desired: ' +  str(zscoring_rawdata_desired))
print_status('baseline_correction_needed: ' +  str(baseline_correction_needed))
 
  

for roi in ROIs_hetero:   # For each heterogeneous condition
#for roi in ['V1_45']:	
	os.chdir(data_path_task)
	print_status('Dealing with ROI: ' + roi)
	
	os.chdir(roi)
	
	ani_list=sorted(os.listdir(),key=extract_number)
	print(ani_list)
	noa=len(ani_list)
	
	print_status('No. of animals in ' + roi + ' is ' + str(noa)) 
	
	# load the neuron preferences and pvalue for all animals 
	# The csv file is created in R (pls. look at the script PrefDirection.R in the scripts folder)
	B = pd.read_csv(os.path.join(pval_pref_path, paradigm,roi+'.csv'))
	
	
	#if roi.startswith('V1'):
	#	list_for_sorting=V1_task
	#else:
	#	list_for_sorting=PPC_task 
	
	#ani_list=[[file for file in ani_list if file.lower().endswith(suffix.lower())] for suffix in list_for_sorting]
	#ani_list= [ani[0] for ani in ani_list] 
	 
	for p in range(len(ani_list)):   # For each animal
		# load that animal data in homo and hetero group  
		# This is not sorted yet using p-value
		df_homo=B[(B['Sub']=='Animal.'+str(p+1))  & (B['Group'] =='homo')]
		df_hetero=B[(B['Sub']=='Animal.'+str(p+1))  & (B['Group'] =='hetero')] 
		
		# resetting the index is needed after subsetting the data
		df_homo.reset_index(drop=True, inplace=True)
		df_hetero.reset_index(drop=True, inplace=True)  # reset the indices
		 
		# arrange accoring to the p-value first
		#idx_homo=np.argsort(-df_homo['Pvalue'])  # decending order of sort
		#df_homo_sorted=df_homo.loc[idx_homo]
		#df_hetero_sorted=df_hetero.loc[idx_homo]  
				 
		print_status('Dealing with the animal ' + ani_list[p])
		
		A=loadmat(ani_list[p])	 # Load the data
		
		# homo and hetero data
		homo_data=A['sample_data_homo']	   #orig.data:	 trials X units X time-pts
		homo_data=homo_data.transpose(1,0,2)   #for decoding:  units X trials X time-pts  
		
		hetero_data=A['sample_data_hetero']
		hetero_data =hetero_data.transpose(1,0,2)  
		
		if zscoring_rawdata_desired: 
			homo_data = zscore(homo_data,axis=2)
			hetero_data = zscore(hetero_data,axis=2)
 
		
		if baseline_correction_needed:
			homo_data -= np.mean(homo_data[:, :, base_idx], axis=-1, keepdims=True)
			hetero_data -= np.mean(hetero_data[:, :, base_idx], axis=-1, keepdims=True)
 
		# homo and hetero data labels
		homo_labels=np.squeeze(A['dirIdx_homo'])  
		hetero_labels=np.squeeze(A['dirIdx_hetero'])  
		
		del A
		
		# Shuffle labels and data before decoding 
		# Not necessary for population decoding though
		idx=np.random.permutation(len(homo_labels))
		homo_labels=homo_labels[idx]
		homo_data=homo_data[:,idx,:]
		
		idx=np.random.permutation(len(hetero_labels))
		hetero_labels=hetero_labels[idx]
		hetero_data=hetero_data[:,idx,:] 
		
		for pp in pvals_threshold: 
		
			path_name=os.path.join(decoding_res_data_path,paradigm,roi,str(pp))
			create_dir(path_name)
			fname=path_name+'/Mouse'+str(p+1) + '.npy'  
			#fname=path_name+'/Mouse'+str(p+1) + '.mat'  			
			if not os.path.exists(fname):
                
                # to get a data frame that is tuned in both the condition
        		# finding indices that match both conditions
				final=pd.merge(df_homo[df_homo['Pvalue'] >= pp], df_hetero[df_hetero['Pvalue'] >= pp], left_index=True, right_index=True)

        		#final=pd.merge(df_homo_sorted[df_homo_sorted['Pvalue'] >= pp], df_hetero_sorted[df_hetero_sorted['Pvalue'] >= pp], left_index=True, right_index=True)
				NeuronIndices=final.index.to_numpy() 
        		
        		
        		# pref_df in only for tuned-tuned
				PrefDir = final[['Preference_x','Preference_y']].copy()

        		#pref_df = pref_df.copy()  # to avoid the copy warnings
        		#pref_df['Neuron'] = final.index.to_list()
				PrefDir=PrefDir.rename(columns={'Preference_x': 'Pref.Homo', 'Preference_y': 'Pref.Hetero'})
        		#pref_df.reset_index(drop=True, inplace=True) 
				PrefDir.loc[:,'Neuron'] = NeuronIndices

				# Get the indices of the other neurons and use them for decoding
				# Only homo data is enough to sort out the neurons
				# as we always choose homo tuned neurons and choose the same in the 
				# hetero condition
				#if pp!=0:  # Other than tuned on both homo and hetero 
					 
					 
				# first remove the tuned tuned one from the original df
				#df_homo_pp=df_homo[df_homo['Pvalue']>pp] 
				#df_hetero_pp=df_hetero[df_homo['Pvalue']>pp]
                
				#homo_indices=df_homo_pp.index.to_numpy()
				#hetero_indices=df_hetero_pp['Unnamed: 0'].to_numpy()					
					# additional neurons
					#nper_cells = int(np.ceil((pp/100)*df_homo_sorted_pp.shape[0])) # get the number of cells 
					
					#additional_neurons=df_homo_sorted_pp[:nper_cells].index.to_numpy()
					 
					
					# Neurons used for decoding
					#NeuronIndices=np.concatenate((TunedIndices,additional_neurons))  
					
					#update_df= pd.DataFrame()
					#update_df['Pref.Homo']=df_homo_sorted_pp[:nper_cells]['Preference']
					#update_df['Pref.Hetero']=df_hetero_sorted_pp[:nper_cells]['Preference']
					
					#update_df['Neuron']=df_homo_sorted_pp[:nper_cells].index
					
					# update the  pref_df (to include the new additional neurons)
					#PrefDirInfo=pd.concat([pref_df, update_df], ignore_index=True)
					
				#else:
				#	NeuronIndices=TunedIndices
				#	PrefDirInfo=pd.DataFrame()
				#	PrefDirInfo=pref_df
				
				homo_data_p=homo_data[NeuronIndices,:,:]  # arranging homo data accroding to sorted pvalue 
				hetero_data_p=hetero_data[NeuronIndices,:,:]  # arranging hetero data accroding to sorted pvalue 
				
				print_status('Homo. shape before decoding is ' + str(homo_data_p.shape))
				print_status('Hetero. shape before decoding is ' + str(hetero_data_p.shape))   
 
				# preference direction computation from data 
				for neu_idx in PrefDir['Neuron'].to_numpy().astype('int32'): 
					homo_data[neu_idx,:,:]=np.nan_to_num(homo_data[neu_idx,:,:],nan=0)
					hetero_data[neu_idx,:,:]=np.nan_to_num(hetero_data[neu_idx,:,:],nan=0)
					
					#neuron_index = PrefDir.loc[PrefDir['Neuron'] == neu_idx].index[0]
                    
					PrefDir.at[neu_idx,'Pref.Homo'] =get_preference(homo_data[neu_idx,:,:],homo_labels)[0]			   
					PrefDir.at[neu_idx,'Pref.Hetero']=get_preference(hetero_data[neu_idx,:,:],hetero_labels)[0]
				
                # Parallel decoding begings here
				st = time.time() 
				A=run_parallel_the_pop_decoding(homo_data_p,hetero_data_p,homo_labels,hetero_labels,PrefDir,nt)  # (Homo result, Hetero. result)
				ed = time.time()
					
				elapsed_time = ed - st
				print_status('Execution time: ' + str(elapsed_time) + ' for animal ' + str(p))   
					 
				print('Shape of A is ', A.shape)
				print_status('Saving the tuning curves') 
				path_name=os.path.join(decoding_res_data_path,paradigm,roi,str(pp))
				create_dir(path_name) 
				print(fname)
							
				np.save(fname,np.squeeze(A))
				#savemat(fname,{'PRC' : np.squeeze(A)})
				del A 
					
				print_status('Done with ' + str(p+1) + '/' + str(noa) + ' for the percentage '+ str(pp),'') 
					 
			else:
				print_status('Already done with ' + str(p+1) + '/' + str(noa) + ' for the percentage '+ str(pp),'') 

 
# Zero-centering and computing slopes for task case
 
##  Slope dynamics computation and storing the results
for roi in ROIs_hetero:   # For each condition  
	for pp in pvals_threshold: # For each percentage of data 
		
		os.chdir(os.path.join(decoding_res_data_path, paradigm))
		print_status('Computing slopes for  ROI' + roi +' for the percentage  ',str(pp)) 
	
		os.chdir(os.path.join(roi,str(pp)))

		ani_list=os.listdir()
		noa=len(ani_list)

		print_status('No. of animals in ' + roi + ' is ' + str(noa))

		st_roi = time.time()

		slopes=np.zeros((noa,nt,2))  # 2 conditions

		for p in range(len(ani_list)):   # For each animal
		
			save_dir=os.path.join(decoding_res_slopes_path,paradigm,roi,str(pp))
			create_dir(save_dir) 
			fname=os.path.join(save_dir,'slopes.npy')
			
			if not os.path.exists(fname): 

				A=np.load(ani_list[p])  # Load the tuning curve data  
				B=zero_center_the_tcs(A,shift=wrap_around)
				C=B.transpose((0,2,1))
				D=avg_across_zero_centered_tcs(C, shift=wrap_around)

				if Gaussian_smoothening_needed:
					print_status('Gaussian smoothening desired!') 
					S=Gaussian_smoothener(D,sig=sig,trun=trun ,ts=ts) 
				else:
					print_status('NO Gaussian smoothening desired!') 
					S=D  

				## Estimate the slope 
				# Homo-case
				print_status('Computing slopes..') 
				for time_pts in range(nt): 
					slopes[p,time_pts,0]=esti_slope(slope_angles,S[:,time_pts,0],intercept=True, standardise=False) 

				# Hetero-case 
				for time_pts in range(nt): 
					slopes[p,time_pts,1]=esti_slope(slope_angles,S[:,time_pts,1],intercept=True, standardise=False)   
				 
				
			else:
				print_status('Already done with slope computations') 
				
		np.save(fname,slopes)   
 		
############## passive data decoding   ##################  

paradigm='passive' 
# mouse name order
#V1_passive=['1L.mat','1R.mat','2L.mat','Ylce.mat','Ylcf.mat','Ylcg.mat']
#PPC_passive=['3L.mat','3R.mat','Ylcb.mat', 'Ylch.mat','Ylci.mat','YLck.mat','Ylcl.mat','Ylcm.mat']

for roi in ROIs_hetero:   # For each heterogeneous condition
#for roi in ['V1_45']:	
	os.chdir(data_path_passive)
	print_status('Dealing with ROI: ' + roi)
	
	os.chdir(roi)
	
	ani_list=sorted(os.listdir(),key=extract_number)
	print(ani_list)
	noa=len(ani_list)
	
	print_status('No. of animals in ' + roi + ' is ' + str(noa)) 
	
	# load the neuron preferences and pvalue for all animals 
	# The csv file is created in R (pls. look at the script PrefDirection.R in the scripts folder)
	B = pd.read_csv(os.path.join(pval_pref_path, paradigm,roi+'.csv'))
	
	
	#if roi.startswith('V1'):
	#	list_for_sorting=V1_task
	#else:
	#	list_for_sorting=PPC_task 
	
	#ani_list=[[file for file in ani_list if file.lower().endswith(suffix.lower())] for suffix in list_for_sorting]
	#ani_list= [ani[0] for ani in ani_list] 
	 
	for p in range(len(ani_list)):   # For each animal
		# load that animal data in homo and hetero group  
		# This is not sorted yet using p-value
		df_homo=B[(B['Sub']=='Animal.'+str(p+1))  & (B['Group'] =='homo')]
		df_hetero=B[(B['Sub']=='Animal.'+str(p+1))  & (B['Group'] =='hetero')] 
		
		# resetting the index is needed after subsetting the data
		df_homo.reset_index(drop=True, inplace=True)
		df_hetero.reset_index(drop=True, inplace=True)  # reset the indices
		 
		# arrange accoring to the p-value first
		#idx_homo=np.argsort(-df_homo['Pvalue'])  # decending order of sort
		#df_homo_sorted=df_homo.loc[idx_homo]
		#df_hetero_sorted=df_hetero.loc[idx_homo]  
				 
		print_status('Dealing with the animal ' + ani_list[p])
		
		A=loadmat(ani_list[p])	 # Load the data
		
		# homo and hetero data
		homo_data=A['sample_data_homo']	   #orig.data:	 trials X units X time-pts
		homo_data=homo_data.transpose(1,0,2)   #for decoding:  units X trials X time-pts  
		
		hetero_data=A['sample_data_hetero']
		hetero_data =hetero_data.transpose(1,0,2)  
		
		if zscoring_rawdata_desired: 
			homo_data = zscore(homo_data,axis=2)
			hetero_data = zscore(hetero_data,axis=2)
 
		
		if baseline_correction_needed:
			homo_data -= np.mean(homo_data[:, :, base_idx], axis=-1, keepdims=True)
			hetero_data -= np.mean(hetero_data[:, :, base_idx], axis=-1, keepdims=True)
 
		# homo and hetero data labels
		homo_labels=np.squeeze(A['dirIdx_homo'])  
		hetero_labels=np.squeeze(A['dirIdx_hetero'])  
		
		del A
		
		# Shuffle labels and data before decoding 
		# Not necessary for population decoding though
		idx=np.random.permutation(len(homo_labels))
		homo_labels=homo_labels[idx]
		homo_data=homo_data[:,idx,:]
		
		idx=np.random.permutation(len(hetero_labels))
		hetero_labels=hetero_labels[idx]
		hetero_data=hetero_data[:,idx,:] 
		
		for pp in pvals_threshold: 
		
			path_name=os.path.join(decoding_res_data_path,paradigm,roi,str(pp))
			create_dir(path_name)
			fname=path_name+'/Mouse'+str(p+1) + '.npy'  
			#fname=path_name+'/Mouse'+str(p+1) + '.mat' 
			
			if not os.path.exists(fname):
                
                # to get a data frame that is tuned in both the condition
        		# finding indices that match both conditions
				final=pd.merge(df_homo[df_homo['Pvalue'] >= pp], df_hetero[df_hetero['Pvalue'] >= pp], left_index=True, right_index=True)

        		#final=pd.merge(df_homo_sorted[df_homo_sorted['Pvalue'] >= pp], df_hetero_sorted[df_hetero_sorted['Pvalue'] >= pp], left_index=True, right_index=True)
				NeuronIndices=final.index.to_numpy() 
        		
        		
        		# pref_df in only for tuned-tuned
				PrefDir = final[['Preference_x','Preference_y']].copy()

        		#pref_df = pref_df.copy()  # to avoid the copy warnings
        		#pref_df['Neuron'] = final.index.to_list()
				PrefDir=PrefDir.rename(columns={'Preference_x': 'Pref.Homo', 'Preference_y': 'Pref.Hetero'})
        		#pref_df.reset_index(drop=True, inplace=True) 
				PrefDir.loc[:,'Neuron'] = NeuronIndices

				# Get the indices of the other neurons and use them for decoding
				# Only homo data is enough to sort out the neurons
				# as we always choose homo tuned neurons and choose the same in the 
				# hetero condition
				#if pp!=0:  # Other than tuned on both homo and hetero 
					 
					 
				# first remove the tuned tuned one from the original df
				#df_homo_pp=df_homo[df_homo['Pvalue']>pp] 
				#df_hetero_pp=df_hetero[df_homo['Pvalue']>pp]
                
				#homo_indices=df_homo_pp.index.to_numpy()
				#hetero_indices=df_hetero_pp['Unnamed: 0'].to_numpy()					
					# additional neurons
					#nper_cells = int(np.ceil((pp/100)*df_homo_sorted_pp.shape[0])) # get the number of cells 
					
					#additional_neurons=df_homo_sorted_pp[:nper_cells].index.to_numpy()
					 
					
					# Neurons used for decoding
					#NeuronIndices=np.concatenate((TunedIndices,additional_neurons))  
					
					#update_df= pd.DataFrame()
					#update_df['Pref.Homo']=df_homo_sorted_pp[:nper_cells]['Preference']
					#update_df['Pref.Hetero']=df_hetero_sorted_pp[:nper_cells]['Preference']
					
					#update_df['Neuron']=df_homo_sorted_pp[:nper_cells].index
					
					# update the  pref_df (to include the new additional neurons)
					#PrefDirInfo=pd.concat([pref_df, update_df], ignore_index=True)
					
				#else:
				#	NeuronIndices=TunedIndices
				#	PrefDirInfo=pd.DataFrame()
				#	PrefDirInfo=pref_df
				
				homo_data_p=homo_data[NeuronIndices,:,:]  # arranging homo data accroding to sorted pvalue 
				hetero_data_p=hetero_data[NeuronIndices,:,:]  # arranging hetero data accroding to sorted pvalue 
				
				print_status('Homo. shape before decoding is ' + str(homo_data_p.shape))
				print_status('Hetero. shape before decoding is ' + str(hetero_data_p.shape))   
 
				# preference direction computation from data 
				for neu_idx in PrefDir['Neuron'].to_numpy().astype('int64'): 
					homo_data[neu_idx,:,:]=np.nan_to_num(homo_data[neu_idx,:,:],nan=0)
					hetero_data[neu_idx,:,:]=np.nan_to_num(hetero_data[neu_idx,:,:],nan=0)
					
					#neuron_index = PrefDir.loc[PrefDir['Neuron'] == neu_idx].index[0]
                    
					PrefDir.at[neu_idx,'Pref.Homo'] =get_preference(homo_data[neu_idx,:,:],homo_labels)[0]			   
					PrefDir.at[neu_idx,'Pref.Hetero']=get_preference(hetero_data[neu_idx,:,:],hetero_labels)[0]
				
                # Parallel decoding begings here
				st = time.time() 
				A=run_parallel_the_pop_decoding(homo_data_p,hetero_data_p,homo_labels,hetero_labels,PrefDir,nt)  # (Homo result, Hetero. result)
				ed = time.time()
					
				elapsed_time = ed - st
				print_status('Execution time: ' + str(elapsed_time) + ' for animal ' + str(p))   
					 
				print('Shape of A is ', A.shape)
				print_status('Saving the tuning curves') 
				path_name=os.path.join(decoding_res_data_path,paradigm,roi,str(pp))
				create_dir(path_name) 
				print(fname)
							
				np.save(fname,np.squeeze(A))
				#savemat(fname,{'PRC' : np.squeeze(A)})
				del A 
					
				print_status('Done with ' + str(p+1) + '/' + str(noa) + ' for the percentage '+ str(pp),'') 
					 
			else:
				print_status('Already done with ' + str(p+1) + '/' + str(noa) + ' for the percentage '+ str(pp),'') 
 
 
# Centering and computing slopes for passive case  
# Slope dynamics computation and storing the results
for roi in ROIs_hetero:   # For each condition  
	for pp in pvals_threshold: # For each percentage of data 
		
		os.chdir(os.path.join(decoding_res_data_path, paradigm))
		print_status('Computing slopes for  ROI' + roi +' for the percentage  ',str(pp)) 
	
		os.chdir(os.path.join(roi,str(pp)))

		ani_list=os.listdir()
		noa=len(ani_list)

		print_status('No. of animals in ' + roi + ' is ' + str(noa))

		st_roi = time.time()

		slopes=np.zeros((noa,nt,2))  # 2 conditions

		for p in range(len(ani_list)):   # For each animal
		
			save_dir=os.path.join(decoding_res_slopes_path,paradigm,roi,str(pp))
			create_dir(save_dir) 
			fname=os.path.join(save_dir,'slopes.npy')
			
			if not os.path.exists(fname): 

				A=np.load(ani_list[p])  # Load the tuning curve data  
				B=zero_center_the_tcs(A,shift=wrap_around)
				C=B.transpose((0,2,1))
				D=avg_across_zero_centered_tcs(C, shift=wrap_around)

				if Gaussian_smoothening_needed:
					print_status('Gaussian smoothening desired!') 
					S=Gaussian_smoothener(D,sig=sig,trun=trun ,ts=ts) 
				else:
					print_status('NO Gaussian smoothening desired!') 
					S=D  

				## Estimate the slope 
				# Homo-case
				print_status('Computing slopes..') 
				for time_pts in range(nt): 
					slopes[p,time_pts,0]=esti_slope(slope_angles,S[:,time_pts,0],intercept=True, standardise=False) 

				# Hetero-case 
				for time_pts in range(nt): 
					slopes[p,time_pts,1]=esti_slope(slope_angles,S[:,time_pts,1],intercept=True, standardise=False)   
				  
			else:
				print_status('Already done with slope computations')  
		np.save(fname,slopes)   
 
##Plotting and montaging  (with task as solid  line and passive as dashed line) 

lwd=3
plt_lwd=3
alp=0.2

# for task
first_sig_task=0.3
second_sig_task=0.29
diff_sig_task=0.31   # black one
# for passive
first_sig_passive=0.26
second_sig_passive=0.25
diff_sig_passive=0.27  # black one


cols=['Paradigm', 'Roi', 'Condition', 'Percentage', 'Cluster p-value','Significant time points']
df = pd.DataFrame(columns=cols)

for roi in ROIs_hetero:  # for each roi
	for pp in pvals_threshold: # for each p-value percentage   
		
		paradigm='task'
		# For task (solid lines) 
		A_task=np.load(os.path.join(decoding_res_slopes_path,paradigm,roi,str(pp),'slopes.npy'))
		
		#A_task[:,:,0]=A_task[:,:,0]-np.mean(A_task[:,17:19,0])
		#A_task[:,:,1]=A_task[:,:,1]-np.mean(A_task[:,17:19,1])
		
		sig1=A_task[:,:,0]  # homo  # no. mouse X no. time points X (homo or hetero)
		sig2=A_task[:,:,1]  # hetero   
		
		#sig1=sig1-np.mean(sig1[:,17:19])
		#sig2=sig2-np.mean(sig2[:,17:19])
		
		
		fig, ax = plt.subplots(1,1,figsize=(7,7))	
		# plot the mean and error bar
		ax.plot(tt,np.mean(sig1,0),'r-',tt,np.mean(sig2,0),'b-',linewidth=plt_lwd) 
		ax.fill_between(tt,np.mean(sig1,0)- (np.std(sig1,0)/np.sqrt(sig1.shape[0])),  
						   np.mean(sig1,0)+ (np.std(sig1,0)/np.sqrt(sig1.shape[0])),alpha=alp,color='r')
		ax.fill_between(tt,np.mean(sig2,0)- (np.std(sig2,0)/np.sqrt(sig2.shape[0])),  
						   np.mean(sig2,0)+ (np.std(sig2,0)/np.sqrt(sig2.shape[0])),alpha=alp,color='b') 

		
		# permutation clustering for homo
		T_obs, clusters, cluster_p_values, H0 = permutation_cluster_1samp_test(sig1,
																			   tail=1,
																			   n_permutations=nperm)
		clus=[]
		p=0
		for k in cluster_p_values:
			if k<(cluster_alpha):
				clus.extend(clusters[p])
				p=p+1
			else:
				p=p+1   
	 
		
		if len(clus):	
			clus=np.concatenate(clus)
			clus=clus[(clus>20) & (clus<101)] 
			
			if len(clus): 
				sig_tt=tt[clus]  # Significant time points
				ax.plot(sig_tt,np.repeat(first_sig_task,len(sig_tt)),'r-', linewidth=lwd) 
				slope_sig1=np.mean(A_task[:,clus,0],1)  
			else:
				cluster_p_values=1
				sig_tt=[None,None]
				print_status('No Significant clusters in ' + roi  + '  ' + str(pp) +' (homo) case')
				slope_sig1=np.zeros(sig1.shape[0])
				
		else:
			cluster_p_values=1
			sig_tt=[None,None]
			print_status('No Significant clusters in ' + roi  + '  ' + str(pp) +' (homo) case')
			slope_sig1=np.zeros(sig1.shape[0])
			
		sig_tt_info = {'Paradigm': paradigm  ,'Roi': roi, 'Condition': 'homo' , 'Percentage': pp ,'Cluster p-value': np.min(cluster_p_values) ,'Significant time points' :  [sig_tt[0],sig_tt[-1]]}
		df.loc[len(df)] = sig_tt_info
		
		#print('Slope sig 1', slope_sig1)
		# permutation clustering for hetero
		T_obs, clusters, cluster_p_values, H0 = permutation_cluster_1samp_test(sig2,
																			   tail=1,
																		   n_permutations=nperm)
		
		clus=[]
		p=0
		for k in cluster_p_values:
			if k<(cluster_alpha):
				clus.extend(clusters[p])
				p=p+1
			else:
				p=p+1   
				
		if len(clus): 
			clus=np.concatenate(clus)
			clus=clus[(clus>20) & (clus<101)]
			
			if len(clus):  
				sig_tt=tt[clus]  # Significant time points
				ax.plot(sig_tt,np.repeat(second_sig_task,len(sig_tt)),'b-', linewidth=lwd) 
				slope_sig2=np.mean(A_task[:,clus,1],1)
				sig_tt=sig_tt[sig_tt<=4.05] 
			else:
				cluster_p_values=1
				sig_tt=[None,None]
				print_status('No Significant clusters in ' + roi + '  ' + str(pp) +' (hetero) case')
				slope_sig2=np.zeros(sig2.shape[0]) 
		else:
			cluster_p_values=1
			sig_tt=[None,None]
			print_status('No Significant clusters in ' + roi + '  ' + str(pp) +' (hetero) case')
			slope_sig2=np.zeros(sig2.shape[0])
		
		sig_tt_info = {'Paradigm': paradigm  ,'Roi': roi, 'Condition': 'hetero' , 'Percentage': pp ,'Cluster p-value': np.min(cluster_p_values) ,'Significant time points' :  [sig_tt[0],sig_tt[-1]]}
		df.loc[len(df)] = sig_tt_info
		

		res=np.column_stack((slope_sig1,slope_sig2))
		np.savetxt(os.path.join('/media/olive/Research/oliver/pop_slopes/',paradigm,roi+'_'+str(pp)+'.csv'),res,delimiter=',')

		# permuation clustering for homo-hetero
		T_obs, clusters, cluster_p_values, H0 = permutation_cluster_1samp_test(sig1-sig2,
																			   tail=0,
																			   n_permutations=nperm)
		clus=[]
		p=0
		for k in cluster_p_values:
			if k<(cluster_alpha):
				clus.extend(clusters[p])
				p=p+1
			else: 
				p=p+1  
 
				
		if len(clus): 
			clus=np.concatenate(clus)
			clus=clus[(clus>20) & (clus<101)]
			
			if len(clus):
				sig_tt=tt[clus]  # Significant time points
				ax.plot(sig_tt,np.repeat(diff_sig_task,len(sig_tt)),'k-', linewidth=lwd) 
				#sig_tt=sig_tt[sig_tt<=4.05]
			else:
				cluster_p_values=1
				sig_tt=[None,None]
		else:
			cluster_p_values=1
			sig_tt=[None,None]
			
		sig_tt_info = {'Paradigm': paradigm  ,'Roi': roi, 'Condition': 'diff' , 'Percentage': pp ,'Cluster p-value': np.min(cluster_p_values) ,'Significant time points' :  [sig_tt[0],sig_tt[-1]]}
		df.loc[len(df)] = sig_tt_info
		
			
		## plotting for passive 
		paradigm='passive'
		A_passive=np.load(os.path.join(decoding_res_slopes_path,paradigm,roi,str(pp),'slopes.npy'))
        
		A_passive=np.nan_to_num(A_passive,nan=0)
		
		#A_passive[:,:,0]=A_passive[:,:,0]-np.mean(A_passive[:,17:19,0])
		#A_passive[:,:,1]=A_passive[:,:,1]-np.mean(A_passive[:,17:19,1])
		
		sig1=A_passive[:,:,0]  # homo  # no. mouse X no. time points X (homo or hetero)
		sig2=A_passive[:,:,1]  # hetero	
 
		#sig1=sig1-np.mean(sig1[:,17:19])
		#sig2=sig2-np.mean(sig2[:,17:19])
		 
		# plot the mean and error bar
		ax.plot(tt,np.mean(sig1,0),'r--',tt,np.mean(sig2,0),'b--',linewidth=plt_lwd) 
		ax.fill_between(tt,np.mean(sig1,0)- (np.std(sig1,0)/np.sqrt(sig1.shape[0])),  
						   np.mean(sig1,0)+ (np.std(sig1,0)/np.sqrt(sig1.shape[0])),alpha=alp,color='r')
		ax.fill_between(tt,np.mean(sig2,0)- (np.std(sig2,0)/np.sqrt(sig2.shape[0])),  
						   np.mean(sig2,0)+ (np.std(sig2,0)/np.sqrt(sig2.shape[0])),alpha=alp,color='b') 

		
		# permutation clustering for homo
		T_obs, clusters, cluster_p_values, H0 = permutation_cluster_1samp_test(sig1,
																			   tail=1,
																			   n_permutations=nperm)
		clus=[]
		p=0
		for k in cluster_p_values:
			if k<(cluster_alpha):
				clus.extend(clusters[p])
				p=p+1
			else:
				p=p+1 
				
		

		if len(clus): 
			
			clus=np.concatenate(clus)
			clus=clus[(clus>20) & (clus<101)] 
			
			if len(clus):
			
				sig_tt=tt[clus]  # Significant time points
				ax.plot(sig_tt,np.repeat(first_sig_passive,len(sig_tt)),'r--', linewidth=lwd) 
				slope_sig1=np.mean(A_passive[:,clus,0],1)	
				sig_tt=sig_tt[sig_tt<=4.05] 
			else:
				cluster_p_values=1
				sig_tt=[None,None]
				print_status('No Significant clusters in ' + roi + '  ' + str(pp) +' (homo) case')
				slope_sig1=np.zeros(sig1.shape[0])
				
		else:
			cluster_p_values=1
			sig_tt=[None,None]
			print_status('No Significant clusters in ' + roi + '  ' + str(pp) +' (homo) case')
			slope_sig1=np.zeros(sig1.shape[0])
		
		sig_tt_info = {'Paradigm': paradigm  ,'Roi': roi, 'Condition': 'homo' , 'Percentage': pp ,'Cluster p-value': np.min(cluster_p_values) ,'Significant time points' :  [sig_tt[0],sig_tt[-1]]}
		df.loc[len(df)] = sig_tt_info
		
		
		#print('Slope sig 1', slope_sig1)
		# permutation clustering for hetero
		T_obs, clusters, cluster_p_values, H0 = permutation_cluster_1samp_test(sig2,
																			   tail=1,
																		   n_permutations=nperm)
		
		clus=[]
		p=0
		for k in cluster_p_values:
			if k<(cluster_alpha):
				clus.extend(clusters[p])
				p=p+1
			else:
				p=p+1  
		
				
		if len(clus): 
			clus=np.concatenate(clus)
			clus=clus[(clus>20) & (clus<101)] 
			
			if len(clus):
				
				sig_tt=tt[clus]  # Significant time points
				ax.plot(sig_tt,np.repeat(second_sig_passive,len(sig_tt)),'b--', linewidth=lwd) 
			
				slope_sig2=np.mean(A_passive[:,clus,1],1)
				sig_tt=sig_tt[sig_tt<=4.05] 
			else:
				cluster_p_values=1
				sig_tt=[None,None]
				print_status('No Significant clusters in ' + roi + '  ' + str(pp) +' (hetero) case')
				slope_sig2=np.zeros(sig2.shape[0])
				
		else:
			cluster_p_values=1
			sig_tt=[None,None]
			print_status('No Significant clusters in ' + roi + '  ' + str(pp) +' (hetero) case')
			slope_sig2=np.zeros(sig2.shape[0])
		
		
		sig_tt_info = {'Paradigm': paradigm  ,'Roi': roi, 'Condition': 'hetero' , 'Percentage': pp ,'Cluster p-value': np.min(cluster_p_values) ,'Significant time points' :  [sig_tt[0],sig_tt[-1]]}
		df.loc[len(df)] = sig_tt_info
		

		res=np.column_stack((slope_sig1,slope_sig2))
		np.savetxt(os.path.join('/media/olive/Research/oliver/pop_slopes/',paradigm,roi+'_'+str(pp)+'.csv'),res,delimiter=',')

		# permuation clustering for homo-hetero
		T_obs, clusters, cluster_p_values, H0 = permutation_cluster_1samp_test(sig1-sig2,
																			   tail=0,
																			   n_permutations=nperm)
		clus=[]
		p=0
		for k in cluster_p_values:
			if k<(cluster_alpha):
				clus.extend(clusters[p])
				p=p+1
			else:
				p=p+1  
				
		if len(clus):
			clus=np.concatenate(clus)
			clus=clus[(clus>20) & (clus<101)]
			
			if len(clus):  
				sig_tt=tt[clus]  # Significant time points
				ax.plot(sig_tt,np.repeat(diff_sig_passive,len(sig_tt)),'k--', linewidth=lwd) 
				sig_tt=sig_tt[sig_tt<=4.05] 
			else:
				cluster_p_values=1
				sig_tt=[None,None]
		else:
			cluster_p_values=1
			sig_tt=[None,None]
			
		sig_tt_info = {'Paradigm': paradigm  ,'Roi': roi, 'Condition': 'diff' , 'Percentage': pp ,'Cluster p-value': np.min(cluster_p_values) ,'Significant time points' :  [sig_tt[0],sig_tt[-1]]}
		df.loc[len(df)] = sig_tt_info
		
		
		ax.set_xticks([0, 1,2,3,4,5]) 
		ax.set_yticks([-0.05,0.1,0.2,0.3]) 
		ax.set_xlim(xmin, xmax) 
		ax.set_ylim(-0.05, 0.32)
		ax.axvline(x=0,color='k',linestyle='--',lw=1)  
		ax.axhline(y=0,color='k',linestyle='--',lw=1)  
		#ax.set_ylim(round(ymin,2), round(ymax,2))  
		ax.spines[['top','right']].set_visible(False) 
		ax.spines[['bottom','left']].set_linewidth(3)
 
		#ax.text(2.5,0.01, '(N=8)',fontsize=32) 
		#ax.text(-1.3,-0.077, '$-0.5$',fontsize=24) 
		ax.tick_params(axis='both', which='major', labelsize=24,length=10,width=2) 
		
		fig.tight_layout(pad=2)   
		#plt.show() 
		save_file_name=  roi + '_'+str(pp)+'.png'
		#save_file_name=   roi + '_'+str(pp)+'.pdf'
		fig.savefig(os.path.join(decoding_res_fig_path,save_file_name)) 
		os.chdir('..') 
		
	

# Extract the significant time points
df.to_excel("/media/olive/Research/oliver/IEMdecodingForCalciumData/neuron_counts/Significant.xlsx", index=False)
	
 
# montaging (this will work only if your system is Linux and montage installed))
if is_montage_installed():
	os.chdir(decoding_res_fig_path)
	create_dir('montages') 
	
	############### Plotting  summary of slopes ############
	fname='montages/V1_vert.png'
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in pvals_threshold for k in ROIs_hetero[:3]]) + '  -tile 3x' + str(len(pvals_threshold))  + '  -geometry +1+1 ' + fname
	status=os.system(strn) 
 
	fname='montages/PPC_vert.png'
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in pvals_threshold for k in ROIs_hetero[3:]]) + '  -tile 3x' + str(len(pvals_threshold))  + '  -geometry +1+1 ' + fname
	status=os.system(strn)  
	
	fname='montages/V1_45.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in pvals_threshold for k in ROIs_hetero[:1]]) + ' -tile ' + str(len(pvals_threshold))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 

	fname='montages/V1_90.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in pvals_threshold for k in ROIs_hetero[1:2]]) + ' -tile ' + str(len(pvals_threshold))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
    
	fname='montages/V1_135.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in pvals_threshold for k in ROIs_hetero[2:3]]) + ' -tile ' + str(len(pvals_threshold))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
			
	fname='montages/PPC_45.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in pvals_threshold for k in ROIs_hetero[3:4]]) + ' -tile ' + str(len(pvals_threshold))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
    
	fname='montages/PPC_90.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in pvals_threshold for k in ROIs_hetero[4:5]]) + ' -tile ' + str(len(pvals_threshold))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
    
	fname='montages/PPC_135.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in pvals_threshold for k in ROIs_hetero[5:6]]) + ' -tile ' + str(len(pvals_threshold))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
	
	os.chdir('montages')
	
	fname='V1_hori.png'
	status=os.system('montage V1_45.png  V1_90.png  V1_135.png  -tile 1x3  -geometry +1+1 ' + fname)

	fname='PPC_hori.png'
	status=os.system('montage PPC_45.png  PPC_90.png  PPC_135.png  -tile 1x3  -geometry +1+1 ' + fname)


else:
	print_status('Montage NOT installed in your computer. Skipping...') 
 
 
## Make the Final.csv in R first (in slopes_summary.R)
  
# Plot the summary of the slopes 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

final=pd.read_csv("/media/olive/Research/oliver/IEMdecodingForCalciumData/neuron_counts/final.csv")
conds=list(np.unique(final['Condition'])) 

#xx = list(range(len(pvals_threshold)))
xx=pvals_threshold[::-1]
percent_data_str=[str(k) for k in xx]

for cond in conds:
	df=final[(final['Condition'].str.contains(cond))] 
	
	df_task_homo=df[(df['Paradigm']=='Task') & (df['variable']=='Homo')]
	df_task_homo.reset_index(drop=True, inplace=True)
	df_task_hetero=df[(df['Paradigm']=='Task') & (df['variable']=='Hetero')]
	df_task_hetero.reset_index(drop=True, inplace=True)
	df_passive_homo=df[(df['Paradigm']=='Passive') & (df['variable']=='Homo')]
	df_passive_homo.reset_index(drop=True, inplace=True)
	df_passive_hetero=df[(df['Paradigm']=='Passive') & (df['variable']=='Hetero')]
	df_passive_hetero.reset_index(drop=True, inplace=True)
	
	fig, ax = plt.subplots(1,1,figsize=(7,7))
	
	# Plot lines 
	ax.plot(xx,df_task_homo['mean_value'],'r-') 
	ax.plot(xx,df_task_homo['mean_value'],'ro')
	ax.errorbar(xx,df_task_homo['mean_value'], yerr=df_task_homo['se_value'],
				fmt='none', capsize=5,color='r')
	
	ax.plot(xx,df_passive_homo['mean_value'],'r--') 
	ax.plot(xx,df_passive_homo['mean_value'],'ro')
	ax.errorbar(xx,df_passive_homo['mean_value'], yerr=df_passive_homo['se_value'],
				fmt='none', capsize=5,color='r') 
	
	ax.plot(xx,df_task_hetero['mean_value'],'b-') 
	ax.plot(xx,df_task_hetero['mean_value'],'bo')
	ax.errorbar(xx,df_task_hetero['mean_value'], yerr=df_task_hetero['se_value'],
				fmt='none', capsize=5,color='b')
	
	ax.plot(xx,df_passive_hetero['mean_value'],'b--') 
	ax.plot(xx,df_passive_hetero['mean_value'],'bo')
	ax.errorbar(xx,df_passive_hetero['mean_value'], yerr=df_passive_hetero['se_value'],
				fmt='none', capsize=5,color='b') 
 
	ax.set_xticks(xx,percent_data_str)
	ax.set_xticklabels( ["", "0.05", "0.5", "0.9"],  rotation=0, fontsize=10)  # to omit the zero
	ax.invert_xaxis()  # for plotting in the reverse order
	#ax.set_xlim(0.9, 0)
	ax.spines[['top','right']].set_visible(False) 
	ax.spines[['bottom','left']].set_linewidth(3) 
	ax.tick_params(axis='both', which='major', labelsize=20,length=5,width=2) 
	#ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
	ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))
	ax.set_ylim(-0.005, 0.20)
	ax.set_yticks([0, 0.1, 0.2])
	#ax.set_yticklabels(["0", "0.1", "0.2"])

		
	fig.tight_layout(pad=2)   
	#plt.show() 
	#save_file_name='Summary_' + cond.replace(' ','_') +'.pdf'
	save_file_name='Summary_' + cond.replace(' ','_') +'.png'
	fig.savefig(os.path.join(decoding_res_fig_path,save_file_name))  
	plt.show()
    
 

# Montage the summary files 
if is_montage_installed():
	os.chdir(decoding_res_fig_path)
	create_dir('montages') 
	
	fname='montages/Summary_V1.png'
	status=os.system('montage Summary_V1_45.png  Summary_V1_90.png  Summary_V1_135.png -tile 3x1  -geometry +1+1 ' + fname) 
 
	fname='montages/Summary_PPC.png'
	status=os.system('montage Summary_PPC_45.png  Summary_PPC_90.png  Summary_PPC_135.png -tile 3x1  -geometry +1+1 ' + fname)  
	
    


    
# For vertical 
os.chdir(decoding_res_fig_path)
strn= 'montage V1_45_0.9.png V1_90_0.9.png V1_135_0.9.png   V1_45_0.5.png V1_90_0.5.png V1_135_0.5.png  V1_45_0.05.png V1_90_0.05.png V1_135_0.05.png  V1_45_0.png V1_90_0.png V1_135_0.png  Summary_V1_45.png  Summary_V1_90.png Summary_V1_135.png  -tile 3x5  -geometry +1+1 montages/V1_vert_summary.png'
os.system(strn) 

strn= 'montage PPC_45_0.9.png PPC_90_0.9.png PPC_135_0.9.png   PPC_45_0.5.png PPC_90_0.5.png PPC_135_0.5.png  PPC_45_0.05.png PPC_90_0.05.png PPC_135_0.05.png  PPC_45_0.png PPC_90_0.png PPC_135_0.png  Summary_PPC_45.png  Summary_PPC_90.png Summary_PPC_135.png  -tile 3x5  -geometry +1+1 montages/PPC_vert_summary.png'
os.system(strn) 
 
"""  
# montaging (this will work only if your system is Linux and montage installed))
if is_montage_installed():
	os.chdir(decoding_res_fig_path)
	create_dir('montages') 
	
	############### Plotting  summary of slopes ############
	fname='montages/V1_vert.png'
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in percent_data for k in ROIs_hetero[:3]]) + '  -tile 3x' + str(len(percent_data_str))  + '  -geometry +1+1 ' + fname
	status=os.system(strn) 
 
	fname='montages/PPC_vert.png'
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in percent_data for k in ROIs_hetero[3:]]) + '  -tile 3x' + str(len(percent_data))  + '  -geometry +1+1 ' + fname
	status=os.system(strn)  
	
	fname='montages/V1_45.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in percent_data for k in ROIs_hetero[:1]]) + ' -tile ' + str(len(percent_data))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 

	fname='montages/V1_90.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in percent_data for k in ROIs_hetero[1:2]]) + ' -tile ' + str(len(percent_data))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
    
	fname='montages/V1_135.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in percent_data for k in ROIs_hetero[2:3]]) + ' -tile ' + str(len(percent_data))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
			
	fname='montages/PPC_45.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in percent_data for k in ROIs_hetero[3:4]]) + ' -tile ' + str(len(percent_data))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
    
	fname='montages/PPC_90.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in percent_data for k in ROIs_hetero[4:5]]) + ' -tile ' + str(len(percent_data))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
    
	fname='montages/PPC_135.png' 
	strn='montage ' + " ".join([k+'_'+str(l)+'.png'  for l in percent_data for k in ROIs_hetero[5:6]]) + ' -tile ' + str(len(percent_data))  + 'x1  -geometry +1+1 ' + fname
	status=os.system(strn) 
	
	os.chdir('montages')
	
	fname='V1_hori.png'
	status=os.system('montage V1_45.png  V1_90.png  V1_135.png  -tile 1x3  -geometry +1+1 ' + fname)

	fname='PPC_hori.png'
	status=os.system('montage PPC_45.png  PPC_90.png  PPC_135.png  -tile 1x3  -geometry +1+1 ' + fname)


else:
	print_status('Montage NOT installed in your computer. Skipping...')  
	
"""

