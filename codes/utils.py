import os
import sys
import numpy as np
import subprocess
import multiprocessing as mp  # For parallel processing

from sklearn import *
from sklearn.model_selection import LeaveOneOut 
from sklearn.utils.validation import check_is_fitted 

from scipy.ndimage import gaussian_filter1d
from scipy import interpolate
from scipy.io import loadmat,savemat   
from scipy.ndimage import gaussian_filter1d
from scipy.stats import zscore 

angles=np.arange(22.5,360,45) # Stimulus angles 
center_around=5  # Center the tuning curves around this angle 
ns=len(angles)


def create_dir(p): 
	"""
	Creates a directory given a path 'p'
	Examples:
	---------
	create_dir('test') -> Creates a folder test in the present working directory
	create_dir('/home/user/test/) -> Creates a folder test in home/user/ directory
	"""
	isExist = os.path.exists(p)  
	if not isExist:  
		os.makedirs(p)
		print('The directory ' +  p   + ' is created!')
	else:
		print('The directory ' +  p   + ' already exists. Nothing created!')
	return p


def print_status(msg,kind='short',line=50):
	
	if kind=='short':
		print('++ '+msg)
	else:
		print('++ '+line*'-')
		print('++ '+msg)
		print('++ '+line*'-')

def avg_across_zero_centered_tcs(data,shift): 
	print_status('Averaging equi-distant orientations')	 
	ns,nt,noc=data.shape
	unsmoothed_tc=np.zeros((shift,nt,noc)) 
  
	for p in range(noc): 
		unsmoothed_tc[:,:,p]=average_zero_centered_curves(data[:,:,p],shift) 

	return unsmoothed_tc  

def average_zero_centered_curves(A,p): 
	# input A is assumed to be 2D matrix
	# First dimension is number of stimulus

	B_left=np.zeros_like(A[:p,:])  # Allocate shape
	B_left=A[:p,:]					# slice upto the zero-centered value
	B_right=np.flip(A[p:,:],axis=0)   # flip the right tail of the tuning curve to the left side
	lr=B_right.shape[0] 
	B_left[(p-1-lr):(p-1),:]= (B_left[(p-1-lr):(p-1),:]+B_right)/2 # average the equi-distant stimulus orientation
	
	return B_left 


def esti_slope(angles,y,intercept=False,standardise=False):
 
	if isinstance(angles,list):
		angles=np.array(angles) 
		
	angles=zscore(angles)	
		
	if intercept:	
		# Design matrix
		X=np.column_stack((np.ones((len(angles),1)),angles)) 
	else:
		X=angles[:,np.newaxis] 
 
	idx =intercept*1
	
	# if zscore for y is needed
	if standardise:
		sha=y.shape
		if len(sha)==1:
			y=zscore(y)
		else:
			for k in range(sha[1]):
				y[:,k]=zscore(y[:,k])
	
	invW=np.linalg.pinv(X)
	return np.matmul(invW,y)[idx]  # return the slope coefficient


def is_montage_installed():
	try:
		subprocess.run(["montage", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		return True
	except subprocess.CalledProcessError:
		return False
	
	
def run_parallel_the_pop_decoding(ho,he,ho_la,he_la,Info,nt):
	
	n_cpu = mp.cpu_count()  # Total number of CPU
	pool = mp.Pool(n_cpu) 
	
	time_step_results = [pool.apply_async(pop_decode_at_a_single_timept,args=(ho[:,:,tr],he[:,:,tr],ho_la,he_la,Info)) for tr in range(nt)]   
	pool.close()
	pool.join()

	results = [r.get() for r in time_step_results]   
 
	return np.stack(results,axis=-1)  # neurons X presented stimuls X (mean, std) X time	

 
def pop_decode_at_a_single_timept(ho,he,ho_la,he_la,Info):
	
	"""
	Population decoding of neural data at single time point for all stimulus values
	
	Given a stimulus, what is the response of neurons that prefers different directions
	(population curve construction)
	"""
	
	# initialization for mean and standard devivation
	res_ho_mean=np.zeros((ns,ns))
	res_he_mean=np.zeros_like(res_ho_mean) 
	res_ho_std=np.zeros_like(res_ho_mean)
	res_he_std=np.zeros_like(res_ho_mean)
	
	# In a given trial a stimuls is presented  
	for k in range(1,ns+1):  # Given a stimulus (directions are coded from 1 to 8)
		
		# presented trials
		homo_trials=np.where(ho_la==k)[0]  # cannot be empty
		hetero_trials=np.where(he_la==k)[0] # cannot be empty
	
		# subset the data based on trials
		ho_subset_1=ho[:,homo_trials]
		he_subset_1=he[:,hetero_trials]
		
		# For each tuning neuron (subset the data based on neurons)
		for l in range(1,ns+1): 
			idx_homo=np.where(Info['Pref.Homo']==l)[0]
			idx_hetero=np.where(Info['Pref.Hetero']==l)[0] 
			
			# subsetted data based on the tuned neurons (if they exists)
			if len(idx_homo)!=0:
				ho_subset2=ho_subset_1[idx_homo,:]
				ho_mean=np.mean(ho_subset2.flatten())
				ho_std=np.std(ho_subset2.flatten())
			else: # if neuron group does not exist
				ho_mean=0
				ho_std=0
			if len(idx_hetero)!=0:
				he_subset2=he_subset_1[idx_hetero,:]
				he_mean=np.mean(he_subset2.flatten())
				he_std=np.std(he_subset2.flatten())
			else:
				he_mean=0
				he_std=0
				 
			res_ho_mean[l-1,k-1]=ho_mean
			res_he_mean[l-1,k-1]=he_mean

			res_ho_std[l-1,k-1]=ho_std
			res_he_std[l-1,k-1]=he_std  
			
	#eturn np.stack((res_ho_mean,res_he_mean,res_ho_std,res_he_std),axis=-1) 
	return np.stack((res_ho_mean,res_he_mean),axis=-1)  


def center_all_and_mean(a,shift=5):  

    res=np.zeros_like(a)
    for k in range(a.shape[0]):
        res[:,k]=np.roll(a[:,k],shift-k-1) 
    
    # return the centered and mean of all that centered
    return np.mean(res,1)

def zero_center_the_tcs(data,shift=5): 
    # this function includes the conditions as well 
    ns,ns, nconds, nt=data.shape
    
    res=np.zeros((ns,nconds,nt))
    
    for t in range(nt):
        for c in range(nconds):
            res[:,c,t]=center_all_and_mean(data[:,:,c,t],shift) 
    return res

def get_tuning_curve(data,labels,dur_from=40,dur_to=120): 
    return np.array([np.mean(data[labels==k,dur_from:dur_to]) for k in range(1,ns+1)])
     
def get_preference(data,labels,dur_from=40,dur_to=120):
    
    tc=get_tuning_curve(data,labels,dur_from=dur_from,dur_to=dur_to)
    return np.where(tc==np.max(tc))[0]+1  # index from 1 to 8  
