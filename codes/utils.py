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


def replace_consecutive_above_threshold(array, threshold=0.3):
	result = array.copy()

	for i in range(array.shape[0]):  # Iterate over rows
		j = 0
		while j < array.shape[1]:
			if array[i, j] > threshold:
				# Find the consecutive region above the threshold
				start_idx = j
				while j < array.shape[1] and array[i, j] > threshold:
					j += 1
				end_idx = j - 1

				# Replace the consecutive values with the average of points before and after the region
				prev_value = array[i, start_idx-1] if start_idx > 0 else array[i, start_idx]
				next_value = array[i, end_idx+1] if end_idx < array.shape[1] - 1 else array[i, end_idx]

				avg_value = (prev_value + next_value) / 2

				# Update the consecutive region with the average value
				result[i, start_idx:end_idx+1] = avg_value
			else:
				j += 1

	return result


def replace_consecutive_below_threshold(array, threshold=-0.3):
	result = array.copy()

	for i in range(array.shape[0]):  # Iterate over rows
		j = 0
		while j < array.shape[1]:
			if array[i, j] < threshold:
				# Find the consecutive region below the threshold
				start_idx = j
				while j < array.shape[1] and array[i, j] < threshold:
					j += 1
				end_idx = j - 1

				# Replace the consecutive values with the average of points before and after the region
				prev_value = array[i, start_idx-1] if start_idx > 0 else array[i, start_idx]
				next_value = array[i, end_idx+1] if end_idx < array.shape[1] - 1 else array[i, end_idx]

				avg_value = (prev_value + next_value) / 2

				# Update the consecutive region with the average value
				result[i, start_idx:end_idx+1] = avg_value
			else:
				j += 1

	return result



def replace_above_threshold(array, threshold=0.3):
	result = array.copy()
	for i in range(array.shape[0]):  # Iterate over animals
		for j in range(1, array.shape[1] - 1):  # Iterate over columns, excluding the first and last columns
			if array[i, j] > threshold:
				# Replace the value with the average of the previous and next time points in the row
				result[i, j] = (array[i, j-1] + array[i, j+1]) / 2
	
	return result
	
def replace_below_threshold(array, threshold=-0.3):
	result = array.copy()
	for i in range(array.shape[0]):  # Iterate over animals
		for j in range(1, array.shape[1] - 1):  # Iterate over columns, excluding the first and last columns
			if array[i, j] < threshold:
				# Replace the value with the average of the previous and next time points in the row
				result[i, j] = (array[i, j-1] + array[i, j+1]) / 2
	
	return result



class InvertedEncoding:
	"""
	Creates an object of Inverted encoding model (IEM) based decoder
	""" 
	 
	def __init__(self,s,tc_type='cosine',p=5,kappa=1,ref=4):
		
		"""
		Initializes Forward Encoding model object

		Parameters
		----------
		tc_type: type of average tuning-curve: 'cosine' (default), 'vonMises' (circular Gaussian) and 'wrappedGaussian' Wrapped Guassian
		p: power of tuning function for cosine: 5 (default)
		s: an array of stimulus values; for orientation it is the angles in deg. Internally radian values are computed 
		ref: reference stimuli at which all tuning curves are centered to
	
		"""	

		#print('++ Initializing Inverted Encoding Decoder with ', tc_type, ' average tuning curves')
		
		if isinstance(s,list):
			s=np.array(s)
		
		if tc_type=='cosine':  
			self.params={'tc_type': tc_type,'p':p} 
		elif tc_type=='vonMises':
			self.params={'vonMises': tc_type,'kappa':kappa}
		else:
			raise NotImplementedError(tc_type + ' not implemented yet') 
			
		self.params['Stimulus values']=s
		self.params['sr']= self.params['Stimulus values']*np.pi/180
		self.params['ref']=ref
		self.params['nos']=s.shape[0]
	
		#print(self.params) 
		

	def get_params(self):
		"""
		Returns the parameters of the model set during initialization
		"""
		return self.params

	
	def TuningCurve(self,sc):
		"""
		Method implements three different tuning curves
		"""
		
		if self.params['tc_type']=='cosine':
			self.tuning_curve=pow((0.5 + (0.5 *np.cos(self.params['sr']-sc))),self.params['p'])
		elif self.params['tc_type']=='vonMises':
			self.tuning_curve=(0.0005 + (0.0095 *np.exp(self.params['kappa']*np.cos(self.params['sr']-sc))))
		elif self.params['tc_type']=='wrappedGaussian':
			pass
		else:
			raise NotImplementedError(self.params['tc_type'] + ' not implemented yet') 
			
		return self.tuning_curve  
	
	def plot_a_tuning_curve(self,sc,tc_type='cosine'):
		if tc_type=='cosine':
			curve=pow((0.5 + (0.5 *np.cos(self.params['sr']-sc))),self.params['p'])
		
		plt.plot(range(self.params['nos']),curve,'bo-',linewidth=2)
		plt.xticks(range(self.params['nos']),self.params['Stimulus values'])
		plt.xlabel("Target angle (deg)")
		plt.ylabel("Tuning curve amplitude")
		plt.show()
		
		#Example
		#from Decoders import InvertedEncodingModel
		#model=InvertedEncodingModel(angles)
		#model.plot_a_tuning_curve(75*np.pi/180)
		
	def make_tuning_matrix(self,y):
		"""
		Tuning matrix

		Parameters
		----------
		y : (n X 1) array of floats
			Angle in radians.

		Returns
		-------
		a tuning matrix (no. stimuls X n)

		""" 
		self.base_tuning=self.TuningCurve(self.params['sr'][self.params['ref']-1])
		
		self.C1=np.zeros((len(self.params['sr']),len(y)))
		for k in range(y.shape[0]): 
			self.C1[:,k]=np.roll(self.base_tuning,y[k]-self.params['ref'])

	def fit(self,X,y,sys_type="under"):
		"""
		Fit the model according to the given training data.
		
		Parameters
		----------
		X : {array-like, sparse matrix} of shape (n_samples, n_features)
			Training vector, where `n_samples` is the number of samples and
			`n_features` is the number of features.
		y : array-like of shape (n_samples,)
			Target vector relative to X. 
		
		Returns
		-------
		self
			Fitted estimator.
		Notes
		-----
		Supports numpy array
		""" 

		# Check that X and y have correct shape 
		#X, y = check_X_y(X.T, y)
		#print(type(y)) 
		
		self.sys_type=sys_type
		
		# Encoding matrix
		self.make_tuning_matrix(y)   # self.C1 is made here and stored in the object
		
		
		# Weight estimator
		if self.sys_type=="under":
			self.esti_w=np.matmul(X,np.linalg.pinv(self.C1))
			#print("Min. Weight:", np.min(self.esti_w),"Max. Weight:", np.min(self.esti_w))
			#fitted values 
			self.fitted_values=np.matmul(self.esti_w,self.C1)
			# residuals of the fit
			self.fit_residuals=X-self.fitted_values
		elif self.sys_type=="over":
			X=X.T
			self.C1=self.C1.T
			self.esti_w=np.matmul(np.linalg.pinv(self.C1),X)
			
			#fitted values 
			self.fitted_values=np.matmul(self.C1,self.esti_w)
			# residuals of the fit
			self.fit_residuals=X-self.fitted_values
		else:
			raise NotImplementedError(self.sys_type + ' not implemented yet') 
		

	def predict(self,X):
		
		# Check if model is fitted first
		check_is_fitted(self,["esti_w"]) 
		#if check_is_fitted(self,["esti_w"]) is None:	   
		#	print('++ Model is fitted already. Predicting the labels for test data') 
		
		if self.sys_type=="under":  # if underdetermined solver is desired 
			s=X.shape
			if len(s) ==1 or len(s)==2:  # 1D or 2D test data
				self.predicted_values = np.matmul(np.linalg.pinv(self.esti_w),X)
				#self.predict_residuals = X-np.matmul(self.esti_w,self.predicted_values) 
			elif len(s)==3:		   # 3D test data; rep 2d testing along third axis
				#pred=np.zeros((self.params['nos'],s[1],s[2]))
				#for k in len(s[2]):
				#	pred[:,:,k] =np.matmul(np.linalg.pinv(self.esti_w),X[:,:,k])
				#self.predicted_values = pred
				invW= np.linalg.pinv(self.esti_w)  # compute inverse only once
				self.predicted_values =np.stack([np.matmul(invW,X[:,:,k]) for k in range(s[2])],-1)
				#del pred # Allocated different memory 
			else:
				raise NotImplementedError('++ Underdetermined solver not implemented yet for more than 3D test data sets') 
				
		elif self.sys_type=="over": # if over-determined solver is desired
			self.esti_w=self.esti_w.T 
			self.predicted_values = np.matmul(np.linalg.pinv(self.esti_w),X)
			#self.predict_residuals = X-np.matmul(self.esti_w,self.predicted_values) 
		else:
			raise NotImplementedError(self.sys_type + ' not implemented yet')

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


def roll_all(X,y,ref): 
	
	s=X.shape
	Ls=X.ndim  # dimension of the data X   
	
	if Ls==1:	#single-trial test set (1D) # This is typical of CV=1 in diagonal decoding 
		res=np.roll(X,ref-y)
	elif Ls==2:  # Muliple trials test set (2D); This is typical of CV=1 in CTG decoding or CV>1 in diagonal decoding
		if s[1]!=len(y):	  # Make sure the no. of cols of X and size of y are the same
			raise ValueError('-- No. of test samples and no. of ground truth labels do not match')
		else:
			res=np.zeros_like(X)
			for k in range(s[1]):  # s[1] can be trial or time axis; if time then all element of y will the same
				#res[:,k]=roll_tuning_curve(X[:,k],y[k],ref)
				res[:,k]=np.roll(X[:,k],ref-y[k])
			 
	elif Ls==3: # Multiple trials test set along time (3D); This is typical of CV>1; in CTG decoding
		if s[1]!=len(y):	  # Trial dimension of X (2nd dim) and trial dimention of y (1st dim) should be the same
			raise ValueError('-- No. of test samples and no. of ground truth labels do not match')
		else:
			res=np.zeros_like(X)
			for t in range(s[2]): # Along time
				for k in range(s[1]): # Along trials
					#res[:,k,t]=roll_tuning_curve(X[:,k,t],y[k,t],ref)
					#print("ref-" ,ref-y[k]) 
					res[:,k,t]=np.roll(X[:,k,t],ref-y[k]) # irrepective of time t and same label 
	else:  # not implemented for more than three dimensions
		raise NotImplementedError('-- roll_all is not implemented for more than 3D :-( ')
		
	return res 
		
def run_parallel_the_time(D1,D2,d1,d2,nt):
 
	n_cpu = mp.cpu_count()  # Total number of CPU
	pool = mp.Pool(n_cpu) 
	
	time_step_results = [pool.apply_async(process_time_step,args=(D1[:,:,tr],D2[:,:,tr],d1,d2)) for tr in range(nt)]   
	pool.close()
	pool.join()

	results = [r.get() for r in time_step_results]   
 
	return np.stack(results,axis=-1).transpose(0,2,1)   # no X time X (homo,hetero)

def process_time_step(Xtra,Xte,ytra,yte):

	# LOO for homogeneous 
	cv=LeaveOneOut()
	cv_results = []
	
	for train, test in cv.split(Xtra.T):  # cross-validation split over trials
		cv_results.append(process_cv(Xtra[:, train], ytra[train], Xtra[:, test], ytra[test]))  
 
	final_result = np.stack(cv_results, axis=-1)   
	cv_res_final=np.mean(final_result,2)   # mean across trials of zero-centered) tuning curves (homo case)
 
	# Homo train- Hetero test 
	model=InvertedEncoding(angles, p=6, ref=center_around) 
	model.fit(Xtra,ytra,sys_type="under") 
 
	model.predict(Xte)   
	predicted=model.predicted_values

	# centering the tuning curve based on the presented orientation
	ypred=roll_all(predicted,yte,center_around)
	ypred=np.mean(ypred,1) # mean across trials of zero-centered) tuning curvesv (hetero case)   
	
	#print(cv_res_final.shape)
	#print(ypred[:,np.newaxis].shape) 
	res=np.concatenate((cv_res_final,ypred[:,np.newaxis]),axis=-1)  
	return res 
 

def process_cv(Xtrain,ytrain,Xtest,ytest): 
 
	model=InvertedEncoding(angles,p=6, ref=center_around)   
	model.fit(Xtrain, ytrain,sys_type="under")  
	
	model.predict(Xtest)	
	predicted=model.predicted_values 
 
	y_pred=roll_all(predicted,ytest,center_around) 
 
	return y_pred  

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