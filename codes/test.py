# Testing the shapes of old and new data sets

import numpy as np
import scipy 
import os



paradigms=['task','passive']
conds=['V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135']
 
# Custom sorting key
def extract_number(animal):
	return int(animal[5:-4]) # Assuming the format is 'animalX'


data_path='/media/olive/Research/oliver/data_down_p/'


for paradigm in paradigms: 
    for cond in conds:  
        print( " Condition :  ",cond)
        os.chdir(os.path.join(data_path,paradigm,cond))
        
        ani_list=sorted(os.listdir(),key=extract_number)
        #print(ani_list)
        noa=len(ani_list)
    	
        #print('No. of animals in ' + roi + ' is ' + str(noa)) 
    	
        for idx,k in enumerate(ani_list):
            A=scipy.io.loadmat(k) 

            homo_data=A['sample_data_homo'] 
            hetero_data=A['sample_data_hetero'] 

            print(k + " Homo shape: ",homo_data.shape)