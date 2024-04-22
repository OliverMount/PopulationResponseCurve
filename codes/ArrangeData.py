## Arranging data for population response curve model 

import os
import scipy
import numpy as np


raw_data_path='/media/olive/Research/oliver/data_down/given_base';
data_save_path='/media/olive/Research/oliver/data_down/';
pval_save_path='/media/olive/Research/oliver/pvals/';
paradigms=['task','passive']
conds=['V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135']
 
 

for paradigm in paradigms:
    os.chdir(os.path.join(raw_data_path,paradigm))
    for cond in conds:
        
        homo_pval=[]
        hetero_pval=[] 
        
        A=scipy.io.loadmat(cond+'.mat')
        
        
        # deal with the p-value
        homo= A['pVal_homo']
        hetero= A['pVal_hetero']
        
        homo_data=A['data_homo_all']
        hetero_data=A['data_hetero_all']
        
        homo_labels=A['dir_homo_all']
        hetero_labels=A['dir_hetero_all']
        
        
        count=1
        for k in range(len(homo)):
            
            if homo[k][0][0].size:  # if there is a meaningful data
            
                mouse_name='Mouse'+str(count)+'.mat'
                
                # Removing the empty list p-values
                homo_pval.append(np.squeeze(homo[k][0][0][0]))
                hetero_pval.append(np.squeeze(hetero[k][0][0][0]))  
                
                
                
                
                # Removing the empty data sets
                #saving the data values
                scipy.io.savemat(os.path.join(data_save_path,paradigm,cond,mouse_name), {'sample_data_homo' : homo_data[k][0],'sample_data_hetero' : hetero_data[k][0],'dirIdx_homo' : homo_labels[k][0][0] ,'dirIdx_hetero': hetero_labels[k][0][0] })
                
                count=count+1
                
                
        # Saving the p-values            
        scipy.io.savemat(os.path.join(pval_save_path,paradigm,cond+'.mat'), { 'homo': homo_pval,'hetero' : hetero_pval})
 