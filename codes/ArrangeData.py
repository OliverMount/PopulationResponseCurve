## Arranging data for population response curve model 

import os
import scipy
import numpy as np


raw_data_path='/media/olive/Research/oliver/data_down/given_raw';
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
        
        
        for k in range(len(homo)):
            if homo[k][0][0].size:  # if there is a meaningful data
                homo_pval.append(np.squeeze(homo[k][0][0][0]))
                hetero_pval.append(np.squeeze(hetero[k][0][0][0]))   
                
    # Saving the p-values            
    scipy.io.savemat(os.path.join(pval_save_path,paradigm,cond+'.mat'), {'homo' : homo_pval,'hetero' : hetero_pval})
        
     
    #saving the data values
        
        




%################# TASK #####################
load('data_V1_task_45_raw_final.mat');
mat_file_name = "/media/olive/Research/oliver/pvals/task/V1_45.mat";

%load('data_V1_task_90_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/task/V1_90.mat";

%load('data_V1_task_135_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/task/V1_135.mat";

%load('data_PPC_task_45_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/task/PPC_45.mat";

%load('data_PPC_task_90_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/task/PPC_90.mat";

%load('data_PPC_task_135_raw_final.mat'); 
%mat_file_name = "/media/olive/Research/oliver/pvals/task/PPC_135.mat";

%#################### PASSIVE ####################### 
%load('data_V1_passive_45_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/passive/V1_45.mat";

%load('data_V1_passive_90_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/passive/V1_90.mat";

%load('data_V1_passive_135_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/passive/V1_135.mat";

%load('data_PPC_passive_45_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/passive/PPC_45.mat";

%load('data_PPC_passive_90_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/passive/PPC_90.mat";

%load('data_PPC_passive_135_raw_final.mat');
%mat_file_name = "/media/olive/Research/oliver/pvals/passive/PPC_135.mat";

% ############### Remove empty cells  #################
non_empty_cells_homo = data_homo_all(~cellfun('isempty', data_homo_all));
non_empty_cells_hetero = data_hetero_all(~cellfun('isempty', data_hetero_all));

labels_homo=dir_homo_all(~cellfun('isempty', dir_homo_all));
labels_hetero=dir_hetero_all(~cellfun('isempty', dir_hetero_all));


homo=pVal_homo(~cellfun('isempty', pVal_homo));
hetero=pVal_hetero(~cellfun('isempty',pVal_hetero ));


save(mat_file_name, 'homo','hetero');
printf("%s \n" , mat_file_name)

% Save each non-empty cell content as a mat file
for i = 1:numel(non_empty_cells_homo)
    mat_file_name = sprintf('Mouse%d.mat', i);
    sample_data_homo = non_empty_cells_homo{i};
    sample_data_hetero = non_empty_cells_hetero{i};
    dirIdx_homo=labels_homo{i};
    dirIdx_hetero=labels_hetero{i};
    save(mat_file_name, 'sample_data_homo','sample_data_hetero','dirIdx_homo','dirIdx_hetero');
    
end 