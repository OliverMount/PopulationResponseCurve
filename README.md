# Population Tuning Curve

###  Population Tuning Curve based Decoding of Mean Direction using Mouse Calcium Imaging Data

This repository contains codes Python codes for obtaining population tuning curves from mouse V1 and PPC calcium imaging data (preprocessed). This repository does not include 

### Instruction for running the code

1. git clone the repository to a local folder in your computer

2. cd in to the scripts folder and set the following paths in the  decoding.py
	a. data_path: This is the path where the preprocessed calcium data are stored
	b. pval_path: pvalues of the tuned and untuned neurons are stored here

	The two levels of the data directory will look like this

	```

	├── passive
	│   ├── PPC_135
	│   ├── PPC_45
	│   ├── PPC_90
	│   ├── V1_135
	│   ├── V1_45
	│   └── V1_90
	├── pvals
	│   ├── task_PPC_135.mat
	│   ├── task_PPC_45.mat
	│   ├── task_PPC_90.mat
	│   ├── task_V1_135.mat
	│   ├── task_V1_45.mat
	│   └── task_V1_90.mat
	└── task
	    ├── PPC_135
	    ├── PPC_45
	    ├── PPC_90
	    ├── V1_135
	    ├── V1_45
	    └── V1_90

	```

3. Run the pop_decoding.py it would create a decoding folder (this is the results folder) in the level as of the scripts folder. The structure  of the decoding folder would be the following

```
.
├── plots
│   ├── montages
│   ├── task_V1_45_10.png
│   ├── task_V1_45_20.png
│   ├──        . 
│   ├──        .  
│   └── task_PPC_135_100.png
├── slopes
│   └── task
│   └── passive
└── tuning_curves
    ├── passive
    └── task
```

Note: Individual figures for each condition is stored inside the plots directory. Montage folder (that cotains sewed figure files) will be recrated only if the platform is Lixux (if montage is installed). 
