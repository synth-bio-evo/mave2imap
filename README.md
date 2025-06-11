mave2imap
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/mave2imap/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/mave2imap/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/mave2imap/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/mave2imap/branch/main)

### Description  
This code is intended for 3D mapping of interface hotspots based on results from MAVE (Multiplexed Assays of Variant Effects).


### Install (Linux)  
$ conda env create -f https://github.com/synth-bio-evo/maveimap_test.git
(https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/environment.yml) 


### Testing  
#### *1) :construction: Create a folder for testing and download testing files*    
$ mkdir /tmp/test  
$ cd /tmp/test  


#### *2) :computer: Run mave2imap pipeline for each targeted region.*   
Exemple:  
$ cd Asf1B+IP3/Asf1_N-Ter  
$ mave2imap -i Asf1_N-ter.ini
$ cd ../Asf1_C-Ter
$ mave2imap -i Asf1_C-ter.ini  

 This will produce the data required for analysis and visualization using the proposed jupyter notebook.   

:microscope: The information copiled in the file, "result_thresh3_2_2_compare_conditions.out", is probably the most relevant to a classical user.


#### *3)  :mag_right: Analyze results using jupyter notebook(s).*  
- enter appropriate folder and launch jupyter-lab  

$ cd ../analysis  
$ jupyter-lab  
- Choose mave2imap kernel  
- Click in "Run" (menu) => "Restart Kernel and Run All Cells"  

*The most perturbed positions should be indicated  below the last cell based on the defined threshold and you should be able to visualized/manipulated the 3D interactive complex (most perturbed regions are indicated by red gradient color)*



### Copyright

Copyright (c) 2025, Raphaël Guérois (CEA-Saclay, DRF/Joliot/I2BC_saclay/SB2SM/LBSR), Oscar H.P. Ramos (CEA-Saclay, DRF/Joliot/MTS/SIMoS/LICB/SBE)


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.11.
