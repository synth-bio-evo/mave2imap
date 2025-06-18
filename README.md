# mave2imap
![mave2imap!](mave2imap.png "mave2imap : graphical abstract")
==============================
### Table of contents

- [mave2imap](#mave2imap)
- [](#)
    - [Table of contents](#table-of-contents)
    - [Description](#description)
    - [Install mave2imap conda environment](#install-mave2imap-conda-environment)
    - [Testing](#testing)
      - [*1) Create a folder to download required data and run the test* :construction:](#1-create-a-folder-to-download-required-data-and-run-the-test-construction)
      - [*2) Run mave2imap pipeline for each targeted region.* :computer:](#2-run-mave2imap-pipeline-for-each-targeted-region-computer)
      - [*3)  Analyze results using jupyter notebook(s).* :mag\_right:](#3--analyze-results-using-jupyter-notebooks-mag_right)
    - [Citing mave2imap](#citing-mave2imap)
    - [Copyright](#copyright)
      - [Acknowledgements](#acknowledgements)

---
### Description  
This code is intended for 3D mapping of interface hotspots based on the most perturbed positions inferred from MAVE (Multiplexed Assays of Variant Effects) results. ([See publication](#citing-mave2imap))  

  
---  
### Install mave2imap conda environment  
$ conda env create -f https://raw.githubusercontent.com/synth-bio-evo/mave2imap/main/mave2imap.yml  

---
### Testing  
*Requires about >= 64 Gb RAM to process the full dataset.  
If you do not dispose of this amount of RAM you can create smaller .fastq.gz files by using the following command:*  

>gunzip -cd \<file>.fastq.gz | head -n 1600000 | gzip > <file_400k_reads>.fastq.gz  
 
- *Replace "\<file>" by your filename*  
- *It will extract and compress 1,6x10⁶ lines from "\<file>.fastq.gz", corresponding to 4x10⁵ reads, and create  "<file_400k_reads>.fastq.gz"*  

<br>  

#### *1) Create a folder to download required data and run the test* :construction:    
>mkdir /tmp/test  
>cd /tmp/test  

Data, notebooks used to produce the manuscript figures, results, and datasets are available at Zenodo (doi: 10.5281/zenodo.15690360; https://zenodo.org/records/15690360). Below you will find an exemple about how to use them.  

If you have aria2c installed (faster):  

>aria2c -j 16 https://zenodo.org/records/15690361/files/Asf1B+IP3.tar.gz?download=1  

else,  

>wget https://zenodo.org/records/15690361/files/Asf1B+IP3.tar.gz?download=1

Uncompress the .tar.gz file 

>tar -xvzf Asf1B+IP3.tar.gz  

<br>  

#### *2) Run mave2imap pipeline for each targeted region.* :computer:   
Example:  
>cd Asf1B+IP3/Asf1_N-Ter  
>mave2imap -i Asf1_N-ter.ini  
>cd ../Asf1_C-Ter  
>mave2imap -i Asf1_C-ter.ini  

 This will produce the data required for analysis and visualization using the proposed jupyter notebook.   

:microscope: *The information available in the output file, "result_thresh3_2_2_compare_conditions.out", is probably the most relevant to a classical user.*  <br>  <br>

#### *3)  Analyze results using jupyter notebook(s).* :mag_right:   

* Enter main folder and launch jupyter-lab  
> cd /tmp/test/Asf1B+IP3  
> jupyter-lab
<br>

Interface mapping (imap) and fitness assessment notebooks are available in the corresponding folders  <br> 
- Open the notebook
- Choose mave2imap kernel  
- If required edit the code according to your specific case (not required for the testing dataset)  
- Click in "Run" (menu) => "Restart Kernel and Run All Cells"  
  - If you are running the imap notebook, the most perturbed positions should be indicated below the last cell based on the defined threshold, and you should be able to visualized/manipulated the 3D interactive complex (most perturbed regions are indicated by reddish gradient). 
  - If you are running the fitness notebook, by mutation results will be outputted to a CSV file.

---
### Citing mave2imap 
"Publication is coming ..."

<br>
<br> </br>  

---


  

### Copyright

Copyright (c) 2025, Raphaël Guérois (CEA-Saclay/DRF/Joliot/I2BC/SB2SM/LBSR), Oscar H.P. Ramos (CEA-Saclay/DRF/Joliot/MTS/SIMoS/LICB/SBE)


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.11.
