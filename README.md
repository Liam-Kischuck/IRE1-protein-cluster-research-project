# IRE1 Protein Cluster Dynamics

This project simulated and analysed the dynamics of IRE1 proteins clusters on the tubular geometry of the endoplasmic reticulum. There are still several other sections of the project, which will be posted soon!

## Description

The results from this project are published in the journal Soft Matter, and can be found here: https://pubs.rsc.org/en/content/articlelanding/2023/sm/d3sm00694h/unauth

### Background

IRE1 proteins reside on the tubular surface of the endoplasmic reticulum inside biological cells. The endoplasmic reticulum is a network of tubes within the cell, and plays a pivital role in the proper delivery of proteins within the cell to their correct location. 
The IRE1 proteins are able to detect and cluster together when the endoplasmic reticulum is overloaded with proteins, so IRE1 act as a protein that can signal increased stress levels. 

In this project, we used simulations to model the IRE1 protein clustering, and found that our simulated clusters can undergo interesting conformational changes that can impact how quickly the IRE1 clusters can grow. We also found that the cluster shape can have a effect on the long-term stability of the IRE1 clusters.

## How to use and build your own simulations

* The code that is used for simulating the IRE1 clusters can be found in the MATLAB files in each section. These are the .m files that typically begin with a the words 'KMC_'. Each MATLAB file is adjusted for the particular simulations that we carried out in that section.
* To get an adequate amount of samples of data for the project, many 100s of simulations had to be performed. For those who would like to replicate these results, the PowerShell files were used to copy the MATLAB files that were then run on a computer cluster. 
* The analysis of the data were done using Python. Each project section has several python files that were used for the data analysis.

This project is closed to contributions, but feel free to download an try building more sophisticated versions of IRE1 clustering, or general protein clustering simulations on the endoplasmic reticulum!
There is a lot of research out there that is starting to look at the exciting properties of proteins on intracellular organelles.
