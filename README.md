# MFPBehavioralAnalysis
This repository contains the open-source code accompanying our research paper, "Identifying behavioral links to neural dynamics of multifiber photometry recordings in a mouse social behavior network." Yibo Chen, Jonathan Chien, Bing Dai, Dayu Lin, Zhe Sage Chen

bioRxiv link：https://www.biorxiv.org/content/10.1101/2023.12.25.573308v1

In this research work, we utilized two models, namely HSMM and PSID. The corresponding code can be found in folders named 'HSMM' and 'PSID' respectively. We have only provided example code; if you wish to access code for other processing methods described in the paper, please feel free to contact the authors.

## Install
### HSMM
Dependencies:
make sure you have following packages installed:
+1. Anaconda2 or Miniconda2
    note: for Miniconda2, these packages are needed:
            pip install numpy scipy matplotlib cython nose (JMC also add: future, requests; also, may be preferrable use "conda install")
+2. gcc-4.8 g++-4.8
    sudo apt install gcc-4.8 g++-4.8 (JMC: this installation command no longer works, see below; also note that gcc 9.4 should be old enough)
+3. install pybasicbayes
    pip install pybasicbayes

Then you need to get three package dependencies: "hips-lib", "pyhsmm", and "pyhsmm_spiketrains".
https://github.com/jonathan-chien/hips-lib
https://github.com/mattjj/pyhsmm
https://github.com/slinderman/pyhsmm_spiketrains

### PSID
This part is entirely Matlab programs. Please add 'Source' to the Matlab path. 

The programs in 'Source' are referenced from https://github.com/ShanechiLab/PSID

## Usage

### HSMM
There is an intruction in folder 'HSMM'.

### PSID
+PSID_AcrossSession.m is used to show the prediction results of all sessions of a animal
+male_tracjory is used to plot neural trajectory of "Attack" behavior
+female_tracjory is used to plot neural trajectory of "Mount" and "Thrust" behavior

## Data
Below is the example data link.
https://www.dropbox.com/t/oRC52RkhYzvuX0oS
