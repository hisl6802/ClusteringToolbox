# Clustering toolbox

## Table of Contents
- [Design](#design)
- [Ensemble Clustering combined with Clustering Optimization (ECCO)](#ensemble-clustering-combined-with-clustering-optimization-ecco)
- [Data Pre-processing notes](#data-pre-processing-notes)
- [Installation and set-up](#installation-and-set-up)
  - [After Installation of python](#after-installation-of-python)
- [Example Files](#example-files)
- [Troubleshooting](#troubleshooting)

## Design
The Clustering Toolbox available here was designed as a user interface that allows for clustering analysis of untargeted metabolomics data. This user interface contains 23 functionalities spanning from hierarchical clustering to clustering comparision to KEGG pathways analysis. Most functionalities within this UI contain user options allowing for customization of analyses. Furthermore, we provide a start-up window that allows the user to define the number of threads (e.g., processes they would like to use during analysis allowing for rapid clustering analysis). We recommend the user use at least one less thread than available on the computer, with the optimal being half of the available threads. 

## Ensemble Clustering combined with Clustering Optimization (ECCO)
Ensemble clustering combined with Clustering Optimization or ECCO is the main functionality within the clustering toolbox. This form of ensemble clustering is discussed in our pre-print [pre-print](https://www.biorxiv.org/content/10.1101/2022.11.03.515009v1.abstract) and can be adapted for your preferences. ECCO is built for ensemble clustering solutions of agglomerative hierarcical clustering algorithms. 

## Data Pre-processing notes
The Clustering Toolbox was built to mimic the pre-processing steps taken during untargeted metabolomics data analysis. We acknowledge that we do not provide all of the available pre-processing steps and recommend pre-processing prior to submission to UI and selecting 'None' and 'None' when prompted to select a transformation and scaling for your data. 

## Installation and set-up
1. **Unzip Release**
   - `unzip ClusteringToolbox_latest.zip`
  
3. **For Mac/Unix Users**: Run the shell script in your terminal
   - `chmod +x RunMeFirst.sh`
   - `./RunMeFirst.sh`

## Example Files
The ExampleFiles directory contains **two example input files**. Please see further documentation for other input files. 

## Troubleshooting
Please submit an Issue in the Issues tab, and I will address as quickly as possible.
