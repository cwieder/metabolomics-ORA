# Pathway analysis in metabolomics: Pitfalls and best practice for the use of Over-representation Analysis

Cecilia Wieder <sup>1</sup>, Clément Frainay <sup>3</sup>, Nathalie Poupin <sup>3</sup>, Pablo Rodríguez-Mier <sup>3</sup>,
Florence Vinson <sup>3</sup>, Juliette Cooke <sup>3</sup>, Rachel PJ Lai <sup>2</sup>, Jacob G Bundy <sup>1</sup>, Fabien Jourdan <sup>3</sup>, Timothy Ebbels <sup>1</sup>

<sup>1</sup> Department of Metabolism, Digestion, and Reproduction, Faculty of Medicine, Imperial College London, London SW7 2AZ, UK

<sup>2</sup> Department of Infectious Disease, Faculty of Medicine, Imperial College London, London SW7 2AZ, UK

<sup>3</sup> INRA, Toulouse University, INP, UMR 1331, Toxalim, Research Centre in Food Toxicology, 180 chemin de Tournefeuille, Toulouse, France


This repository contains the code to run the simulations presented in the study. The Python code to generate the results 
is contained within the Jupyter notebook **src/reproducible_simulations.ipynb**. Users may adapt the code in the notebook to perform the simulations on their own data. 
All code has been tested using Python 3.8 on MacOS (v11.2.3) with standard hardware (16GB RAM). 

<h2>Getting started: local installation</h2>
Clone the repository

```
git clone https://github.com/cwieder/metabolomics-ORA.git
```

Install the required packages

```
cd metabolomics-ORA/src
pip3 install -r requirements.txt
```
Cloning the repository and installing the dependencies should take less than 10 minutes on a standard desktop computer. 
<h2>Usage</h2>
Launch the reproducible_simulations.ipynb Jupyter notebook and run the code cells

```
cd metabolomics-ORA/src
jupyter-notebook reproducible_simulations.ipynb
```
<h2>Run in any browser using Google Colaboratory</h2>
As an alternative to local installation, the Jupyter notebook can now be run 
as an analogous Google <a href="https://research.google.com/colaboratory/faq.html"> Colab </a> version. This does not require any local installation of code, packages, or dependencies, and all
simulations are run in a browser window. All code is run on one of Google's virtual machines and 
therefore does not require access to the user's hardware. 


To get started, open the <a href="https://colab.research.google.com/drive/1Ga_PasVyIXOQYwrlYZGEU9iNjqFnOrRh">Colab notebook</a> 
and run the cells.

<h2>License</h2>
MIT

<h2>Contact</h2>
<a href="mailto:cw2019@ic.ac.uk">cw2019@ic.ac.uk</a>
