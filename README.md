# Pathway analysis in metabolomics: Pitfalls and best practice for the use of Over-representation Analysis

Cecilia Wieder <sup>1</sup>, Clément Frainay <sup>3</sup>, Nathalie Poupin <sup>3</sup>, Pablo Rodríguez-Mier <sup>3</sup>,
Florence Vinson <sup>3</sup>, Juliette Cooke <sup>3</sup>, Rachel Lai <sup>2</sup>, Jake Bundy <sup>1</sup>, Fabien Jourdan <sup>3</sup>, Timothy Ebbels <sup>1</sup>

<sup>1</sup> Department of Metabolism, Digestion, and Reproduction, Faculty of Medicine, Imperial College London, London SW7 2AZ, UK

<sup>2</sup> Department of Infectious Disease, Faculty of Medicine, Imperial College London, London SW7 2AZ, UK

<sup>3</sup> INRA, Toulouse University, INP, UMR 1331, Toxalim, Research Centre in Food Toxicology, 180 chemin de Tournefeuille, Toulouse, France


This repository contains the code to run the simulations presented in the study. The Python code to generate the results 
is contained within the Jupyter notebook **src/reproducible_simulations.ipynb**. Users may adapt the code in the notebook to perform the simulations on their own data. 
All code has been tested using Python 3.8. 

<h2>Getting started</h2>
Clone the repository

```
git clone https://github.com/cwieder/metabolomics-ORA.git
```

Install the required packages

```
cd metabolomics-ORA/src
pip3 install -r requirements.txt
```

<h2>Usage</h2>
Launch the reproducible_simulations.ipynb Jupyter notebook and run the code cells

```
cd metabolomics-ORA/src
jupyter-notebook reproducible_simulations.ipynb
```

<h2>License</h2>
MIT

<h2>Contact</h2>
<a href="mailto:cw2019@ic.ac.uk">cw2019@ic.ac.uk</a>
