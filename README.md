<div align="center">    
 
# QFold     

[![arXiv](http://img.shields.io/badge/arXiv-2101.10279-B31B1B.svg)](https://arxiv.org/pdf/2101.10279.pdf)
[![Journal](http://img.shields.io/badge/Quantum_Science_and_Technology-2022-4b44ce.svg)](https://iopscience.iop.org/article/10.1088/2058-9565/ac4f2f)
<!--
[![Conference](http://img.shields.io/badge/ICLR-2019-4b44ce.svg)](https://papers.nips.cc/book/advances-in-neural-information-processing-systems-31-2018)
[![Conference](http://img.shields.io/badge/AnyConference-year-4b44ce.svg)](https://papers.nips.cc/book/advances-in-neural-information-processing-systems-31-2018)  

ARXIV   
[![Paper](http://img.shields.io/badge/arxiv-quant.ph:arXiv:2101.10279-B31B1B.svg)](https://arxiv.org/pdf/2101.10279.pdf)
-->


<!--  
Conference   
-->   
</div>
 
## Description   
Software to solve the folding protein problem using Quantum Computing and Machine Learning  

## Installation  
First, install dependencies   
```bash
# clone project   
git clone https://github.com/roberCO/QFold/
```

Then we create a conda environment
```
conda create -n qfold python=3.6
conda activate qfold
```

Then we must install a few packages
```
pip install numpy
pip install scipy
pip install tensorflow
pip install keras
pip install qiskit
pip install matplotlib
pip install bokeh
pip install functools
pip install progressbar

```

## How to run
Next, run it.   
```
python main.py [peptide_name] [# aminoacids] [# rotation bits] [initialization: random/minifold] [mode: simulation/experiment]

```
Example:
```
python3 main.py glycylglycine GG 2 minifold simulation
```
## Configuration
It is possible to configure different parameters during QFold execution. In the config/config.json it is possible to modify the value of the parameters.

The most important parameter is the psi4 library path. The variable "psi4_path" containts a path where the psi4 binary file is stored.

For example: /home/user/installations/psi4conda/bin/psi4. The binary execution file of psi4 can be downloaded from https://psicode.org/installs/v15/


### Citation   
```
@article{QFold,
  title={QFold: Quantum Walks and Deep Learning to Solve Protein Folding},
  author={Casares, PAM and Campos, Roberto and Martin-Delgado, MA},
  journal={arXiv preprint arXiv:2101.10279},
  year={2021}
}
```   
