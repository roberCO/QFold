<div align="center">    
 
# QFold     

[![arXiv](http://img.shields.io/badge/arXiv-2101.10279-B31B1B.svg)](https://arxiv.org/abs/2101.10279)
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
See also configuration for the installation of `psi4`

## Configuration
It is possible to configure different parameters during QFold execution. In the config/config.json it is possible to modify the value of the parameters.

The most important parameter is the psi4 library path. The variable "psi4_path" containts a path where the psi4 binary file is stored.

For example: /home/user/installations/psi4conda/bin/psi4. The binary execution file of psi4 can be downloaded from https://psicode.org/installs/v15/

This repository should work with the following qiskit versions
```
qiskit                    0.29.0                   pypi_0    pypi
qiskit-aer                0.8.2                    pypi_0    pypi
qiskit-aqua               0.9.4                    pypi_0    pypi
qiskit-ibmq-provider      0.16.0                   pypi_0    pypi
qiskit-ignis              0.6.0                    pypi_0    pypi
qiskit-terra              0.18.1                   pypi_0    pypi
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

### Citation   
```
@article{casares2022qfold,
  title={QFold: quantum walks and deep learning to solve protein folding},
  author={Casares, Pablo Antonio Moreno and Campos, Roberto and Martin-Delgado, Miguel Angel},
  journal={Quantum Science and Technology},
  volume={7},
  number={2},
  pages={025013},
  year={2022},
  publisher={IOP Publishing}
}
```   
