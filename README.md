# Action potential model experimental design

Voltage clamp (VC) and current clamp (CC) protocol optimisation for action potential (AP) models.
Consider we have a pool of AP models, and that we will have to pick the one that best fits to the VC or CC recordings by fitting only the conductance values.
Can we design a voltage protocol that can allow us to distinguish between different AP models?


### Requirements

The code requires Python (3.5+) and the following dependencies:
[PINTS](https://github.com/pints-team/pints#installing-pints),
[Myokit](http://myokit.org/install/).

To setup, navigate to the path where you downloaded this repo and run
```console
$ python3 -m venv env
$ source env/bin/activate
$ pip install -r requirements.txt
```


### Outline

1. [design](design): Contains all the optimal experimental design.
2. [study](study): Contains all the analysis of the results in [design](design).


### Supporting files

- [method](method): Contains all the Python helper modules, classes and functions.
- [mmt](./mmt): [Myokit](http://myokit.org/) model files, contains all the AP models.
- [cellml](./cellml): CellML versions of those in [mmt](./mmt).
- [test](test): For debug purposes, can be removed later.
