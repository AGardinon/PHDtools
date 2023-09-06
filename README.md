# PHDtools

PHDtools is a Python library that contains a collection of generic tools that I developed to aid my computational work during my PhD.

__Disclaimer__ : the package is generic and experimental, and contains tools that were useful to me during my PhD.  
Handle with care!  


## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install **PHDtools**.

```bash
pip install phdtools
```


## Structure

__Version__ : 1.0.1

Package oranization:
------------

    ├── README.md               <- README file of the project.    
    |
    ├── example
    |   ├── data                <- data used in the tutorial examples.  
    |   └── notebook            <- tutorial python notebooks.  
    |
    ├── phdtools
    |   |
    |   ├── ASEtools            <- implementation of python tools 
    |   |                          based on the ase.atoms library
    |   |
    |   ├── computes            <- generic tools for scientific data analysis
    |   |                          and MD simulations. 
    |   |
    |   ├── descriptors         <- implementation of python tools 
    |   |   |                      based on the quip descriptor library.
    |   |   └── quip
    |   |
    |   └── plots               <- generic tools for aided python image plotting.
    |
    └── setup.py                <- package manager installer file.
------------


## Roadmap

- [x] Add MD simulation tools
    - [x] trajectory unwrapper (v1.0.1)
- [x] Add [ase.atoms](https://wiki.fysik.dtu.dk/ase/index.html) support (v1.0.1)
- [x] Add [quip](https://libatoms.github.io/GAP/index.html) SOAP/turbo-SOAP descriptors support (v1.0.1)
- [x] Add FES computes and plot support (v1.0.1)
- [ ] Add project manager for ASEtools
- [ ] Add MD simulation tools
    - [ ] RDF
    - [ ] MDAnalysis tools
- [ ] Add kinetic models / MSM model support
- [ ] Add ML tools


## Support

Any kind of feedback is welcome, just be patient.


## Authors and acknowledgment

All the materials have been developed by me (this is probably the reason why is not perfect), nonetheless I grateful for all the usefull discussion I had suring my PhD, which inspired me to create this repository.  
In particular I would like to acknowledge:  
- Dr. Riccardo Capelli
- Dr. Giovanni Doni
- Dr. Ioan-Bogdan Magdau

