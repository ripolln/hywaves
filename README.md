# HyWaves

A hybrid metamodel (numerical-statistical) approach to optimize complex waves calculations.

## Table of contents
1. [Description](#desc)
2. [Main Contents](#mc)
3. [Documentation](#doc)
4. [Schemes](#sch)
5. [Install](#ins)
    1. [Install from sources](#ins_src)
    2. [Install SWAN numerical model](#ins_swn)
6. [Examples](#exp)
7. [Contributors](#ctr)
8. [License](#lic)


<a name="desc"></a>
## Description

Solving open source waves numerical models executions could require considerable computational resources and time.
This could be a problem when trying to solve a big dataset of numerical cases.

This toolbox includes a hybrid (metamodel) methodology to optimize this process:

- Using MaxDiss classification a small but representative subset of the input dataset is selected.
- Numerical cases are solved with the customized python wrappers given with the source code.
- a Radial Basis Function interpolator is feed with the selected subset and its corresponding Numerical cases output.
- Now this Radial Basis function interpolator can be used to solve the entire dataset with almost no computer resources.
 

<a name="mc"></a>
## Main contents

SWAN Numerical model wrapper

hywaves include our SWAN numerical model python wrapper, detailed information, documentation and examples can be found at <https://gitlab.com/geoocean/bluemath/numerical-models-wrappers/wrap_swan>

[wswan](./hywaves/wswan/): SWAN numerical model toolbox 
- [io](./hywaves/wswan/io.py): SWAN numerical model input/output operations
- [wrap](./hywaves/wswan/wrap.py): SWAN numerical model python wrap 
- [geo](./hywaves/wswan/geo.py): azimuth distance function
- [storms](./hywaves/wswan/storms.py): storm parameters function 
- [stopmotion](./hywaves/wswan/stopmotion.py): stopmotion library 
- [vortex](./hywaves/wswan/vortex.py): vortex TC model 
- [plots](./hywaves/wswan/plots/): plotting module 

Statistical modules 

- [mda](./hywaves/statistical/mda.py): MaxDiss classification module 
- [rbf](./hywaves/statistical/rbf.py): Radial Basis Function module
- [plots](./hywaves/statistical/plots.py): statistical related plots 


<a name="doc"></a>
## Documentation

Camus et al. (2011). A hybrid efficient method to downscale wave climate to coastal areas. Coastal Engineering 58 (2011), 851-862. <http://doi.org/10.1016/j.coastaleng.2011.05.007>

SWAN numerical model detailed documentation can be found at: <http://swanmodel.sourceforge.net/> 

- [SWAN install/compile manual](http://swanmodel.sourceforge.net/download/download.htm)
- [SWAN user manual](http://swanmodel.sourceforge.net/online_doc/swanuse/)


<a name="sch"></a>
## Schemes

Statistical Metamodel
![picture](docs/img/metamodel.svg)

HyWaves - SWAN STATIONARY Methodology
![picture](docs/img/mdaswanrbf.svg)


<a name="ins"></a>
## Install

Source code is currently privately hosted on GitLab at:  <https://gitlab.com/geoocean/bluemath/hybrid-models/hywaves/tree/master> 


<a name="ins_src"></a>
### Install from sources

Use Python3.7 for full compatibility.

Navigate to the base root of [hywaves](./)

Using a Python virtual environment is recommended

```sh
# install virtualenv package 
python3.7 -m pip install virtualenv

# create a new virtual environment for installation
python3.7 -m virtualenv venv

# now activate the virtual environment
source venv/bin/activate
```

Install requirements.

```bash 
pip install -r requirements.txt
```

Then install hywaves:

```bash
python setup.py install
```

<a name="ins_swn"></a>
### Install SWAN numerical model 

Download and Compile SWAN numerical model:

```bash
# you may need to install a fortran compiler
sudo apt install gfortran

# download and unpack
wget http://swanmodel.sourceforge.net/download/zip/swan4131.tar.gz
tar -zxvf swan4131.tar.gz

# compile numerical model
cd swan4131/
make config
make ser
```

Copy SWAN binary file to module resources

```bash
# Launch a python interpreter
$ python

Python 3.6.9 (default, Apr 18 2020, 01:56:04) 
[GCC 8.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
  
>>> from hywaves import wswan
>>> wswan.set_swan_binary_file('swan.exe')
```


<a name="exp"></a>
## Examples:

### HySWAN metamodel: MDA -> SWAN -> RBF

- [notebook - HyWaves](./notebooks/01_MDA_STATIONARY_RBF.ipynb): Waves dataset propagation to a point using MDA classfication, SWAN-STATIONARY numerical simulation and RBF reconstruction. 


<a name="ctr"></a>
## Contributors:

Nicolas Ripoll Cabarga (ripolln@unican.es)\
Alba Ricondo Cueva (ricondoa@unican.es)\
Sara Ortega Van Vloten (sara.ortegav@unican.es)\
Fernando Mendez Incera (fernando.mendez@unican.es)


<a name="lic"></a>
## License

This project is licensed under the MIT License - see the [license](./LICENSE.txt) file for details

