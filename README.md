# Code to analyze Nodal Seismic Data collected along the Parks Highway in Feb-March of 2019
## University of Alaska Fairbanks

## Data

Between Febuary 11th and March 26th of 2019 a set of 303 Fairfield Nodal 3C 5Hz sensors were deployed along the Parks Highway in south-central Alaska between the towns of Nenana (north) and Trapper Creek (south). A map of these can be found from the [FDSN network page](http://ds.iris.edu/gmap/#network=ZE&maxlat=64.8752&maxlon=-147.5002&minlat=62.227&minlon=-151.5871&drawingmode=box&planet=earth). 

## Code Organization

Most of the code is contained inside jupyter notebooks aside from prelude.py which contains shared functions usefull for looking at the data.

## Installation
To download the parkshwynodal project code, type this:

```git clone https://github.com/uafgeotools/parkshwynodal.git```


Next, to enter the repository type:

```cd parkshwynodal```

To create the conda environment and install dependencies type:

```conda env create environment.yml```

To enter this environment, type: 

```conda activate denalinodal```

Once this is all done try running the sample script by typing:

```python sample_inversion_script.py```

## Example Data product
![stations](/input/sample_image.png)



