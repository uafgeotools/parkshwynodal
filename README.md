# Code to analyze Nodal Seismic Data collected along the Parks Highway in Feb-March of 2019
## University of Alaska Fairbanks

## Data

Between Febuary 11th and March 26th of 2019 a set of 303 Fairfield Nodal 3C 5Hz sensors were deployed along the Parks Highway in south-central Alaska between the towns of Nenana (north) and Trapper Creek (south). A map of these can be found from the [FDSN network page](http://ds.iris.edu/gmap/#network=ZE&maxlat=64.8752&maxlon=-147.5002&minlat=62.227&minlon=-151.5871&drawingmode=box&planet=earth). 

## Code Organization

Most of the code is scripts that read in functions from doppler_funcs.py, containing tools to run the inverisons, and main_inv_fig_functions.py, containing tools to plot and display flight crossings in the nodal seismic data. The input folder contains our database of plane that come within a 2km horizontal distance of one of the nodes along with other information used to classify aircrafts in the data. The output folder contains the inversion results. Dtabase information can also be found on [Zenodo](https://zenodo.org/records/16997158).

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



