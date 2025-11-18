# Analysis of Nodal Seismic Data collected along the Parks Highway in Feb-March of 2019
## University of Alaska Fairbanks

## Reference

Seppi, I., C. Tape, and D. Fee, 2025, [Classification of aircraft types using seismic data in Alaska](https://pubs.geoscienceworld.org/ssa/tsr/article/5/4/330/688104/Classification-of-Aircraft-Types-Using-Seismic), *The Seismic Record*.

## Data

Between February 11th and March 26th of 2019 a set of 303 Fairfield Nodal 3C 5Hz sensors were deployed along the Parks Highway in south-central Alaska between the towns of Nenana (north) and Trapper Creek (south). A map of these can be found from the [FDSN network page](http://ds.iris.edu/gmap/#network=ZE&maxlat=64.8752&maxlon=-147.5002&minlat=62.227&minlon=-151.5871&drawingmode=box&planet=earth).


## Code Organization

Most of the code is scripts that read in functions from doppler_funcs.py, containing tools to run the inversions, and from main_inv_fig_functions.py, containing tools to plot and display flight crossings in the nodal data. The input folder contains our database of aircraft that come within a 2km horizontal distance of one of the nodes along with other information used to classify aircrafts in the data. The output folder contains the inversion results. Database information can also be found on [Zenodo](https://zenodo.org/records/16997158).

## Installation

To download the parkshwynodal project code from your terminal, type this:

```git clone https://github.com/uafgeotools/parkshwynodal.git```

Next, to enter the repository:

```cd parkshwynodal```

To create the conda environment and install dependencies:

```conda env create environment.yml```

To enter this environment:

```conda activate denalinodal```

Once this is all done try running the sample script from your terminal:

```python sample_inversion_script.py```

A spectrogram should pop up and you can follow instructions printed in the terminal to select data points on the figure. You will need to confirm you want to keep the points you picked after closing the image or you will be prompted to repick them. Three images will pop up for you to pick data on and then a final example image should appear. After picking data the inversion results for each iteration will be printed in the terminal followed by misfit values at the end of each round of iteration. The final plot should look similar to the plot on the left side of the example data product below. The plot header and doppler estimation (blue curves) will vary depending on how the data was picked.

## Example Data product
![stations](/input/sample_image.png)



