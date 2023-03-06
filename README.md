# Modal Anaylsis with poly-reference least squares comples frequency method (p-LSCF)

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

Implementation of poly-reference least squares complex frequency method in Matlab.

## Table of contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Background

Operational Modal Analysis (OMA) identifies the modal parameters (natural frequency, damping ratio and eigenform) of a structure from experimentally determined measured data.
The special feature of the OMA is that the excitation of the structure is not measured and is thus unknown.
Given the structure excitation by white noise, the power spectral density of the measured system response contains the complete information on the modal parameters.

The poly-reference least squares complex frequency method (p-LSCF) is the current industry standard of OMA methods in the frequency domain.
In a first least squares step, the stabilization diagram is constructed based on the parametric model to identify the stable poles of the structure.
The eigenforms are determined in a second least-square step.
This method can also identify closely spaced eigenmodes and provides clear, easily interpretable stabilization diagrams.

The p-LSCF method implemented as a fully automated program in MATLAB in this repository was built in a Bachelor's thesis at KIT, Germany (see [here](#reference) and verified by several data sets.

## Install

This project uses [Matlab](https://www.mathworks.com/products/matlab.html), so make sure you have installed it before using this software.
The following toolboxes are required:
- Signal Processing Toolbox

The script has been tested with the following Matlab versions (other versions are likely to work too):
- R2020b
- R2022a
- R2022b

You can just clone this repository or download the code to your PC for installation.

```bash
git clone https://github.com/KIT-FAST/modalAnalysis-OMA-EMA.git
```

## Usage

The code is contained within the folder `modalAnalysis`.
You can open the file `modalAnalysis_pLSCF_main.m` in your Matlab GUI and run the script.

All the configuration values for the user are placed in section `2. User Input`.
Just change this values before executing the script.

The script is provided with some example data so that it is possible to run and test the script directly after downloading.
If you want to use your own data, you can do as follows for OMA or EMA.
The required format is described in more detail in the [pdf](https://github.com/KIT-FAST/modalAnalysis-OMA-EMA/blob/master/Implementierung%20des%20p-LSCF-Algorithmus%20zur%20Operational%20Modal%20Analysis.pdf) on page 46, unfortunately it is written in German.

Further information can be found [here](#reference).


### Input data format for EMA

For EMA the expected data is a matrix `FRF` with the frequencies in the 1st dimension and the sensors in the 2nd dimension, so each column represents a sensor. Moreover, a column vector `frequencyBand` with the corresponding frequencies is expected.

`size(FRF) = n_frequencies, n_sensors`
`size(frequencyBand) = n_frequencies`

![ema_trim](https://user-images.githubusercontent.com/13416487/86216540-7b78a900-bb7e-11ea-9b8d-5858e035ef32.png)

### Input data format for OMA (single measurement)

For OMA with a single measurement, a matrix `records` with the time in the 1st dimension and the sensors in the 2nd dimension is expected, so each column represents a sensor.

`size(records) = n_time, n_sensors`

![oma_single_trim](https://user-images.githubusercontent.com/13416487/86216624-9f3bef00-bb7e-11ea-9fd4-7723b0bdaddc.png)

### Input data format for OMA (multiple measurements)

For OMA with multiple measurements, the matrix `records` becomes 3-dimensional with the 3rd dimension representing each measurement. The reference sensors must always be the first columns in the matrix.

`size(records) = n_time, n_sensors, n_measurements`

![oma_multiple_trim](https://user-images.githubusercontent.com/13416487/86216634-a236df80-bb7e-11ea-8536-2248995bd00e.png)


## Contributing

If you have questions about the usage of this repository, then head over to the issues.
First search if you questions has already been asked and answered (see also the closed issues!).
If your questions has not been asked, just create a new issue.
Please provide a detailed description, this makes it easier to answer your question.

As this repository is no longer further developed by us, there is no policy on how you can contribute improvements to this repository.
If you want to improve this repository, just raise an issue and we can discuss how this can be done.

## License

[GNU General Public License v3.0 only" (GPL-3.0)](LICENSE.txt) © [raphajaner](https://github.com/raphajaner)

## Reference
This repository was built in following Bachelor's thesis at KIT, Germany:
```git 
@thesis{Trumpp2017_1000156492,
    author       = {Trumpp, Raphael},
    year         = {2017},
    title        = {Implementierung des poly-reference least square complex frequency (p-LSCF) Algorithmus zur Operational Modal Analysis},
    type         = {Abschlussarbeit - Bachelor},
    school       = {Karlsruher Institut für Technologie (KIT)},
    language     = {german}
}
```
