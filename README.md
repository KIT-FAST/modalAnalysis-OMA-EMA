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

The p-LSCF method implemented as a fully automated program in MATLAB in this repository was built in a Bachelor thesis and verified by several data sets.

## Install

## Usage

## Contributing

If you have questions about the usage of this repository, then head over to the issues.
First search if you questions has already been asked and answered (see also the closed issues!).
If your questions has not been asked, just create a new issue.
Please provide a detailed description, this makes it easier to answer your question.

As this repository is not developed further from our site, there is no policy on how you can contribute improvements to this repository.
If you want to improve this repository, just raise an issue and we can discuss how this can be done.

## License

[GNU General Public License v3.0 only" (GPL-3.0)](LICENSE.txt) © [raphajaner](https://github.com/raphajaner)
