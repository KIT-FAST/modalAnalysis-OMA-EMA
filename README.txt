# Modal Anaylsis with poly-reference least squares comples frequency method (p-LSCF)

Operational Modal Analysis (OMA) identifies the modal parameters (natural frequency, damping ratio and eigenform) of a structure from experimentally determined measured data.
The special feature of the OMA is that the excitation of the structure is not measured and is thus unknown.
Given the structure excitation by white noise, the power spectral density of the measured system response contains the complete information on the modal parameters.

The poly-reference least squares complex frequency method (p-LSCF) is the current industry standard of OMA methods in the frequency domain.
In a first least squares step, the stabilization diagram is constructed based on the parametric model to identify the stable poles of the structure.
The eigenforms are determined in a second least-square step.
This method can also identify closely spaced eigenmodes and provides clear, easily interpretable stabilization diagrams.

The p-LSCF method implemented as a fully automated program in MATLAB in this repository was built in a Bachelor thesis and verified by several data sets.
