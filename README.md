MethodOfMomentsExample
======================

A toy implementation of a 1D method of moments alongside a fit.
For debugging and learning purposes.

What it does
------------
This code uses the RooFit framework to generate a toy dataset whose shape is a simple polynomial (although you can modify this to any RooFit pdf you like).
It then runs the method of moments (MoM) over this toy and calculates the chi2 / ndf and p-values for: the original pdf (as a control); the MoM pdf; and an optional fitted pdf.
If you run this over many toy datasets, it will produce a distribution of the chi2's and p-values.

Requirements
------------
* root.cern.ch > 6.34/18
* gcc >= 4.8.1

If you're an LHCb collaborator, I made this in a `Urania v2r0p1` environment.
So for an LHCb environment do
```bash
SetupProject Urania v2r0p1
```
before the instructions in 'Quick Start' below.

Quick start
-----------
```bash
git clone https://github.com/green0eggs/MethodOfMomentsExample.git # grab this code
cd MethodOfMomentsExample
make
./bin/run --help # read the options
```


```bash
./bin/run # does one toy by default
# look at ./plots/toy0.pdf with your favourite viewer
```


```bash
./bin/run --ntoys=10
# look at ./plots/{chi2s,pvals}.pdf with your favourite viewer.
```
