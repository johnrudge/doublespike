# The double spike toolbox

Error propagation and data reduction for the double spike technique in mass spectrometry.

The workings of the toolbox are described in:

Rudge J.F., Reynolds B.C., Bourdon B. The double spike toolbox (2009) Chem. Geol. 265:420-431
https://dx.doi.org/10.1016/j.chemgeo.2009.05.010

and at:

https://johnrudge.com/doublespike

This is the MATLAB version of the code. A python version is also available at:
https://github.com/johnrudge/pyspike

Requirements: MATLAB 7.1 (R14SP3) or higher with the Statistics and Optimization Toolboxes.

Quick start instructions:

1, Extract the contents of the zip file to a folder on the hard drive.
2, Add the doublespike folder to your MATLAB path. (From MATLAB's menu: File-> Set Path->Add Folder)
3, Type 'doc doublespike' for a brief guide to the toolbox or 'help doublespike' to get a list of commands.

"maininput.csv" contains the isotope numbers, atomic weights, standard compositions, and ORNL spike compositions for all the elements. These can be modified as desired (e.g. with more precise standard compositions or different commercial spike compositions). 

Copyright (c) John F. Rudge 2009-2021



