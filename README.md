# cs-kg4vasp
slight change of kg4vasp so that it handles collinear spin cases for DC cases ONLY

Original kg4vasp package is here: https://github.com/conodipaola/kg4vasp

use NISS=1 for non-spin-polarized cases and NISS=2 for spin-polarized cases in INPUT

tested using vasp5.4.4

compile the two .f90 files in places of the original files (without the _2 suffix)

everything is dropped except for DC conductivity.
