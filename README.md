# iron-snow-fortran
Fortran model for computing the thermal evolution of a planetary core with an iron snow zone


Code is contained within the `code` folder. Inside, the makefile will compile to a `thermal_history` executable which requires an input file.
The input file for the reference case is contained in the `reference_case` folder, with the following parameters

- rho        = 7211 kg/m3
- P_cmb      = 21 GPa
- k          = 40 W/m/K
- r_cmb      = 1627 km
- T_cmb(t=0) = 2400 K

as per the reference case described in Davies and Pommier (2018).

The code on the main branch has been modified slightly from the original code. The original code is available on the `original` branch. The modification is in the gravitational energy released associated with the moving boundary of the snow zone, which has been changed to be proportional to the jump in Sulphur mass fraction (c^s-c^l), rather than just the Sulphur mass fraction (c^s) (see Davies and Pommier (2018) equation 17).


