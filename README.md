========================================================================
   PW_DBEAM: Double-Beam Plasma Dispersion Solver 
========================================================================

[Authors]      Yasuhito Narita (1), Uwe Motschmann (1), Horia Comisel (2), and Daniel Schmid (3)

[Contact]     y.narita@tu-braunschweig.de

[Affiliation] (1) Institute of Theoretical Physics, Technical University of Braunschweig, Braunschweig, Germany. 

(2) Institute for Space Sciences, Bucharest-Magurele, Romania

(3) Space Research Institute, Austrian Academy of Sciences, Graz, Austria

[License]     MIT License

[DOI]         https://doi.org/10.3847/1538-4357/adc1bc

------------------------------------------------------------------------
1. OVERVIEW
------------------------------------------------------------------------
This repository hosts a self-contained scientific script for analyzing 
plasma instabilities.
The code **`pw_dbeam.py`** solves the complex dispersion relations 
for right-hand modes propagating parallel to the background magnetic 
field in a double-beam electron-ion plasma system (e-i-b-d). 
The term "pw" in the dispersion solver stands for the plasma waves.
The code operates within the displacement-current-free approximation, 
valid for low-frequency waves under dilute plasma (Alfven speed must
be sufficiently lower than the speed of light).
The code tracks growth rates gamma for unstable modes 
excited by two beam species.

The solver outputs figures in PDF format and serializes 
raw multidimensional calculation matrices into NumPy binary 
archive tracks (`.npz`) for post-processing and analysis.

Details of the algorithm are presented in the following journal article:

    Narita, Y., Motschmann, U., Comi\c{s}el, H., and Schmid, D.:
    Double beam instability for the Mercury upstream waves,
    Astrophys. J., 983, 125 (2026.
    https://doi.org/10.3847/1538-4357/adc1bc

------------------------------------------------------------------------
2. REQUIREMENTS
------------------------------------------------------------------------
- Python 3.x
- NumPy
- Matplotlib (optional, for plotting)

------------------------------------------------------------------------
3. INSTALLATION
------------------------------------------------------------------------
The code pw_dbeam.py can run on the terminal.

------------------------------------------------------------------------
4. BASIC USAGE
------------------------------------------------------------------------
Set the following parameters: 

 (1) Particle mass parameters.
 The electron mass is set to unity. The bulk ion mass
 is set by the parameter omega_ce (with the minus sign),
 and the beam ion masses (beam b and beam d) by omega_cb
 and omega_cd, respectively. 
 The default is the proton mass for both the bulk ions and the beam ions,
  omega_ce = -1836, omega_ci = 1, omega_cb = 1, omega_cd = 1.

 (2) Particle density parameters.
 The bulk ion density (number density) is set by 
  the ion plasma frequency omega_pi (default 3000),
  and the beam densities by omega_pb (default 100)
  and omega_pd (default 100).
  The plasma frequency is proportional to the square density
  of the respective species. The electron density omega_pe is
  automatically determined by the charge neutrality condition.

 (3) Particle velocity parameters.
  The frequencies are solved in the frame of vanishing
  total momentum. The beam speeds are given by
  v_b (the first beam species, default 7 Alfven speed)
  and v_d (the second beam species, defaul 5 Alfven speed).
  The bulk electron and ion velocities are chosen such that 
  the total momentum and the total current vanish.

 (4) Wavenumber array (parallel wavenumbers).
  Give the starting and ending wavenumbers 
  in units of the ion inertial wavenumer (k c /omega_pi)
  and the number of wavenumber arrays. Default is 
  kcw_min = 0.0, kcw_max = 0.5, and kcw_ticks = 1000.


Run the code on the terminal by typing, e.g.,

\>\> python3 pw_dbeam.py

The outputs are the PDF file showing the dispersion relation
and the growth rate, and the NPZ binary data file containing
(1) the wavenumbers, (2) the real part of frequencies,
(3) the wavenumbers of the unstable modes (growth rate gamma/omega_ci 
larger than 1.0e-4, (4) the real part of frequencies of
the unstable modes, and (5) the growth rate (imaginary
part of frequencies) of the unstable modes.


------------------------------------------------------------------------
5. CITATION
------------------------------------------------------------------------
If you use this code in your research, please cite it as:

    Narita, Y., Motschmann, U., Comi\c{s}el, H., and Schmid, D.:
    Double beam instability for the Mercury upstream waves,
    Astrophys. J., 983, 125 (2026.
    https://doi.org/10.3847/1538-4357/adc1bc

------------------------------------------------------------------------
6. LICENSE
------------------------------------------------------------------------
This project is licensed under the MIT License.


