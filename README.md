# The Tidal.ninja model
This code performs the 'core' calculations of the Tidal.ninja model, which is described in a pre-print manuscript.

It delivers an hourly time series of capacity factors for tidal energy plants for a given location, modelling either tidal stream or tidal range.  It includes an R wrapper for the [OSU TPXO Tide Models](https://www.tpxo.net/) and soft-links this with models for tidal turbine operation.



## Setup
You must first download a copy of the [OTPS software](https://www.tpxo.net/otps), the version which reads binary files rather than OTPSnc.  You need to compile `predict_tide.f90` on your system, so that you have an executable file which can be called from R's `system()` command.

You must then obtain a copy of the [TPXO 10 Atlas (v2)](https://www.tpxo.net/global/tpxo10-atlas) in the OTPS binary format.  The TPXO Atlas requires registration to access, but is freely available for non-commercial purposes.

Tidal.ninja requires `R` version 4+, with the `terra` library.  Tested on Windows, should also work on Linux and macOS.


## Usage
The calculation code lives inside `tidal_ninja.r`. In that file, you just need to change the path to say where you have downloaded all the files.

Two examples of how to simulate tidal plants at a single location are given in `example_tidal_stream.r` and `example_tidal_range.r`.


## Contact
Contact [Iain Staffell](i.staffell@imperial.ac.uk) if you have questions about Tidal.ninja.


## Citation
If you use Tidal.ninja or code derived from it in academic work, please cite:

I. Staffell, A. Pusey, J.J. Xie, and N. Johnson (2025). _The global potential and hourly capacity factors of tidal range and tidal stream technologies_. Pre-print.


If you use the OTPS software which Tidal.ninja relies upon, please cite:

G.D. Egbert and S.Y. Erofeeva (2002). [_Efficient Inverse Modeling of Barotropic Ocean Tides_](https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2). Journal of Atmospheric and Oceanic Technology, 19, 2, 183–204.


## License
The Tidal.ninja code is available under the BSD-3-Clause license.  The OTPS Software and TPXO Atlas must be obtained separately, and are copyright Oregon State University.  The OSU terms around non-commercial usage (without license) must be respected.

