LSE
===

This is the Legged State Estimator (LSE) library.
It contains an Observability Constrained Extended Kalman Filter and a routine for time delay calibration between different sensor modalities.
The Manager class represents the interface. It handles the different measurements.

INSTALLATION:
- go to top folder (the folder this file is located in)
- execute cmake: "cmake ."
- either simply compile the library: "make", or call the installation: "sudo make install" (this will copy the library and header to your system folder and create a corresponding FindLSE.cmake file)

DEPENDENCIES:
- Standard library
- Eigen3
- The required tinyxml library is provided together with this installation
