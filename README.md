# MRXImage
This program is to reconstruct the dipole sources that are responsible for the magnetic field detected by the SQUID detectors in the MRX.  

Program files:
    MRXImage.m: Main file with optimization for source reconstruction, using CVX.  User input chooses between SeDuMi and SDPT3 algorithms.
    header.m: header file that defines physical parameters such as FOV, detector positions, which detectors to use, etc.
    headerVariables.mat: a mat file that stores the variables defined in the header.m code, to be loaded into subroutines.
    makeBField.m: Calculates the initial B field vector for a user-defined source.
    plotSourceScatter.m: Plots output from MRXImage.m in a scatter plot
    gaussianDistance.m: Calculates the gaussian distance metric between a reconstructed source and the true source
    forwardproblem.m: Calculates the full B field using the Biot Savart Law.
