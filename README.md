# exact_Compton
These routines compute the exact redistribution function due to Compton scattering

Included is a driving program (driveSRF.f) which calls all pertinent routines.

The main output of the code is what we calle the "super redistribution function
(SRF)" for Compton scatterting. For a given gas temperature T and final photon
energy Ef, the SRF is defined for a set of initial photon energies Ei as:

         SRF(T,Ef,Ei) = IRF(Ef,Ei)/N(Ei)*skn(Ei)*dEi/Ei
     
where IRF(Ef,Ei) is in fact the inverse redistribution function for the Compton
scattering of a photon from initial energy Ei to final energy Ef; N(Ei) is the
normalization to ensure photon number conservation; and skn(Ei) is the
Klein-Nishina cross section.  This routine implements the exact Compton RF from
Madej et al. (2017).  The SRF contains all the information needed for the
convolution of a given spectrum to account for the Compton scattering at the
given temperature.  Only significant values of the SRF are actually written,
i.e., when RF > limit.

The output of this code is used by the XILLVER model.

## Input file structure
Five parameters can be read in by the program by creating the file <I>drive_params.dat</I>.  The parameters are:
* nmaxp (number of energy points in the grid)
* itrans (number of temperatures in the grid)
* mgi (number of angles to evaluate at)
* limit (the fractional limit of the maximum value of the SRF, below which it's assumed to be 0.)
* avangle (whether to return the averaged/integrated SRF over angles)

The file should be formatted as so (these are the default parameters within the code)


      nmaxp 500
      itrans 70
      mgi 3000
      limit 0.001
      avangle true
## Output file structure
The fits file produced contains five HDUs.  It will be called <I>angles.fits</I> for the angular resolved code, else <I>table.fits</I> for the angular integrated code.
### 1.  Primary HDU
This doesn't contain any infomation, and exists simply for the datastructure of the fits file.
### 2.  Parameters
This contains two or four columns, and one row.  The last two columns only exist when the code produces an angular resolved SRF.
1. TEMP (K) - array of temperature points used.
2. ENERGIES (eV) - array energy points used.
3. cos_ANGLES - The cosine of the angle of the Guassian quadratures used.
4. WEIGHTS - The weighting that needs to be applied to the SRF at each angle.

### 3.  KN_cross
One column, with one row.  This contains the Klein–Nishina cross-section for every photon energy and electron temperature.

Each element in the array corresponds to a temperature/energy grid pair.  For N energy points, and T temperature points, the first N elements correspond to the first temperature grid point; the elements (N+1)...(2N) correspond to the second temperature point, and so on.
### 4.  SRF_Pointers
This contains two columns, with a row for each angle (only one row for the case of an angular integrated SRF)
1. IND - the starting indices of the energy array which the first value of the saved iSRF array corresponds to.
Each element here is an array, containing values for each row of the iSRF.
2. LEN - the length of the saved iSRF array. Each element here is an array, containing values for each row of the iSRF.
### 5.  iSRF
This contains one column per angle (with only one column for the angular integrated SRF).  Each row corresponds to a final energy & temperature pair, in a similar way as the elements of the KN cross-section (table 3.)

Each element in each column is a variable length array, having trimmed out all the elements below the cut-off limit.  Use the IND and LEN columns from table 4 to reconstruct the full iSRF.
## Code info
     Version: 0.5.0 - Tue Feb 28th 14:35:00 PDT 2023

     Authors: Javier Garcia (javier@caltech.edu)
              Ekaterina Sokolova-Lapa (ekaterina.sokolova-lapa@fau.de)
                (see Garcia et al. 2020 in prep)
              Jameson Dong (jdong2@caltech.edu)
              Isabel Franco (francog@caltech.edu)
              Guglielmo Mastroserio (guglielmo.mastroserio@inaf.it)
              Edward Nathan (enathan@caltech.edu)
     With routines provided by J. Madej and A. Rozanska (see Madej et al. 2017).
