# DKIST_analysis_package

OWNER
    Cole A. Tamburri
    cole.tamburri@colorado.edu
    
INSTITUTION
    University of Colorado Boulder
    National Solar Observatory
    Laboratory for Atmospheric and Space Physics
    
MOST RECENT UPDATE
    3 December 2023
    
QUICK SUMMARY
    Analysis package for DKIST data, particularly ViSP and VBI.  Includes
    intensity calibration, line fitting, co-alignment routines.  Intention is 
    applicability to all DKIST proposals, but original version written to 
    properly analyze ViSP and VBI spectra for pid_1_84 ("flare patrol"), which
    leverages high resolution of DKIST instruments to study chromospheric 
    diagnostics during impulsive flare phase, in particular the spatial and
    temporal differences in evolution between Ca II H 396.8 nm and H-epsilon
    397.0 nm lines. File headers should include author information for any
    added functionality.
    
FILE SUMMARIES
    dkistpkg_ct.py: Main working functions used by C.A.T.
    19_August_2022_DKIST_analysis.py: Initial intensity calibration, line 
        evolution (incl. widths, strengths, etc.) for pid_1_84 observations of
        decay phase of C6.7 class flare beginning 20:42 UT on 19 August 2022.
        Observations are limited to the decay phase, and seeing conditions
        are variable, making tracking of chromospheric lines difficult in a
        meaningful way, but functionality of package is demonstrated.
    19_Aug_2022_co_align.py: Demonstration of co-alignment routines for flare
        observations referenced immediately above.  L1 header information for
        pointing of DKIST instrument is inaccurate, so co-alignment necessary
        for accurate representation of solar spatial coordinates.  First ViSP 
        is co-aligned with VBI such that the two share a reference frame, and 
        then VBI is co-aligned with SDO instruments (AIA or HMI; flexibility in
        this functionality currently limited) and both ViSP and VBI placed in
        that reference frame.
        
ACKNOWLEDGEMENT, COLLABORATION
    As the availiability and applicability of DKIST observations grows, we 
    welcome collaboration.  Please contact cole.tamburri@colorado.edu with 
    request to contribute; see https://github.com/coletamburri2112/DKIST_analysis_package
    for current source code.
    
    Please reference package owner, author of specific code, and GitHub 
    repository when using source code found here.
