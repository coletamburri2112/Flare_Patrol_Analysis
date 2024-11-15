# Flare Patrol Analysis 
# Converted DKIST working code from original DKIST_Analysis_Package

OWNER
    Cole A. Tamburri
    cole.tamburri@colorado.edu
    
INSTITUTION
    University of Colorado Boulder
    National Solar Observatory
    Laboratory for Atmospheric and Space Physics
    
MOST RECENT UPDATE
    15 November 2024
    
QUICK SUMMARY
    Analysis package for DKIST data, particularly ViSP and VBI. 
    
    Includes intensity calibration, line fitting, co-alignment routines.  Original intention was 
    flexibility for application to many DKIST proposals and datasets, but original version written to 
    properly analyze ViSP and VBI spectra for pid_1_84 ("flare patrol"), which
    leverages high resolution of DKIST instruments to study chromospheric 
    diagnostics during impulsive flare phase, in particular the spatial and
    temporal differences in evolution between Ca II H 396.8 nm and H-epsilon
    397.0 nm lines, with H-beta data added by PID_2_11 in September 2024. 
    File headers should include author information for any added functionality.
    
    November 2024: Addition of RADYN+RH output processing for application to DKIST spectral profiles
    specifically for H-epsilon and Ca II H, but intention to add H-beta.
    
    
FILE SUMMARIES
    WORKING_SOURCE folder includes .py files (mainly, dkistpkg_ct.py) handling much of the
    co-alignment and intensity calibration.  Also other, more direct working files to handle specific
    tasks unrelated directly to data processing.
    
    RHD_COMPARE folder includes code to compare RADYN+RH RHD simulations (originally 
    using runs from the F-CHROMA grid (Carlsson et al. 2023) of RADYN models) to DKIST data 
    directly, for publication in Tamburri et al. 2025
    
    PID_X_YY folders include specific working scripts for data collected in execution of e.g. PID_1_84
    and PID_2_11 where X.YY Is the proposal ID corresponding to DKIST OCP observing window X.
    
    MISC_PROPOSALS folder includes miscellaneous scripts used to investigate data from other
    publicly-available flare data.
    
    NOTEBOOKS folder includes .ipynb file formats used to test data processing code before formalziation
    in the WORKING_SOURCE folder or specific PID_X_YY folder.
        
ACKNOWLEDGEMENT, COLLABORATION
    As the availability and applicability of DKIST observations grows, we 
    welcome collaboration.  Please contact cole.tamburri@colorado.edu with 
    request to contribute; see https://github.com/coletamburri2112/Flare_Patrol_Analysis
    for current source code.
    
    Please reference package owner, author of specific code, and GitHub 
    repository when using source code found here.
