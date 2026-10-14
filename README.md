# Flare Patrol Analysis 
### Converted DKIST working code from original DKIST_Analysis_Package

OWNER
    
    Cole A. Tamburri
    cole.tamburri@colorado.edu
    
INSTITUTION
    
    University of Colorado Boulder
    National Solar Observatory
    Laboratory for Atmospheric and Space Physics
    
MOST RECENT UPDATE
    
    14 October 2026
    
QUICK SUMMARY
    
    Analysis codebase for DKIST flare data data, especially ViSP with some VBI.
    
    Includes intensity calibration, line fitting, co-alignment routines.  Original intention was 
    flexibility for application to many DKIST proposals and datasets, but original version was written to 
    properly analyze ViSP and VBI spectra for pid_1_84 ("flare patrol"), which
    leverages high resolution of DKIST instruments to study chromospheric 
    diagnostics during impulsive flare phase, in particular the spatial and
    temporal differences in evolution between Ca II H 396.8 nm and H-epsilon
    397.0 nm lines.  H-beta data added by PID_2_11 in September 2024. 
    File headers should include author information for any added functionality.
    
    November 2024: Addition of RADYN+RH output processing for application to DKIST spectral profiles
    specifically for H-epsilon and Ca II H, but intention to add H-beta.

    October 2026: In the past two years, have added analysis of X-class flare loop widths
    and k-means clustering of ViSP spectra, now with a .py file that creates a
    k-means clustering object with flexibility for many k values, performance metrics, etc.
    Much of the intensity calibration has also been improved, and specific co-alignment
    for the ViSP arms and VBI images in pid 2.11 added.  Co-alignment is done on a case-by-case
    basis and shouldn't be trusted for anything other than these two proposals; updates planned
    Fall 2026.
    
    
FILE SUMMARIES
    
    WORKING_SOURCE folder includes .py files (mainly, dkistpkg_ct.py) handling much of the
    co-alignment and intensity calibration.  Also other, more direct working files to handle specific
    tasks unrelated directly to data processing. k-means clustering files include the original script
    "flare_kmeans_script.py" and more formal object-oriented version "Clusterer.py" which is called
    in specific flare_kmeans_script_[PID].py" files, applied to PID 4.17 and PID 2.11 as of Oct 2026.
    
    RHD_COMPARE folder includes code to compare RADYN+RH RHD simulations (originally 
    using runs from the F-CHROMA grid (Carlsson et al. 2023) of RADYN models) to DKIST data 
    directly, for publication in Tamburri et al. 2026a.
    
    PID_X_YY folders include specific working scripts for data collected in execution of e.g. PID_1_84
    and PID_2_11 where X.YY Is the proposal ID corresponding to DKIST OCP observing window X.
    
    MISC_PROPOSALS folder includes miscellaneous scripts used to investigate data from other
    publicly-available flare data.
    
    NOTEBOOKS folder includes .ipynb file formats used to test data processing code before formalziation
    in the WORKING_SOURCE folder or specific PID_X_YY folder. There are also many notebooks used to generate 
    the figures and some extra analysis for Tamburri et al. 2026b.
        
ACKNOWLEDGEMENT, COLLABORATION
    
    As the availability and applicability of DKIST observations grows, we 
    welcome collaboration.  Please contact c.tamburri@qub.ac.uk with 
    request to contribute; see https://github.com/coletamburri2112/Flare_Patrol_Analysis
    for current source code.
    
    Please reference package owner, author of specific code, and GitHub 
    repository when using source code found here.
