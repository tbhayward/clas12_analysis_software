# CLAS12 analysis framework developed by Timothy B. Hayward

* repository for various analysis codes (EventBuilder, fitters, etc.) used for analyzing CLAS12 data at Jefferson Lab. Primarily SIDIS focused. Modern iteration of [my previous analysis software](https://github.com/tbhayward/clas_analysis_code)
* README last updated April 9, 2023

## Table of Contents
1. [QA Information](#info)
2. [Processing CLAS12 data](#processing)
    1. [Inclusive](#inclusive)
    2. [Single Hadron](#single-hadron)
    3. [Dihadron](#dihadron)
    4. [Trihadron](#trihadron)
3. [Relevant Files](#files)
    1. [/processing_classes/src/extended_kinematic_fitters/analysis_fitter.java](#analysis-fitter)
    2. [/processing_classes/src/analyzers/Inclusive.java](#inclusive-class)
    3. [/processing_classes/src/analyzers/Hadron.java](#hadron-class)
    4. [/processing_classes/src/analyzers/Dihadrons.java](#dihadrons-class)
    5. [/processing_classes/src/analyzers/Trihadron.java](#trihadron-class)
    6. [/processing_classes/dist/processing_classes.jar](#processing-classes-jar)
    7. [processing_single_hadrons.groovy](#single-hadrons-groovy)
    8. [processing_dihadrons.groovy](#dihadrons-groovy)

<a name="info"></a>
## QA Information
Processing scripts make use of the [clasqa database](https://github.com/JeffersonLab/clasqaDB) developed by C. Dilks for the collaboration. By default the "OKForAsymmetry" criteria is enforced. Other options are described in the clasqaDB repository, including a list of "Golden Runs". qadb is sourced automatically when processing scripts are run without the need for user input.

<a name="processing"></a>
## Processing CLAS12 data
Users can process CLAS12 data with the "processing_scripts/processing.csh" shell script and providing the required input arguments. 

### Inclusive
```processing_scripts/processing.csh inclusive [arg2] [arg3] [arg4] ...```

### Single Hadron
```processing_scripts/processing.csh single_hadron [arg2] [arg3] [arg4] ...```

### Dihadron
```processing_scripts/processing.csh dihadrons [arg2] [arg3] [arg4] ...```

### Trihadron
```processing_scripts/processing.csh trihadrons [arg2] [arg3] [arg4] ...```

<a name="files"></a>
## Relevant Files
### /processing_classes/src/extended_kinematic_fitters/analysis_fitter.java
This is the class for the kinematic fitter I use to build events (takes the detector responses to assign particle ID to tracks and adds them to the event). The idea is to take the CLAS12 EventBuilder as a basis and enhance the PID on top of that. Loops through all particles in REC::Particle bank and
