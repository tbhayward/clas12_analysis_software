# CLAS12 analysis framework developed by Timothy B. Hayward
* repository for various analysis codes (EventBuilder, fitters, etc.) used for analyzing CLAS12 data at Jefferson Lab. Primarily SIDIS focused. Modern iteration of [my previous analysis software](https://github.com/tbhayward/clas_analysis_code).
* README last updated April 9, 2023

### Table of Contents
1. [QA Information](#info)
2. [Processing CLAS12 data](#processing)
3. [Analyzing CLAS12 data](#analyzing)
4. [Relevant Files](#files)

<a name="info"></a>
# QA Information
Processing scripts make use of the [clasqa database](https://github.com/JeffersonLab/clasqaDB) developed by C. Dilks for the collaboration. By default the "OKForAsymmetry" criteria is enforced. Other options are described in the clasqaDB repository, including a list of "Golden Runs". qadb is sourced automatically when processing scripts are run without the need for user input.

<a name="processing"></a>
# Processing CLAS12 data
Users can process CLAS12 data with the "processing_scripts/processing.csh" shell script and providing the required input arguments. 

```processing_scripts/processing.csh [processing script] [arg2] [arg3] [arg4] ...```

where [processing script] corresponds to one of the groovy processing scripts also included in the processing_scripts directory (e.g. the inclusive script to process e'X events, the single hadron script to process e'hX events, etc.) and [arg2], [arg3], [arg4], ... are the necessary arguments to execute those scripts (defined below). Generally the first argument provides a directory with hipo (CLAS12) data files you wish to analyze, the next arguments give the [PDG PID](https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf) numbers for the analyzed particles (e.g. pi+ = 211), the next argument gives an output filename and the final argument defines how many files to analyze. Only the processing_script argument is strictly necessary and the shell script will run with assumptions about the rest of your arguments (and provide appropriate warnings), however, it is highly recommended that you provide all arguments. 

The output of the processing scripts will be text data files that place the relevant event-by-event kinematics and physics observables (e.g. particle 4-vectors, Q^2, W, x, y, z, etc.) in each of the columns.  

### Examples
1. Process the inclusive electron sample 

```processing_scripts/processing.csh processing_scripts/processing_inclusive.groovy /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/ /scratch/thayward/inclusive.txt```

2. Processing the semi-inclusive pi+ sample 

```processing_scripts/processing.csh processing_scripts/processing_single_hadrons.groovy /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/ 211 /scratch/thayward/single_hadrons.txt```

3. Processing the semi-inclusive k+ pi- sample 

```processing_scripts/processing.csh processing_scripts/processing_dihadronss.groovy /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/ 321 -211 /scratch/thayward/dihadrons.txt```

<a name="analyzing"></a>
# Analyzing CLAS12 data

Section under development. See [BSA fits](https://github.com/tbhayward/clas12_analysis_software/blob/main/analysis_scripts/BSA_fits.C) for a current iteration of my ROOT script to perform chi2 and maximum-likelihood fits for the ongoing CLAS12 SIDIS back-to-back dihadron analysis. Will update with further scripts and options soon.

--------

<a name="files"></a>
# Relevant Files
I. **/processing_classes/src/extended_kinematic_fitters/analysis_fitter.java**  
&nbsp;&nbsp;&nbsp;This is the class for the kinematic fitter I use to build events (takes the detector responses to assign particle ID to tracks and adds them to the event). The idea is to take the CLAS12 EventBuilder as a basis and enhance the PID on top of that. Loops through all particles in REC::Particle bank and sees if they pass the enhanced particle PID cuts (e.g. tightened sampling fraction, fiducial cuts, chi2pid cuts for hadron identification etc.) Start reading around line 700,  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"public PhysicsEvent getPhysicsEvent(DataEvent event) {"

IIa. **/processing_classes/src/analyzers/Inclusive.java**  
&nbsp;&nbsp;&nbsp;This is the class used to calculate relevant kinematic variables (Q2, W, Mx, xF, PT, phi_trento, etc.)  for ep -> e'X events (inclusive DIS). 
  
IIb. **/processing_classes/src/analyzers/Hadron.java**  
&nbsp;&nbsp;&nbsp;This is the class used to calculate relevant kinematic variables (Q2, W, Mx, xF, PT, phi_trento, etc.)  for ep -> e'hX events (single hadron SIDIS). 

Columns: &nbsp;&nbsp;&nbsp; 1:runnum, 2:evnum, 3:helicity, 4:e_p, 5:e_theta, 6:e_phi, 7:vz_e, 8:p_p, 9:p_theta, 10:p_phi, 11:vz_p, 12:Q2, 13:W, 14:Mx, 15: Mx2, 16:x, 17:y, 18:z, 19:xF, 20:pT, 21:zeta, 22:eta, 23:phi (trento), 24:DepA, 25:DepB, 26:DepC, 27:DepV, 28:DepW

IIc. **/processing_classes/src/analyzers/Dihadrons.java**
  Extension (technically written first) of the single hadron case to two hadrons, ep -> e' h1 h2 X, includes additional variables for each hadron, e.g. z1, z2, PT1, PT2, etc.
  
Columns: &nbsp;&nbsp;&nbsp; 1:runnum, 2:evnum, 3:helicity, 4:e_p, 5:e_theta, 6:p1_p, 7:p1_theta, 8:p2_p, 9:p2_theta, 10:Q2, 11:W, 12:Mx, 13:Mx1, 14:Mx2, 15:x, 16:y, 17:z, 18:z1, 19:z2, 20:Mh, 21:xF, 22:xF1, 23:xF2, 24:pT, 25:pT1, 26:pT2, 27:pTpT 28:zeta, 29:eta, 30:eta1, 31:eta2, 32:Delta_eta, 33:eta1_gN, 34:eta2_gN, 35:phi1, 36:phi2, 37:Delta_phi, 38:phih, 39:phiR, 40:theta, 41:DepA, 42:DepB, 43:DepC, 44:DepV, 45:DepW, 46: vz_e, 47: vz_p1, 48: vz_p2.
  
IId. **/processing_classes/src/analyzers/Trihadron.java**
  Extension of the single hadron case to three hadrons, ep -> e' h1 h2 h3 X, includes additional variables for each hadron, e.g. z1, z2, z3, PT1, PT2, PT3 etc. as well as all permutations of hadron combinations, i.e. Mh12, Mh13, Mh23

&nbsp;&nbsp;&nbsp;The above physics classes start with a "channel_test" function that allows for cuts on Q2, W, xF, y, Mx, etc. By default Q2 > 1, W > 2, y < 0.8* are enabled and the others are commented out. If you are sure of the cuts you desire, more can be turned on to significantly increase processing speed. *I'd recommend double checking these.

III. **/processing_classes/dist/processing_classes.jar**
&nbsp;&nbsp;&nbsp;The distribution version of all my classes if you want to include the package in your own build.

IV. **processing_single_hadrons.groovy**
&nbsp;&nbsp;&nbsp;I prefer to process the clas12 hipo4 files and create text outputs that I can then import into Mathematica, ROOT, etc. It's not the most efficient but it is very nicely universal. Includes a section (commented out by default to increase computation speed) for the printing of revelant RICH variables for hadron PID studies. If no track is present in the RICH (which only exists in one sector; second sector installed for RGC forward) then it prints PID = 0. This script accepts 4 input arguments:  
  1. hipo file directory  
  a directory such as, /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/, that contains hipo4 data files.  
  2. pid for the hadron i.e. 211 for pi+  
  3. output text file name  
  4. number of files to process in the directory
  5. 
IVb. **processing_dihadrons.groovy**
&nbsp;&nbsp;&nbsp;dihadron version of the above. This script accepts 5 input arguments:  
  1. hipo file directory  
  a directory such as, /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/, that contains hipo4 data files.  
  2. pid for p1 (hadron 1), i.e. 211 for pi+  
  3. pid for p2 (hadron 2), i.e. -211 for pi-
  4. output text file name  
  5. number of files to process in the directory
