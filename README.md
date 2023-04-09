# analysis codes for CLAS12 SIDIS analyses
## Tmothy B. Hayward, README last updated April 9, 2023

repository for various analysis codes (EventBuilder, fitters, etc.) used for analyzing CLAS12 data at Jefferson Lab. Primarily SIDIS focused. Modern iteration of [my previous analysis software](https://github.com/tbhayward/clas_analysis_code)

--------

included files:  
I. **/processing_classes/src/extended_kinematic_fitters/analysis_fitter.java**  
&nbsp;&nbsp;&nbsp;This is the class for the kinematic fitter I use to build events (takes the detector responses to assign particle ID to tracks and adds them to the event). The idea is to take the CLAS12 EventBuilder as a basis and enhance the PID on top of that. Loops through all particles in REC::Particle bank and sees if they pass the enhanced particle PID cuts (e.g. tightened sampling fraction, fiducial cuts, chi2pid cuts for hadron identification etc.) Start reading around line 700,  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"public PhysicsEvent getPhysicsEvent(DataEvent event) {"

IIa. **/processing_classes/src/analyzers/Inclusive.java**  
&nbsp;&nbsp;&nbsp;This is the class used to calculate relevant kinematic variables (Q2, W, Mx, xF, PT, phi_trento, etc.)  for ep -> e'X events (inclusive DIS). 
  
IIb. **/processing_classes/src/analyzers/Hadron.java**  
&nbsp;&nbsp;&nbsp;This is the class used to calculate relevant kinematic variables (Q2, W, Mx, xF, PT, phi_trento, etc.)  for ep -> e'hX events (single hadron SIDIS). 

IIc. **/processing_classes/src/analyzers/Dihadrons.java**
  Extension (technically written first) of the single hadron case to two hadrons, ep -> e' h1 h2 X, includes additional variables for each hadron, e.g. z1, z2, PT1, PT2, etc.
  
IId. **/processing_classes/src/analyzers/Trihadron.java**
  Extension of the single hadron case to three hadrons, ep -> e' h1 h2 h3 X, includes additional variables for each hadron, e.g. z1, z2, z3, PT1, PT2, PT3 etc. as well as all permutations of hadron combinations, i.e. Mh12, Mh13, Mh23

&nbsp;&nbsp;&nbsp;The above physics classes start with a "channel_test" function that allows for cuts on Q2, W, xF, y, Mx, etc. By default Q2 > 1, W > 2, y < 0.8* are enabled and the others are commented out. If you are sure of the cuts you desire, more can be turned on to significantly increase processing speed. *I'd recommend double checking these.
