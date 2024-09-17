The code in this directory is used for estimating a systematic uncertainty on the accumulated FC charge, based on <https://github.com/JeffersonLab/clas12-qadb/issues/48>.

### Procedure

- edit the scripts such that the specific list of runs is used
- `run-groovy calculateGapCharge.groovy`
- `root analyzeGapCharge.C`
- the mean of the resulting fractional charge-difference distribution is the systematic uncertainty
