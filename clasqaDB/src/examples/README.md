# Examples reading HIPO files using Coatjava

Here are some examples of reading HIPO files in an event loop, using
[`coatjava`](https://github.com/JeffersonLab/clas12-offline-software)
If you do not use `coatjava` in your analysis, you can still use these
examples to see how to implement `QADB.groovy` in your analysis. See
also the C++ implementation.


## Instructions
* make sure environment variable `$COATJAVA` is set to your `coatjava` installation
* use `${COATJAVA}/bin/run-groovy` to execute these examples
  * it is useful to have `${COATJAVA}/bin` in your `$PATH`
* most of the example scripts require a HIPO file as an argument
* example scripts named `test*.groovy` are actually developer tools for
  testing the QADB, but may also be useful as additional examples
