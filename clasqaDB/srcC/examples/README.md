# Examples using Clas12Root

Here are some examples of reading HIPO files in an event loop, using
[`clas12root`](https://github.com/JeffersonLab/clas12root)
If you do not use `clas12root` in your analysis, you can still use these
examples to see how to implement `QADB.h` in your analysis. See also the
groovy implementation.

## Instructions
* make sure environment variable `$CLAS12ROOT` is set to your `clas12root` installation
* be sure `$LD_LIBRARY_PATH` includes `${CLAS12ROOT}/lib`
* compile examples with `make`
* executables `*.exe` will appear; most of them require a HIPO file as an
  argument
* example programs named `test*.cpp` are actually developer tools for
  testing the QADB, but may also be useful as additional examples
