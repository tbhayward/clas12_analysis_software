#! /usr/bin/env python

from __future__ import print_function
import platform, sysconfig, os, sys
from glob import glob


## Build dirs
srcdir = os.path.abspath("/home/thayward/clas12_analysis_software/LHAPDF/LHAPDF-6.5.4/src")
libdir = os.path.abspath("/home/thayward/clas12_analysis_software/LHAPDF/LHAPDF-6.5.4/src/.libs")
incdirs = [os.path.abspath("/home/thayward/clas12_analysis_software/LHAPDF/LHAPDF-6.5.4/include"),
           os.path.abspath("/home/thayward/clas12_analysis_software/LHAPDF/LHAPDF-6.5.4/include")]

## The extension source file (built by separate Cython call)
srcname = "lhapdf.cpp"
srcpath = os.path.join("/home/thayward/clas12_analysis_software/LHAPDF/LHAPDF-6.5.4/wrappers/python", srcname)
if not os.path.isfile(srcpath): # distcheck has it in srcdir
    srcpath = os.path.join("/home/thayward/clas12_analysis_software/LHAPDF/LHAPDF-6.5.4/wrappers/python", srcname)

## Include args
incargs = " ".join("-I{}".format(d) for d in incdirs)
incargs += " -I/u/home/thayward/clas12_analysis_software/LHAPDF/include"
incargs += "  -I/u/home/thayward/clas12_analysis_software/LHAPDF/include"

## Compile args
cmpargs = ""  #"@PYEXT_CXXFLAGS@"
cmpargs += " -O3"

## Link args -- base on install-prefix, or on local lib dirs for pre-install build
linkargs = " -L/home/thayward/clas12_analysis_software/LHAPDF/LHAPDF-6.5.4/src/.libs" if "LHAPDF_LOCAL" in os.environ else "-L/u/home/thayward/clas12_analysis_software/LHAPDF/lib"
linkargs += "  -L/u/home/thayward/clas12_analysis_software/LHAPDF/lib"

## Library args
libraries = ["stdc++", "LHAPDF"]
libargs = " ".join("-l{}".format(l) for l in libraries)

## Python compile/link args
pyargs = "-I" + sysconfig.get_config_var("INCLUDEPY")
libpys = [os.path.join(sysconfig.get_config_var(ld), sysconfig.get_config_var("LDLIBRARY")) for ld in ["LIBPL", "LIBDIR"]]
libpys.extend( glob(os.path.join(sysconfig.get_config_var("LIBPL"), "libpython*.*")) )
libpys.extend( glob(os.path.join(sysconfig.get_config_var("LIBDIR"), "libpython*.*")) )
libpy = None
for lp in libpys:
    if os.path.exists(lp):
        libpy = lp
        break
if libpy is None:
    print("No libpython found in expected location exiting")
    print("Considered locations were:", libpys)
    sys.exit(1)
pyargs += " " + libpy
pyargs += " " + sysconfig.get_config_var("LIBS")
pyargs += " " + sysconfig.get_config_var("LIBM")
#pyargs += " " + sysconfig.get_config_var("LINKFORSHARED")


## Assemble the compile & link command
compile_cmd = "  ".join([os.environ.get("CXX", "g++"), "-shared -fPIC",
                         "-o", srcname.replace(".cpp", ".so"),
                         srcpath, incargs, cmpargs, linkargs, libargs, pyargs])
print("Build command =", compile_cmd)


## (Re)make the build/rivet dir and copy in Python sources
import shutil
try:
    shutil.rmtree("build/lhapdf")
except:
    pass
try:
    os.makedirs("build/lhapdf")
except FileExistsError:
    pass
for pyfile in glob("./*.py"):
    if "build.py" not in pyfile:
        shutil.copy(pyfile, "build/lhapdf/")


## Run the extension compilation in the build dir
import subprocess
sys.exit(subprocess.call(compile_cmd.split(), cwd="/home/thayward/clas12_analysis_software/LHAPDF/LHAPDF-6.5.4/wrappers/python/build/lhapdf"))
