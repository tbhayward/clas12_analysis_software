#!/bin/tcsh

# Use the current working directory as the script directory
set this_dir = $cwd
setenv QADB $this_dir

set src_dir=${QADB}/src/
if (! $?JYPATH) then
  setenv JYPATH $src_dir
else
  if ("$JYPATH" == "") then
    setenv JYPATH $src_dir
  else
    setenv JYPATH $src_dir\:$JYPATH
  endif
endif