#!/bin/tcsh
set this_cmd=($_)
set this_dir=`dirname ${this_cmd[2]}`
setenv QADB `cd ${this_dir} && pwd -P`

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
