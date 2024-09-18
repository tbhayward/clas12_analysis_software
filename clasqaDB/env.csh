#!/bin/tcsh

setenv QADB `pwd -P`
setenv JYPATH "${QADB}/src/"
env|grep --color -w JYPATH
