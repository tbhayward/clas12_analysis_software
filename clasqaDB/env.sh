#!/bin/bash

export QADB=$(dirname $(realpath $0))

JYPATH="${JYPATH}:${QADB}/src/"
export JYPATH=$(echo $JYPATH | sed 's/^://')

env|grep --color -w QADB
env|grep --color -w JYPATH
