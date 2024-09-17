#!/bin/bash

if [ -z "${BASH_SOURCE[0]}" ]; then
  export QADB=$(dirname $(realpath $0))
else
  export QADB=$(dirname $(realpath ${BASH_SOURCE[0]}))
fi

# class path for groovy
export JYPATH="${QADB}/src/${JYPATH:+:${JYPATH}}"
