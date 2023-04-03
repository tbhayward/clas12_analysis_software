#!/bin/bash
# - makes a latex table for the analysis note
# - runs makeLatexTable2.groovy
# - columns:
#   - run number
#   - fraction of files that pass QA cuts
#   - list of files that fail QA cuts
# - only runs which have some files that pass QA cuts are
#   included in the list
# - cf. makeLatexTables.sh

if [ -z "$QADB" ]; then
  echo "ERROR: you must source env.sh first"; exit
fi
pushd $QADB/qadb > /dev/null

# set run period
if [ $# -ne 1 ]; then
  echo """
  USAGE: $0 [run-period-number]
  - see case statement in script
  - edit the script to add your own run period
    (a workaround, since QADB does not (yet) store run period name)
  """
  exit 2
fi
periodNum=$1
case $periodNum in
  0)
    period="RGA Fall 2018, inbending"
    runrange="5032 5419"
    ;;
  1)
    period="RGA Fall 2018, outbending"
    runrange="5422 5666"
    ;;
  2)
    period="RGA Spring 2019, inbending"
    runrange="6616 6783"
    ;;
  3)
    period="RGB Spring 2019, inbending"
    runrange="6156 6603"
    ;;
  4)
    period="RGB Fall 2019, outbending"
    runrange="11093 11283"
    ;;
  5)
    period="RGB Spring 2020, inbending"
    runrange="11323 11571"
    ;;
  6)
    period="RGK 7.5 GeV"
    runrange="5674 5870"
    ;;
  7)
    period="RGK 6.5 GeV"
    runrange="5875 6000"
    ;;
  *)
    echo "ERROR: unknown run-period-number (edit the script!)"
    exit 1
    ;;
esac

run-groovy ../util/makeLatexTable2.groovy $runrange

mv qaTable.tex{,.tmp}

function app { echo "$1" >> qaTable.tex; }

app "\\textbf{Run Period: $period}"
app '\begin{center}'
app '\begin{longtable}{|p{0.15\linewidth}|p{0.2\linewidth}|p{0.6\linewidth}|}'
app '\hline {\bf Run} & {\bf Pass Fraction}& {\bf Excluded Files}\\\hline\hline\endhead'

while read line; do
  line=$(echo $line | sed 's/\.00//g')
  app "  $line"
done < qaTable.tex.tmp
app "\\caption{Run list for $period period.}"
app "\\label{tab:runlist_${periodNum}}"
app '\end{longtable}'
app '\end{center}'

rm qaTable.tex.tmp
mv qaTable{,_$periodNum}.tex
echo "done: see ${QADB}/qadb/qaTable_${periodNum}.tex"
popd > /dev/null
