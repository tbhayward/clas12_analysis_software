#!/bin/bash
# - makes a latex table for the analysis note
# - this is meant to be run after `makeTables.sh`
# - columns:
#   - run number
#   - fraction of golden files
#   - fraction of files with defect bit 0
#   - fraction of files with defect bit 1
#   - fraction of files with defect bit 2
#   - fraction of files with defect bit 3
#   - fraction of files with defect bit 4
#   - fraction of files with defect bit 5
# - cf. makeLatexTables2.sh
if [ -z "$QADB" ]; then
  echo "ERROR: you must source env.sh first"; exit
fi
pushd $QADB/qadb > /dev/null
run-groovy ../util/makeLatexTable.groovy

nrow=45 # number of rows per table

mv qaTable.tex{,.tmp}

function app { echo "$1" >> qaTable.tex; }
function header {
  app '  \begin{tabular}{|c|c|c|c|c|c|c|c|}'
  app '    \hline {\bf Run} & {\bf Golden Fraction} & {\bf Bit 0}& {\bf Bit 1}& {\bf Bit 2}& {\bf Bit 3}& {\bf Bit 4}& {\bf Bit 5}\\\hline'
}

app '\begin{table}'
app '  \begin{center}'
header

rcnt=0
tcnt=1
while read line; do
  line=$(echo $line | sed 's/\.00//g')

  #echo $line
  app "    $line"
  ((rcnt++))
  if [ $rcnt -eq $nrow ]; then
    app '  \end{tabular}'
    if [ $tcnt -eq 1 ]; then app "  \\caption{Run QA table.}"
    else app "  \\caption{Continuation of run QA table.}"; fi
    app "  \\label{tab_qa_result_${tcnt}}"
    app '  \end{center}'
    app '\end{table}'
    app '\begin{table}'
    app '  \begin{center}'
    header
    ((tcnt++))
    rcnt=0
  fi
done < qaTable.tex.tmp
app '  \end{tabular}'
app '  \caption{Continuation of run QA table.}'
app "  \\label{tab_qa_result_${tcnt}}"
app '  \end{center}'
app '\end{table}'

rm qaTable.tex.tmp

echo "done: see ${QADB}/qadb/qaTable.tex"

popd > /dev/null
