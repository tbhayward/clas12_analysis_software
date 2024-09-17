#!/bin/bash

dataTable=$1 # data_table.dat file

QADB=QADB.sqlite
> $QADB

# create tables
sqlite3 $QADB << EOF
CREATE TABLE charge(
  runnum INTEGER,
  filenum INTEGER,
  fcStart REAL,
  fcStop REAL,
  ufcStart REAL,
  ufcStop REAL
);
EOF
sqlite3 $QADB << EOF
CREATE TABLE yield(
  runnum INTEGER,
  filenum INTEGER,
  sector INTEGER,
  nElec REAL,
  nElecFT REAL
);
EOF
sqlite3 $QADB << EOF
CREATE TABLE qa(
  runnum INTEGER,
  filenum INTEGER,
  evnumMin INTEGER,
  evnumMax INTEGER,
  defect INTEGER
);
CREATE UNIQUE INDEX dstId
ON qa (runnum,filenum);
EOF


# loop through dataTable, filling tables
while IFS=" " read -r\
  runnum\
  filenum\
  evnumMin\
  evnumMax\
  sector\
  nElec\
  nElecFT\
  fcStart\
  fcStop\
  ufcStart\
  ufcStop
do
  if [ $filenum -eq 0 -a $sector -eq 1 ]; then
    echo "read run $runnum"
  fi

  sqlite3 $QADB << EOF
  INSERT INTO yield(
    runnum,
    filenum,
    sector,
    nElec,
    nElecFT
  ) VALUES (
    $runnum,
    $filenum,
    $sector,
    $nElec,
    $nElecFT
  );
EOF

  # sector-independent tables
  if [ $sector -eq 1 ]; then
    sqlite3 $QADB << EOF
    INSERT INTO charge(
      runnum,
      filenum,
      fcStart,
      fcStop,
      ufcStart,
      ufcStop
    ) VALUES (
      $runnum,
      $filenum,
      $fcStart,
      $fcStop,
      $ufcStart,
      $ufcStop
    );
EOF
    sqlite3 $QADB << EOF
    INSERT INTO qa(
      runnum,
      filenum,
      evnumMin,
      evnumMax,
      defect
    ) VALUES (
      $runnum,
      $filenum,
      $evnumMin,
      $evnumMax,
      0
    );
EOF
  fi # eo sector-independent tables
done < $dataTable
