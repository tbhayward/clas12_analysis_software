![QADB](/doc/logo.png)

# CLAS12 Quality Assurance Database
Provides storage of and access to the QA monitoring results for the 
CLAS12 experiment at Jefferson Lab

### Table of Contents
1. [QA Information](#info)
1. [Database Access](#access)
1. [Data storage](#storage)
1. [Faraday Cup Charge Access](#charge)
1. [Database Management](#dev)
1. [Contributions](#contributions)

<a name="info"></a>
# QA Information

> [!CAUTION]
> The QADB for older data sets may have some issues. It is **HIGHLY recommended** to [check the known important issues](/doc/issues.md) to see if any issues impact your analysis.

## Available Data Sets
The following tables describe the available datasets in the QADB. The columns are:
- **Pass**: the Pass number of the data set (higher is newer)
- **Data Set Name**: a unique name for the data-taking period; click it to see the corresponding QA timelines
  - Typically `[RUN_GROUP]_[RUN_PERIOD]`
  - `[RUN_PERIOD]` follows the convention `[SEASON(sp/su/fa/wi)]_[YEAR]`, and sometimes includes an additional keyword
- **Run range**: the run numbers in this data set
- **Status**:
  - _Up-to-Date_: this is the most recent Pass of these data, and the QADB has been updated for it
  - _Deprecated_: a newer Pass exists for these data, but the QADB for this version is still preserved
  - _TO DO_: the Pass for these data exist, but the QADB has not yet been updated for it
- **Data files**: the DST files used for the QA

### Run Group A

| Pass | Data Set Name and Timelines Link                                                           | Run Range   | Status       | Data Files                                                                  |
| ---  | ---                                                                                        | ---         | ---          | ---                                                                         |
| 2    | `rga_sp19`                                                                                 | 6616 - 6783 | _TO DO_      | `/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/recon`    |
| 1    | [`rga_fa18_inbending`](https://clas12mon.jlab.org/rga/pass1/qa/fa18_inbending/tlsummary)   | 5032 - 5419 | _Up-to-Date_ | `/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/recon`   |
| 1    | [`rga_fa18_outbending`](https://clas12mon.jlab.org/rga/pass1/qa/fa18_outbending/tlsummary) | 5422 - 5666 | _Up-to-Date_ | `/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v0/dst/recon`   |
| 1    | [`rga_sp19`](https://clas12mon.jlab.org/rga/pass1/qa/sp19/tlsummary)                       | 6616 - 6783 | _Up-to-Date_ | `/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v0/dst/recon` |

### Run Group B

| Pass | Data Set Name and Timelines Link                                     | Run Range     | Status       | Data Files                                                                   |
| ---  | ---                                                                  | ---           | ---          | ---                                                                          |
| 2    | `rgb_sp19`                                                           | 6156 - 6603   | _TO DO_      | `/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/recon/` |
| 1    | [`rgb_sp19`](https://clas12mon.jlab.org/rgb/pass1/qa/sp19/tlsummary) | 6156 - 6603   | _Up-to-Date_ | `/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/recon`  |
| 1    | [`rgb_fa19`](https://clas12mon.jlab.org/rgb/pass1/qa/fa19/tlsummary) | 11093 - 11300 | _Up-to-Date_ | `/cache/clas12/rg-b/production/recon/fall2019/torus+1/pass1/v1/dst/recon`    |
| 1    | [`rgb_wi20`](https://clas12mon.jlab.org/rgb/pass1/qa/wi20/tlsummary) | 11323 - 11571 | _Up-to-Date_ | `/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass1/v1/dst/recon`  |

### Run Group C

| Pass | Data Set Name and Timelines Link                                                              | Run Range     | Status       | Data Files                                     |
| ---  | ---                                                                                           | ---           | ---          | ---                                            |
| 1    | [`rgc_su22`](https://clas12mon.jlab.org/rgc/Summer2022/qa-physics/pass1-sidisdvcs/tlsummary/) | 16042 - 16771 | _Up-to-Date_ | `/cache/clas12/rg-c/production/summer22/pass1` |

### Run Group F

| Pass | Data Set Name and Timelines Link | Run Range     | Status  | Data Files                                                                                  |
| ---  | ---                              | ---           | ---     | ---                                                                                         |
| 1    | `rgf_sp20_torusM1`               | 12210 - 12329 | _TO DO_ | `/cache/clas12/rg-f/production/recon/spring2020/torus-1_solenoid-0.8/pass1v0/dst/recon`     |
| 1    | `rgf_su20_torusPh`               | 12389 - 12434 | _TO DO_ | `/cache/clas12/rg-f/production/recon/summer2020/torus+0.5_solenoid-0.745/pass1v0/dst/recon` |
| 1    | `rgf_su20_torusMh`               | 12436 - 12443 | _TO DO_ | `/cache/clas12/rg-f/production/recon/summer2020/torus-0.5_solenoid-0.745/pass1v0/dst/recon` |
| 1    | `rgf_su20_torusM1`               | 12447 - 12951 | _TO DO_ | `/cache/clas12/rg-f/production/recon/summer2020/torus-1_solenoid-0.745/pass1v0/dst/recon`   |

### Run Group K

| Pass | Data Set Name and Timelines Link                                                   | Run Range   | Status       | Data Files                                                                        |
| ---  | ---                                                                                | ---         | ---          | ---                                                                               |
| 1    | [`rgk_fa18_7.5GeV`](https://clas12mon.jlab.org/rgk/pass1/qa/fa18_7.5GeV/tlsummary) | 5674 - 5870 | _Up-to-Date_ | `/cache/clas12/rg-k/production/recon/fall2018/torus+1/7546MeV/pass1/v0/dst/recon` |
| 1    | [`rgk_fa18_6.5GeV`](https://clas12mon.jlab.org/rgk/pass1/qa/fa18_6.5GeV/tlsummary) | 5875 - 6000 | _Up-to-Date_ | `/cache/clas12/rg-k/production/recon/fall2018/torus+1/6535MeV/pass1/v0/dst/recon` |

### Run Group M

| Pass | Data Set Name and Timelines Link                                                     | Run Range     | Status       | Data Files                                                  |
| ---  | ---                                                                                  | ---           | ---          | ---                                                         |
| 1    | [`rgm_fa21`](https://clas12mon.jlab.org/rgm/pass1_finalqadb/rgm_fall2021/tlsummary/) | 15019 - 15884 | _Up-to-Date_ | `/cache/clas12/rg-m/production/pass1/allData_forTimelines/` |


## Defect Bit Definitions

* QA information is stored for each DST file, in the form of "defect bits"
  * the user needs only the run number and event number to query the QADB
  * queries will find the DST file associated with the event, and are only
    performed "as needed"
  * full dumps of the QADB are also possible, for browsing
* N/F is defined as the electron yield N, normalized by the Faraday Cup charge F; the
  electron yield is for Forward Detector electrons with `status<0`, unless specified otherwise
  * The QA checks for outliers of N/F, along with several other miscellaneous criteria
  * The term "golden" means that a file has ***no*** defects
* The table below lists the defect bits 
  (Example: `defect=0b11000` has defects `SectorLoss` and `LowLiveTime`)

| Bit | Name                | Description                                                                           |
| --- | ---                 | ---                                                                                   |
| 0   | `TotalOutlier`      | outlier N/F, but not terminal, marginal, or sector loss, for FD electron              |
| 1   | `TerminalOutlier`   | outlier N/F of first or last file of run, not marginal, for FD electron               |
| 2   | `MarginalOutlier`   | marginal outlier N/F, within one standard deviation of cut line, for FD electron      |
| 3   | `SectorLoss`        | N/F diminished within a FD sector for several consecutive files                       |
| 4   | `LowLiveTime`       | live time < 0.9                                                                       |
| 5   | `Misc`              | miscellaneous defect, documented as comment                                           |
| 6   | `TotalOutlierFT`    | outlier N/F, but not terminal, marginal, or `LossFT`, FT electron                     |
| 7   | `TerminalOutlierFT` | outlier N/F of first or last file of run, not marginal, FT electron                   |
| 8   | `MarginalOutlierFT` | marginal outlier N/F, within one standard deviation of cut line, FT electron          |
| 9   | `LossFT`            | N/F diminished within FT for several consecutive files                                |
| 10  | `BSAWrong`          | Beam Spin Asymmetry is the wrong sign                                                 |
| 11  | `BSAUnknown`        | Beam Spin Asymmetry is unknown, likely because of low statistics                      |
| 12  | `TSAWrong`          | __[not yet used]__ Target Spin Asymmetry is the wrong sign                            |
| 13  | `TSAUnknown`        | __[not yet used]__ Target Spin Asymmetry is unknown, likely because of low statistics |
| 14  | `DSAWrong`          | __[not yet used]__ Double Spin Asymmetry is the wrong sign                            |
| 15  | `DSAUnknown`        | __[not yet used]__ Double Spin Asymmetry is unknown, likely because of low statistics |
| 16  | `ChargeHigh`        | FC Charge is abnormally high                                                          |
| 17  | `ChargeNegative`    | FC Charge is negative                                                                 |
| 18  | `ChargeUnknown`     | FC Charge is unknown; the first and last time bins always have this defect            |
| 19  | `PossiblyNoBeam`    | Both N and F are low, indicating the beam was possibly off                            |

<a name="access"></a>
# Database Access

## Text Access
  * this only provides human-readable access; see below for access with
    common programming languages and software used at CLAS
    * using the Groovy or C++ access is the preferred method to apply QA cuts
  * the human-readable tables are stored in `qadb/*/qaTree.json.table`; see
    the section *QA data storage, Table files* below for details for how
    to read these files
  * QADB JSON files are stored in `qadb/*/qaTree.json`
    * the JSON files are the QADB files
      * for now we use JSON out of convenience, although it's not a proper
        database format; future development plans include considering more
        efficient formats, such as `SQLlite`
  * there are also some text files stored in `text/`:
    * `text/listOfGoldenRuns.txt`: list of runs, each classified as one of the following:
      * `gold`: all files have no defects. Note that this is **very strict**,
        so not many runs are `gold`, since most runs have at least one file
        with a defect; in practice it is better to apply QA cuts per file,
        using the QADB software
      * `silver`: the only defects are terminal outliers (first or last file is
        an outlier); note that `gold` runs are, by definition, also `silver`.
        This is also **very strict**: so far, only about half the runs are
        `silver`
      * `defect`: not `gold` or `silver`
    * `text/listOfGoldenFiles.txt`: list of files with no defects
    * `text/summary.txt`: summary table, where for each file the QA criteria result
      (`Golden`, `OkForAsymmetry`, etc.) is provided
  * it is also possible to produce `latex` tables; see
    `util/makeLatexTables.sh` and `util/makeLatexTables2.sh`

## Software Access

Classes in both C++ and Groovy are provided, for access to the QADB within analysis code

### Groovy Access
* first set environment variables by running `source environ.sh`
  * `bash` is recommended, though if you choose to use `tcsh`, run
    instead `source environ.csh`
* then proceed following `src/README.md`

### C++ Access
* __NOTE:__ [`clas12root`](https://github.com/JeffersonLab/clas12root) now provides
  access to the QADB
* needs [`rapidjson`](https://github.com/Tencent/rapidjson/) library; 
  it is a submodule of this repository and can be obtained by
  ```
  git clone --recurse-submodules https://github.com/JeffersonLab/clas12-qadb.git
  ```
* first set environment variables by running `source environ.sh`
  * alternatively, set environment variable `$QADB` to the path to this
    `clas12-qadb` repository
  * `bash` is recommended, though if you choose to use `tcsh`, run
    instead `source environ.csh`
* then proceed following `srcC/README.md`

<a name="storage"></a>
# Data Storage

## Table files
Human-readable format of QA result, stored in `qadb/*/*/qaTree.json.table`
* each run begins with the keyword `RUN:`; lines below are for each of that 
  run's file and its QA result, with the following syntax:
  * `run number` `file number`  `defect bits` `comment`
  * the `defect bits` are listed by name, and the numbers in the `[brackets]`
    indicate which sectors have that defect
  * if a comment is included, it will be printed after the defect bits, following the
    `::` delimiter
* these table files can be generated from the JSON files using `bin/makeTables.sh`

## JSON files

### `qaTree.json`
* The QADB itself is stored as JSON files in `qadb/*/*/qaTree.json`
* the format is a tree (nested maps):
```
qaTree.json ─┬─ run number 1
             ├─ run number 2 ─┬─ file number 1
             │                ├─ file number 2
             │                ├─ file number 3 ─┬─ evnumMin
             │                │                 ├─ evnumMax
             │                │                 ├─ sectorDefects
             │                │                 ├─ defect
             │                │                 └─ comment
             │                ├─ file number 4
             │                └─ file number 5
             ├─ run number 3
             └─ run number 4
```
* for each file, the following variables are defined:
  * `evnumMin` and `evnumMax` represent the range of event numbers associated
    to this file; use this to map a particular event number to a file number
  * `sectorDefects` is a map with sector number keys paired with lists of associated
    defect bits
  * `defect` is a decimal representation of the `OR` of each sector's defect bits, for
    example, `11=0b1011` means the `OR` of the defect bit lists is `[0,1,3]`
  * `comment` stores an optional comment regarding the QA result

### `chargeTree.json`
* the charge is also stored in JSON files in `qadb/*/*/chargeTree.json`, with
  a similar format:
```
chargeTree.json ─┬─ run number 1
                 ├─ run number 2 ─┬─ file number 1
                 │                ├─ file number 2
                 │                ├─ file number 3 ─┬─ fcChargeMin
                 │                │                 ├─ fcChargeMax
                 │                │                 ├─ ufcChargeMin
                 │                │                 ├─ ufcChargeMax
                 │                │                 └─ nElec ─┬─ sector 1
                 │                │                           ├─ sector 2
                 │                │                           ├─ sector 3
                 │                │                           ├─ sector 4
                 │                │                           ├─ sector 5
                 │                │                           └─ sector 6
                 │                ├─ file number 4
                 │                └─ file number 5
                 ├─ run number 3
                 └─ run number 4
```
* for each file, the following variables are defined:
  * `fcChargeMin` and `fcChargeMax` represent the minimum and maximum DAQ-gated
    Faraday cup charge, in nC
  * `ufcChargeMin` and `ufcChargeMax` represent the minimum and maximum FC charge,
    but not gated by the DAQ
  * the difference between the maximum and minimum charge is the accumulated charge
    in that file
  * `nElec` lists the number of electrons from each sector


<a name="charge"></a>
# Faraday Cup Charge Access
* the charge is stored in the QADB for each DST file, so that it is possible to
  determine the amount of accumulated charge for data that satisfy your
  specified QA criteria.
* see `src/examples/chargeSum.groovy` or `srcC/examples/chargeSum.cpp` for
  usage example in an analysis event loop
  * call `QADB::AccumulateCharge()` within your event loop, after your QA cuts
    are satisfied; the QADB instance will keep track of the accumulated charge
    you analyzed (accumulation performed per DST file)
  * at the end of your event loop, the total accumulated charge you analyzed is
    given by `QADB::getAccumulatedCharge()`
* note: for Pass 1 QA results, we find some evidence that the charge from file to file may slightly overlap,
  or there may be gaps in the accumulated charge between each file; the former leads to
  a slight over-counting and the latter leads to a slight under-counting
  * for RGK, we find the correction to this issue would be very small
    (no more than the order of 0.1%)
  * corrections of this issue are therefore not applied
  * if you require higher precision of the accumulated charge than what is
    provided, contact the developers to discuss an implementation of the
    corrections


<a name="dev"></a>
# QADB Management

Documentation for QADB maintenance and revision

## Adding to or revising the QADB
* the QADB files are produced by [`clasqa` timeline-production code](https://github.com/c-dilks/clasqa);
  if you have produced QA results for a new data set, and would like to add
  them to the QADB, or if you would like to update results for an existing
  dataset, follow the following procedure:
  * `mkdir qadb/pass${pass}/${dataset}/`, then copy the final `qaTree.json` and
    `chargeTree.json` to that directory
  * add/update a symlink to this dataset in `qadb/latest`, if this is a new Pass
  * run `source environ.sh`
  * run `bin/makeTables.sh`
  * run `bin/makeTextFiles.sh`
  * update customized QA criteria sets, such as `OkForAsymmetry`
  * update the above table of data sets
  * use `git status` and `git diff` to review changes, then add and commit to
    git, and push to the remote branch

## Adding new defect bits
* defect bits must be added in the following places:
  * Groovy:
    * `src/clasqa/Tools.groovy` (copy from `clasqa` repository version)
    * `src/clasqa/QADB.groovy`
    * `src/examples/dumpQADB.groovy` (optional)
  * C++:
    * `srcC/include/QADB.h`
    * `srcC/examples/dumpQADB.cpp` (optional)
  * Documentation:
    * bits table in `README.md`


<a name="contributions"></a>
# Contributions

All contributions are welcome, whether to the code, examples, documentation, or the QADB itself. You are welcome to open an issue and/or a pull request. If the maintainer(s) do not respond in a reasonable time, send them an email.
