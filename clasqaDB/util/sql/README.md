# SQLite QADB

This is an idea for an alternative storage option for the QADB, as a SQLite
database. The script `buildSQL.sh` reads the sample `data_table.dat` file and
produces `QADB.sqlite`; a sample query is provided in `testQuery.sql`, which
can be executed as `sqlite3 QADB.sqlite < testQuery.sql`. Inserting QA information
from `qaTree.json` still needs to be implemented.
