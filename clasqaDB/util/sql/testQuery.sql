SELECT runnum, filenum, sector
FROM qa
INNER JOIN yield USING(runnum,filenum);
