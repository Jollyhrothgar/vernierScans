This Version of HPSS is added to the vernierScans repository for completeness, as the
analysis used the files grabbed here (ppg PRDFFS, i.e. scaler-events only)  The code is
actually maintained as a separate package, which will be later added somehere on github
and called something like "PHENIX_TOOLS/HPSS" or someting. The main branch of HPSS
currently resides on a private CVS repository in Brookhaven National Lab.

(1) Modify the file.list file to contain the HPSS directory of the files you want.
(1.a.) To browse the HPSS directories, use the "hsi" command from the command line 
(2) Modify hpssAll.sh to reflect your user name and the directory to which you copy the
    HPSS files
(3) Run hpssAll.sh, a file called "hpssAll.txt" is created.
(4) Modify get.sh to the analysis username you use (phnxspin, in this case)
(5) run get.sh
(6) Wait
(7) Keep waiting
(8) After about 15-30 minutes you can check the status with (a modified) mysqlCheck.sh to
    be your user name.
