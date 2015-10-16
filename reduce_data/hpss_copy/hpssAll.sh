#!/bin/sh
rm -f  hpssAll.txt
touch  hpssAll.txt
for i in `cat file.list`
do
  echo $i pftp://beaumim@rcas2066/direct/phenix+spin2/beaumim/vernierScans/prdf_analysis/HPSS/files/`basename $i` >>  hpssAll.txt
done
