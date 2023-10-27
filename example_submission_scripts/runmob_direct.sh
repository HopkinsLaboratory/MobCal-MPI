#!/bin/sh
 ls *.mfj > mobfit.in 
 input="mobfit.in"
 suffix=".mfj"
 cnt=0
 one=1
 while IFS= read -r var
 do
  inp=${var%$suffix}
  echo "$inp" > temp.in
  tail -n +2 $var >> temp.in
  cp temp.in $var
  rm temp.in
  echo "$inp.mfj" > mobcal.run
  echo "$inp.mout" >> mobcal.run
  mobinp="mobcal.run_$cnt"
  mv mobcal.run $mobinp
  echo "running $inp.mfj"
  mpirun -np 4 ./MobCal_MPI_201.exe $mobinp &
  sleep 0.5s
  cnt=$(($cnt+$one))
 done < "$input"
 rm $input
