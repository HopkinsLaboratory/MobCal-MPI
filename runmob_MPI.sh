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
  cp mobcal.run $mobinp
  runfile="runit_$cnt"
  echo "#!/bin/bash" > $runfile
  echo "#SBATCH --account=rrg-shopkins-ab" >> $runfile
  echo "#SBATCH --nodes=1" >> $runfile
  echo "#SBATCH --ntasks-per-node=8" >> $runfile
  echo "#SBATCH --mem=8000M" >> $runfile
  echo "#SBATCH --time=0-03:00" >> $runfile
  echo " srun MobCal_MPI.exe $mobinp" >> $runfile
  sbatch $runfile
  sleep 0.5s
  cnt=$(($cnt+$one))
 done < "$input"
