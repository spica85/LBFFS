grep "^Time" log | awk '{print $3}' > time.dat

grep CyIBM: log | awk '{print $4}' | sed 's/,//g' > CyIBM.dat
grep CxIBM: log | awk '{print $2}' | sed 's/,//g' > CxIBM.dat
paste time.dat CxIBM.dat > t_CxIBM.dat
paste time.dat CyIBM.dat > t_CyIBM.dat
sed '1s/^/t\ty\n/' t_CxIBM.dat > tmp; mv tmp t_CxIBM
sed '1s/^/t\ty\n/' t_CyIBM.dat > tmp; mv tmp t_CyIBM