grep "^Time" log | awk '{print $3}' > time.dat

grep Cy: log | awk '{print $4}' | sed 's/,//g' > Cy.dat
grep Cx: log | awk '{print $2}' | sed 's/,//g' > Cx.dat
paste time.dat Cx.dat > t_Cx.dat
paste time.dat Cy.dat > t_Cy.dat
sed '1s/^/t\ty\n/' t_Cx.dat > tmp; mv tmp t_Cx
sed '1s/^/t\ty\n/' t_Cy.dat > tmp; mv tmp t_Cy
