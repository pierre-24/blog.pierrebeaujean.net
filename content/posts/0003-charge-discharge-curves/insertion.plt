set term pngcairo enhanced font 'Helvetica,12' size 600,700 dashed
set output "insertion.png"

set multiplot layout 2,1


set samples 500
unset key

set xrange [0:1]

unset ytics
set ylabel "Free energy, G_{A_xBX}"

f(x) = -.5*x + x*log(x)+(1-x)*log(1-x)
df(x) = -.5 + log(x/(1-x))

xint=0.55

set lmargin 7
set rmargin 7
set arrow from xint,f(xint) to xint,-1.9 nohead dashtype 2

set linestyle 1 black

plot f(x) ls 1 lw 2, (df(xint)*x + (f(xint) - df(xint)*xint)) ls 1 dashtype 4


unset arrow
set ylabel "Potential, E^0"
set xlabel "x"
set yrange [-5:5]

set label at xint,-df(xint) "•" center font ",16"

set label at graph 0,-.12 "BX" center font ",16"
set label at graph 1,-.12 "ABX" center font ",16"

plot -df(x) ls 1 lw 2
