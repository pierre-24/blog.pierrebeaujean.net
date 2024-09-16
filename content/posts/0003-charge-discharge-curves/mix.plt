set term pngcairo enhanced font 'Helvetica,20' size 600,700 dashed
set output "mix.png"

set samples 500
unset key

set xrange [0:1]
set yrange [-.7:0.6]

unset ytics
set xlabel "Molar fraction in A, x_A"

f(x) = .5*x + x*log(x)+(1-x)*log(1-x)
df(x) = .5 + log(x/(1-x))

xint=0.4

set lmargin 7
set rmargin 7

set linestyle 1 black

set arrow from -.02,0 to  -.02,(f(xint) - df(xint)*xint) filled
set arrow from 1.02,0.5 to  1.02,(df(xint) + f(xint) - df(xint)*xint) filled
set arrow from xint,f(xint) to xint,-.7 nohead dashtype 2
set arrow from xint,.5*xint to xint,f(xint) filled

set label at -.24,0 "µ^@0_A=G_A"
set label at -.125,(f(xint) - df(xint)*xint) "µ^@{/Symbol a}_A"
set label at -.075,-.4 "+ RT ln(a^@{/Symbol a}_A)" rotate by 90
set label at 1.05,0.5 "µ^@0_B=G_B"
set label at 1.05,(df(xint) + f(xint) - df(xint)*xint) "µ^@{/Symbol a}_B"
set label at 1.075,-.15 "+ RT ln(a^@{/Symbol a}_B)" rotate by 90
set label at xint,-.1 " {/Symbol D}G_{mix}"
set label at .7,f(.7) " G^@{/Symbol a}_{AB}"
set label at .5,.5*.5+.075 "x_A G_A + x_B G_B" rotate by 30

plot f(x) ls 1 lw 2, .5*x ls 1 dashtype 2, (df(xint)*x + (f(xint) - df(xint)*xint)) ls 1 dashtype 4
