set term pngcairo enhanced font 'Helvetica,12' size 600,700 dashed
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

set lmargin 10
set rmargin 10

set linestyle 1 black

set arrow from -.02,0 to  -.02,(f(xint) - df(xint)*xint) filled
set arrow from 1.02,0.5 to  1.02,(df(xint) + f(xint) - df(xint)*xint) filled
set arrow from xint,f(xint) to xint,-.7 nohead dashtype 2
set arrow from xint,.5*xint to xint,f(xint) filled

set label at -.15,0 "µ^@0_B=G_B" font ",16"
set label at -.1,(f(xint) - df(xint)*xint) "µ^@{/Symbol a}_B" font ",16"
set label at -.075,-.4 "+ RT ln(a^@{/Symbol a}_B)" rotate by 90 font ",16"
set label at 1.05,0.5 "µ^@0_A=G_A" font ",16"
set label at 1.05,(df(xint) + f(xint) - df(xint)*xint) "µ^@{/Symbol a}_A" font ",16"
set label at 1.075,-.15 "+ RT ln(a^@{/Symbol a}_A)" rotate by 90 font ",16"
set label at xint,-.1 " {/Symbol D}G_{mix}" font ",16"
set label at .7,f(.7) " G^@{/Symbol a}_{AB}" font ",16"
set label at .5,.5*.5+.075 "x_A G_A + x_B G_B" rotate by 30 font ",16"

set label at graph 0,-.07 "B" center font ",16"
set label at graph 1,-.07 "A" center font ",16"

plot f(x) ls 1 lw 2, .5*x ls 1 dashtype 2, (df(xint)*x + (f(xint) - df(xint)*xint)) ls 1 dashtype 4
