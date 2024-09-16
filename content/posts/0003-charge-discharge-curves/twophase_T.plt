set term pngcairo enhanced font 'Helvetica,20' size 600,1000 dashed
set output "twophase_T.png"

set multiplot layout 3,1

set samples 500
unset key

set xrange [0:1]

unset ytics
set format x ""

dGt(T) = .7-.25*T
falpha(x,T) = .5*x + T*(x*log(x)+(1-x)*log(1-x))
fbeta(x,T) = dGt(T) - .5*x + T*(x*log(x)+(1-x)*log(1-x))

set lmargin 7
set rmargin 7

set tmargin 3
set bmargin 0

set linestyle 1 black lw 2

set label 1 at .2,falpha(.2,10) "  G^{/Symbol a}"
set label 2 at .8,fbeta(.8,10) "  G^{/Symbol b}" tc "red"

T=10

set label 3 at -.15,0 "µ^@{0,{/Symbol a}}_A"
set label 4 at 1.03,.5 "µ^@{0,{/Symbol a}}_B"
set label 5 at -.15,dGt(T) "µ^@{0,{/Symbol b}}_A" tc "red"
set label 6 at 1.03,-.5+dGt(T) "µ^@{0,{/Symbol b}}_B" tc "red"
set label 9 at graph .5,.8 "T_1"

plot falpha(x,T) ls 1, fbeta(x,T) ls 1 lc "red"

unset label 1
unset label 2

T=1
set yrange[-.65:.5]

set bmargin 1.5
set tmargin 1.5

set label 5 at -.15,dGt(T)
set label 6 at 1.03,-.5+dGt(T)
set label 9 "T_2"

set arrow 3 from .2,falpha(.3,T)+.01 to .7,fbeta(.6,T)-.03 nohead

set label 10 at .35,falpha(.35,T) "•\nP_1" center
set label 11 at .55,fbeta(.55,T) "•\nQ_1" center

plot falpha(x,T) ls 1, fbeta(x,T) ls 1 lc "red"

T=0.5
set xlabel "x_A"
set bmargin 3
set tmargin 0
set format x "%.1f"
set yrange[-.2:.6]

set arrow 1 from -.02,0 to -.02,dGt(T) filled
set arrow 2 from 1.02,.5 to 1.02,-.5+dGt(T) filled
set label 5 at -.15,dGt(T)
set label 6 at 1.03,-.5+dGt(T)

set arrow 3 from .2,falpha(.3,T)-.02 to .9,fbeta(.8,T)+.01 nohead

set label 10 at .3,falpha(.3,T) "•\nP_2" center
set label 11 at .8,fbeta(.8,T) "•\nQ_2" center


set label 7 at -.075,.15 "{/Symbol D}G^@{{/Symbol a}→{/Symbol b}}_A" rotate by 90
set label 8 at 1.075,.42 "{/Symbol D}G^@{{/Symbol a}→{/Symbol b}}_B" rotate by -90
set label 9 "T_3"

plot falpha(x,T) ls 1, fbeta(x,T) ls 1 lc "red"
