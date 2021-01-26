clear all
sysuse auto, clear
gen s_1=rnormal(0,1)
gen s_2=rnormal(0,1)
reg mpg weight length, robust
scpc 
