// Set working directory
if "`c(username)'" == "alvaro" cd ~/Repos/robttest
else cd "c:/dropbox/mystuff/heavymean/2020/MATA"

clear all
do code/heavymeanmata.do
do robttest.ado
sysuse auto, clear
//  drop if _n > 20
reg mpg weight length, robust
robttest
