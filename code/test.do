// Set working directory
if "`c(username)'" == "alvaro" cd ~/Repos/heavymean
else cd "c:/dropbox/mystuff/heavymean/2020/MATA"

do code/heavymeanmata.do
do code/robttest.ado
sysuse auto, clear
reg mpg weight length, robust
robttest
