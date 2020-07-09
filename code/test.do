// Set working directory
if "`c(username)'" == "alvaro" cd ~/Repos/robttest
else cd "c:/dropbox/mystuff/heavymean/2020/MATA"

clear all
sysuse auto, clear
reg mpg weight length, robust
robttest
