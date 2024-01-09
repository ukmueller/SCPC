{smcl}
{* *! version 1.1.0 Jan 2022}{...}
{title:Title}

{pstd}
{hi:scpc} {hline 2} Spatial Correlation Robust Inference via SCPC.


{title:Syntax}

{p 8 16 2}
{cmd:scpc} [{cmd:,} {it:options}]
{p_end}

{synoptset 11}{...}
{synopthdr}
{synoptline}
{synopt :{opt avc(#)}}overrides default value of 0.03 for the maximal average pairwise correlation (must be between 0.001 and 0.99){p_end}
{synopt :{opt uncond}}overrides default of computing critical values conditional on regressors to compute "unconditional" SCPC inference as described in {help scpc##mainpaper:Müller and Watson (2021)}{p_end}
{synopt :{opt latlong}}if present, expects s_1 to contain latitude and s_2 to be longitude and uses implied great-circle distance; if not present, computes Euclidian distance between coordinates defines by s_* variables{p_end}
{synopt :{opt k(#)}}overrides default value of 10 for the number of coefficients SCPC t-statistics are computed for{p_end}
{synopt :{opt cvs}}prints a table of two-sided critical values of SCPC t-statistic of level 32%, 10%, 5% and 1%{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
This Stata package implements the Spatial Correlation Principal Components (SCPC) method described in {help scpc##mainpaper:Müller and Watson (2022, 2023)} for the construction of confidence intervals that account for many forms of spatial correlation.
The {it:scpc} command expects the locations of the observations to be stored in the variables s_*. If the option {it:latlong} is present, then s_1 is interpreted as the latitude, and s_2 as the longitude of the observation, and distances between observations are computed with the great-circle formula. If the option {it:latlong} is not present, then all {it:p} variables beginning with s_ are interpreted as points in {it:p}-dimensional space, and distances are computed with the Euclidian norm. 
scpc is implemented as a postestimation command that can be used after the Stata commands {manhelp regress R:regress}, {manhelp ivregress R:ivregress}, {manhelp areg R:areg}, {manhelp logit R:logit} or {manhelp probit R:probit},
as long as these are used with the standard error option {it:robust} or {it:cluster} (see {manhelp vce_option R:vce option}). If the estimation uses the {it:cluster} option, 
then {it:scpc} corrects for spatial correlations between clusters, assuming that all observations within a cluster share the same location.
Also by default, scpc computes conditional critical values described in Müller and Watson (2023). These critical values are adjusted based on the realized values of the regressors, and is implemented only for {manhelp regress R:regress} and {manhelp ivregress R:ivregress}. Standard (unconditional) SCPC inference can be obtained via the {it:uncond} option, which is necessary when using scpc after {manhelp areg R:areg}, {manhelp logit R:logit} or {manhelp probit R:probit}.
By default, {it:scpc} conducts spatially robust inference under the assumption that the largest average pairwise correlation between the observations / clusters is no larger than 0.03. This default can be overridden by the
{it:avc} option. Note that computation times increase for smaller values of {it:avc}. To make the algorithm faster, SCPC inference is computed only for the first 10 coefficients. This can be changed by the {it:k} option. 
The underlying algorithm can also handle large datasets; internally, a different algorithm is used when the number of observations / clusters exceeds 4500. Computing time is appproximately linear in the number of observations, and is roughtly one minute for 5000 observations.
Note: The implementation generates and deletes variables with names starting with scpc_*, so variables of that name in the original data set are corrupted/deleted.
 

{marker options}{...}
{title:Options}

{phang}
{opt avc(#)} overrides default value of 0.03 for the maximal average pairwise correlation; see {help scpc##mainpaper:Müller and Watson (2021)} for details.

{phang}
{opt avc(#)} overrides default of computing conditional critical value in regression; see {help scpc##mainpaper:Müller and Watson (2021)} for details.

{phang}
{opt latlong} if present, expects s_1 to contain latitude and s_2 to be longitude and uses implied great-circle distance; if not present, computes Euclidian distance between coordinates defines by variables beginning with "s_".

{phang} 
{opt cvs} prints a table of one- and two-sided critical values of SCPC t-statistic of level 32%, 10%, 5% and 1%.

{phang} 
{opt k(#)} overrides default value of 10 for the number of coefficients for which SCPC inference is computed.

{marker examples}{...}
{title:Examples}

{phang}{cmd:. sysuse auto}{p_end}
{phang}{cmd:. gen s_1=rnormal(0,1)}{p_end}
{phang}{cmd:. gen s_2=rnormal(0,1)}{p_end}
{phang}{cmd:. regress mpg weight length, robust}{p_end}
{phang}{cmd:. scpc}{p_end}
{phang}{cmd:. scpc ,avc(0.03)}{space 6}(equivalent to above command){p_end}
{phang}{cmd:. scpc ,k(4)}(equivalent to above command){p_end}
{phang}{cmd:. scpc ,avc(0.01) cvs k(1)}{p_end}
{phang}{cmd:. gen clust=round(rnormal(0,10),1)}{p_end}
{phang}{cmd:. regress mpg weight length, cluster(clust)}{p_end}
{phang}{cmd:. scpc}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:scpc} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(scpcstats)}}Matrix of SCPC inference results (same as printed table){p_end}
{synopt:{cmd:e(scpccvs)}}Matrix of SCPC critical values (same as printed table, only computed under {it:cvs} option){p_end}
{p2colreset}{...}

{pstd}
Additionally, {cmd:scpc} preserves all macros and scalars of the estimated model in memory.



{marker authors}{...}
{title:Author}

{pstd}
Ulrich K. Müller{break}
Princeton University{break}


{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
This software is provided "as is", without warranty of any kind.
If you have suggestions or want to report problems, please create a new issue in the {browse "https://github.com/ukmueller/SCPC/issues":project repository} or contact the project maintainer.


{marker references}{...}
{title:References}

{marker mainpaper}{...}
{phang}Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference", Econometrica 90 (2022), 2901–2935. {browse "https://www.princeton.edu/~umueller/SHAR.pdf"}.

{marker followpaper}{...}
{phang}Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference in Linear Regression and Panel Models", Journal of Business and Economic Statistics 41 (2023), 1050–1064. {browse "https://www.princeton.edu/~umueller/SpatialRegression.pdf"}.


