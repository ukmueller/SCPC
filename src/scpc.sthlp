{smcl}
{* *! version 1.0.0 Jul 2020}{...}
{title:Title}

{pstd}
{hi:scpc} {hline 2} Spatial Correlation Robust Inference via SCPC.


{title:Syntax}

{p 8 16 2}
{cmd:scpc} {it:modelspec} [{cmd:,} {it:options}]
{p_end}

{phang}
{it:modelspec} specifies a model estimated by {manhelp regress R:regress},
{manhelp ivregress R:ivregress}, {manhelp areg R:areg}, {manhelp logit R:logit},
or {manhelp probit R:probit}, estimated with the standard error option {it:robust} or {it:cluster} (see {manhelp vce_option R:vce option}).
{it:modelspec} is

                {it:name}{c |}{cmd:.}{c |} {cmd:(}{it:namelist}{cmd:)}

{pmore}
{it:name} is the name under which estimation results were stored using
{helpb estimates store:estimates store}, and "{cmd:.}" refers to the last
estimation results, whether or not these were already stored.
If {it:modelspec}
is not specified, the last estimation result is used; this is equivalent to
specifying {it:modelspec} as "{cmd:.}".
{p_end}


{synoptset 11}{...}
{synopthdr}
{synoptline}
{synopt :{opt avc(#)}}overrides default value of 0.05 for maximal average pairwise correlation {p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
This Stata package implements the Spatial Correlation Principal Components (SCPC) method described in {help SCPC##mainpaper:Müller and Watson (2021)} for the construction of confidence intervals that account for many forms of spatial correlation. 
The scpc command expects the locations of the observations to be stored in the variables s_*. For instance,
if s_1 and s_2 are the only variables whose name begins with "s_", then the method uses 2-dimensional locations.
It is implemented as a postestimation command that can be used after the Stata commands {manhelp regress R:regress}, {manhelp ivregress R:ivregress}, {manhelp areg R:areg}, {manhelp logit R:logit} or {manhelp probit R:probit},
as long as these are used with the standard error option {it:robust} or {it:cluster} (see {manhelp vce_option R:vce option}). If the estimation uses the {it:cluster} option, then scpc corrects for spatial correlations between clusters, assuming that all observations within a cluster share the same location.
 

{marker options}{...}
{title:Options}

{phang}
{opt k(#)} overrides default value of 0.05 for maximal average pairwise correlation; see {help scpc##mainpaper:Müller and Watson (2021)} for details.



{marker examples}{...}
{title:Examples}

{phang}{cmd:. sysuse auto}{p_end}
{phang}{cmd:. gen s_1=rnormal(0,1)}{p_end}
{phang}{cmd:. gen s_2=rnormal(0,1)}{p_end}
{phang}{cmd:. regress mpg weight length, robust}{p_end}
{phang}{cmd:. estimates store A}{p_end}
{phang}{cmd:. scpc}{p_end}
{phang}{cmd:. scpc .}{space 6}(equivalent to above command){p_end}
{phang}{cmd:. scpc A}{space 6}(equivalent to above command){p_end}
{phang}{cmd:. scpc, avc(0.01)}{p_end}


{marker results}{...}
{title:Stored results}

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
{phang}Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference" Working Paper. February 2021. {browse "https://www.princeton.edu/~umueller/SHAR.pdf"}.

