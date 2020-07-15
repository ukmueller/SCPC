{smcl}
{* *! version 1.0.0 Jul 2020}{...}
{title:Title}

{pstd}
{hi:robttest} {hline 2} Robust t-test for potentially heavy-tailed observations/scores


{title:Syntax}

{p 8 16 2}
{cmd:robttest} {it:modelspec} [{cmd:,} {it:options}]
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
{synopt :{opt k(#)}}override default value for {it:k}th order statistic{p_end}
{synopt :{opt v:erbose}}display the 2{it:k} extreme terms, in multiples of standard deviation of sum of middle {it:n}-2{it:k} terms{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
Standard inference about a scalar parameter estimated via GMM amounts to 
applying a t-test to a particular set of observations. If the number of 
observations is not very large, then moderately heavy tails can lead to poor
behaviour of the t-test.
This package implements the method described in 
{help robttest##mainpaper:Müller (2020)}, which combines extreme value theory for the smallest and
largest observations with a normal approximation for the average of the 
remaining observations to construct a more robust alternative to the t-test.
This new test controls size more succesfully in small samples compared to 
existing methods.
It is implemented as a postestimation command that can be used after the Stata commands {manhelp regress R:regress}, {manhelp ivregress R:ivregress}, {manhelp areg R:areg}, {manhelp logit R:logit} or {manhelp probit R:probit},
as long as these are used with the standard error option {it:robust} or {it:cluster} (see {manhelp vce_option R:vce option}).
P-values take values on the grid 0.002,0.004,0.006,0.008,0.01,0.02,0.030,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.12,0.14,0.16,0.18,0.2,0.250,0.3,0.4,0.5,1.0.
The reported p-value is the largest value on this grid for which all tests of level as small or smaller reject. Thus, the reported p-value is to be interpreted as an upper bound, that is if it was computed with a finer grid of levels, it would be smaller than the one that is reported, but not as small as the next smallest value in this grid. In particular, a value of 0.002 means that tests of all considered levels reject, and a p-value of 1.0 means that the test of level 0.5 does not reject.


{marker options}{...}
{title:Options}

{phang}
{opt k(#)} overrides default value for {it:k}th order statistic; see {help robttest##mainpaper:Müller (2020)} for details.

{phang}
{opt verbose} displays the 2{it:k} extreme terms, in multiples of standard deviation of sum of middle {it:n}-2{it:k} terms.


{marker examples}{...}
{title:Examples}

{phang}{cmd:. sysuse auto}{p_end}
{phang}{cmd:. regress mpg weight length, robust}{p_end}
{phang}{cmd:. estimates store A}{p_end}
{phang}{cmd:. robttest}{p_end}
{phang}{cmd:. robttest .}{space 6}(equivalent to above command){p_end}
{phang}{cmd:. robttest A}{space 6}(equivalent to above command){p_end}
{phang}{cmd:. robttest, k(6) verbose}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:robttest} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(k)}}value of {it:k} used in computations{p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(rob_pvals)}}robust p-values{p_end}
{synopt:{cmd:e(rob_CIs)}}robust confidence intervals{p_end}
{synopt:{cmd:e(rob_WL)}}{it:k} extreme values (left tail){p_end}
{synopt:{cmd:e(rob_WR)}}{it:k} extreme values (right tail){p_end}
{p2colreset}{...}

{pstd}
Additionally, {cmd:robttest} preserves all macros and scalars of the estimated model in memory.


{marker authors}{...}
{title:Authors}

{pstd}
Álvaro Carril (maintainer){break}
Princeton University{break}
acarril@princeton.edu

{pstd}
Ulrich K. Müller{break}
Princeton University{break}


{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
This software is provided "as is", without warranty of any kind.
If you have suggestions or want to report problems, please create a new issue in the {browse "https://github.com/acarril/robttest/issues":project repository} or contact the project maintainer.
All remaining errors are our own.


{marker references}{...}
{title:References}

{marker mainpaper}{...}
{phang}Müller, Ulrich K. "A More Robust t-Test" Working Paper. July 2020. {browse "https://www.princeton.edu/~umueller/heavymean.pdf"}.

