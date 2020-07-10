{smcl}
{* *! version 1.0.0 Jul 2020}{...}
{title:Title}

{pstd}
{hi:robttest} {hline 2} Robust t-test alternative for small samples


{title:Syntax}

{p 8 16 2}
{cmd:rddsga} {depvar} {it:assignvar} [{indepvars}] {ifin}
{cmd:,} {it:options}
{p_end}

{phang}
{it:depvar} is the outcome variable, {it:assignvar} is the assignment variable for which there is a known cutoff and {it:indepvars} are a set of control variables.
{p_end}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Subgroup and RDD}
{p2coldent:* {opth sg:roup(varname)}}subgroup indicator variable{p_end}
{synopt :{opth f:uzzy(varname)}}indicator for actual treatment status; if not specified, a sharp RDD is assumed{p_end}
{synopt :{opt c:utoff(real)}}specifies the cutoff value in {it:assignvar}; default is 0{p_end}
{p2coldent:* {opt bw:idth(real)}}specifies the bandwidth around the cutoff{p_end}
{synopt :{opt k:ernel(kernelfn)}}specifices the kernel function used to construct the local-polynomial estimator(s); default is {opt k:ernel(uniform)}

{syntab :Balance}
{p2coldent:+ {opth bal:ance(varlist)}}variables that enter the propensity score estimation; default is {indepvars}{p_end}
{synopt :{opt probit}}predict propensity score using a {manhelp probit R:probit} model; default is {manhelp logit R:logit}{p_end}
{synopt :{opt com:sup}} restrict sample to area of common propensity score support and a variable will be generated {p_end}

{syntab :Model}
{synopt :{opt first:stage}}estimate the discontinuity in the treatment probability using OLS{p_end}
{synopt :{opt reduced:form}}estimate the reduced form effect using OLS{p_end}
{synopt :{opt iv:regress}}estimate the treatment effect using instrumental variable regression; requires that a treatment variable is specified in {opt treatment(varname)}{p_end}
{synopt :{opt p:(int)}}specifies the order of the local polynomial; default is p(1) {p_end}
{synopt :{opth vce(vcetype)}}{it:vcetype} may be {opt un:adjusted},
   {opt r:obust}, {opt cl:uster} {it:clustvar}, {opt boot:strap},
   {opt jack:knife}, or {opt hac} {help ivregress##kernel:{it:kernel}}; default is {opt bootstrap}{p_end}
{synopt :{opt bsr:eps(#)}}perform # bootstrap replications; default is {opt bsreps(50)}{p_end}
{synopt :{opt norm:al}}compute and report normal approximation p-values and CIs; default is percentile{p_end}
{synopt :{opt noboot:strap}}do not compute bootstrap standard errors for RD estimates{p_end}
{synopt :{opt noipsw}}do not use inverse propensity score weighting{p_end}
{synopt :{opt fixed:bootstrap}}compute bootstrap with first stage fixed when instrumental variable regression is used. {p_end}
{synopt :{opt block:btrp(varlist)}}specifies bootstrap samples that are selected within each stratum. {p_end}
{synopt :{opt weights(weightsvar)}}is the variable used for optional weighting of the estimation procedure. The unit-specific weights multiply the kernel function. {p_end}


{syntab :Reporting and Output}
{synopt :{opt dibal:ance}}display original balance and propensity score-weighted balance tables and statistics{p_end}
{synopt :{opth ipsw:eight(newvar)}}name of new variable for the inverse propensity score weight for each observation; if not specified, no variable will be generated{p_end}
{synopt :{opth psc:ore(newvar)}}name of new variable for the estimated propensity score; if not specified, no variable will be generated{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* These options must be specified.{p_end}
{p 4 6 2}+ This option must be specified if {it:indepvars} is empty.{p_end}
{p 4 6 2}{it:indepvars} may contain factor variables; see {help fvvarlist}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
Standard inference about a scalar parameter estimated via GMM amounts to 
applying a t-test to a particular set of observations. If the number of 
observations is not very large, then moderately heavy tails can lead to poor
behaviour of the t-test. This package implements the method described in 
{help robttest##mainpaper:Müller (2020)}, which combines extreme value theory for the smallest and
largest observations with a normal approximation for the average of the 
remaining observations to construct a more robust alternative to the t-test.
This new test controls size more succesfully in small samples compared to 
existing methods. Analytical results in the canonical inference for the mean 
problem demonstrate that the new test provides a refinement over the full
sample t-test under more than two but less than three moments, while the
bootstrapped t-test does not.

{pstd}
{cmd:DESCRIBE SYNTAX} allows to conduct a binary subgroup analysis in RDD settings based on inverse propensity score weights (IPSW).
Observations in each subgroup are weighted by the inverse of their conditional probabilities to belong to that subgroup, given a set of covariates.
Analyzing the differential treatment effect in the reweighted sample helps isolating the difference due to the subgroup characteristic of interest from other observable dimensions.

{pstd}
{cmd:DESCRIBE SYNTAX} computes IPSW based on the vector of covariates in {it:indepvars}, which are also used as control variables in the model.
A separate set of variables to compute IPSW may be optionally specified with {opt balance(varlist)}.
A new variable with IPSW may be generated using {opt psweight(newvar)}.
In order to assess the statistical significance of the difference in means for each covariate, {cmd:rddsga} uses a t-test on the equality of means and reports the resulting p-value, as well as the (weighted) standardized mean difference.
Joint significance is assessed using an F-test.
The resulting balance tables are stored as matrices (see {help rddsga##results:stored results} below), and may also be displayed using {opt dibalance}.

{pstd}
Options {opt firststage}, {opt reducedform} and {opt ivregress} may be used to estimate the first stage, reduced form and treatment effects using OLS and IV methods, respectively.
Standard variance estimator options may be used for these estimators with {opt vce(vcetype)}.
The estimated coefficients for the interaction of indicator variables for each subgroup and a treatment indicator are reported, along with bootstrap standard errors (unless other {it:vcetype} is specified).

{pstd}
Additional details regarding the methodology implemented by {cmd: rddsga} can be found in the project's {browse "https://gitlab.com/acarril/rddsga/wikis/home":repository}.
An application can be found in {help rddsga##mainpaper:Gerardino, Litschig, Olken and Pomeranz (2017)}.


{marker options}{...}
{title:Options}

{dlgtab:RD design}

{phang}
{opt sgroup(varname)} specifies a subgroup indicator variable.
This variable must be a {it:dummy} (values 0 or 1).
This option must be specified.

{phang}
{opt fuzzy(varname)} specifies the treatment status variable used to implement fuzzy RD estimation. Default is Sharp RD design and hence this option is not used.

{phang}
{opt cutoff(real)} specifies the cutoff value in {it:assignvar}; default is 0 (assuming normalized {it:assignvar}).

{phang}
{opt bwidth(real)} specifies a symmetrical the bandwidth around the cutoff.
This option must be specified.

{phang}
{opt kernel(kernelfn)} specifies the kernel function used to construct the
        local-polynomial estimator(s). Options are: triangular, epanechnikov, and
        uniform.  Default is {opt kernel(uniform)}.


{dlgtab:Balance}

{phang}
{opt balance(varlist)} specifies the variables that enter the propensity score estimation.
If not specified, variables in {it:indepvars} are used.
This option is useful if one wants to balance a different set of covariates than the ones used as controls in the model.
This option must be specified if {it:indepvars} is empty and {it:ipsw} is used.

{phang}
{opt probit} indicates that the propensity score is computed after a fitting a  {manhelp probit R:probit} model; default is {manhelp logit R:logit}.

{phang}
{opt comsup} indicates that the sample is to be restricted to the area of common support and the respective variable will be generated.

{marker options_model}{...}
{dlgtab:Model}

{phang}
{opt firststage} estimates the discontinuity in the treatment probability using OLS.

{phang}
{opt reducedform} estimates the reduced form effect using OLS.

{phang}
{opt ivregress} estimates the treatment effect using instrumental variable regression.
If specified, it requires that a treatment variable is also specified in {opt treatment(varname)}.
See {help ivregress} for additional information.

{phang}
{opt p(int)} specifies the order of the local polynomial used to construct the point estimator.  Default is p(1) (local linear regression).

{phang}
{opt vce(vcetype)} specifies the variance estimators options.
{it:vcetype} may be unadjusted, robust, cluster clustvar, bootstrap, jackknife, or hac kernel (see {help  vce_option}).
Bootstrap is used by default, unless {opt nobootstrap} is used (see below), or a different {it:vcetype} is specified.

{pmore}
{cmd:rddsga} computes the bootstrap estimate of the variance-covariance matrix by taking a random sample (with replacement) of size {it:N} from the observations and computing the estimate as usual.
This process is repeated the number of times specied in {opt bsreps(#)} (see below), which yields a bootstrap distribution of estimates.
See {it:Methods and formulas} in {manlink R bootstrap} for further details.

{pmore}
Confidence intervals and p-values are estimated using the percentile method by default. See the {opt normal} option below for further details and other options.

{phang}
{opt bsreps(#)} perform # bootstrap replications; default is {opt bsreps(50)}.

{phang}
{opt normal} computes and report normal-based p-values and confidence intervals.
By default {cmd:rddsga} uses the so-called percentile method, which means we compute the empirical p-value, i.e. the proportion of bootstrap estimates further in the tails than the (absolute value of the) original estimate.
The percentile method yields nonparametric confidence intervals that correspond to the {it:p}th quantile of the bootstrap distribution.
See {it:Methods and formulas} in {manlink R bootstrap} for further details.

{pmore}
If {opt normal} is specified, p-values are obtained using the normality assumption, which yields confidence intervals that correspond to the {it:p}th quantile of the standard normal distribution.

{phang}
{opt nobootstrap} prevents computing bootstrap standard errors, which is the default.

{phang}
{opt noipsw} prevents employing inverse propensity score weighting for the estimations.

{phang}
{opt fixedbootstrap} if RDD is fuzzy, bootstrap is computed only in the reduced form and keeping the first stage is fixed; default is bootstrap computed in both first stage and reduced form. 

{phang}
{opt blockbtrp(varlist)} specifies the variables identifying strata. If blockbtrp() is specified, bootstrap samples are selected within each stratum.

{marker options_reporting}{...}
{dlgtab:Reporting}

{phang}
{opt dibalance} display original balance and propensity score weighting balance tables and statistics.
This balance is computed for each covariate in {it:indepvars}, unless {opt balance(varlist)} is specified.

{phang}
{opt psweight(newvar)} specifies a name for a new variable with the propensity score weighting. If not specified, no variable will be generated.

{phang}
{opt pscore(newvar)} specifies a name for a new variable with the propensity score. If not specified, no variable will be generated.


{marker examples}{...}
{title:Examples}

{pstd}Setup (click {browse "https://github.com/acarril/rddsga#installation":here} for details on getting ancillary files){p_end}
{phang2}{cmd:. use rddsga_synth}{p_end}

{pstd}Assess covariate imbalance using one covariate{p_end}
{phang2}{cmd:. rddsga Y Z, balance(X1) sgroup(G) bwidth(10) dibal}{p_end}

{pstd}Silently store computed IPSW based on X1 and X2{p_end}
{phang2}{cmd:. rddsga Y Z, balance(X1 X2) sgroup(G) bwidth(10) ipsweight(ipsw)}{p_end}

{pstd}Fit reduced form model{p_end}
{phang2}{cmd:. rddsga Y Z, balance(X1 X2) sgroup(G) bwidth(10) reduced}{p_end}

{pstd}Estimate treatment effect in fuzzy RDD using instrumental variables regression and 200 bootstrap replications, comparing output without and with IPSW{p_end}
{phang2}{cmd:. rddsga Y Z X1 X2, sgroup(G) bwidth(6) ivreg bsreps(100) treatment(T) noipsw}{p_end}
{phang2}{cmd:. rddsga Y Z X1 X2, sgroup(G) bwidth(6) ivreg bsreps(100) treatment(T)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:rddsga} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(ipsw_N_G0)}}number of observations in subgroup 0 (IPSW){p_end}
{synopt:{cmd:e(ipsw_N_G1)}}number of observations in subgroup 1 (IPSW){p_end}
{synopt:{cmd:e(ipsw_Fstat)}}F-statistic (IPSW){p_end}
{synopt:{cmd:e(ipsw_pvalue)}}F-statistic p-value (IPSW){p_end}
{synopt:{cmd:e(ipsw_avgdiff)}}Average of absolute values of standardized differences (IPSW){p_end}

{synopt:{cmd:e(unw_N_G0)}}number of observations in subgroup 0 (unweighted){p_end}
{synopt:{cmd:e(unw_N_G1)}}number of observations in subgroup 1 (unweighted){p_end}
{synopt:{cmd:e(unw_Fstat)}}F-statistic (unweighted){p_end}
{synopt:{cmd:e(unw_pvalue)}}F-statistic p-value (unweighted){p_end}
{synopt:{cmd:e(unw_avgdiff)}}Average of absolute values of standardized differences (unweighted){p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}subgroup estimators coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synopt:{cmd:e(ipsw)}}balance table matrix (IPSW){p_end}
{synopt:{cmd:e(unw)}}balance table matrix (unweighted){p_end}
{p2colreset}{...}

{pstd}
Additionally, {cmd:rddsga} stores all macros and scalars of the estimated model ({opt firststage}, {opt firststage} or {opt ivregress}; see {help rddsga##options_model:Model options} above).
Refer to the {help regress##results:Stored Results} section of {manhelp regress R:regress} for the full list of macros and scalars additionally stored in {cmd:e()}.


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
If you have suggestions or want to report problems, please create a new issue in the {browse "https://github.com/acarril/rddsga/issues":project repository} or contact the project maintainer.
All remaining errors are our own.


{marker references}{...}
{title:References}

{marker mainpaper}{...}
{phang}Müller, Ulrich K. "A More Robust t-Test" Working Paper. July 2020. {browse "https://www.princeton.edu/~umueller/heavymean.pdf"}.

