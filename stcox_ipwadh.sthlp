{smcl}
{* *! version 0.1.0  10oct2026}{...}
{vieweralsosee "ipwadh" "help ipwadh"}{...}
{vieweralsosee "ipwadh_eval" "help ipwadh_eval"}{...}
{vieweralsosee "[ST] stcox" "help stcox"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:stcox_ipwadh} {hline 2}}Optional: standard errors of the weighted Cox model that account for the estimation of the weights{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:stcox_ipwadh} {varlist}[{cmd:,} {opth str:ata(varlist)} {opt nohr} {opt l:evel(#)} {opt eval:at(end|mid|start)} {opt pre:fix(name)} {opt gen:erate(name)} {opth cl:uster(varname)} {opt tol:erance(#)} {opt det:ail}]

{pstd}{it:varlist} may contain factor variables.{p_end}


{title:Description}

{pstd}
Runs {cmd:stcox} {it:varlist}{cmd:, breslow vce(cluster} {it:id}{cmd:)} on data that were {cmd:stset} with the {cmd:ipwadh} weights and replaces the variance by that of the stacked
estimating equations: the weighted Cox score plus the scores of all Poisson rate models that produced the weights (adherence, death, loss to follow-up; numerator and denominator models; all strata).
The coefficients are those of {cmd:stcox}.

{pstd}
By default the output is the {cmd:stcox} table with the adjusted standard errors and one line with the ratios adjusted / cluster-robust; {opt detail} adds the comparison table and the reproduction checks.

{pstd}
{bf:Use it as a supplement.} The plain cluster-robust standard error (weights treated as known) remains the recommended default. In simulations the correction was small (about +/-4%) and
went in both directions; it did not repair the under-coverage that heavy-tailed weights cause. {cmd:stcox_ipwadh} is useful for a sensitivity analysis and to show that the conclusion does not depend on the variance estimator.

{pstd}
{bf:Requirements} (checked):

{phang2}1. {cmd:ipwadhc} (version 0.3.1 or later) was run with {opt parts} on the unsplit episodes; it stores the model specification in the data. Rows and covariates must not be changed afterwards.{p_end}
{phang2}2. {cmd:ipwadh_split}, then {cmd:ipwadh_eval} (same {opt evalat()}).{p_end}
{phang2}3. {cmd:stset t1 [pw=_ipw_psw} (or {cmd:_ipw_pswt}{cmd:] if adh==1, time0(t0) failure(...)}: pweights, no {opt id()}, no {opt origin()}/{opt scale()}.{p_end}
{phang2}4. The model contains the arm and exactly the {opt stcoxvars()} given to {cmd:ipwadhc}; {opt strata()} is allowed. Not allowed: {cmd:efron}/{cmd:exactp} ties, {opt tvc()}, {opt shared()}, {opt offset()}.{p_end}

{pstd}
The rate models are refitted from the stored specification. The recomputed weights must equal the stored weights (relative difference below {opt tolerance()}, default 1e-6) and the recomputed cluster-robust variance must equal
{cmd:e(V)} of {cmd:stcox}; otherwise the command stops.

{pstd}
Truncated pieces ({cmd:ipwadh_eval, truncate()}) have no weight derivative: the percentile bounds are treated as fixed. Truncation therefore also removes part of the correction.


{title:Stored results}

{pstd}
All {cmd:stcox} results, with {cmd:e(V)} replaced by the adjusted variance, and

{synoptset 28 tabbed}{...}
{synopt:{cmd:e(V_naive)}}cluster-robust variance of {cmd:stcox}{p_end}
{synopt:{cmd:e(ipwadh_adjusted)}}1{p_end}
{synopt:{cmd:e(ipwadh_maxdiff_parts)}}reproduction check: largest difference of the recomputed weight parts{p_end}
{synopt:{cmd:e(ipwadh_maxdiff_w)}}reproduction check: largest relative difference of the recomputed weights{p_end}
{synopt:{cmd:e(ipwadh_maxdiff_V)}}reproduction check: largest difference of the recomputed cluster-robust variance{p_end}
{synopt:{cmd:e(ipwadh_maxscore)}}largest absolute score at the estimate{p_end}
{synopt:{cmd:e(ipwadh_nclipped)}}pieces with truncated weights in the Cox sample{p_end}
{synopt:{cmd:e(ipwadh_vce)}}{cmd:adjusted}{p_end}
{p2colreset}{...}


{marker ai}{...}
{title:AI disclosure}

{pstd}
This command was developed using substantial work from Claude AI (Anthropic), with human supervision.
{p_end}

{marker author}{...}
{title:Author}

{pstd}
Barna Szabó-Söderberg{break}
{it:Karolinska Institutet}{break}
{ul:barna.szabo@ki.se}
{p_end}
