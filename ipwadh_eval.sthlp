{smcl}
{* *! version 0.3.0  10oct2026}{...}
{vieweralsosee "ipwadh" "help ipwadh"}{...}
{vieweralsosee "ipwadhc" "help ipwadhc"}{...}
{vieweralsosee "ipwadh_split" "help ipwadh_split"}{...}
{vieweralsosee "stcox_ipwadh" "help stcox_ipwadh"}{...}
{title:Title}

{p2colset 5 21 23 2}{...}
{p2col :{cmd:ipwadh_eval} {hline 2}}Weights of the pieces, diagnostics, truncation and weight warning{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:ipwadh_eval} {ifin}{cmd:,} {opth time0(varname)} {opth time1(varname)} [{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opth adh:erence(varname)}}restrict diagnostics and truncation to adherent pieces{p_end}
{synopt :{opth arm(varname)}}diagnostics by arm as well{p_end}
{synopt :{opt trunc:ate(# #)}}clip at person-time-weighted percentiles, for example {cmd:truncate(1 99)}; creates {cmd:_ipw_pswt}{p_end}
{synopt :{opt eval:at(end|mid|start)}}point of the piece; default {cmd:end} (exact for the exact split){p_end}
{synopt :{opt esscaution(#)}}ESS percent below which a caution is printed; default 20{p_end}
{synopt :{opt esssevere(#)}}ESS percent below which a severe warning is printed; default 5{p_end}
{synopt :{opt maxwarn(#)}}largest weight above which a note is printed; default 100{p_end}
{synopt :{opt nowarn}}suppress the warning text (the values are still returned){p_end}
{synopt :{opt det:ail}}full diagnostics table (by arm too) and the analysis hints (default: three summary lines){p_end}
{synopt :{opt pre:fix(name)}}prefix of the variables from {cmd:ipwadhc}/{cmd:ipwadh_split}; default {cmd:_ipw}{p_end}
{synopt :{opt gen:erate(name)}}prefix of the new variables; default = prefix{p_end}
{synopt :{opt replace}}replace existing weight variables{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
Requires {cmd:ipwadhc, parts} and {cmd:ipwadh_split}. Computes for every piece {cmd:log W(t) = E0 - dA (t - s) + J0 + Ec0 - dC (t - s)} at the chosen point t, where s is the start
of the original row. Creates {cmd:_ipw_pswa} (adherence), {cmd:_ipw_pswc} (censoring), {cmd:_ipw_psw} (product) and, with {opt truncate()}, {cmd:_ipw_pswt}. Zero-length
rows have missing weights.

{pstd}
Diagnostics are person-time weighted (each piece counts with its length): mean, minimum, p1, median, p99, maximum and the effective sample size
ESS = (sum w t)^2 / (sum t sum w^2 t), as a share of the person-time, overall and by arm. Truncation bounds are the person-time weighted percentiles, so they do not depend on how finely the rows were split.

{pstd}
{bf:Weight warning.} The smallest ESS over the overall and by-arm values (untruncated weight) is compared with {opt esscaution()} and {opt esssevere()}. Simulations: ESS 40-55% was stable, about 19% acceptable after truncation,
6% and 1% gave unstable estimates (SD about 3 times the unweighted SD) with too small standard errors. The thresholds are heuristics, not tests, and cannot detect unmeasured confounding.
When a warning appears, simplify the adherence model, use truncation and report unweighted, weighted untruncated and weighted truncated estimates side by side.

{pstd}
Then analyse with {cmd:stset} {it:time1} {cmd:[pw=_ipw_psw} (or {cmd:_ipw_pswt}{cmd:] if} {it:adherence}{cmd:==1, time0(}{it:time0}{cmd:) failure(}{it:failure}{cmd:)} (no {opt id()}) and {cmd:stcox ..., breslow vce(cluster} {it:id}{cmd:)}.


{title:Stored results}

{pstd}
{cmd:r(N)}, {cmd:r(n_missing)}, {cmd:r(n_used)}, {cmd:r(pt)}, {cmd:r(mean)}, {cmd:r(min)}, {cmd:r(p1)}, {cmd:r(p50)}, {cmd:r(p99)}, {cmd:r(max)}, {cmd:r(ess)}, {cmd:r(ess_min)}, {cmd:r(warn)}, {cmd:r(warn_max)};
with {opt arm()}: {cmd:r(n_levels)}, {cmd:r(level_#)}, {cmd:r(ess_#)} and others; with {opt truncate()}: {cmd:r(trunc_lo)}, {cmd:r(trunc_hi)}, {cmd:r(share_lo)}, {cmd:r(share_hi)}, {cmd:r(ess_trunc)};
macros {cmd:r(evalat)} and {cmd:r(weightvar)} (the variable to use in {cmd:stset}).


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
