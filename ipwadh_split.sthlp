{smcl}
{* *! version 0.2.1  10oct2026}{...}
{vieweralsosee "ipwadh" "help ipwadh"}{...}
{vieweralsosee "ipwadhc" "help ipwadhc"}{...}
{vieweralsosee "[ST] stsplit" "help stsplit"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:ipwadh_split} {hline 2}}Exact split of adherent episodes at the failure times{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:ipwadh_split} {ifin}{cmd:,} {opth id(varname)} {opth time0(varname)} {opth time1(varname)} {opth adh:erence(varname)} {opth fail:ure(varname)} [{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opth str:ata(varlist)}}cut each row only at failure times of its own stratum (strata constant within a row){p_end}
{synopt :{opth last:only(varlist)}}indicators (for example {cmd:dead ltfu}) that must be 1 only on the last piece of a row; {opt failure()} is always included{p_end}
{synopt :{opt band:points(# | auto)}}approximate bands for very large data; default 1 = exact{p_end}
{synopt :{opt maxgb(#)}}memory limit in GB; default 8{p_end}
{synopt :{opt maxr:ows(#)}}row limit{p_end}
{synopt :{opt force}}attempt the split although the memory estimate is above the limits (never overrides the observation limit of the Stata edition){p_end}
{synopt :{opt allr:ows}}cut all rows, not only adherent rows{p_end}
{synopt :{opt check}}report the size only; do not split{p_end}
{synopt :{opt det:ail}}print the size table also for the split itself (default: one summary line; the table is always shown for banded splits and memory warnings){p_end}
{synopt :{opt gen:erate(name)}}prefix of the linking variables; default {cmd:_ipw}{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
Cuts rows at failure times strictly inside the row, like {cmd:stsplit, at(failures)}, but only where the on-treatment analysis needs it: only adherent rows, only at
failure times of adherent rows (and of the row's own stratum with {opt strata()}). The default {cmd:bandpoints(1)} is exact: the Breslow estimates and robust standard
errors equal those of the fully split data. The command prints a size check before it splits and refuses (r(950)) if the estimated memory exceeds {opt maxgb()}
unless {opt force} is given.

{pstd}
{opt bandpoints(#)} (approximation, only if the exact split is too large) puts one boundary per stratum in every group of # consecutive cut points; {cmd:bandpoints(auto)}
chooses the smallest # that fits the budget. Evaluate the weights in banded data at the midpoint ({cmd:ipwadh_eval, evalat(mid)}) and compare with another # as a sensitivity analysis.

{pstd}
Creates {cmd:_ipw_orig} (original observation number) and {cmd:_ipw_piece} (piece within the row). The failure and {opt lastonly()} indicators are set to 0 on all but the last piece.
Rows outside {cmd:if}/{cmd:in} are untouched. Split once per analysis. With multiply imputed data the structure is the same in every imputation, so {helpb ipwadh_mi} splits once and recomputes only the weights in each imputation.


{title:Stored results}

{pstd}
{cmd:r(exact)}, {cmd:r(m)}, {cmd:r(N_before)}, {cmd:r(N_after)}, {cmd:r(N_exact)}, {cmd:r(n_cutpoints)}, {cmd:r(n_boundaries)}, {cmd:r(max_bandpoints)}, {cmd:r(est_gb)}, {cmd:r(max_N)}, {cmd:r(ok)}, {cmd:r(forced)}.


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
