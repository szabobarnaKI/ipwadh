{smcl}
{* *! version 0.3.1  10oct2026}{...}
{vieweralsosee "ipwadh" "help ipwadh"}{...}
{vieweralsosee "ipwadh_split" "help ipwadh_split"}{...}
{vieweralsosee "ipwadh_eval" "help ipwadh_eval"}{...}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:ipwadhc} {hline 2}}Continuous-time stabilized weights for time-varying non-adherence, death and loss to follow-up{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:ipwadhc} {ifin}{cmd:,} {opth id(varname)} {opth time0(varname)} {opth time1(varname)} {opth adh:erence(varname)} {opth arm(varname)} [{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Plain-language covariate lists (recommended)}
{synopt :{opt stcox:vars(string)}}covariates of your {cmd:stcox} besides the arm; enter both the numerator and the denominator models{p_end}
{synopt :{opt adhp:redictors(string)}}additional predictors of stopping / restarting treatment (denominator only){p_end}
{synopt :{opt censorp:redictors(string)}}additional predictors of death / loss to follow-up (denominator only){p_end}
{syntab:Censoring}
{synopt :{opth dead(varname)}}death indicator{p_end}
{synopt :{opth ltfu(varname)}}loss-to-follow-up indicator{p_end}
{synopt :{opt term:inal(censor|none)}}death as informative censoring (default) or as independent censoring{p_end}
{syntab:Model}
{synopt :{opt byarm}}fit all rate models separately by arm{p_end}
{synopt :{opt pois:sonopts(string)}}options for every Poisson fit; replaces the default tight tolerances{p_end}
{syntab:Output}
{synopt :{opt parts}}keep the pieces of the weight formula; needed by {helpb ipwadh_split} / {helpb ipwadh_eval} and {helpb stcox_ipwadh}{p_end}
{synopt :{opt trunc:ate(# #)}}clip row-end weights at row-weighted percentiles (quick approximation only){p_end}
{synopt :{opt gen:erate(name)}}prefix of the new variables; default {cmd:_ipw}{p_end}
{synopt :{opt replace}}replace existing variables{p_end}
{synopt :{opt noi:sily}}show the Poisson fits{p_end}
{synopt :{opt det:ail}}full output: covariate block, table of all weight components, usage advice (default: a short summary){p_end}
{syntab:Advanced (do not mix with the plain-language lists)}
{synopt :{opt num:vars(string)}}explicit numerator covariates of the adherence model{p_end}
{synopt :{opt den:vars(string)}}explicit denominator covariates of the adherence model{p_end}
{synopt :{opt cnum:vars(string)}}explicit numerator covariates of the censoring model{p_end}
{synopt :{opt cden:vars(string)}}explicit denominator covariates of the censoring model{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
{cmd:ipwadhc} needs {it:unsplit} episode data: one row per period with constant covariates (a new row at every adherence change and at every change of a time-varying
covariate); {opt adherence()} is 0/1. Adherence changes (stopping, restarting), death and loss to follow-up are modelled by piecewise-constant rates (Poisson regression
with the row length as exposure), with a numerator and a denominator model for each. The result is the stabilized weight

{p 8 8 2}{cmd:_ipw_swa} (adherence) {c 215} {cmd:_ipw_swc} (censoring by death and loss to follow-up) = {cmd:_ipw_sw}

{pstd}
evaluated at the {it:end of each row}. Weights are left-continuous: the jump at an adherence change enters from the next row on. With {opt parts} the pieces
{cmd:_ipw_E0 _ipw_J0 _ipw_dA _ipw_Ec0 _ipw_dC} are kept, so that the weight at any time t inside a row is
{cmd:log W(t) = E0 - dA(t-time0) + J0 + Ec0 - dC(t-time0)}; {helpb ipwadh_eval} uses them after the rows have been split.

{pstd}
The Cox model that uses the weights must contain the arm and exactly the covariates in {opt stcoxvars()}. The command prints the line to use. Time-fixed or
time-varying covariates in {opt adhpredictors()} are kept out of the comparison: they explain adherence but are not in the outcome model.

{pstd}
Ignored rows (weight missing, no error): rows outside {cmd:if}/{cmd:in}, rows with missing id / time / adherence / arm, and zero-length rows ({it:time1} = {it:time0}), which
{cmd:stset} ignores as well. Rows with {it:time1} < {it:time0} are an error. Baseline non-adherence is not modelled: everyone starts in the state of the first usable row.

{pstd}
Poisson fits use tight tolerances by default and fall back to Stata's defaults (with a note) if a tight fit does not converge. A warning is shown when a stratum
loses rows to perfect prediction.


{title:Terminal events}

{phang}{cmd:terminal(censor)} (default): death is informative censoring and weighted; needs {opt dead()}. Estimand: rate ratio had nobody died.{p_end}
{phang}{cmd:terminal(none)}: no death weights; {opt dead()} is only checked. Interpretation: infections per unit time alive.{p_end}
{phang}Composite event: {cmd:gen evt = inf | dead}, {cmd:failure(evt)} in {helpb ipwadh_split} and {cmd:stset}, and {cmd:terminal(none)}.{p_end}
{phang}{opt ltfu()}: weighted (IPCW) under every strategy; without it loss to follow-up is independent censoring.{p_end}


{title:Stored results}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}usable rows{p_end}
{synopt:{cmd:r(n_zero)}}zero-length rows ignored{p_end}
{synopt:{cmd:r(n_changes)}}adherence changes{p_end}
{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:r(weightvar)}}name of the row-end weight variable{p_end}
{synopt:{cmd:r(censoring)}}censoring causes weighted{p_end}
{synopt:{cmd:r(terminal)}}strategy used{p_end}
{p2colreset}{...}

{pstd}
Variables: {cmd:_ipw_swa _ipw_swc _ipw_sw} (and with {opt truncate()} {cmd:_ipw_swt}); with {opt parts} also {cmd:_ipw_E0 _ipw_J0 _ipw_dA _ipw_Ec0 _ipw_dC}.
The model specification is stored in data characteristics {cmd:_dta[ipwadh_*]} for {helpb stcox_ipwadh}.


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
