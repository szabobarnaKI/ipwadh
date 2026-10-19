{smcl}
{* *! version 0.2.1  11oct2026}{...}
{vieweralsosee "ipwadhc" "help ipwadhc"}{...}
{vieweralsosee "ipwadh_split" "help ipwadh_split"}{...}
{vieweralsosee "ipwadh_eval" "help ipwadh_eval"}{...}
{vieweralsosee "stcox_ipwadh" "help stcox_ipwadh"}{...}
{vieweralsosee "ipwadh_mi" "help ipwadh_mi"}{...}
{vieweralsosee "[ST] stcox" "help stcox"}{...}
{viewerjumpto "Description" "ipwadh##description"}{...}
{viewerjumpto "Main route" "ipwadh##workflow"}{...}
{viewerjumpto "Syntax" "ipwadh##syntax"}{...}
{viewerjumpto "Options" "ipwadh##options"}{...}
{viewerjumpto "Multiple imputation" "ipwadh##mi"}{...}
{viewerjumpto "Advice" "ipwadh##advice"}{...}
{viewerjumpto "Terminal events" "ipwadh##terminal"}{...}
{viewerjumpto "Examples" "ipwadh##examples"}{...}
{viewerjumpto "Stored results" "ipwadh##results"}{...}
{viewerjumpto "AI disclosure" "ipwadh##ai"}{...}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:ipwadh} {hline 2}}Inverse-probability weights for treatment non-adherence (and death / loss to follow-up) for an on-treatment Cox or Lin-Wei-Yang-Ying analysis{p_end}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
The {cmd:ipwadh} suite estimates {it:stabilized} inverse-probability weights for time-varying non-adherence in a randomized trial with recurrent events,
so that the {it:on-treatment} (per-protocol-like) comparison, which uses only the adherent person-time of each patient, can be analysed with a weighted
Cox model (Lin-Wei-Yang-Ying marginal means model for recurrent events, or an ordinary Cox model). Adherence changes, death and loss to
follow-up are modelled as piecewise-constant rates (Poisson regression, one row of the episode data = one period with constant covariates). Death and loss to
follow-up can be handled as informative censoring (inverse-probability-of-censoring weights).

{pstd}
The recommended route has two steps: (1) {cmd:ipwadh, weightsonly} computes the weights in one call and leaves the data {cmd:stset} with them; (2) you fit the
weighted Cox / LWYY model yourself with {cmd:stcox} (or {helpb stcox_ipwadh}). Because the model is yours, you can compare it with unweighted and untruncated versions,
change covariates, strata or ties handling, and reuse the weights.
{cmd:ipwadh} without {opt weightsonly} is a {it:shortcut}: it also fits the standard Cox model and restores your original data (see {help ipwadh##syntax:Syntax of the shortcut}).

{pstd}
For {it:multiply imputed} covariates see {help ipwadh##mi:Multiple imputation} and {helpb ipwadh_mi}: with {opt frame()} the same command returns the split data and the weights of every imputation as an mi dataset.

{pstd}
{cmd:ipwadh} runs four stand-alone commands in this order; each can also be used alone when you want to inspect or change one step:

{p2colset 7 24 26 2}{...}
{p2col :{helpb ipwadhc}}rate models and stabilized weights; with option {opt parts} it keeps the pieces needed to evaluate the weight inside a row{p_end}
{p2col :{helpb ipwadh_split}}cuts the adherent rows exactly at the failure times (as {cmd:stsplit, at(failures)} would, but only where needed){p_end}
{p2col :{helpb ipwadh_eval}}weight of every piece, diagnostics, truncation, weight-instability warning{p_end}
{p2col :{helpb stcox_ipwadh}}{it:optional}: standard errors that also account for the estimation of the weights{p_end}
{p2colreset}{...}

{marker workflow}{...}
{title:Main route: weights first, then your own stcox}

{phang2}{cmd:. ipwadh, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) failure(inf) dead(dead) ltfu(ltfu) ///}{p_end}
{phang2}{cmd:      stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) weightsonly}{p_end}
{phang2}{cmd:. stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}

{pstd}
{opt weightsonly} runs the rate models, the exact split, the piece weights (truncated at p1/p99 unless {opt notruncate}), prints the weight diagnostics and warning, and
leaves the {it:split} data in memory, {cmd:stset} with {cmd:[pw=_ipw_pswt]} (or {cmd:_ipw_psw} with {opt notruncate}) for the adherent pieces. The original data are not restored,
so reload them before a new run. The last line of the output repeats the {cmd:stset} used and the Cox model to fit.

{pstd}
{bf:Keep your own data untouched.} With {opt frame(name)} the same call runs on a copy and leaves the weighted, {cmd:stset} split data in the frame {it:name} (add {opt change} to switch to it); your data stay as they are.

{pstd}
{bf:Which variables go into stcox afterwards?} The treatment arm and {it:exactly} the variables you gave in {opt stcoxvars()}, nothing more and nothing less
(the example above: {opt stcoxvars(c.age i.male)} gives {cmd:stcox i.arm c.age i.male}). Reason: the numerator models of the stabilized weights contain the arm and the
{opt stcoxvars()}, so these are the covariates the weights are valid for. Variables in {opt adhpredictors()} or {opt censorpredictors()} explain adherence or censoring;
they stay out of {cmd:stcox}. With {opt strata()} use the same variables in {cmd:stcox, strata()} and leave them out of the covariate list. Adding other covariates, or leaving
one out, changes the model and the weights no longer match it.

{pstd}
{bf:Recommended: run three versions} (see {help ipwadh##advice:Advice}): the same {cmd:stcox} on the truncated weights (the default above), on the untruncated weights
({cmd:ipwadh, ... notruncate weightsonly} on freshly loaded data), and unweighted ({cmd:stset t1 if adh==1, time0(t0) failure(inf)} on the same data, then {cmd:stcox}).

{pstd}
{bf:The same in single steps} (equivalent; use when you want to change or inspect one step):

{phang2}{cmd:. ipwadhc, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu) ///}{p_end}
{phang2}{cmd:      stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) parts replace}{p_end}
{phang2}{cmd:. ipwadh_split, id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)}{p_end}
{phang2}{cmd:. ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace}{p_end}
{phang2}{cmd:. stset t1 [pw=_ipw_pswt] if adh==1, time0(t0) failure(inf)}{p_end}
{phang2}{cmd:. stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}

{pstd}
Rules that make the weights valid for the Cox model:

{phang2}1. The Cox model contains the treatment arm and {it:exactly} the covariates given in {opt stcoxvars()} of {cmd:ipwadhc} (the numerator models of the stabilized
weights contain the same variables). Strata variables may be used in {cmd:stcox, strata()}; keep them in {opt stcoxvars()} as well.{p_end}
{phang2}2. {cmd:ipwadhc} must be run on the {it:unsplit} episode data. Split afterwards.{p_end}
{phang2}3. {cmd:stset} with {opt pweight}s, {it:without} {opt id()} and without {opt exit()} (Stata refuses time-varying weights together with {opt id()}). Use {opt breslow}
ties and {cmd:vce(cluster} {it:id}{cmd:)}.{p_end}
{phang2}4. Do not use a {opt tvc()} or {opt offset()} specification with the weights.{p_end}


{marker syntax}{...}
{title:Syntax of the shortcut (and of weightsonly)}

{p 8 16 2}
{cmd:ipwadh} {ifin}{cmd:,} {opth id(varname)} {opth time0(varname)} {opth time1(varname)} {opth adh:erence(varname)} {opth arm(varname)} {opth fail:ure(varname)} [{it:options}]

{synoptset 32 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt :{opth dead(varname)}}death indicator (1 on the last row of a patient who died){p_end}
{synopt :{opth ltfu(varname)}}loss-to-follow-up indicator; weighted (IPCW) when given{p_end}
{synopt :{opt term:inal(censor|none)}}handling of death; default {cmd:censor}{p_end}
{synopt :{opt stcox:vars(string)}}covariates of the Cox model besides the arm{p_end}
{synopt :{opt adhp:redictors(string)}}extra predictors of stopping / restarting treatment{p_end}
{synopt :{opt censorp:redictors(string)}}extra predictors of death / loss to follow-up{p_end}
{synopt :{opt byarm}}fit the rate models separately by arm{p_end}
{synopt :{opt pois:sonopts(string)}}options passed to every Poisson fit{p_end}
{synopt :{opth str:ata(varlist)}}stratified split and stratified Cox model{p_end}
{syntab:Weights}
{synopt :{opt trunc:ate(# #)}}truncation percentiles; default {cmd:truncate(1 99)}{p_end}
{synopt :{opt notrun:cate}}use the untruncated weights{p_end}
{synopt :{opt band:points(# | auto)}}banded split for very large data; default exact split{p_end}
{synopt :{opt maxgb(#)}}memory limit for the split; default 8{p_end}
{synopt :{opt force}}attempt the split although the memory estimate exceeds the limits{p_end}
{synopt :{opt eval:at(end|mid|start)}}point of the piece where the weight is evaluated; default {cmd:end}{p_end}
{synopt :{opt esscaution(#)}}ESS percent below which a caution is printed; default 20{p_end}
{synopt :{opt esssevere(#)}}ESS percent below which a severe warning is printed; default 5{p_end}
{synopt :{opt maxwarn(#)}}largest weight above which a note is printed; default 100{p_end}
{syntab:Output}
{synopt :{opt adj:ustse}}standard errors that account for the estimated weights ({helpb stcox_ipwadh}){p_end}
{synopt :{opt nohr}}report coefficients instead of hazard ratios{p_end}
{synopt :{opt l:evel(#)}}confidence level{p_end}
{synopt :{opt keepdata}}leave the split data, the weights and the {cmd:stset} in memory{p_end}
{synopt :{opt weightsonly}}stop after {cmd:stset}; fit your own model{p_end}
{synopt :{opt verb:ose}}show the output of the individual steps{p_end}
{syntab:Result in a frame}
{synopt :{opt fr:ame(name [, replace])}}leave the split data, the weights and the {cmd:stset} in a new frame, your data untouched; {it:required} for mi data (the result is an mi dataset){p_end}
{synopt :{opt ch:ange}}switch to the result frame (needs {opt frame()}){p_end}
{syntab:Multiply imputed (mi set) data only}
{synopt :{opt mistyle(flong|wide)}}mi style of the result; default {cmd:flong}{p_end}
{synopt :{opt mistset(psw|pswt)}}run {cmd:mi stset} with this weight at the end (flong){p_end}
{synopt :{opt mik:eep(varlist)}}keep only these variables besides the ones needed (memory){p_end}
{synoptline}
{p2colreset}{...}


{marker options}{...}
{title:What the command does and options}

{pstd}
In this order: {helpb ipwadhc} {cmd:..., parts replace}; {helpb ipwadh_split} {cmd:..., lastonly(}{it:dead ltfu}{cmd:)}; {helpb ipwadh_eval};
{cmd:stset} {it:time1} {cmd:[pw=}{it:weight}{cmd:] if} {it:adherence}{cmd:==1,} {cmd:time0(}{it:time0}{cmd:) failure(}{it:failure}{cmd:)};
{cmd:stcox i.}{it:arm} {it:stcoxvars}{cmd:, breslow vce(cluster} {it:id}{cmd:)}. The data must be the unsplit episode data (variables {cmd:_ipw_orig} / {cmd:_ipw_piece} must not exist).

{phang}{opt truncate(# #)} clips the piece weights at the person-time-weighted percentiles; default {cmd:truncate(1 99)}. {opt notruncate} uses the untruncated weight.
Results with and without truncation should both be looked at (see {help ipwadh##advice:Advice}).

{phang}{opt strata(varlist)} cuts each row only at failure times of its own stratum and fits {cmd:stcox, strata()}. The strata variables are removed from the covariates of
{cmd:stcox} but stay in the weight models if they appear in {opt stcoxvars()}. The stratified split is exact only together with a Cox model stratified by the
same variables ({cmd:stcox ..., strata()}); an unstratified Cox model should use the split without {opt strata()}.

{phang}{opt adjustse} replaces the cluster-robust variance by the variance that accounts for the estimation of the weights ({helpb stcox_ipwadh}). It is supplementary; see {help ipwadh##advice:Advice}.

{phang}{opt keepdata} leaves the split data in memory together with {cmd:_ipw_psw}, {cmd:_ipw_pswt} and the {cmd:stset}. Without it, the original unsplit data are restored and only
the estimation results remain. {opt weightsonly} implies {opt keepdata} and fits no model.

{phang}{opt verbose} shows the full output ({opt detail}) of the individual steps. By default only a short header, the compact {helpb ipwadh_eval} report with the weight warning, and the Cox table are printed.

{phang}{opt bandpoints()}, {opt maxgb()}, {opt force}, {opt evalat()}, {opt esscaution()}, {opt esssevere()}, {opt maxwarn()}, {opt poissonopts()}, {opt byarm}: as in
{helpb ipwadhc}, {helpb ipwadh_split} and {helpb ipwadh_eval}.

{phang}{opt frame(name [, replace])} runs the pipeline on a {it:copy} of your data in a new {help frame:frame} and leaves the split data, the weights and the {cmd:stset} (as with {opt weightsonly}) there;
the data in your current frame are not changed and no model is fitted ({opt adjustse} is not allowed). An existing frame is only overwritten with {opt replace}; the current frame cannot be the result frame; if the command fails, no frame is left behind.
The data exist twice while the command runs. Then {cmd:frame change} {it:name} (or the option {opt change}) and fit your model there. Needs Stata 16 or later. With mi data the result is an mi dataset, see {help ipwadh##mi:Multiple imputation}.

{phang}{opt change} switches to the result frame at the end.

{phang}{opt mistyle()}, {opt mistset()}, {opt mikeep()}: for mi data only (see {help ipwadh##mi:Multiple imputation}); an error on ordinary data.


{marker mi}{...}
{title:Multiple imputation of covariates}

{pstd}
If the data in memory are {cmd:mi set} (covariates such as a biomarker imputed), {cmd:ipwadh} needs {opt frame(name [, replace])} and hands the call to {helpb ipwadh_mi} with the same options.
Your mi data are not changed. The structure of the episode data (id, times, adherence, events, strata, arm) must be the same in every imputation, so the rows are cut {it:once}; the rate models and the weights
are computed {it:in every imputation}, and the result frame holds the split data and the weights as an mi dataset (default {cmd:flong}). You then call {helpb mi estimate} yourself:

{phang2}{cmd:. ipwadh, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) failure(inf) dead(dead) ltfu(ltfu) ///}{p_end}
{phang2}{cmd:      stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) frame(res) mistset(psw) change}{p_end}
{phang2}{cmd:. mi estimate: stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}
{phang2}{cmd:. mi estimate: mean age biomarker [pw=_ipw_ptw] if adh==1, over(arm)}{space 4}// weighted baseline table{p_end}

{pstd}
No model is fitted by {cmd:ipwadh} itself, {opt adjustse} is not available, and {opt keepdata} / {opt weightsonly} are implied. {cmd:mi estimate} cannot wrap the whole pipeline (the split changes the number
of observations); the result frame is the way around that. What the weights, the person-time weights and the estimation steps are, which limits apply (imputed events, times or adherence are not supported;
the strata variables of {opt strata()} may not be imputed; memory) and how to read the output is described in {helpb ipwadh_mi}.
{cmd:mi estimate: stcox} with {opt strata()} is certified: use the same strata variables in the split and in {cmd:stcox, strata()}.


{marker advice}{...}
{title:Advice: standard errors, truncation and how to judge the weights}

{pstd}
{bf:Default standard errors.} Use the cluster-robust standard errors ({cmd:vce(cluster} {it:id}{cmd:)}, which is what {cmd:ipwadh} and the workflow above give you).
They treat the weights as known. Accounting for the estimation of the weights ({helpb stcox_ipwadh}, option {opt adjustse}) changes the standard errors by only
a few percent in either direction in our simulations, and it does {it:not} repair the main problem of weighted analyses, namely unstable weights. {cmd:stcox_ipwadh}
is therefore optional and supplementary.

{pstd}
{bf:Always compare three analyses} of the same model: (1) unweighted (on-treatment person-time without weights), (2) weighted, untruncated, (3) weighted, truncated
({cmd:truncate(1 99)}). If the three agree, the weights are not driving the result. If the weighted estimate differs strongly from the unweighted one,
report all three and discuss why. Truncation lowers the variance but moves the estimate towards the unweighted one, so neither the truncated nor the untruncated analysis should be
presented alone when weights are heavy-tailed.

{pstd}
{bf:Read the weight warning} that {helpb ipwadh_eval} (and {cmd:ipwadh}) prints. The effective sample size (ESS, as share of the person-time) summarises how much of the data
the weights effectively use: about 40-55% was stable in the simulations; below 20% (default {opt esscaution()}) use truncation and compare; below 5% (default {opt esssevere()})
the weighted estimate was unstable and its standard error too small. A very large single weight ({opt maxwarn()}) deserves a look at the pieces concerned. Extreme weights typically
arise when adherence depends strongly on a persistent covariate (near violation of positivity): simplify the adherence model rather than trusting the standard error.
These thresholds are heuristics, not tests, and cannot detect unmeasured confounding of adherence.


{marker terminal}{...}
{title:Terminal events (death)}

{pstd}
The estimand depends on how death is handled; choose it with {opt dead()} and {opt terminal()} (see {helpb ipwadhc}).

{phang2}{cmd:terminal(censor)} (default). Death is informative censoring and is weighted (IPCW). Estimand: the rate ratio had nobody died (hypothetical). Needs {opt dead()}.{p_end}
{phang2}{cmd:terminal(none)}. No death weights; death ends follow-up as independent censoring. Interpretation: infections per unit time alive (while-alive rate ratio).{p_end}
{phang2}Composite event (infection or death): create {cmd:gen evt = inf | dead}, use {cmd:failure(evt)} and {cmd:terminal(none)}.{p_end}
{phang2}{opt ltfu()}: loss to follow-up is independent censoring unless the option is given; then it is weighted under every strategy.{p_end}

{pstd}
Ghosh-Lin-type analyses (dead persons remain in the risk set) are not implemented.


{marker examples}{...}
{title:Examples}

{pstd}Main route, weights only, then your own models:{p_end}
{phang2}{cmd:. ipwadh, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) failure(inf) dead(dead) ltfu(ltfu) ///}{p_end}
{phang2}{cmd:      stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) weightsonly}{p_end}
{phang2}{cmd:. stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}
{phang2}{cmd:. stcox_ipwadh i.arm c.age i.male}{space 8}// optional: SE that accounts for the estimated weights{p_end}

{pstd}Unweighted comparison on the same data:{p_end}
{phang2}{cmd:. stset t1 if adh==1, time0(t0) failure(inf)}{p_end}
{phang2}{cmd:. stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}

{pstd}Untruncated weights (reload the original data first):{p_end}
{phang2}{cmd:. ipwadh, ... notruncate weightsonly}{p_end}

{pstd}Stratified analysis:{p_end}
{phang2}{cmd:. ipwadh, ... stcoxvars(c.age i.male) strata(male) weightsonly}{p_end}
{phang2}{cmd:. stcox i.arm c.age, breslow vce(cluster id) strata(male)}{p_end}

{pstd}Shortcut: weights and the standard Cox model in one call (original data restored; defaults truncation p1/p99, cluster-robust SE):{p_end}
{phang2}{cmd:. ipwadh, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) failure(inf) dead(dead) ltfu(ltfu) ///}{p_end}
{phang2}{cmd:      stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)}{p_end}
{phang2}{cmd:. ipwadh, ... notruncate}{space 8}// untruncated version{p_end}
{phang2}{cmd:. ipwadh, ... adjustse}{space 11}// adjusted SE{p_end}
{phang2}{cmd:. ipwadh, ... keepdata}{space 11}// leave split data, weights and stset in memory{p_end}

{pstd}Leave your own data untouched; the weighted split data go to the frame {cmd:w}:{p_end}
{phang2}{cmd:. ipwadh, ... frame(w) change}{p_end}
{phang2}{cmd:. stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}
{phang2}{cmd:. frame change default}{p_end}

{pstd}Multiply imputed covariates (the data are {cmd:mi set}); the result is an mi dataset in frame {cmd:res}:{p_end}
{phang2}{cmd:. ipwadh, ... frame(res) mistset(psw) change}{p_end}
{phang2}{cmd:. mi estimate: stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}
{phang2}{cmd:. mi estimate: mean age [pw=_ipw_ptw] if adh==1, over(arm)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ipwadh} stores the results of {helpb stcox} in {cmd:e()} and, in addition:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(ipwadh_ess)}}effective sample size, untruncated weight (share of person-time){p_end}
{synopt:{cmd:e(ipwadh_ess_min)}}smallest ESS over the arms{p_end}
{synopt:{cmd:e(ipwadh_ess_trunc)}}ESS of the truncated weight (missing without truncation){p_end}
{synopt:{cmd:e(ipwadh_wmax)}}largest weight{p_end}
{synopt:{cmd:e(ipwadh_warn)}}0 / 1 caution / 2 severe weight warning{p_end}
{synopt:{cmd:e(ipwadh_n_pieces)}}number of pieces after the split{p_end}
{synopt:{cmd:e(ipwadh_n_rows)}}number of usable episode rows{p_end}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(ipwadh_weight)}}name of the weight variable used{p_end}
{synopt:{cmd:e(ipwadh_se)}}{cmd:cluster-robust} or {cmd:adjusted}{p_end}
{synopt:{cmd:e(ipwadh_cmd)}}Cox covariates used{p_end}
{p2colreset}{...}

{pstd}
Without {opt keepdata} the data in memory are the original data, so {cmd:e(sample)} does not refer to them and {cmd:predict} is not available.
With {opt frame()} no model is fitted and nothing is stored in {cmd:e()}. With mi data the model specification is in the characteristics {cmd:_dta[ipwadhmi_*]} of the result frame (see {helpb ipwadh_mi}).


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
