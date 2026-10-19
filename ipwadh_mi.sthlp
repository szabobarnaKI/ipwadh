{smcl}
{* *! version 0.2.1  11oct2026}{...}
{vieweralsosee "ipwadh" "help ipwadh"}{...}
{vieweralsosee "ipwadhc" "help ipwadhc"}{...}
{vieweralsosee "ipwadh_split" "help ipwadh_split"}{...}
{vieweralsosee "ipwadh_eval" "help ipwadh_eval"}{...}
{vieweralsosee "[MI] mi estimate" "help mi estimate"}{...}
{vieweralsosee "[MI] mi stset" "help mi stset"}{...}
{viewerjumpto "Description" "ipwadh_mi##description"}{...}
{viewerjumpto "Syntax" "ipwadh_mi##syntax"}{...}
{viewerjumpto "Options" "ipwadh_mi##options"}{...}
{viewerjumpto "The result frame" "ipwadh_mi##result"}{...}
{viewerjumpto "Using the result" "ipwadh_mi##using"}{...}
{viewerjumpto "Limits" "ipwadh_mi##limits"}{...}
{viewerjumpto "Examples" "ipwadh_mi##examples"}{...}
{viewerjumpto "AI disclosure" "ipwadh_mi##ai"}{...}
{title:Title}

{p2colset 5 19 21 2}{...}
{p2col :{cmd:ipwadh_mi} {hline 2}}Inverse-probability weights for non-adherence on multiply imputed data; the weighted exact-split data are returned as an mi dataset{p_end}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ipwadh_mi} does what {helpb ipwadh} does, but for {it:mi set} data in which covariates (for example a biomarker) were multiply imputed. The weights depend on the
covariates, so they differ between the imputations, while the {it:structure} of the episode data (who, which rows, times, adherence, events) is the same in every imputation.
{cmd:ipwadh_mi} uses this: it cuts the rows {it:once} (exact split at the failure times, as {helpb ipwadh_split}), then computes the rate models and the weights {it:in every imputation}
({helpb ipwadhc}, {helpb ipwadh_eval}) and stores them in an {cmd:mi} dataset in a new {help frame:frame}.

{pstd}
Your mi data are never changed. The result frame can be used directly with {helpb mi estimate}: {cmd:mi estimate: stcox} for the weighted hazard ratio, and
{cmd:mi estimate: mean} or {cmd:mi estimate: proportion} for weighted summaries of any variable (a baseline or characteristics table). Nothing is fitted by {cmd:ipwadh_mi} itself;
the point is that the split data and the weights stay in memory so that you can fit and summarise whatever you need, pooled by Rubin's rules.

{pstd}
{cmd:ipwadh} calls {cmd:ipwadh_mi} when the data in memory are mi set and {opt frame()} is given; both commands have the same options. ({cmd:ipwadh, frame()} on ordinary data does the analogous thing for ordinary data: the split data and the weights go to a new frame.) Requires Stata 16 or later (frames).


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:ipwadh_mi} {ifin}{cmd:,} {opth id(varname)} {opth time0(varname)} {opth time1(varname)} {opth adh:erence(varname)} {opth arm(varname)} {opth fail:ure(varname)}
{opt fr:ame(name [, replace])} [{it:options}]

{synoptset 32 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model, weights}
{synopt :{opth dead(varname)}}death indicator (1 on the last row of a patient who died){p_end}
{synopt :{opth ltfu(varname)}}loss-to-follow-up indicator; weighted (IPCW) when given{p_end}
{synopt :{opt term:inal(censor|none)}}handling of death; default {cmd:censor}{p_end}
{synopt :{opt stcox:vars(string)}}covariates of the Cox model besides the arm{p_end}
{synopt :{opt adhp:redictors(string)}}extra predictors of stopping / restarting treatment{p_end}
{synopt :{opt censorp:redictors(string)}}extra predictors of death / loss to follow-up{p_end}
{synopt :{opt byarm}}fit the rate models separately by arm{p_end}
{synopt :{opt pois:sonopts(string)}}options passed to every Poisson fit{p_end}
{synopt :{opth str:ata(varlist)}}stratified split{p_end}
{synopt :{opt trunc:ate(# #)}}truncation percentiles; default {cmd:truncate(1 99)}{p_end}
{synopt :{opt notrun:cate}}use the untruncated weights{p_end}
{synopt :{opt band:points(# | auto)}}banded split for very large data; default exact split{p_end}
{synopt :{opt maxgb(#)}}memory limit for the split; default 8{p_end}
{synopt :{opt force}}attempt the split although the memory estimate exceeds the limits{p_end}
{synopt :{opt eval:at(end|mid|start)}}point of the piece where the weight is evaluated; default {cmd:end}{p_end}
{synopt :{opt esscaution(#)}}ESS percent below which a caution is printed; default 20{p_end}
{synopt :{opt esssevere(#)}}ESS percent below which a severe warning is printed; default 5{p_end}
{synopt :{opt maxwarn(#)}}largest weight above which a note is printed; default 100{p_end}
{syntab:Result}
{synopt :{opt fr:ame(name [, replace])}}{it:required}: name of the result frame{p_end}
{synopt :{opt mistyle(flong|wide)}}mi style of the result; default {cmd:flong}{p_end}
{synopt :{opt mistset(psw|pswt)}}run {cmd:mi stset} with that weight at the end (needs {cmd:flong}){p_end}
{synopt :{opt keep(varlist)}}keep only these variables besides the ones ipwadh needs (saves memory){p_end}
{synopt :{opt ch:ange}}switch to the result frame at the end{p_end}
{synopt :{opt verb:ose}}show the output of the individual steps{p_end}
{synoptline}
{p2colreset}{...}
{pstd}
{opt if} and {opt in} are evaluated on the {it:m}=0 data; do not use imputed variables in them.
In {cmd:ipwadh} the option {opt keep()} is called {opt mikeep()} (because {opt keep()} would be ambiguous with {opt keepdata}).


{marker options}{...}
{title:Options}

{phang}
{opt frame(name [, replace])} names the result frame. An existing frame is only overwritten with {opt replace}; the current frame cannot be the result frame. If the command fails, no result frame is left behind.

{phang}
{opt mistyle(flong|wide)} chooses the mi style of the result.

{phang2}
{cmd:flong} (default). The result has (M+1) x (number of pieces) rows. The weights {cmd:_ipw_psw}, {cmd:_ipw_pswt}, {cmd:_ipw_ptw}, {cmd:_ipw_ptwt} are {it:unregistered}
("super-varying") variables: the {it:m}=0 rows hold the value 1 (and {cmd:_ipw_ptw} = {it:time1}-{it:time0}), the rows of imputation {it:m} hold the weights of imputation {it:m}. This is the only form that
{helpb mi stset} and {cmd:mi estimate: stcox} accept natively with weights that differ between imputations ({cmd:mi stset} refuses weights that are registered as imputed or passive).{p_end}
{phang2}
{cmd:wide}. The result has one row per piece and one extra column per imputation for each varying variable. The weights are registered {it:imputed}, with the {it:m}=0 value missing.
Much cheaper in memory, and {cmd:mi estimate: mean} / {cmd:proportion} with {cmd:[pw=_ipw_ptw]} work, but {cmd:mi stset} and {cmd:mi estimate: stcox} do not (see {help ipwadh_mi##limits:Limits}).{p_end}

{phang}
{opt mistset(psw|pswt)} runs {cmd:mi stset} {it:time1} {cmd:[pw=_ipw_psw]} (or {cmd:_ipw_pswt}) {cmd:if} {it:adherence}{cmd:==1, time0(}{it:time0}{cmd:) failure(}{it:failure}{cmd:)} on the result, so that
{cmd:mi estimate: stcox} can be called at once. {cmd:pswt} needs truncated weights. Without the option the result is {it:not} stset; use the line printed at the end of the output.

{phang}
{opt keep(varlist)} drops all variables of the result except the structure variables, the imputed and passive variables, the weights, the internal keys and {it:varlist}. Use it for large data: the flong
result stores every kept variable (M+1) times.

{phang}
{opt strata(varlist)} cuts each row only at the failure times of its own stratum (as in {helpb ipwadh_split}); the exact result then needs a Cox model stratified by the same variables
(see {help ipwadh_mi##using:Using the result}). The strata variables must not be registered imputed or passive.

{phang}
{opt change} switches to the result frame when the command is done (otherwise the current frame is unchanged and the output tells you how to switch).

{phang}
All other options have the meaning described in {helpb ipwadh}.


{marker result}{...}
{title:The result frame}

{pstd}
The result frame holds the exact-split pieces. Besides your variables it contains

{p2colset 7 28 30 2}{...}
{p2col :{cmd:_ipw_orig}}row of the unsplit {it:m}=0 data from which the piece comes{p_end}
{p2col :{cmd:_ipw_piece}}running number of the piece within its row{p_end}
{p2col :{cmd:_ipw_psw}}stabilized weight at the end of the piece (the weight at the failure time), untruncated{p_end}
{p2col :{cmd:_ipw_pswt}}the same, truncated (not with {opt notruncate}){p_end}
{p2col :{cmd:_ipw_ptw}, {cmd:_ipw_ptwt}}person-time weights: weight x ({it:time1} - {it:time0}) of the piece{p_end}
{p2col :{cmd:_ipw_row}, {cmd:_ipw_sel}}row key and selection indicator used internally{p_end}
{p2colreset}{...}

{pstd}
The imputed variables of your data are registered imputed in the result and carry the imputed values of every imputation; passive variables are carried as passive, with their values {it:copied}, not re-derived.
Everything else is registered regular. The model specification is stored in the characteristics {cmd:_dta[ipwadhmi_*]} of the result frame (id, times, adherence, arm, failure, M, mi style, style of the source, truncation).

{pstd}
In the weights, pieces of zero length ({it:time0} >= {it:time1} in the data) have a missing weight in every imputation. {cmd:mi stset} reports them as entries on or after exit; they carry no person-time.

{pstd}
At the end a table shows, for every imputation, the effective sample size of the weights (as a share of the person-time), the smallest ESS over the arms, the largest weight and the weight warning of {helpb ipwadh_eval}.
Look at it: a weight problem in a single imputation matters for the pooled result.


{marker using}{...}
{title:Using the result}

{pstd}
{bf:Weighted hazard ratio / LWYY model} (result in {cmd:flong}):

{phang2}{cmd:. frame change res}{p_end}
{phang2}{cmd:. mi stset t1 [pw=_ipw_psw] if adh==1, time0(t0) failure(inf)}{space 4}// not needed after mistset(psw){p_end}
{phang2}{cmd:. mi estimate: stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}

{pstd}
The Cox model is the same as for {cmd:ipwadh}: the arm and exactly the variables of {opt stcoxvars()}. In the tests, the pooled coefficients equal the mean of the per-imputation {cmd:stset}/{cmd:stcox} fits
(difference about 1e-16), including the within-imputation variance. Use {cmd:_ipw_pswt} for the truncated weights. Run the three versions recommended in {help ipwadh##advice:Advice}
(unweighted, untruncated, truncated).

{pstd}
{bf:Stratified Cox model.} {cmd:mi estimate: stcox ..., strata()} works on the result (certified: the pooled coefficients and the total variance equal the per-imputation {cmd:stset}/{cmd:stcox, strata()} fits combined by Rubin's rules).
Use the same strata variables in the split and in {cmd:stcox}:

{phang2}{cmd:. ipwadh_mi, ... stcoxvars(c.age i.male) strata(male) frame(res) mistset(psw) change}{p_end}
{phang2}{cmd:. mi estimate: stcox i.arm c.age, strata(male) breslow vce(cluster id)}{p_end}

{pstd}
The stratified split has fewer pieces and gives {it:exactly} the same stratified estimate as the unstratified exact split. It is {it:not} exact for a Cox model that is not stratified by the same variables:
for an unstratified {cmd:stcox} use the split without {opt strata()}. A time-varying strata variable that is constant within a row (such as time since infection) and several strata variables work.
The strata variables of the {it:split} may not be imputed; {cmd:stcox, strata()} may use an imputed variable when the split is not stratified.

{pstd}
{bf:Weighted baseline / characteristics table.} The natural weight of a summary over the on-treatment follow-up is the {it:person-time weight}: each adherent piece counts with its weight times its length.

{phang2}{cmd:. mi estimate: mean age biomarker [pw=_ipw_ptw] if adh==1, over(arm)}{p_end}
{phang2}{cmd:. mi estimate: proportion male [pw=_ipw_ptw] if adh==1, over(arm)}{p_end}

{pstd}
This works in {cmd:flong} and {cmd:wide} results and does not need {cmd:mi stset}. Use {cmd:_ipw_ptwt} for the truncated weights. In the tests the pooled means equal the mean over the imputations of the per-imputation
weighted means. Standard deviations and medians are not combined by Rubin's rules; compute them per imputation (for example with {cmd:summarize} {it:x} {cmd:[aw=_ipw_ptw]} after {cmd:mi extract}) and report what you did.
Note that at the {it:first} row of a patient the stabilized weights are 1, so a table at baseline would equal the unweighted table; the person-time weighted table describes the weighted follow-up.

{pstd}
{bf:Messages that are not errors.} {cmd:mi varying} lists the weights as "unregistered super varying": a report. {cmd:mi estimate} prints "Number of observations in a subpopulation varies among imputations":
the weighted number of observations of the {cmd:if} subpopulation differs between imputations. Both are expected here.


{marker limits}{...}
{title:Limits}

{phang}1. {it:The structure may not vary.} The variables {opt id()}, {opt time0()}, {opt time1()}, {opt adherence()}, {opt failure()}, {opt dead()}, {opt ltfu()}, {opt strata()} and {opt arm()} must not be registered imputed or passive,
and the rows must be the same in every imputation (this is checked; the command stops otherwise). With imputed events, times or adherence use the per-imputation loop with the full pipeline inside each imputation.{p_end}
{phang}2. At least two imputations; imputed and passive variables must be numeric; source data in style wide, flong or mlong ({cmd:mi extract} keeps the row order for wide and flong but not for mlong; {cmd:ipwadh_mi} sorts by an explicit row key, so all three work).{p_end}
{phang}3. {opt adjustse} ({helpb stcox_ipwadh}) is not available with mi data.{p_end}
{phang}4. {cmd:mistyle(wide)}: {cmd:mi stset} refuses weights that are registered imputed or passive, and in a wide dataset an unregistered weight cannot differ between imputations ({cmd:mi update} resets the values of {it:m}>0 to the {it:m}=0 value).
Hazard ratios from a {cmd:wide} result need a loop over {cmd:mi extract} {it:m} and Rubin's rules.{p_end}
{phang}5. Memory: the data are copied once per imputation while the command runs, and the flong result stores (M+1) x pieces x kept variables. Use {opt keep()}. A full-size run (M=10) has not been timed yet.{p_end}
{phang}6. {cmd:mi estimate: stcox} with {opt strata()} is certified (see {help ipwadh_mi##using:Using the result}). A stratified split is exact only together with a Cox model stratified by the same variables.{p_end}
{phang}7. Wrapping the whole pipeline in {cmd:mi estimate} does not work (the split changes the number of observations); this is why the result is an mi dataset of the split data instead.{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Impute, compute the weights, switch to the result frame, estimate:{p_end}
{phang2}{cmd:. mi set wide}{p_end}
{phang2}{cmd:. mi register imputed biomarker abx}{p_end}
{phang2}{cmd:. mi impute chained (regress) biomarker (logit, augment) abx = arm age male adh i.tsi inf dead ltfu t1, add(10) rseed(1)}{p_end}
{phang2}{cmd:. ipwadh_mi, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) failure(inf) dead(dead) ltfu(ltfu) ///}{p_end}
{phang2}{cmd:      stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) frame(res) mistset(psw) keep(age male) change}{p_end}
{phang2}{cmd:. mi estimate: stcox i.arm c.age i.male, breslow vce(cluster id)}{p_end}
{phang2}{cmd:. mi estimate: mean age biomarker [pw=_ipw_ptw] if adh==1, over(arm)}{p_end}
{phang2}{cmd:. frame change default}{p_end}

{pstd}The same through {cmd:ipwadh} (the options are the same; {opt keep()} is {opt mikeep()} there):{p_end}
{phang2}{cmd:. ipwadh, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) failure(inf) dead(dead) ltfu(ltfu) ///}{p_end}
{phang2}{cmd:      stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) frame(res) mistset(psw) change}{p_end}

{pstd}Stratified by sex (split and Cox model):{p_end}
{phang2}{cmd:. ipwadh_mi, ... stcoxvars(c.age i.male) strata(male) frame(res) mistset(psw) change}{p_end}
{phang2}{cmd:. mi estimate: stcox i.arm c.age, strata(male) breslow vce(cluster id)}{p_end}

{pstd}Cheaper result for the table only:{p_end}
{phang2}{cmd:. ipwadh_mi, ... frame(tab, replace) mistyle(wide)}{p_end}
{phang2}{cmd:. frame change tab}{p_end}
{phang2}{cmd:. mi estimate: mean age [pw=_ipw_ptw] if adh==1, over(arm)}{p_end}


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
