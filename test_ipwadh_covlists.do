* test_ipwadh_covlists.do -- run in Stata 19 (folder with ipwadh_covlists.ado and data/).
* v2: group 1 passed in Stata (2026-10-09); removed a stray comma before -quiet- in group 2 (syntax error, not a command bug).
clear all
use "data/S3_N200/episodes.dta", clear

* 1. stcox list goes to both models; predictor that duplicates a stcox term is not repeated; c.age == age
ipwadh_covlists, arm(arm) adherence(adh) stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker age)
assert "`r(numvars)'" == "c.age i.male"
assert "`r(denvars)'" == "i.tsi i.abx c.biomarker"
assert "`r(cnumvars)'" == "c.age i.male" & "`r(cdenvars)'" == "i.tsi i.abx c.biomarker"
assert "`r(stcox_hint)'" == "i.arm c.age i.male"
assert "`r(ignored_dup)'" == "age"

* 2. arm typed in any notation is ignored (arm is automatic); arm interaction kept without byarm
ipwadh_covlists, arm(arm) adherence(adh) stcoxvars(i.arm c.age ib1.arm) adhpredictors(c.arm) quiet
assert "`r(numvars)'" == "c.age" & "`r(denvars)'" == ""
assert "`r(ignored_arm)'" == "i.arm ib1.arm c.arm"
ipwadh_covlists, arm(arm) adherence(adh) stcoxvars(i.arm#c.age c.age) quiet
assert "`r(numvars)'" == "i.arm#c.age c.age"

* 3. byarm: every term involving arm is ignored (the models are fitted separately by arm)
ipwadh_covlists, arm(arm) adherence(adh) stcoxvars(i.arm#c.age c.age) byarm quiet
assert "`r(numvars)'" == "c.age" & "`r(ignored_arm)'" == "i.arm#c.age"

* 4. extra death/loss predictors only enter the censoring full model
ipwadh_covlists, arm(arm) adherence(adh) stcoxvars(c.age) adhpredictors(i.tsi) censorpredictors(i.abx) quiet
assert "`r(denvars)'" == "i.tsi" & "`r(cdenvars)'" == "i.tsi i.abx"

* 5. errors: unknown variable (198), missing values (416) unless allowmissing
capture ipwadh_covlists, arm(arm) adherence(adh) stcoxvars(nosuchvar) quiet
assert _rc == 198
use "data/S8_N200/episodes.dta", clear                 // S8 has missing biomarker / abx
capture ipwadh_covlists, arm(arm) adherence(adh) stcoxvars(c.age) adhpredictors(c.biomarker i.abx) quiet
assert _rc == 416
ipwadh_covlists, arm(arm) adherence(adh) stcoxvars(c.age) adhpredictors(c.biomarker) allowmissing quiet
assert "`r(denvars)'" == "c.biomarker"
di as res "test_ipwadh_covlists.do: all assertions passed"
