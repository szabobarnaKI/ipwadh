* ipwadh_check_simdata.do  -- Stata-side sanity checks of a Phase 0a dataset.  NOT YET RUN IN STATA (no Stata in the
* build container). The authoritative checks are in ipwadh_check_simdata.py; this file confirms the .dta files load
* with the right storage types and that the key structural properties hold in Stata itself.
* Usage: do ipwadh_check_simdata.do data/S3_N2000
args dir
use "`dir'/episodes.dta", clear
foreach v in t0 t1 age biomarker {
    local ty : type `v'
    assert "`ty'" == "double"
}
sort id rowno
by id: assert t0 == t1[_n-1] if _n > 1
by id: assert t0 == 0 if _n == 1
assert t1 >= t0
count if t1 == t0
di as txt "zero-length rows: " r(N)
by id: assert dead == 0 if _n < _N
by id: assert ltfu == 0 if _n < _N
assert inf + dead + ltfu <= 1
by id: gen long cuminf = sum(inf)
by id: assert ninf == cuminf - inf
* NOTE: (id, t0) is NOT unique for persons with a zero-length first episode; rowno is the row key.
duplicates report id t0
isid id rowno
di as res "ipwadh_check_simdata.do: all assertions passed for `dir'"
