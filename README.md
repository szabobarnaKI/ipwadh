# ipwadh
A Stata command to calculate inverse probability of treatment adherence weights, to allow per-protocol analysis of randomized trials excluding treatment non-adherence periods.

ipwadh calculates inverse probability of adherence weights to compensate for the imbalance caused by the exclusion of non-adherent follow-up periods. It allows stratified analysis, time-varying covariates and both single event and repeated events data. The command can handle mi imputed data as well. After the calculation of IPW, the per-protocol model can be fit using stcox. The help file provides detailed description on how to use ipwadh as well as the correct setup of the stcox command, and model assumptions to check.

To install it, run the following command in Stata:
##
    net from "https://szabobarnaki.github.io/ipwadh/"
