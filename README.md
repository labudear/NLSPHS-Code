# NLSPHS-Code
use "C:\Users\Lava\Dropbox\Data\NLSPHS\Analysis\data\NLSPHSNACCHOARFAll4Waves_wts_peer_rurb7.dta", clear

/*
Impute missing population for HI in 2010 using:
http://health.hawaii.gov/kauai/files/2013/07/KAUAI-CHNA-July-2013.pdf
http://www.census.gov/quickfacts/table/PST045215/
*/

list nacchoid yearsurvey pop if state2014=="HI"

replace pop=185079 if pop==. & nacchoid=="HI002" & yearsurvey==2012 /*Hawaii (county) district health office*/
replace pop=154924 if pop==. & nacchoid=="HI005" & yearsurvey==2012 /*Maui county*/
replace pop=171191 if pop==. & nacchoid=="HI002" & yearsurvey==2006 /*Hawaii (county) district health office:  http://www.hawaii.edu/hivandaids/USA/HI/PH/State_of_Hawaii_Data_Book_2006.pdf*/
replace pop=141440 if pop==. & nacchoid=="HI005" & yearsurvey==2006 /*Maui county*/

gen nlsphssamp=.
replace nlsphssamp=1 if survresp!=.
replace nlsphssamp=0 if survresp==.
tab nlsphssamp

/*Checking the distribution of some variables and using natural log transformed variables in the GLLAMM model*/

log using "C:\Users\LAVA\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMIterations_LNTRANS3.smcl", replace 

hist bedpcap
hist mdpcap
hist pop
hist expcap
hist popdens
hist incpcap
hist fte 

sum bedpcap mdpcap pop expcap popdens incpcap fte

gen lnsqrtbedcap=log(sqrt(bedpcap))
gen lnmdpcap=log(mdpcap+0.0001)
gen lnpop=log(pop)
gen lnexpcap=log(expcap+0.0001)
gen lnpopdens=log(popdens+0.0001)
gen ftepcp=fte/pop*1000
gen lnftepcp=log(ftepcp)

hist lnsqrtbedcap
hist lnmdpcap
hist lnpop
hist lnexpcap
hist lnpopdens
hist lnftepcp

sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemp incpcap pct65
log close
/*GLLAMM Procedure using 12 integration points, Newton-Rhapson iterations and adaptive quadrature*/
/*

Newton-Raphson iterations to maximize the likelihood for a given number of mass-points.
Estimation of the conventional model assuming a normal exposure distribution required a large number of quadrature points.  
This is likely to be due to the integrands having sharp peaks that can easily fall between adjacent quadrature locations.  
Adaptive quadrature can overcome these problems and has been implemented in gllamm.

source: http://www.gllamm.org/statmod_03.pdf

*/
gen nacyear=.
replace nacyear=1997 if yearsurvey==1998
replace nacyear=2005 if yearsurvey==2006
replace nacyear=2010 if yearsurvey==2012
replace nacyear=2013 if yearsurvey==2014

tab nacyear
tab comp yearsurvey, col
tab clusttot yearsurvey, col
tab phtyp yearsurvey, col

/*Missing data imputation */
log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\GLAMM Accessories\GLLAMMIterations_misschk.smcl", replace 
misschk boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 if nacyear!=.
misschk boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 if nacyear!=. & nacyear!=1997
log close

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\GLAMM Accessories\GLLAMMIterations_imputation.smcl", replace 
ice boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 if nacyear!=. , ///
cmd(boh central:ologit) saving("C:\GLLAMM Accessories\GLLAMMimputed_dataset_naccho.dta", replace) m(100) seed(533265)
log close

use "C:\GLLAMM Accessories\GLLAMMimputed_dataset_naccho.dta", clear

sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemprate incpcap pct65
sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemprate incpcap pct65 if _mj==0
sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemprate incpcap pct65 if _mj==1
sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemprate incpcap pct65 if _mj==2

tab _mj comp

save fullIMP, replace

keep if _mj==1

save gllammMJ1, replace
keep state
duplicates drop state, force
count

gen statenum=_n
br
sort state
save statenum, replace


use gllammMJ1, clear
sort state
merge m:1 state using statenum
save gllammMJ1, replace
save "C:\GLLAMM Accessories\gllammMJ1.dta", replace

log using "C:\GLLAMM Accessories\GLLAMMIterations_LNTRANS_FTE_MI2_nip25.smcl", replace 

/*Random effect for state and fixed effect for year*/
xi:gllamm comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 i.nacyear, link(logit) fam(binom) i(statenum) nip(25) adapt 

gllapred prob, mu fsample /*We will be using this posterior probability in our further calculations*/ 


/*Generating confidence intervals of the predicted probabilities*/
ci_marg_mu LCLProb95 UCLProb95, level(95) reps(1000) dots

*********************************
***********************************

sort  nacyear statenum prob
qui by nacyear statenum: gen f=_n==1 /*Creates a dummy for a state by year (1st Obs)*/

/*Verify if the sorting is correct so that f was created correctly*/
br comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 nacyear state prob f

/*Mean values and SD for each state by year*/
bysort nacyear statenum: egen sumCOMP=sum(comp) 

list state statenum nacyear sumCOMP prob if f==1 & nacyear!=., clean

gen POPCOMP=prob*pop

bysort nacyear statenum : egen sumPOPCOMP=sum(POPCOMP) if prob!=.
bysort nacyear statenum : egen meanPOPCOMP=mean(POPCOMP) if prob!=.
bysort nacyear statenum : egen sumpop=sum(pop) if prob!=.
bysort nacyear statenum : egen meanpop=mean(pop) if prob!=.

gen sumPCTPOPCOMP=sumPOPCOMP/sumpop
gen meanPCTPOPCOMP=meanPOPCOMP/meanpop

br prob state statenum nacyear sumPCTPOPCOMP if f==1
list prob state statenum nacyear sumPCTPOPCOMP if f==1 & nacyear!=.

/*Confidence intervals of sumPCTPOPCOMP*/
gen POPCOMPLC95=LCLProb95*pop
bysort nacyear statenum : egen sumPOPCOMPLC95=sum(POPCOMPLC95)
gen sumPCTPOPCOMPLC95=sumPOPCOMPLC95/sumpop

gen POPCOMPUC95=UCLProb95*pop
bysort nacyear statenum : egen sumPOPCOMPUC95=sum(POPCOMPUC95)
gen sumPCTPOPCOMPUC95=sumPOPCOMPUC95/sumpop


/*Prediction at National Levels*/

xi:gllamm comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 i.nacyear, link(logit) fam(binom) i(statenum) nip(25) adapt 
gllapred probUSA, mu fsample /*We will be using this posterior probability in our further calculations*/ 

/*Generating confidence intervals of the predicted probabilities*/

ci_marg_mu LCLProbUSA95 UCLProbUSA95, level(95) reps(1000) dots

sort nacyear probUSA
qui by nacyear: gen f1=_n==1 /*Creates a dummy for a each year(1st Obs)*/

/*Verify if the sorting is correct so that f1 was created correctly*/
br comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 nacyear state probUSA f1

/*Mean values and SD for the nation by year*/
bysort nacyear: egen sumCOMPUSA=sum(comp) 

list state statenum nacyear sumCOMPUSA probUSA if f1==1 , clean

gen POPCOMPUSA=probUSA*pop

bysort nacyear: egen sumPOPCOMPUSA=sum(POPCOMPUSA) if probUSA!=.
bysort nacyear: egen sumpopUSA=sum(pop) if probUSA!=.

gen sumPCTPOPCOMPUSA=sumPOPCOMPUSA/sumpopUSA

br nacyear sumPCTPOPCOMPUSA if f1==1

/*Confidence intervals of sumPCTPOPCOMPUSA*/
gen POPCOMPLC95USA=LCLProbUSA95*pop
bysort nacyear: egen sumPOPCOMPLC95USA=sum(POPCOMPLC95USA)
gen sumPCTPOPCOMPUSALC95=sumPOPCOMPLC95USA/sumpopUSA

gen POPCOMPUC95USA=UCLProbUSA95*pop
bysort nacyear: egen sumPOPCOMPUC95USA=sum(POPCOMPUC95USA)
gen sumPCTPOPCOMPUSAUC95=sumPOPCOMPUC95USA/sumpopUSA

br nacyear sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 if f1==1

/*Significant differences between state and national average of % of population served by comprehensive PH system*/
gen sig=.
replace sig=1 if (sumPCTPOPCOMPLC95>sumPCTPOPCOMPUSAUC95 & sumPCTPOPCOMPUC95>sumPCTPOPCOMPUSAUC95)
replace sig=1 if (sumPCTPOPCOMPLC95<sumPCTPOPCOMPUSALC95 & sumPCTPOPCOMPUC95<sumPCTPOPCOMPUSALC95)
replace sig=1 if (sumPCTPOPCOMPLC95>sumPCTPOPCOMPUSALC95 & sumPCTPOPCOMPLC95>sumPCTPOPCOMPUSAUC95)
replace sig=1 if (sumPCTPOPCOMPUC95<sumPCTPOPCOMPUSALC95 & sumPCTPOPCOMPUC95<sumPCTPOPCOMPUSALC95)
replace sig=0 if sig==. & sumPCTPOPCOMP!=. & sumPCTPOPCOMPUSA!=.

gen higherLHD=.
replace higherLHD=1 if sumPCTPOPCOMP>sumPCTPOPCOMPUSA
replace higherLHD=0 if sumPCTPOPCOMP<=sumPCTPOPCOMPUSA
replace higherLHD=. if prob==.

log close

log using "C:\GLLAMM Accessories\sumPCTPOPCOMP_lntrans_fte_mi2_nip25.smcl", replace 
br nacyear yearsurvey state sumPCTPOPCOMP sumPCTPOPCOMPLC95 sumPCTPOPCOMPUC95 sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 sig higherLHD if f==1
list nacyear yearsurvey state sumPCTPOPCOMP sumPCTPOPCOMPLC95 sumPCTPOPCOMPUC95 sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 sig higherLHD f if f==1 & nacyear!=.

list nacyear yearsurvey state sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 f1 if f1==1 & nacyear!=.
log close

save "C:\GLLAMM Accessories\gllammfull_lntrans_fte_mi2_nip25.dta", replace

log using "C:\GLLAMM Accessories\compPROBdesc_lntrans_fte_mi2_nip25.smcl", replace 

/*Diagnostics*/
/*mean, median, variance, quintiles or deciles) for COMP vs predicted comp (prob) by year in sample only*/
bysort yearsurvey: sum comp prob if nlsphssamp==1, detail
 
pctile quintCOMP98 = comp if nlsphssamp==1 & yearsurvey==1998, nq(5) gen(quintcomp98) 
pctile quintCOMP06 = comp if nlsphssamp==1 & yearsurvey==2006, nq(5) gen(quintcomp06) 
pctile quintCOMP12 = comp if nlsphssamp==1 & yearsurvey==2012, nq(5) gen(quintcomp12) 
pctile quintCOMP14 = comp if nlsphssamp==1 & yearsurvey==2014, nq(5) gen(quintcomp14) 

pctile quintPROB98 = prob if nlsphssamp==1 & yearsurvey==1998, nq(5) gen(quintprob98) 
pctile quintPROB06 = prob if nlsphssamp==1 & yearsurvey==2006, nq(5) gen(quintprob06) 
pctile quintPROB12 = prob if nlsphssamp==1 & yearsurvey==2012, nq(5) gen(quintprob12) 
pctile quintPROB14 = prob if nlsphssamp==1 & yearsurvey==2014, nq(5) gen(quintprob14) 

list quint*98 in 1/5
list quint*06 in 1/5
list quint*12 in 1/5
list quint*14 in 1/5

pctile decCOMP98 = comp if nlsphssamp==1 & yearsurvey==1998, nq(10) gen(deccomp98) 
pctile decCOMP06 = comp if nlsphssamp==1 & yearsurvey==2006, nq(10) gen(deccomp06) 
pctile decCOMP12 = comp if nlsphssamp==1 & yearsurvey==2012, nq(10) gen(deccomp12) 
pctile decCOMP14 = comp if nlsphssamp==1 & yearsurvey==2014, nq(10) gen(deccomp14) 

pctile decPROB98 = prob if nlsphssamp==1 & yearsurvey==1998, nq(10) gen(decprob98) 
pctile decPROB06 = prob if nlsphssamp==1 & yearsurvey==2006, nq(10) gen(decprob06) 
pctile decPROB12 = prob if nlsphssamp==1 & yearsurvey==2012, nq(10) gen(decprob12) 
pctile decPROB14 = prob if nlsphssamp==1 & yearsurvey==2014, nq(10) gen(decprob14) 

list dec*98 in 1/10
list dec*06 in 1/10
list dec*12 in 1/10
list dec*14 in 1/10

log close

log using "C:\GLLAMM Accessories\countProbnmis_lntrans_fte_mi2_nip25.smcl", replace 
bysort nacyear state: count if prob!=.
bysort nacyear state: count
log close

log using "C:\GLLAMM Accessories\countvar_lntrans_fte_mi2_nip25.smcl", replace 
bysort nacyear state: mdesc comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65
log close

view "C:\GLLAMM Accessories\countvar_lntrans_fte_mi2_nip25.smcl"


log using "C:\GLLAMM Accessories\countPrednmis_lntrans_fte_mi2_nip25.smcl", replace 
gen missfromod=.
replace missfromod=0 if boh!=. & lnsqrtbedcap!=. & lnmdpcap!=. & ///
 lnpop!=. & lnftepcp!=. & lnpopdens!=. & central!=. & pctnonwh!=. & ///
 unemp!=. & incpcap!=. & pct65!=.
replace missfromod=1 if boh==. | lnsqrtbedcap==. | lnmdpcap==. | ///
 lnpop==. | lnftepcp==. | lnpopdens==. | central==. | pctnonwh==. | ///
 unemp==. | incpcap==. | pct65==.

tab missfromod
tab comp

bysort nacyear state : tab missfromod
bysort nacyear state : tab comp, missing

log close

view "C:\GLLAMM Accessories\countProbnmis_lntrans_fte_mi2_nip25.smcl"

log using "C:\GLLAMM Accessories\countrnmis_lntrans_fte_mi2_nip25.smcl", replace 
 /*findit rmiss*/
egen nmis=rmiss2(boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65)
tab nmis
bysort nacyear state : tab nmis
log close

view "C:\GLLAMM Accessories\countrnmis_lntrans_fte_mi2_nip25.smcl"


br comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 statenum state nacyear prob nmis //if missfromod==1

sort nacyear state nacchoid

/*Exporting the output from GLLAMM*/
keep nacchoid state statenum nacyear arm sumCOMP pop lnpop prob probUSA POPCOMP sumPOPCOMP sumpop sumPCTPOPCOMP meanPOPCOMP meanpop meanPCTPOPCOMP LCLProb95 UCLProb95 sumPCTPOPCOMPLC95 sumPCTPOPCOMPUC95 sumPCTPOPCOMPUSA LCLProbUSA95 UCLProbUSA95 sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 sig higherLHD f f1 nmis
*export delimited using "X:\xDATA\NLSPHS 2014\Analysis\data\Lava\Analytical data files\NLSPHSSNAbyYEAR\GLLAMMOUTPUT_FINAL_LNTRANS.csv", replace
export delimited using "C:\GLLAMM Accessories\GLLAMMOUTPUT_FINAL_LNTRANS_FTE_MI2_NIP25.csv", replace

/*Exporting the output by statenum and nacyear*/
keep state statenum nacyear sumCOMP sumPOPCOMP sumpop prob LCLProb95 UCLProb95 sumPCTPOPCOMP sumPCTPOPCOMPLC95 sumPCTPOPCOMPUC95 probUSA LCLProbUSA95 UCLProbUSA95 sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 sig higherLHD f f1 nmis
keep if f==1 | f1==1
br statenum state nacyear prob sumPCTPOPCOMP probUSA sumPCTPOPCOMPUSA nmis f f1

*export delimited using "X:\xDATA\NLSPHS 2014\Analysis\data\Lava\Analytical data files\NLSPHSSNAbyYEAR\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS.csv", replace
export delimited using "C:\GLLAMM Accessories\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS_FTE_MI2_NIP25.csv", replace


/*Use the following data file for creating state level maps*/
save "C:\GLLAMM Accessories\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS_FTE_MI2_nip25.dta", replace

view "C:\GLLAMM Accessories\GLLAMMIterations_LNTRANS_FTE_MI2_nip25.smcl"
*view "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMIterations_LNTRANS_FTE_MI.smcl"





























/*






/*Missing data imputation */
log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\GLAMM Accessories\GLLAMMIterations_misschk3.smcl", replace 
misschk boh lnsqrtbedcap lnmdpcap lnpop lnexpcap lnftepcp lnpopdens central pctnonwh unemp incpcap pct65 if nacyear!=.
*misschk boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65 if nacyear!=. & nacyear!=1997
*misschk boh lnsqrtbedcap lnmdpcap lnpop lnexpcap lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65 if nacyear==2010 | nacyear==2013
log close

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMIterations_imputation3.smcl", replace 
ice boh lnsqrtbedcap lnmdpcap lnpop lnexpcap lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65 if nacyear!=. , ///
cmd(boh central:ologit) saving("C:\Users\Lava\Documents\Lava\GLLAMMimputed_dataset_naccho3.dta", replace) m(200) seed(533265)
log close

use "C:\Users\Lava\Documents\Lava\GLLAMMimputed_dataset_naccho3.dta", clear

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\Summary_imputation3.smcl", replace 
sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemprate incpcap pct65
sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemprate incpcap pct65 if _mj==0
sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemprate incpcap pct65 if _mj==1
sum lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens lnftepcp pctnonwh unemprate incpcap pct65 if _mj==2
log close

*save fullIMP, replace

keep  if _mj==1

count

save "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\ImputedMJ1.dta", replace

br if state==""
*replace state=stabr_f if state==""

save gllammMJ1, replace
keep state
duplicates drop state, force
count

gen statenum=_n
br
sort state
save statenum, replace

/*
keep state 
sort state
quietly by state:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup
*/


use gllammMJ1, clear
drop _merge statenum
sort state
merge m:1 state using statenum
save gllammMJ1, replace
save "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\ImputedMJ1.dta", replace
save "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\ImputedMJ1.dta", replace
saveold "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\ImputedMJ1_old.dta", replace

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMIterations_LNTRANS_FTE_MI3.smcl", replace 

/*Random effect for state and fixed effect for year*/
xi:gllamm comp boh lnsqrtbedcap lnmdpcap lnpop lnexpcap lnpopdens central pctnonwh unemprate incpcap pct65 i.nacyear, link(logit) fam(binom) i(statenum) nip(12) adapt 

/*
*gllapred c, cooksd /*Influence of top level unit on likelihood*/
/*Raw and standardised residuals*/
gllapred res_, u       /* posterior mean and sd in res_m1 res_s1  */
gllapred stres_, ustd /*also called diagnostic standard errors*/

gllapred p,xb /*the linear predictor for the fixed effects is returned for those who were in NLSPHS sample only*/
gllapred prob2, xb fsample /*predict using the full sample*/

gllapred eb, u /* the posterior means and standard deviation of the latent variables or random effects*/

/*
The expectation of the response (see mu option). By default,
the expectation is with respect to the posterior distribution
of the latent variables, but marginal option gives the expectation 
with respect to the prior distribution
*/
gllapred prob3, mu marginal fsample /*Predict marginal probability that comp=1 (with respect to random effects) for full sample*/

gllapred prob4, mu marginal /*Predict marginal probability that comp=1 (with respect to random effects) for those in NLSPHS sample*/

gllapred prob1, mu
*/
gllapred prob, mu fsample /*We will be using this posterior probability in our further calculations*/ 


/*Generating confidence intervals of the predicted probabilities*/
/*Package need for postgllamm command*/
/*
ssc describe ci_marg_mu
ssc install ci_marg_mu
*/

ci_marg_mu LCLProb95 UCLProb95, level(95) reps(1000) dots

*********************************
***********************************

sort  nacyear statenum prob
qui by nacyear statenum: gen f=_n==1 /*Creates a dummy for a state by year (1st Obs)*/

/*Verify if the sorting is correct so that f was created correctly*/
br comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65 nacyear state prob f

/*Mean values and SD for each state by year*/
bysort nacyear statenum: egen sumCOMP=sum(comp) 

list state statenum nacyear sumCOMP prob if f==1 & nacyear!=., clean

gen POPCOMP=prob*pop

bysort nacyear statenum : egen sumPOPCOMP=sum(POPCOMP) if prob!=.
bysort nacyear statenum : egen meanPOPCOMP=mean(POPCOMP) if prob!=.
bysort nacyear statenum : egen sumpop=sum(pop) if prob!=.
bysort nacyear statenum : egen meanpop=mean(pop) if prob!=.

gen sumPCTPOPCOMP=sumPOPCOMP/sumpop
gen meanPCTPOPCOMP=meanPOPCOMP/meanpop

br prob state statenum nacyear sumPCTPOPCOMP if f==1
list prob state statenum nacyear sumPCTPOPCOMP if f==1 & nacyear!=.

/*Confidence intervals of sumPCTPOPCOMP*/
gen POPCOMPLC95=LCLProb95*pop
bysort nacyear statenum : egen sumPOPCOMPLC95=sum(POPCOMPLC95)
gen sumPCTPOPCOMPLC95=sumPOPCOMPLC95/sumpop

gen POPCOMPUC95=UCLProb95*pop
bysort nacyear statenum : egen sumPOPCOMPUC95=sum(POPCOMPUC95)
gen sumPCTPOPCOMPUC95=sumPOPCOMPUC95/sumpop


/*Prediction at National Levels*/

xi:gllamm comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65 i.nacyear, link(logit) fam(binom) i(statenum) nip(12) adapt 
gllapred probUSA, mu fsample /*We will be using this posterior probability in our further calculations*/ 

/*Generating confidence intervals of the predicted probabilities*/

ci_marg_mu LCLProbUSA95 UCLProbUSA95, level(95) reps(1000) dots

sort nacyear probUSA
qui by nacyear: gen f1=_n==1 /*Creates a dummy for a each year(1st Obs)*/

/*Verify if the sorting is correct so that f1 was created correctly*/
br comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65 nacyear state probUSA f1

/*Mean values and SD for the nation by year*/
bysort nacyear: egen sumCOMPUSA=sum(comp) 

list state statenum nacyear sumCOMPUSA probUSA if f1==1 , clean

gen POPCOMPUSA=probUSA*pop

bysort nacyear: egen sumPOPCOMPUSA=sum(POPCOMPUSA) if probUSA!=.
bysort nacyear: egen sumpopUSA=sum(pop) if probUSA!=.

gen sumPCTPOPCOMPUSA=sumPOPCOMPUSA/sumpopUSA

br nacyear sumPCTPOPCOMPUSA if f1==1

/*Confidence intervals of sumPCTPOPCOMPUSA*/
gen POPCOMPLC95USA=LCLProbUSA95*pop
bysort nacyear: egen sumPOPCOMPLC95USA=sum(POPCOMPLC95USA)
gen sumPCTPOPCOMPUSALC95=sumPOPCOMPLC95USA/sumpopUSA

gen POPCOMPUC95USA=UCLProbUSA95*pop
bysort nacyear: egen sumPOPCOMPUC95USA=sum(POPCOMPUC95USA)
gen sumPCTPOPCOMPUSAUC95=sumPOPCOMPUC95USA/sumpopUSA

br nacyear sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 if f1==1

/*Significant differences between state and national average of % of population served by comprehensive PH system*/
gen sig=.
replace sig=1 if (sumPCTPOPCOMPLC95>sumPCTPOPCOMPUSALC95 & sumPCTPOPCOMPLC95>sumPCTPOPCOMPUSAUC95)
replace sig=1 if (sumPCTPOPCOMPLC95<sumPCTPOPCOMPUSALC95 & sumPCTPOPCOMPLC95<sumPCTPOPCOMPUSAUC95)
replace sig=0 if sig==. & sumPCTPOPCOMP!=. & sumPCTPOPCOMPUSA!=.

gen higherLHD=.
replace higherLHD=1 if sumPCTPOPCOMP>sumPCTPOPCOMPUSA
replace higherLHD=0 if sumPCTPOPCOMP<=sumPCTPOPCOMPUSA
replace higherLHD=. if prob==.

log close

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\sumPCTPOPCOMP_lntrans_fte_mi3.smcl", replace 
br nacyear yearsurvey state sumPCTPOPCOMP sumPCTPOPCOMPLC95 sumPCTPOPCOMPUC95 sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 sig higherLHD if f==1
list nacyear yearsurvey state sumPCTPOPCOMP sumPCTPOPCOMPLC95 sumPCTPOPCOMPUC95 sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 sig higherLHD f if f==1 & nacyear!=.

list nacyear yearsurvey state sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 f1 if f1==1 & nacyear!=.
log close

*save "X:\xDATA\NLSPHS 2014\Analysis\data\Lava\Analytical data files\NLSPHSSNAbyYEAR\gllammfull_lntrans.dta", replace
save "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\gllammfull_lntrans_fte_mi3.dta", replace
save "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\gllammfull_lntrans_fte_mi3.dta", replace
saveold "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\gllammfull_lntrans_fte_mi3_old.dta", replace

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\compPROBdesc_lntrans_fte_mi3.smcl", replace 

/*Diagnostics*/
/*mean, median, variance, quintiles or deciles) for COMP vs predicted comp (prob) by year in sample only*/
bysort yearsurvey: sum comp prob if nlsphssamp==1, detail
 
pctile quintCOMP98 = comp if nlsphssamp==1 & yearsurvey==1998, nq(5) gen(quintcomp98) 
pctile quintCOMP06 = comp if nlsphssamp==1 & yearsurvey==2006, nq(5) gen(quintcomp06) 
pctile quintCOMP12 = comp if nlsphssamp==1 & yearsurvey==2012, nq(5) gen(quintcomp12) 
pctile quintCOMP14 = comp if nlsphssamp==1 & yearsurvey==2014, nq(5) gen(quintcomp14) 

pctile quintPROB98 = prob if nlsphssamp==1 & yearsurvey==1998, nq(5) gen(quintprob98) 
pctile quintPROB06 = prob if nlsphssamp==1 & yearsurvey==2006, nq(5) gen(quintprob06) 
pctile quintPROB12 = prob if nlsphssamp==1 & yearsurvey==2012, nq(5) gen(quintprob12) 
pctile quintPROB14 = prob if nlsphssamp==1 & yearsurvey==2014, nq(5) gen(quintprob14) 

list quint*98 in 1/5
list quint*06 in 1/5
list quint*12 in 1/5
list quint*14 in 1/5

pctile decCOMP98 = comp if nlsphssamp==1 & yearsurvey==1998, nq(10) gen(deccomp98) 
pctile decCOMP06 = comp if nlsphssamp==1 & yearsurvey==2006, nq(10) gen(deccomp06) 
pctile decCOMP12 = comp if nlsphssamp==1 & yearsurvey==2012, nq(10) gen(deccomp12) 
pctile decCOMP14 = comp if nlsphssamp==1 & yearsurvey==2014, nq(10) gen(deccomp14) 

pctile decPROB98 = prob if nlsphssamp==1 & yearsurvey==1998, nq(10) gen(decprob98) 
pctile decPROB06 = prob if nlsphssamp==1 & yearsurvey==2006, nq(10) gen(decprob06) 
pctile decPROB12 = prob if nlsphssamp==1 & yearsurvey==2012, nq(10) gen(decprob12) 
pctile decPROB14 = prob if nlsphssamp==1 & yearsurvey==2014, nq(10) gen(decprob14) 

list dec*98 in 1/10
list dec*06 in 1/10
list dec*12 in 1/10
list dec*14 in 1/10

log close

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\countProbnmis_lntrans_fte_mi3.smcl", replace 
bysort nacyear state: count if prob!=.
bysort nacyear state: count
log close

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\countvar_lntrans_fte_mi3.smcl", replace 
bysort nacyear state: mdesc comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65
log close

view "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\countvar_lntrans_fte_mi3.smcl"


log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\countPrednmis_lntrans_fte_mi3.smcl", replace 
gen missfromod=.
replace missfromod=0 if boh!=. & lnsqrtbedcap!=. & lnmdpcap!=. & ///
 lnpop!=. & lnftepcp!=. & lnpopdens!=. & central!=. & pctnonwh!=. & ///
 unemprate!=. & incpcap!=. & pct65!=.
replace missfromod=1 if boh==. | lnsqrtbedcap==. | lnmdpcap==. | ///
 lnpop==. | lnftepcp==. | lnpopdens==. | central==. | pctnonwh==. | ///
 unemprate==. | incpcap==. | pct65==.

tab missfromod
tab comp

bysort nacyear state : tab missfromod
bysort nacyear state : tab comp, missing

log close

view "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\countProbnmis_lntrans_fte_mi3.smcl"

log using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\countrnmis_lntrans_fte_mi3.smcl", replace 
 /*findit rmiss*/
egen nmis=rmiss2(boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65)
tab nmis
bysort nacyear state : tab nmis
log close

view "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\countrnmis_lntrans_fte_mi3.smcl"


br comp boh lnsqrtbedcap lnmdpcap lnpop lnftepcp lnpopdens central pctnonwh unemprate incpcap pct65 statenum state nacyear prob nmis //if missfromod==1

sort nacyear state nacchoid

/*Exporting the output from GLLAMM*/
keep nacchoid state statenum nacyear arm0 sumCOMP pop lnpop prob probUSA POPCOMP sumPOPCOMP sumpop sumPCTPOPCOMP meanPOPCOMP meanpop meanPCTPOPCOMP LCLProb95 UCLProb95 sumPCTPOPCOMPLC95 sumPCTPOPCOMPUC95 sumPCTPOPCOMPUSA LCLProbUSA95 UCLProbUSA95 sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 sig higherLHD f f1 nmis
*export delimited using "X:\xDATA\NLSPHS 2014\Analysis\data\Lava\Analytical data files\NLSPHSSNAbyYEAR\GLLAMMOUTPUT_FINAL_LNTRANS.csv", replace
export delimited using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMOUTPUT_FINAL_LNTRANS_FTE_MI3.csv", replace
export delimited using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\GLLAMMOUTPUT_FINAL_LNTRANS_FTE_MI3.csv", replace

/*Exporting the output by statenum and nacyear*/
keep state statenum nacyear sumCOMP sumPOPCOMP sumpop prob LCLProb95 UCLProb95 sumPCTPOPCOMP sumPCTPOPCOMPLC95 sumPCTPOPCOMPUC95 probUSA LCLProbUSA95 UCLProbUSA95 sumPCTPOPCOMPUSA sumPCTPOPCOMPUSALC95 sumPCTPOPCOMPUSAUC95 sig higherLHD f f1 nmis
keep if f==1 | f1==1
br statenum state nacyear prob sumPCTPOPCOMP probUSA sumPCTPOPCOMPUSA nmis f f1

*export delimited using "X:\xDATA\NLSPHS 2014\Analysis\data\Lava\Analytical data files\NLSPHSSNAbyYEAR\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS.csv", replace
export delimited using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS_FTE_MI3.csv", replace
export delimited using "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS_FTE_MI3.csv", replace


/*Use the following data file for creating state level maps*/
save "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS_FTE_MI3.dta", replace
save "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS_FTE_MI3.dta", replace
saveold "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\CoH Metrics\GLLAMMOUTPUT_PCTPOPCOMP_LNTRANS_FTE_MI3_old.dta", replace

view "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMIterations_LNTRANS_FTE_MI3.smcl"
view "C:\Users\Lava\Dropbox\NLSPHS 2014\ALL4Waves\GLLAMMIterations_LNTRANS_FTE_MI3.smcl"

