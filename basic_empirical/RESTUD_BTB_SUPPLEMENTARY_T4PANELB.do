
***********************************************************
*****  OBTAIN ALL CITING PATENTS NUMBER OF CITATIONS  *****
***********************************************************

*  MERGE WITH PANEL
use panel_allfile_updatedgyear_restud.dta, replace
merge patent year using citationperyear_all_updatedappyear_notrunc_restud.dta
drop if _merge==2

*  ARRANGE N_CITATIONS
foreach var of varlist TotN_citations_* AV_TotN_citations_*{
replace `var'=0 if _merge==1
}
drop _merge

*  GENERATE BDAY PATENT AND DISCARD PREPATENT HISTORY
gen tempcited_bday=year if dummy==1
bysort patent: egen cited_bday=max(tempcited_bday)
drop if cited_bday>year
drop tempcited_bday

*  GENERATE LAG
gen age_cohort=year-cited_bday

*  DES STATS
sum TotN_citations_15yrs AV_TotN_citations_15yrs if age_cohort==1

*  GENERATE SPLIT
** NO SPLIT
gen all=1
** SPLIT PRIVATE VERSUS PUBLIC (ALL)
gen public=.
replace public=0 if asscode==2 | asscode==3
replace public=1 if asscode==6 | asscode==7

*  KEEP FRANCE
keep if country=="FR"
keep if cited_bday>1974 &  cited_bday<1986

*  WHAT VARIABLES?
keep if age_cohort<11
keep AV_TotN_citations_10yrs AV_TotN_citations_10yrs_SAMESIC1 age_cohort public 


*  SAVE
save temp_france_restud.dta, replace



*****************************************************
*****  GRAPHICAL ANALYSIS - INITIALIZE PROGRAM  *****
*****************************************************

*  INITIALIZE
use temp_france_restud.dta, replace

*  PROGRAM MYDID2
capture program drop mydid2
program define mydid2

* `1' = variable of interest
* `2' = splitting variable
* `3' = time variable
* `4' = start date
* `5' = end date
* `6' = increment

capture drop d1 
capture drop d2 
capture drop d3 
capture drop ul 
capture drop ll 
capture drop k1 
capture drop k2 
capture drop k3 
capture drop ulk 
capture drop llk
capture drop diff 
capture drop diff_up 
capture drop diff_low 
capture drop year

clear matrix
forvalues a=`4'(`6')`5'{
		
		*  DES STAT
		ttest `1' if `3'==`a', by(`2') unequal
		*ranksum `1' if `3'==`a', by(`2') 

		*  OBTAIN AT
		matrix at`a'=`a'
		matrix list at`a'

		*  OBTAIN COEFFICIENT
		matrix ba`a'= r(mu_2)-r(mu_1)
		matrix list ba`a'
		
		*  OBTAIN SE
		matrix sea`a'=r(se) 
				
		*  INTO 1 MATRIX
		matrix da`a'=at`a',ba`a',sea`a'

		*  UPDATE THE FULL MATRIX
		if `a'==`4'{
		matrix C`a'=da`a'
		}
		else{
		local n=`a'-`6'
		matrix C`a' = C`n'\da`a'
}
}
*

*  CONSTRUCT VARIABLES
svmat C`5', names(d)
generate ul = d2 + 1.96*d3
generate ll = d2 - 1.96*d3
rename d2 diff
rename d1 year
rename ul diff_up
rename ll diff_low
end




*  INITIALIZE LOG FILE
capture log close
log using BTB_COOLDOWNSTATSDIFF_ALL_20180920,t replace
use temp_france, replace
** RUN THE MYDID2 PROGRAM
mydid2 AV_TotN_citations_10yrs public age_cohort 1 10 1 
** GRAPH
preserve
keep if year!=.
twoway (connected diff year, sort msize(large) lwidth(medthick) mcolor(red) msymbol(triangle) lcolor(red)) ///
	   (line diff_up year, sort lcolor(red) lwidth(medthin) lpattern(dash_dot_dot)) ///
	   (line diff_low year, sort lcolor(red) lwidth(medthin) lpattern(dash_dot_dot)), ///
	   ytitle(Average 10 Yr Forward Citations of Citing) xtitle(Age Cohort) title(Difference in Citation Patterns of Public and Private Patents) ///
	   legend(order(1 "Public Minus Private Patents")) ///
	   note("Note: French public and private patent between 75 and 85 with 10 year forward citations.")
	   gr export Fig1b_all.eps, replace
restore

*  INITIALIZE LOG FILE
capture log close
log using BTB_COOLDOWNSTATSDIFF_SAMESIC_20180920,t replace
use temp_france, replace
** RUN THE MYDID2 PROGRAM
mydid2 AV_TotN_citations_10yrs_SAMESIC1 public age_cohort 1 10 1 
** GRAPH
preserve
keep if year!=.
twoway (connected diff year, sort msize(large) lwidth(medthick) mcolor(red) msymbol(triangle) lcolor(red)) ///
	   (line diff_up year, sort lcolor(red) lwidth(medthin) lpattern(dash_dot_dot)) ///
	   (line diff_low year, sort lcolor(red) lwidth(medthin) lpattern(dash_dot_dot)), ///
	   ytitle(Average 10 Yr Forward Citations of Citing) xtitle(Age Cohort) title(Difference in Citation Patterns of Public and Private Patents) ///
	   legend(order(1 "Public Minus Private Patents")) ///
	   note("Note: French public and private patent between 75 and 85 with 10 year forward citations.")
	   gr export Fig1b_samesic.eps, replace
restore

