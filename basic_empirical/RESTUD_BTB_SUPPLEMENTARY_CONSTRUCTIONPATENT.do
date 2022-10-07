*       FINAL ORGANIZATION OF NBER PATENT DATA
*       INPUT FILE:     pat76_06_assg / cite76_06 / ipcsic
*       OUTPUT FILE:    PATENT DATA 
*       DO FILE WRITTEN FOR STATA 10
*       DATA MAY BE STORED UNDER DIFFERENT PATH
*  DOWNLOAD ORIGINAL FILES FROM 2 SOURCES:
** MAIN SOURCE: https://sites.google.com/site/patentdataproject/Home
** ADDITIONAL SOURCE: 


*****************************
*****  ORGANIZE ALL CS  *****
*****************************

*  INITIALIZE & ORGANIZE DATA I
use pat76_06_assg.dta, replace
keep patent appyear gyear icl_class
drop if appyear==.
drop if gyear==.
duplicates drop patent, force
sort patent

*  MERGE
sort icl_class
merge m:1 icl_class using ipcsic.dta
drop icl_class
keep patent appyear gyear sic11
capture drop _merge

*  SAVE
sort patent
save temp0_updated_restud.dta, replace


********************************
*****  ORGANIZE GVT PANEL  *****
********************************

*  INITIALIZE & ORGANIZE DATA II
use pat76_06_assg.dta, replace
keep if asscode==6 | asscode==7
keep if appyear>1974
drop if appyear==.
duplicates drop patent, force
keep patent asscode appyear gyear country

*  DROP DUPLICATES
duplicates drop patent, force

*  REFORMAT
gen dummy=1
reshape wide dummy, i(patent) j(gyear)
foreach var of varlist dummy*{
replace `var'=0 if `var'==.
}
*
reshape long dummy, i(patent) j(year)
sort patent year
save panel_gvtfile_updated_restud.dta, replace


***************************************************************
*****  ORGANIZE ALL PANEL: GRANTED YEAR AND APPLIED YEAR  *****
***************************************************************

*  INITIALIZE & ORGANIZE DATA II: GRANT YEAR
use pat76_06_assg.dta, replace
keep if appyear>1974
drop if appyear==.
keep patent asscode appyear gyear country

*  DROP DUPLICATES
duplicates drop patent, force

*  REFORMAT
gen dummy=1
reshape wide dummy, i(patent) j(gyear)
foreach var of varlist dummy*{
replace `var'=0 if `var'==.
}
*
reshape long dummy, i(patent) j(year)
sort patent year
save panel_allfile_updatedgyear_restud.dta, replace

*  INITIALIZE & ORGANIZE DATA II: APPLICATION YEAR
use pat76_06_assg.dta, replace
keep if appyear>1974
drop if appyear==.
keep patent asscode appyear gyear country

*  DROP DUPLICATES
duplicates drop patent, force

*  REFORMAT
gen dummy=1
reshape wide dummy, i(patent) j(appyear)
foreach var of varlist dummy*{
replace `var'=0 if `var'==.
}
*
reshape long dummy, i(patent) j(year)
sort patent year
save panel_allfile_updatedappyear_restud.dta, replace


************************************************************************
*****  OBTAIN CITING PATENTS CITATIONS WITHIN X YEARS: GRANT YEAR  *****
************************************************************************

*  INITIALIZE CITING DATA
use cite76_06.dta, clear
keep cited citing
order cited citing
sort cited citing
rename citing patent
sort patent

*  MERGE ORIGINAL DATA AND CITATIONS TO GET BDAY OF CITING
merge patent using temp0_updated_restud.dta
keep if _merge==3
drop _m
rename gyear year
rename patent citing

*  PREPARE TO RETURN TO PANEL
sort cited year
rename cited patent

*  TOTAL CITATIONS PER YEAR
bysort patent year: gen N_citations=_N

*  GENERATE THE NUMBER OF DISTINCT 1 DIGIT SIC 
forvalues r=0(1)9{
	bysort patent year: gen ncites_mainsic`r'=_N if sic11==`r'
	replace ncites_mainsic`r'=0 if ncites_mainsic`r'==.
}
*

*  REFORMAT & SAVE
duplicates drop patent year, force
drop citing 
sort patent year
save citationperyear_all_updatedgranted_restud, replace

*  MERGE WITH PANEL
use panel_allfile_updatedgyear_restud.dta, replace
merge patent year using citationperyear_all_updatedgranted_restud.dta
drop if _merge==2

*  ARRANGE N_CITATIONS
** TOTAL CITATIONS
replace  N_citations=0 if _merge==1
** CITATIONS FROM MAIN SIC CLASSES
forvalues r=0(1)9{
	replace ncites_mainsic`r'=0 if _merge==1
}
*
drop _merge

*  GENERATE BDAY PATENT AND DISCARD PREPATENT HISTORY
gen tempcited_bday=year if dummy==1
bysort patent: egen cited_bday=max(tempcited_bday)
drop if cited_bday>year
drop tempcited_bday

*  GENERATE LAG
gen age_cohort=year-cited_bday

*  INITIALIZE LOOP ON AGE COHORTS FOR OVERALL CITATIONS
forvalues i=10(5)15 {

		*  INITIALIZE LOOP ON CITATIONS
		foreach var of varlist N_citations {

					*  CONSTRUCT VARIABLE
					bysort patent: egen tempTot`var'_`i'yrs=total(`var') if age_cohort<=`i'
					bysort patent: egen Tot`var'_`i'yrs=max(tempTot`var'_`i'yrs)
					drop tempTot`var'_`i'yrs
}
}
*

*  INITIALIZE LOOP ON AGE COHORTS ACCORDING TO MAIN SIC
forvalues i=10(5)10 {

		*  INITIALIZE LOOP ON CITATIONS
		foreach var of varlist ncites_mainsic0 ncites_mainsic1 ncites_mainsic2 ncites_mainsic3 ncites_mainsic4 ncites_mainsic5 ncites_mainsic6 ncites_mainsic7 ncites_mainsic8 ncites_mainsic9 {

					*  CONSTRUCT VARIABLE
					bysort patent: egen tempTot`var'_`i'yrs=total(`var') if age_cohort<=`i'
					bysort patent: egen Tot`var'_`i'yrs=max(tempTot`var'_`i'yrs)
					drop tempTot`var'_`i'yrs
}
}
*

*  SAVE
*save allcitingcitations_updatedgrantyear_restud.dta, replace
*use allcitingcitations_updatedgrantyear_restud.dta, replace

*  DROP DUPLICATES
duplicates drop patent, force

*  LOOP OVER YEAR HORIZONS
forvalues y=10(5)10{

	*  NUMBER OF DISTINCT SIC CLASSES OF CITATIONS
	forvalues r=0(1)9{
			*  GENERATE THE DUMMY
			gen Dcites_mainsic`r'_`y'yrs=(Totncites_mainsic`r'_`y'yrs!=0)
	}
	*

*  GENERATE THE ROWTOTAL
egen cites_nummainsic_`y'yrs= rowtotal(Dcites_mainsic0_`y'yrs Dcites_mainsic1_`y'yrs Dcites_mainsic2_`y'yrs Dcites_mainsic3_`y'yrs Dcites_mainsic4_`y'yrs Dcites_mainsic5_`y'yrs Dcites_mainsic6_`y'yrs Dcites_mainsic7_`y'yrs Dcites_mainsic8_`y'yrs Dcites_mainsic9_`y'yrs)

}
*

*  DROP VARIABLES
drop Dcites_mainsic*

*  SAVE
sort patent
save allcitingcitations_onerow_updatedgrantyear_restud.dta, replace


******************************************************************************
*****  OBTAIN CITING PATENTS CITATIONS WITHIN X YEARS: APPLICATION YEAR  *****
******************************************************************************

*  INITIALIZE CITING DATA
use cite76_06.dta, clear
keep cited citing
order cited citing
sort cited citing
rename citing patent
sort patent

*  MERGE ORIGINAL DATA AND CITATIONS TO GET BDAY OF CITING
merge patent using temp0_updated_restud.dta
keep if _merge==3
drop _m
rename appyear year
rename patent citing

*  PREPARE TO RETURN TO PANEL
sort cited year
rename cited patent

*  TOTAL CITATIONS PER YEAR
bysort patent year: gen N_citations=_N

*  GENERATE THE NUMBER OF DISTINCT 1 DIGIT SIC 
forvalues r=0(1)9{
	bysort patent year: gen ncites_mainsic`r'=_N if sic11==`r'
	replace ncites_mainsic`r'=0 if ncites_mainsic`r'==.
}
*

*  REFORMAT & SAVE
duplicates drop patent year, force
drop citing 
sort patent year
save citationperyear_all_updated_restud, replace

*  MERGE WITH PANEL
use panel_allfile_updatedappyear_restud.dta, replace
merge patent year using citationperyear_all_updated_restud.dta
drop if _merge==2

*  ARRANGE N_CITATIONS
** TOTAL CITATIONS
replace  N_citations=0 if _merge==1
** CITATIONS FROM MAIN SIC CLASSES
forvalues r=0(1)9{
	replace ncites_mainsic`r'=0 if _merge==1
}
*
drop _merge

*  GENERATE BDAY PATENT AND DISCARD PREPATENT HISTORY
gen tempcited_bday=year if dummy==1
bysort patent: egen cited_bday=max(tempcited_bday)
drop if cited_bday>year
drop tempcited_bday

*  GENERATE LAG
gen age_cohort=year-cited_bday

*  INITIALIZE LOOP ON AGE COHORTS
forvalues i=10(5)15 {

		*  INITIALIZE LOOP ON CITATIONS
		foreach var of varlist N_citations {

					*  CONSTRUCT VARIABLE
					bysort patent: egen tempTot`var'_`i'yrs=total(`var') if age_cohort<=`i'
					bysort patent: egen Tot`var'_`i'yrs=max(tempTot`var'_`i'yrs)
					drop tempTot`var'_`i'yrs
}
}
*

*  INITIALIZE LOOP ON AGE COHORTS
forvalues i=10(5)10 {

		*  INITIALIZE LOOP ON CITATIONS
		foreach var of varlist ncites_mainsic0 ncites_mainsic1 ncites_mainsic2 ncites_mainsic3 ncites_mainsic4 ncites_mainsic5 ncites_mainsic6 ncites_mainsic7 ncites_mainsic8 ncites_mainsic9 {

					*  CONSTRUCT VARIABLE
					bysort patent: egen tempTot`var'_`i'yrs=total(`var') if age_cohort<=`i'
					bysort patent: egen Tot`var'_`i'yrs=max(tempTot`var'_`i'yrs)
					drop tempTot`var'_`i'yrs
}
}
*

*  DROP DUPLICATES
duplicates drop patent, force

*  LOOP OVER YEAR HORIZONS
forvalues y=10(5)10{

	*  NUMBER OF DISTINCT SIC CLASSES OF CITATIONS
	forvalues r=0(1)9{
			*  GENERATE THE DUMMY
			gen Dcites_mainsic`r'_`y'yrs=(Totncites_mainsic`r'_`y'yrs!=0)
	}
	*

*  GENERATE THE ROWTOTAL
egen cites_nummainsic_`y'yrs= rowtotal(Dcites_mainsic0_`y'yrs Dcites_mainsic1_`y'yrs Dcites_mainsic2_`y'yrs Dcites_mainsic3_`y'yrs Dcites_mainsic4_`y'yrs Dcites_mainsic5_`y'yrs Dcites_mainsic6_`y'yrs Dcites_mainsic7_`y'yrs Dcites_mainsic8_`y'yrs Dcites_mainsic9_`y'yrs)

}
*

*  DROP VARIABLES
drop Dcites_mainsic*

*  SAVE
sort patent
save allcitingcitations_onerow_updatedappyear_restud.dta, replace


***************************************************************************************
*****  PLUG BACK TO ORIGINAL CITED TO OBTAIN AVERAGE FORWARD CITATIONS OF CITING  *****
***************************************************************************************

*  INITIALIZE CITING DATA
use cite76_06.dta, clear
order cited citing
sort cited citing
rename citing patent
sort patent

*  MERGE ORIGINAL DATA AND CITATIONS TO GET BDAY OF CITING
merge patent using temp0_updated_restud.dta
keep if _merge==3
drop _m
rename appyear year

*  MERGE TO FWD CITATION INFORMATION OF CITING
sort patent
merge patent using allcitingcitations_onerow_updatedappyear_restud.dta, keep(Tot*)
keep if _merge==3
drop _m

*  PREPARE TO RETURN TO PANEL
rename patent citing
sort cited year
rename cited patent

*  GET THE SIC 1 OF THE CITED
rename sic11 sic11_citing
merge patent using temp0_updated_restud.dta, keep(sic11)
rename sic11 sic11_cited
rename sic11_citing sic11 

*  AVERAGE CITATIONS OF CITING PATENTS
foreach var of varlist TotN_citations_*{
bysort patent year: egen AV_`var'=mean(`var')
bysort patent year: egen AV_`var'_SAME=mean(`var') if sic11==sic11_cited
bysort patent year: egen AV_`var'_SAMESIC1=max(AV_`var'_SAME)
drop AV_`var'_SAME
replace AV_`var'_SAMESIC1=. if AV_`var'==.
}
*

*  REFORMAT & SAVE
duplicates drop patent year, force
drop citing 
sort patent year
drop _merge
sort patent year
save citationperyear_all_updatedappyear_notrunc_restud, replace


