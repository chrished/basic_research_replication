*       FINAL ORGANIZATION & ANALYSIS OF MULTI MARKET ACTIVITY AND BASIC R&D
*       INPUT FILE:     LIFI, EAE, BRN, R&D, PATENT DATA
*       OUTPUT FILE:    FIRM LEVEL DATA FOR ANALYSIS
*       DO FILE WRITTEN FOR STATA 10
*       DATA MAY BE STORED UNDER DIFFERENT PATH

*  INITIALIZE SET UP
capture log close
log using RESTUD_BTB_SUPPLEMENTARY_CONSTRUCTIONFIRM,t replace
di "log file printed on $S_DATE at $S_TIME"
set more off
clear all
macro drop _all
set mem 15g
set matsize 800
version 10


								**************************************************************************
								*****  DATA ORGANIZATION: BALANCE SHEET AND MULTI-INDUSTRY PRESENCE  *****
								**************************************************************************

*  USE AND INITIALIZE LOOP
use lifi_firms_100, replace
forvalues i=100/106 {


**********************************
*****  MULTIDIVISIONAL DATA  *****
**********************************

*  ORGANIZE DATA
use eaebr`i', replace
sort siren 
save eaebr`i', replace


*****************************
*****  OWNERSHIP DATA  *****
*****************************

*  INITIALIZE DATA
use lifi_firms_`i', replace

*  SELECT AGAIN MAJORITY OWNERSHIP
drop if txcontr<5001 | txcontr==.
duplicates tag sirlifi, gen(tag)
tab tag
drop tag

*  FOREIGN BUSINESS GROUPS
** DESCRIPTIVE STATISTICS
tab paystg
** GENERATE INDICATORS ON NATIONALITY OF HEAD OF BG
gen foreign_head=.
replace foreign_head=1 if paystg!="061"
gen US_head=.
replace US_head=1 if paystg=="404"
gen UK_head=.
replace UK_head=1 if paystg=="132"
gen GER_head=.
replace GER_head=1 if paystg=="109"

*  STATE OWNED BG
gen state_owned=0
replace state_owned=1 if natact=="2"

*  RENAME
rename sirlifi siren

*  ORGANIZE & SAVE
sort siren
keep siren sirtg txcontr rstg foreign_head US_head UK_head GER_head state_owned
save temp_restud_bg1_`i', replace


********************************************
*****  OWNERSHIP & BALANCE SHEET DATA  *****
********************************************

*  MERGE OWNERSHIP INFO WITH BALANCE SHEET
use brn`i'new, replace
sort siren
merge siren using temp_restud_bg1_`i'
rename _merge merge_lifibrn
duplicates drop siren, force

*  MERGE NEW DATA WITH MULTIDIVISIONAL DATA
sort siren
merge siren using eaebr`i', keep (s003  apbr)
rename _merge merge_eaebr

*  DROP ALL LIFI FIRMS THAT WE CANNOT FIND IN BALANCE SHEETS
drop if merge_lifibrn==2

*  KEEP MAIN VARIABLES
keep siren merge_lifibrn naf600 apbr merge_eaebr sirtg ebe immob effec q s003 r216 r217 va foreign_head US_head UK_head GER_head state_owned

*  MAIN ACTIVITY
gen main_activity=0
replace main_activity=1 if naf600==apbr
replace main_activity=1 if merge_eaebr==1
sort siren main_activity

*  UNIQUE IDENTIFIER
sort siren main_activity
gen uniq_siren=1
bysort siren: replace uniq_siren=. if siren[_n]==siren[_n+1]
count if uniq_siren==1
distinct siren

*  VIRTUAL BG IDENTIFIER FOR FIRMS NOT AFFILIATED
gen var1="FIRM"
gen siren_alt=var1+siren
gen sirtg_alt=sirtg
replace sirtg_alt=siren_alt if sirtg==""
rename sirtg sirtg_old
rename sirtg_alt sirtg

*  OUTLIERS
** 5X INTERQUARTILE RANGE
gen prof=ebe/immob
sum prof,d
gen out_mean_prof=0 
_pctile prof, percentiles(1 5 25 50 75 95 99)
replace out_mean_prof=1 if abs(prof-r(r4) )> 5*(r(r5)-r(r3))
** DROP EXTREME VALUES
drop if out_mean_prof==1

*  DROP NO EMPLOYMENT FIRMS/BG
drop if effec==0

*  DROP NEGATIVE OR O SALES
drop if q<=0 & q!=.

*  DROP NON CLASSIFIED FIRMS
drop if naf600==";" | naf600=="N"

*  ESTABLISH MARKET CLASSIFICATIONS
** MAIN ACTIVITY CLASSIFICATION
gen marketmain_4SIC=substr(naf600,1,4)
gen marketmain_3SIC=substr(naf600,1,3)
gen marketmain_2SIC=substr(naf600,1,2)
gen marketmain_1SIC=substr(naf600,1,1)
** ALL ACTIVITIES
gen market_4SIC=apbr
replace market_4SIC=naf600 if market_4SIC==""
gen market_3SIC=substr(apbr,1,3)
replace market_3SIC=marketmain_3SIC if market_3SIC==""
gen market_2SIC=substr(apbr,1,2)
replace market_2SIC=marketmain_2SIC if market_2SIC==""
gen market_1SIC=substr(apbr,1,1)
replace market_1SIC=marketmain_1SIC if market_1SIC==""

*  FINANCIAL FIRMS IDENTIFIER
gen financial=0
replace financial=1 if market_2SIC=="65" | market_2SIC=="66" | market_2SIC=="67"
*  WHAT ARE THE EXCLUDED CATEGORIES?
*  65.1A Banque centrale
*  65.1C Banques
*  65.1D Banques mutualistes
*  65.1E Caisses d'épargne
*  65.1F Intermédiations monétaires n.c.a.
*  65.2A Crédit-bail
*  65.2C Distribution de crédit
*  65.2E Organismes de placement en valeurs mobilières
*  65.2F Intermédiations financières diverses
*  66.0A Assurance-vie et capitalisation
*  66.0C Caisses de retraite
*  66.0E Assurance dommages
*  66.0F Réassurance
*  66.0G Assurance relevant du code de la mutualité
*  67.1A Administration de marchés financiers
*  67.1C Gestion de portefeuilles
*  67.1E Autres auxiliaires financiers
*  67.2Z Auxiliaires d'assurance


******************************************************
*****  NUMBER OF PRODUCT MARKETS AND INDUSTRIES  *****
******************************************************

*  START LOOP
forvalues s=1(1)4{
*  1. COUNT UNIQUE MARKETS FOR EACH SIC LEVEL
quietly unique market_`s'SIC, by(sirtg) gen(nber1_`s'SICtemp)
quietly unique market_`s'SIC if s003>9 | s003==., by(sirtg) gen(nber2_`s'SICtemp)
quietly unique marketmain_`s'SIC if uniq_siren==1, by(sirtg) gen(nber3_`s'SICtemp)
quietly unique marketmain_`s'SIC if uniq_siren==1 & financial!=1, by(sirtg) gen(nber4_`s'SICtemp)
*  2. PUT INFORMATION INTO 1 LINE 
egen nber1_`s'SIC=mean(nber1_`s'SICtemp), by(sirtg)
egen nber2_`s'SIC=mean(nber2_`s'SICtemp), by(sirtg)
egen nber3_`s'SIC=mean(nber3_`s'SICtemp), by(sirtg)
egen nber4_`s'SIC=mean(nber4_`s'SICtemp), by(sirtg)
drop nber1_`s'SICtemp nber2_`s'SICtemp nber3_`s'SICtemp nber4_`s'SICtemp
*  3. REPLACE 0 BY 1 IN COUNT OF PRODUCT MARKETS
** 1ST & 3RD MEASURE NO 0
** 2ND MEASURE: HAD TOO LITTLE IN EACH BRANCH
replace nber2_`s'SIC=1 if nber2_`s'SIC==0
** 4TH MEASURE: ALL BANKS
replace nber4_`s'SIC=1 if nber2_`s'SIC==0
}
*END LOOP
save temp_restud_bg2_`i', replace


*****************************************
*****  FIRM AND BG CHARACTERISTICS  *****
*****************************************

*  FINISHED WITH MULTIDIVISIONAL DATA
duplicates drop siren, force

*  FRANCS INTO EUROS
local variables "r216 r217 q va ebe immob"
foreach x of local variables{
		replace `x'=`x'/6.55957 if `i'==100
}
*

*  HQ LOCATION
replace foreign_head=0 if foreign_head==.
replace US_head=0 if US_head==.
replace UK_head=0 if UK_head==.
replace GER_head=0 if GER_head==.

*  CREATE BG CHARACTERISTICS
** EMPLOYMENT
bysort sirtg: egen totalbg_employment=total(effec)
bysort sirtg: egen totalbg2_employment=total(effec) if marketmain_2SIC!="65" & marketmain_2SIC!="66" & marketmain_2SIC!="67"
** SALES
bysort sirtg: egen totalbg_sales=total(q)
bysort sirtg: egen totalbg2_sales=total(q) if marketmain_2SIC!="65" & marketmain_2SIC!="66" & marketmain_2SIC!="67"
** PROFITS
bysort sirtg: egen totalbg_profits=total(ebe)
bysort sirtg: egen totalbg2_profits=total(ebe) if marketmain_2SIC!="65" & marketmain_2SIC!="66" & marketmain_2SIC!="67"
** LABOR COSTS
bysort sirtg: egen totalbg_laborcosts=total(r216) 
bysort sirtg: egen totalbg2_laborcosts=total(r216) if marketmain_2SIC!="65" & marketmain_2SIC!="66" & marketmain_2SIC!="67"
** VA 
bysort sirtg: egen totalbg_va=total(va) 
bysort sirtg: egen totalbg2_va=total(va) if marketmain_2SIC!="65" & marketmain_2SIC!="66" & marketmain_2SIC!="67"
** ASSETS
bysort sirtg: egen totalbg_assets=total(immob) 
bysort sirtg: egen totalbg2_assets=total(immob) if marketmain_2SIC!="65" & marketmain_2SIC!="66" & marketmain_2SIC!="67"

*  SHARE OF FINANCIAL ACTIVITY
** EMPLOYMENT
bysort sirtg: egen totalbg_employmentfin=total(effec) if marketmain_2SIC=="65" | marketmain_2SIC=="66" | marketmain_2SIC=="67"
replace totalbg_employmentfin=0 if totalbg_employmentfin==.
bysort sirtg: gen empshare_temp_fin=totalbg_employmentfin/totalbg_employment
bysort sirtg: egen empshare_fin=max(empshare_temp_fin)
**SALES
bysort sirtg: egen totalbg_salesfin=total(q) if marketmain_2SIC=="65" | marketmain_2SIC=="66" | marketmain_2SIC=="67"
replace totalbg_salesfin=0 if totalbg_salesfin==.
bysort sirtg: gen salesshare_temp_fin=totalbg_salesfin/totalbg_sales
bysort sirtg: egen salesshare_fin=max(salesshare_temp_fin)

*  WEIGHTS OF EACH UNIT IN BG
gen shareinbg_sales=q/totalbg2_sales if market_2SIC!="65" & market_2SIC!="66" & market_2SIC!="67"
gen shareinbg_emp=effec/totalbg2_employment if market_2SIC!="65" & market_2SIC!="66" & market_2SIC!="67"

*  COMPUTE WEIGHTED ROS 
gen ros=ebe/q if out_mean_prof!=1
bysort sirtg: egen bgros_wemp=wtmean(ros), weight(shareinbg_emp)

*  INDUSTRY DUMMIES
tabulate  marketmain_1SIC, gen(pippo_1SIC)
forvalues x = 1(1)10{
bysort sirtg: egen bg_ind1SIC`x'=max(pippo_1SIC`x')
}
drop pippo*

*REALITY CHECK
list naf600 effec totalbg_employment if sirtg=="015450752"
list naf600 q totalbg_sales if sirtg=="015450752"
count
sort siren sirtg

*  LABOR EXPENSE
gen temp_totlaborcosts=r216+r217
gen temp_totlaborcosts_q=temp_totlaborcosts/q
bysort sirtg: egen totalbg_laborcosts_plus=total(temp_totlaborcosts)
capture drop temp_totlaborcosts temp_totlaborcosts_q
** ADD AGE 
merge 1:1 siren using temp_btb_age.dta
drop if _merge==2
drop _merge
** GENERATE CREATION YEAR
gen entry_year=substr(datcr,1,4)
destring entry_year, force replace
replace entry_year=. if entry_year==0
** CONVERT YEAR
gen calendar_year=.
replace calendar_year=2000 if `i'==100
replace calendar_year=2001 if `i'==101
replace calendar_year=2002 if `i'==102
replace calendar_year=2003 if `i'==103
replace calendar_year=2004 if `i'==104
replace calendar_year=2005 if `i'==105
replace calendar_year=2006 if `i'==106
** GENERATE AGE
gen age=calendar_year-entry_year
** GENERATE AVERAGE AGE
bysort sirtg: egen bg_avage_wemp=wtmean(age), weight(shareinbg_emp)

*  SAVE
save temp_restud_bg4_`i', replace
}
*END LOOP


								***************************************************
								*****  DATA ORGANIZATION: PATENT INFORMATION  *****
								***************************************************

* INITIALIZE EPO AND MERGE
use oeb1993_2003, replace
drop at be ch de dk es fr gb gr ie it li lu mc nl pt se fi cy
drop if siren==""
sort siren

*  APLICATION YEAR
gen depot_year=93 if depot==1993
replace depot_year=94 if depot==1994
replace depot_year=95 if depot==1995
replace depot_year=96 if depot==1996
replace depot_year=97 if depot==1997
replace depot_year=98 if depot==1998
replace depot_year=99 if depot==1999
replace depot_year=100 if depot==2000
replace depot_year=101 if depot==2001
replace depot_year=102 if depot==2002
replace depot_year=103 if depot==2003
sort siren

* GENERATE INDICATORS
gen applied=1
gen d_received=1 if delivrance!=. 
replace d_received=0 if d_received!=1

* CONSTRUCT TOTAL INDUSTRY INDICATORS
bysort siren: egen totsiren_app9303=total(applied)
bysort siren: egen totsiren_gra9303=total(d_received)

* GENERATE PRE-SAMPLE INDICATORS
drop if delivrance>2001
bysort siren: egen totsiren_app9300=total(applied)
bysort siren: egen totsiren_gra9300=total(d_received)

* ORGANIZE & CLEAN
keep siren totsiren* depot_year
duplicates drop siren, force
save restud_panelstock, replace

* ADD TO BG DATA
forvalues i=100/106 {
	
	* INITIALIZE DATA
	use temp_restud_bg4_`i', replace
	capture drop totbg_app9303 totbg_gra9303 totbg_app9300 totbg_gra9300 nototbg_app9303 nototbg_gra9303 nototbg_app9300 nototbg_gra9300

	* MERGE WITH PRESAMPLE STOCK
	merge 1:1 siren using restud_panelstock
	drop if _merge==2

	* GENERATE TOTAL BG PATENTS
	bysort sirtg: egen totbg_app9303=total(totsiren_app9303)
	bysort sirtg: egen totbg_gra9303=total(totsiren_gra9303)
	bysort sirtg: egen totbg_app9300=total(totsiren_app9300)
	bysort sirtg: egen totbg_gra9300=total(totsiren_gra9300)
	
	* GENERATE DUMMIES FOR 0
	gen nototbg_app9303=1 if totbg_app9303==0
	gen nototbg_gra9303=1 if totbg_gra9303==0
	gen nototbg_app9300=1 if totbg_app9300==0
	gen nototbg_gra9300=1 if totbg_gra9300==0
	
	*TOTAL NUMBER OF PATENTS
	bysort sirtg depot_year: gen n_patents=_n
	bysort sirtg depot_year: egen total_patents=max(n_patents)

	* ORGANIZE & SAVE 
	drop _merge
	sort siren sirtg
	save temp_restud_bg5_`i', replace

}
*


								*****************************************************************
								*****  DATA ORGANIZATION: HISTORICAL OWNERSHIP INFORMATION  *****
								*****************************************************************

*  INITIALIZE HEAD OF GROUP DATA
** NOTE: PUBLIC CLASSIFICATION CHANGES ACROSS YEARS
** 1985
use tg85, replace
gen year=85
gen public_head_8587=1 if natact=="1"
** 1986
append using tg86
replace year=86 if year==.
replace public_head_8587=1 if natact=="2" & year==86
** 1987
append using tg87
replace year=87 if year==.
replace public_head_8587=1 if natact=="1" & year==87
** PUBLIC DUMMY
replace public_head_8587=0 if public_head_8587==.
** ORGANIZE DATA
keep siren nomgroup public_head_8587
duplicates drop siren, force
** VARIABLES
gen sirtg=siren
** SAVE
sort siren
save restud_iv_headsirenlist_8587, replace

*  MERGE IT 
forvalues i=100/106 {

*  INITIALIZE DATA
use temp_restud_bg5_`i', replace

*  PREPARE IDENTIFIER
rename siren siren_str
gen siren=siren_str
destring siren, replace force

*  MERGE DATASETS
capture drop _merge
sort siren
merge siren using restud_iv_headsirenlist_8587

*  PUT EVERYTHING BACK
drop siren
rename siren_str siren
drop if _merge==2

*  ORGANIZE IVs
** TEMPORARY VARIABLES
bysort sirtg: egen pippo_present_8587=max(_merge)
bysort sirtg: egen pippo_public_8587=max(public_head_8587)
** DEFINITE VARIABLES
***PRESENCE IN 8587
gen d2_present_8587=1 if pippo_present_8587==3
replace d2_present_8587=0 if pippo_present_8587==1
***PUBLIC IN 8587
gen d2_public_8587=1 if pippo_public_8587==1
replace d2_public_8587=0 if pippo_public_8587==0
***DROP VARIABLES
drop pippo_present_8587 pippo_public_8587

*  SAVE
capture drop _merge
save temp_restud_bg5_`i', replace
}
*



								************************************************
								*****  DATA ORGANIZATION: R&D INFORMATION  *****
								************************************************


*  USE AND INITIALIZE LOOP
use lifi_firms_100, replace
forvalues i=100/106 {

*  INITIALIZE & SORT DATA
** SORT FIRM DATA
use r&d_firm_`i', replace
sort siren
d
list in 1/1
save r&d_firm_`i', replace
** START FROM BRANCH
use r&d_branch_`i', replace
sort siren
d
list in 1/1
** NOTE: NEED TO START FROM BRANCH BECAUSE MOST R&D VARIABLES ARE ONLY CONTAINED IN THE BRANCH PART

*  IDENTIFY DUPLICATES
duplicates tag siren, generate (dup)
** NOTE: TWO DIFFERENT LOOPS BECAUSE OF CHANGE IN VARIABLE NAMES
if `i'<105{
		*AGGREGATE EVERYTHING INTO FIRM LEVEL BEFORE 2005
		**SURVEY
		tab qs_nr
		gen survey=1 if qs_nr==""
		replace survey=0 if qs_nr!=""
		**COMPOSITION
		bysort siren: egen basic_rd=total(rech_fond) 
		bysort siren: egen applied_rd=total(rech_app) 
		bysort siren: egen dev_rd=total(devel)
		**OUTSOURCING
		bysort siren: egen f_ent_gr=total(ent_gr)
		bysort siren: egen f_etr_gr=total(etr_gr) 
		bysort siren: egen f_ens_sup=total(ens_sup) 
		bysort siren: egen f_org_pub=total(org_pub)
		**INCOMING
		bysort siren: egen f_r_ent_gr=total(r_ent_gr)
		bysort siren: egen f_r_etr_gr=total(r_etr_gr)
		}
else if `i'>104{
		*AGGREGATE EVERYTHING INTO FIRM LEVEL AFTER 2005
		**SURVEY
		gen survey=1 if nr=="1"
		replace survey=0 if nr!="1"
		**COMPOSITION
		bysort siren: egen basic_rd=total(di_rech_fond)
		bysort siren: egen applied_rd=total(di_rech_app)
		bysort siren: egen dev_rd=total(di_devel)
		**OUTSOURCING
		bysort siren: egen f_ent_gr=total(de_ent_gr)
		bysort siren: egen f_etr_gr=total(de_etr_gr)
		bysort siren: egen f_ens_sup=total(de_ens_sup)
		bysort siren: egen f_org_pub=total(de_org_pub)
		**INCOMING
		bysort siren: egen f_r_ent_gr=total(r_ent_gr)
		bysort siren: egen f_r_etr_gr=total(r_etr_gr)
		**TOTALS
		bysort siren: egen dird=total(di_dird)
		bysort siren: egen derd=total(de_derd)
		bysort siren: egen budgetot=total(bra_budgetot)
		}

*DROP DUPLICATES
drop if dup!=0

*ORGANIZATION
keep entreprise siren dird derd budgetot basic_rd applied_rd dev_rd f_ent_gr f_etr_gr f_ens_sup f_org_pub f_r_ent_gr f_r_etr_gr survey pond 
rename dird internal_rd
rename derd external_rd
rename budgetot total_rd
rename pond weighting
save restud_temp_rd1`i', replace

*  MERGE FIRM AND BRANCH DATA
** START FROM BRANCH AND USE FIRM DATA
sort siren
merge siren using r&d_firm_`i'
drop if _merge!=3
** CONSISTENCY CHECK
gen ok_internal=1 if dird==internal_rd
tab ok_internal
gen ok_external=1 if derd==external_rd
tab ok_external
gen ok_total=1 if budgetot==total_rd
tab ok_total
list budgetot total_rd if ok_total!=1
gen ok_weights=1 if pond==weighting
tab ok_weights

*  UNWEIGHT QS OBSERVATIONS
if `i'<104{ 
local variables "total_rd internal_rd external_rd basic_rd applied_rd dev_rd f_ent_gr f_etr_gr f_ens_sup f_org_pub f_r_ent_gr f_r_etr_gr financ_pub"
foreach x of local variables{
		replace `x'=`x'/weighting
}
}
*

*  FRANCS INTO EUROS
local variables "total_rd internal_rd external_rd basic_rd applied_rd dev_rd f_ent_gr f_etr_gr f_ens_sup f_org_pub f_r_ent_gr f_r_etr_gr financ_pub"
foreach x of local variables{
		replace `x'=`x'/6.55957 if `i'==100
}
*

*PREPARE FOR MERGE WITH BG DATA
keep siren entreprise total_rd internal_rd external_rd basic_rd applied_rd dev_rd f_ent_gr f_etr_gr f_ens_sup f_org_pub f_r_ent_gr f_r_etr_gr survey weighting financ_pub pinformat pbiotech penvir pmater pshs
rename f_ent_gr outsource_allnaf_frbg
rename f_etr_gr outsource_allnaf_rowbg
rename f_ens_sup outsource_fruniv
rename f_org_pub outsource_rowpub
rename f_r_ent_gr incoming_targeted_frbg
rename f_r_etr_gr incoming_targeted_rowbg 
rename financ_pub public_financing
rename pinformat rdsh_software
rename pbiotech rdsh_biotech
rename penvir rdsh_environment
rename pmater rdsh_materials
rename pshs rdsh_socsciences

*SAVE
sort siren
save restud_temp_rd2`i', replace
}
*END LOOP



								***************************************************************************************
								*****  DATA ORGANIZATION: MERGING AND ORGANIZING R&D AND MULTI-INDUSTRY PRESENCE  *****
								***************************************************************************************

								
***********************************
*****  FINAL PUTTING TOGETHER *****
***********************************

*  QUICK SORT
forvalues i=100/106 {
use temp_restud_bg5_`i', replace
sort siren
save temp_restud_bg5_`i', replace
}
*

*  USE AND INITIALIZE LOOP
use lifi_firms_100, replace
forvalues i=100/106 {
*  MERGE
use restud_temp_rd2`i', replace
sort siren
merge siren using temp_restud_bg5_`i'
rename _merge final_merge

*  DROP
drop if final_merge!=3

*  YEAR IDENTIFIER
gen year=`i'

*  ULTIMATE SAVE
save restud_nelson_`i', replace
}
*  END LOOP


*  INTERMEDIATE STEP FOR FORMATTING YEAR 2000 R&D DATA 
use restud_nelson_100, replace
foreach var of varlist rdsh_software rdsh_biotech rdsh_environment rdsh_materials rdsh_socsciences{
		*  NEW VARIABLES
		gen `var'_2digittemp1=substr(`var',1,2)
		gen `var'_1digittemp1=substr(`var',1,1)
		*  IDENTIFIER
		gen tempvar=substr(`var',2,1)
		gen dummycomma=1 if tempvar==","
		*  NEW VARIABLES & IDENTIFIER
		gen `var'_2digittemp2=`var'_2digittemp1 if dummycomma!=1
		replace `var'_2digittemp2=`var'_1digittemp1 if dummycomma==1
		*  DESTRING
		destring `var'_2digittemp2, replace
		replace `var'_2digittemp2=100 if `var'=="100,00"
		drop `var' dummycomma tempvar `var'_1digittemp1 `var'_2digittemp1
		rename `var'_2digittemp2 `var'
		replace `var'=0 if `var'==.
}
save restud_nelson_100, replace
*

*  ORGANIZE DATA
use restud_nelson_100, replace
append using restud_nelson_101
append using restud_nelson_102
append using restud_nelson_103
append using restud_nelson_104
append using restud_nelson_105
append using restud_nelson_106


************************************
*****  VARIABLES ORGANIZATION  *****
************************************

*  OUTCOME VARIABLES
gen mix_applied_rd=basic_rd/applied_rd

*  COVARIATES
** LOG SIZE, EMPLOPYMENT OF ALL THE GROUP
gen ltotbgemployment=log(totalbg_employment)
gen ltotalbg2_employment=log(totalbg2_employment)
** ORGANIZATIONAL FIXED EFFECT FOR STAND ALONE
gen sa1=0
replace sa1=1 if nber1_4SIC==1
gen sa2=0
replace sa2=1 if nber2_4SIC==1
gen sa3=0
replace sa3=1 if nber3_4SIC==1
gen sa4=0
replace sa4=1 if nber4_4SIC==1

*  MARKET FIXED EFFECTS
encode naf600, generate(sect)
** 1SIC Classification
gen sic1=substr(naf600,1,1) 
destring sic1, gen(code) force
drop sic1
rename code sic1
** 2SIC Classification
gen sic2=substr(naf600,1,2)
destring sic2, gen(code) force
drop sic2
rename code sic2
** 3SIC Classification
gen sic3=substr(naf600,1,3) 
destring sic3, gen(code) force
list sic3 code in 1/10
drop sic3
rename code sic3

*  LOG MULTI MARKET MEASURES
forvalues s=1/4 {
	gen lnber1_`s'SIC=log(nber1_`s'SIC)
	gen lnber2_`s'SIC=log(nber2_`s'SIC)
	gen lnber3_`s'SIC=log(nber3_`s'SIC)
	gen lnber4_`s'SIC=log(nber4_`s'SIC)
}
*

*  DROP NEGATIVE VALUES
drop if mix_applied_rd<0

*  GENERATE OUTCOME VARIABLE AT THE BG LEVEL
bysort sirtg year: egen totbg_basic_rd=total(basic_rd)
bysort sirtg year: egen totbg_applied_rd=total(applied_rd)
bysort sirtg year: egen totbg_total_rd=total(total_rd)
list  sirtg siren year basic_rd applied_rd totbg_basic_rd totbg_applied_rd in 1/10
gen mixbg_applied_rd=totbg_basic_rd/totbg_applied_rd
gen mixbg_total_rd=totbg_basic_rd/totbg_total_rd
list  sirtg siren year totbg_basic_rd totbg_applied_rd mixbg_applied_rd in 1/10

*  TAKE ONLY ONE OBSERVATION
bysort year: distinct sirtg
gen uniq_bg=1
bysort sirtg year: replace uniq_bg=. if sirtg[_n]==sirtg[_n-1]
bysort year: count if uniq_bg==1

*  BG OUTLIERS
gen outbg_mixapplied=0 if uniq_bg==1
_pctile mixbg_applied_rd if mixbg_applied_rd!=0 & mixbg_applied_rd!=. & uniq_bg==1, percentiles(1 5 25 50 75 95 99)
replace outbg_mixapplied=1 if abs(mixbg_applied_rd-r(r4) )> 5*(r(r5)-r(r3)) & uniq_bg==1

*  NEW IDENTIFIER
sort sirtg year
gen count=_n
bysort sirtg: egen new_identifier=max(_n)

*  DUMMY FINANCIAL INTERMEDIARIES
gen d_fin=0
replace d_fin=1 if empshare_fin>0 & empshare_fin!=1

*  TREATMENT DUMMY
gen dbasic=1 if mixbg_applied_rd>0 & mixbg_applied_rd!=. & outbg_mixapplied==0 & uniq_bg==1
replace dbasic=0 if mixbg_applied_rd==0 & outbg_mixapplied==0 & uniq_bg==1

*  UNIVERSITY OUTSOURCING
bysort sirtg year: egen outsourcebg_univ=total(outsource_fruniv)
gen dout_univ=1 if outsourcebg_univ!=0 & mixbg_applied_rd!=. & outbg_mixapplied==0 & uniq_bg==1
replace dout_univ=0 if outsourcebg_univ==0 & mixbg_applied_rd!=. & outbg_mixapplied==0 & uniq_bg==1

*  PUBLIC FINANCING
bysort sirtg year: egen totbg_public_rd=total(public_financing)
gen ltotbg_public_rd=log(1+totbg_public_rd)
gen rtotbg_public_rd=totbg_public_rd/totbg_total_rd
gen dtotbg_public_rd=1 if totbg_public_rd>0 & totbg_public_rd!=. 
replace dtotbg_public_rd=0 if dtotbg_public_rd!=1

*  R&D AREA
** WEIGHTS
gen shareinbg_rdtot=total_rd/totbg_total_rd
** WEIGHTED SHARES
foreach var of varlist rdsh_software rdsh_biotech rdsh_environment rdsh_materials rdsh_socsciences{
		*NEW VARIABLES
		bysort sirtg year: egen bg`var'=wtmean(`var'), weight(shareinbg_rdtot)
}
*

*SAVE DATA
save restud_nelson, replace




*****************************************
*****  GROWTH AND EXIT INFORMATION  *****
*****************************************

*DECLARE LOCAL														
local i=100
*DECLARE LOOP (HANDLE BLOCKS OF 3 YEARS)
while `i' <= 104 {														
		*DECLARE LOCALS
		local j=`i'+1
		local k=`i'+2
		di `i'
		di `j'
		di `k'
		
		*INITIALIZE DATA
		use temp_restud_bg4_`i', replace
		gen year=`i'
		append using temp_restud_bg4_`j'
		replace year=`j' if year==. 
		append using temp_restud_bg4_`k'
		replace year=`k' if year==.

		*DESTRING
		destring siren, force gen(code)
		drop if code==.

		*DECLARE PANEL
		tsset code year		
		keep code year effec siren sirtg shareinbg_sales naf600


*************************************************
**GROWTH RATES IN SALES, EMPLOYMENT AND PROFITS**
*************************************************

*ARRANGE DATA
gen emp=effec

*TAKE ONLY ONE OBSERVATION
bysort year: distinct sirtg
gen uniq_bg=1
bysort sirtg year: replace uniq_bg=. if sirtg[_n]==sirtg[_n-1]
bysort year: count if uniq_bg==1

* EMPLOYMENT GROWTH OF BG
**COMPUTE GROWTH RATES IN EMPLOYMENT 
sort code year
gen growth_firmemp=(emp-l.emp)/l.emp
sort code year
bysort code: gen fwd_growth_firmemp=growth_firmemp[_n+1]
**  BEAM IT UP TO THE GROUP LEVEL
bysort sirtg year: egen fwd_bggrowth_wemp=wtmean(fwd_growth_firmemp), weight(shareinbg_sales)


*********************************
**ENTRY AND EXIT RATES OF FIRMS**
*********************************

*  EXIT RATES BY FIRMS
** GENERATE RUNS
gen run=.
sort code year
bysort code: replace run=cond(l.run==., 1, l.run + 1)
by code: egen maxrun=max(run)
** DROP THOSE WITH HOLES
gen obs=1
bysort code: gen totobs=_N
gen tag=1 if totobs==2 & maxrun==1
** DEFINE LAST PERIOD
bysort code: egen lastyear=max(year)
** DEFINE EXIT GENERAL
gen gen_exit=.		
replace gen_exit=1 if totobs==2 & lastyear!=`k' & tag!=1
replace gen_exit=1 if totobs==1 & lastyear!=`k' & tag!=1
**DEFINE EXIT TIMING
gen exit=0 if tag!=1
replace exit=1 if year==lastyear & gen_exit==1 & tag!=1

*  EXIT RATES BY BG
** COUNT OPEN LINES IN BG
bysort sirtg year naf600: gen open_lines=_N
** COUNT CLOSING LINES IN BG
bysort sirtg year naf600: egen closed_lines=total(exit) 
gen bg_openclosed=closed_lines/open_lines
gen bg_exit=1 if bg_openclosed==1
replace bg_exit=0 if bg_openclosed!=1 & bg_openclosed!=.

*SAVE
**GENERAL SAVE
save moments_economypanel`i'`k'_20200630, replace
**SAVE R&D SAMPLE
sort siren year
merge siren year using sirenmomentlist_20101005
drop if _merge!=3
save moments_panel`i'`k'_20200630, replace
**REGULAR SAVE
drop if year==`k'
sort siren year
save moments1_panel`i'`k'_20200630, replace

*RELAUNCH LOOP
local i = `i' + 2
}
*





log close






