*  INITIALIZE DATA
use restud_nelson.dta, replace	

*SET SEED
set seed 14
version 12

*******************************
*****  DATA ORGANIZATION  *****
*******************************

*  CONSTRUCT TOTAL R&D TO SALES AT BG LEVEL AND CORRECT OUTLIERS
gen totrdsales_bg=totbg_total_rd/totalbg_sales 
gen outbg_totrdsales=0 if uniq_bg==1
_pctile totrdsales_bg if uniq_bg==1, percentiles(1 5 25 50 75 95 99)
replace outbg_totrdsales=1 if abs(totrdsales_bg-r(r4) )> 5*(r(r5)-r(r3)) & uniq_bg==1

*  CENSORED INDUSTRY DISTRIBUTION
gen nber1ALT_1SIC=nber1_1SIC
replace nber1ALT_1SIC=8 if nber1_1SIC>8

*  STATE SHAREHOLDER
** STATE PARTICIPATION
gen iv2=1 if d2_public_8587==1
replace iv2=0 if iv2==.
** STILL STATE OWNED TODAY
replace state_owned=0 if state_owned==.

*  PATENTS
replace nototbg_gra9300=0 if totbg_gra9300!=0
gen d_nopat=0
replace d_nopat=1 if  total_patents==. 
replace total_patents=0 if total_patents==.
gen lpat=log(1+total_patents)
bysort sirtg: egen prepatents=mean(total_patents)
gen lprepatents=log(1+prepatents)

*  REPLACE
foreach var of varlist bgrdsh_software bgrdsh_biotech bgrdsh_environment bgrdsh_materials bgrdsh_socsciences {
		replace `var'=`var'/100
}
*
keep 	mixbg_applied_rd lnber1_1SIC ltotbgemployment sa1 lprepatents d_nopat bg_ind1SIC* ///
	bgrdsh_software bgrdsh_biotech bgrdsh_environment bgrdsh_materials bgrdsh_socsciences financial ///
	year outbg_mixapplied uniq_bg outbg_totrdsales new_identifier iv2 state_owned dtotbg_public_rd ///
	dout_univ nber1ALT_1SIC totbg_total_rd totalbg_sales bgros_wemp bg_avage_wemp totalbg_employment dbasic siren totalbg_laborcosts_plus

*  BASELINE REGRESSION SAMPLE	
xi: tobit mixbg_applied_rd lnber1_1SIC ltotbgemployment sa1 i.year if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0,  ll(0) vce(cluster new_identifier)
gen sample=e(sample)
save tempfirm, replace
keep if sample==1
	
*  VARIABLE LABELS
label var mixbg_applied_rd "Basic R\&D Intensity"
label var lnber1_1SIC "Log \# of 1 Digit Product Markets"
label var ltotbgemployment "BG Employment"
label var bgrdsh_software "\% in Software"
label var bgrdsh_biotech "\% in Biotech"
label var bgrdsh_environment "\% in Environment"
label var bgrdsh_materials "\% in Materials"
label var bgrdsh_socsciences "\% in Social Science"
label var dout_univ "Outsourcing to Univ."
label var dtotbg_public_rd "Public R\&D Funding"
label var iv2 "State Participation in 1985"
label var state_owned "SOE Today"
label var nber1ALT_1SIC "\# of Industries"


*********************
*****  TABLE 2  *****
*********************

*  ESTIMATES
** COLUMN 1
xi: tobit mixbg_applied_rd lnber1_1SIC ltotbgemployment sa1 i.year if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0,  ll(0) vce(cluster new_identifier)
margins, dydx(*) predict(ystar(0,.)) at(sa1=(0) _Iyear_101=(0) _Iyear_102=(1) _Iyear_103=(0) _Iyear_104=(0) ///
		_Iyear_105=(0) _Iyear_106=(0)) post vsquish
eststo reg1
** COLUMN 2
xi: tobit mixbg_applied_rd lnber1_1SIC ltotbgemployment sa1 financial i.year bg_ind1SIC* if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0,  ll(0) vce(cluster new_identifier)
margins, dydx(*) predict(ystar(0,.)) at(sa1=(0) _Iyear_101=(0) _Iyear_102=(1) _Iyear_103=(0) _Iyear_104=(0) ///
		_Iyear_105=(0) _Iyear_106=(0)) post vsquish
eststo reg2
** COLUMN 3
xi: tobit mixbg_applied_rd lnber1_1SIC ltotbgemployment sa1 financial i.year bg_ind1SIC* bgrdsh_software bgrdsh_biotech bgrdsh_environment bgrdsh_materials bgrdsh_socsciences if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0,  ll(0) vce(cluster new_identifier)
margins, dydx(*) predict(ystar(0,.)) at(sa1=(0) _Iyear_101=(0) _Iyear_102=(1) _Iyear_103=(0) _Iyear_104=(0) ///
		_Iyear_105=(0) _Iyear_106=(0)) post vsquish
eststo reg3
** COLUMN 4
xi: tobit mixbg_applied_rd lnber1_1SIC ltotbgemployment sa1 lprepatents d_nopat bg_ind1SIC* bgrdsh_software bgrdsh_biotech bgrdsh_environment bgrdsh_materials bgrdsh_socsciences i.year if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0,  ll(0) vce(cluster new_identifier)
margins, dydx(*) predict(ystar(0,.)) at(sa1=(0) _Iyear_101=(0) _Iyear_102=(1) _Iyear_103=(0) _Iyear_104=(0) ///
		_Iyear_105=(0) _Iyear_106=(0)) post vsquish
eststo reg4		
*  TABLE 2: COLUMNS (1) TO (4)
estout reg1 reg2 reg3 reg4 using RESTUD_BTB_T2.tex, replace title(LOOK FILE NAME) ///
		substitute(_ \_) cells(b(star label(Coef.) fmt(%9.3f)) se(par fmt(%9.3f))) starlevels(* 0.10 ** 0.05 *** 0.01) ///
		stats(r2 N, fmt(%9.3f %9.0g) labels(R-squared)) ///
		collabels(, none) nolegend label varlabels(_cons Constant) style(tex) ///
		prehead(\begin{tabular}{l*{@M}{c}}) postfoot(\end{tabular}) 

		
*********************
*****  TABLE 3  *****
*********************

*  ESTIMATES
** COLUMN 1 & 4
version 12
** NOTE LATER VERSIONS OF STATA dydx(_all)
xi: ivtobit mixbg_applied_rd ltotbgemployment sa1 financial i.year ///
			bgrdsh_software bgrdsh_biotech bgrdsh_environment bgrdsh_materials bgrdsh_socsciences ///
			state_owned (lnber1_1SIC = iv2) if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0, ll(0) vce(cluster new_identifier) diff first
			margins, dydx(*) predict(e(0,.)) at(sa1=(0) lnber1_1SIC=(0) _Iyear_101=(0) ///
			_Iyear_102=(1) _Iyear_103=(0) _Iyear_104=(0) _Iyear_105=(0) _Iyear_106=(0) iv2=(1) ///
			state_owned==(0)) post vsquish
eststo reg1
** COLUMN 2 & 5
xi: ivtobit mixbg_applied_rd ltotbgemployment sa1 financial dtotbg_public_rd i.year ///
			bgrdsh_software bgrdsh_biotech bgrdsh_environment bgrdsh_materials bgrdsh_socsciences ///
			state_owned (lnber1_1SIC = iv2) if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0, ll(0) vce(cluster new_identifier) diff first
			margins, dydx(*) predict(e(0,.)) at(sa1=(0) lnber1_1SIC=(0) _Iyear_101=(0) ///
			_Iyear_102=(1) _Iyear_103=(0) _Iyear_104=(0) _Iyear_105=(0) _Iyear_106=(0) iv2=(1) ///
			state_owned==(0)) post vsquish
eststo reg2
** COLUMN 3 & 6
xi: ivtobit mixbg_applied_rd ltotbgemployment sa1 financial dout_univ i.year ///
			bgrdsh_software bgrdsh_biotech bgrdsh_environment bgrdsh_materials bgrdsh_socsciences ///
			state_owned (lnber1_1SIC = iv2) if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0, ll(0) vce(cluster new_identifier) diff first
			margins, dydx(*) predict(e(0,.)) at(sa1=(0) lnber1_1SIC=(0) _Iyear_101=(0) ///
			_Iyear_102=(1) _Iyear_103=(0) _Iyear_104=(0) _Iyear_105=(0) _Iyear_106=(0) iv2=(1) ///
			state_owned==(0)) post vsquish
eststo reg3
*  TABLE 3: COLUMNS (1) TO (6)
estout reg1 reg2 reg3 using RESTUD_BTB_T3.tex, replace title(LOOK FILE NAME) ///
		substitute(_ \_) cells(b(star label(Coef.) fmt(%9.3f)) se(par fmt(%9.3f))) starlevels(* 0.10 ** 0.05 *** 0.01) ///
		stats(r2 N, fmt(%9.3f %9.0g) labels(R-squared)) ///
		collabels(, none) nolegend label varlabels(_cons Constant) style(tex) ///
		prehead(\begin{tabular}{l*{@M}{c}}) postfoot(\end{tabular}) 




******************************************************
*****  FIGURE 4: DISTRIBUTION ACROSS INDUSTRIES  *****
******************************************************

*  FIGURE
twoway  (histogram nber1ALT_1SIC if uniq_bg==1 & sample==1, start(1) discrete percent color(gs13)) ///
		(histogram nber1ALT_1SIC if uniq_bg==1 & sample==1, start(1) discrete percent ///
		fcolor(none) lcolor(black)), ///
			   graphregion(fcolor(white) lcolor(white)) ///
			   xlabel(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6" 7 "7" 8 "8+", labels labcolor(black)) ///
			   legend(off) xtitle(NUMBER OF 1 DIGIT SIC INDUSTRIES) 
graph export histPROD1whitefinal.pdf, replace


*****************************************************************
*****  FIGURE 5: BASIC RESEARCH INTENSITY ACROSS INDUSTRIES  ****
*****************************************************************

*  PARTIAL OUT SIZE & ORGANIZATIONAL FORM
xi: reg mixbg_applied_rd ltotbgemployment sa1 if outbg_mixapplied==0 & uniq_bg==1 & outbg_totrdsales==0, vce(cluster new_identifier)
predict xbnew if e(sample)
gen yhat=mixbg_applied_rd-xbnew+.0480563
				
*  AVERAGE BY INDUSTRY
bysort nber1ALT_1SIC: egen mixind_applied_rd=mean(yhat) if sample==1

*  NEW BEST GRAPH
label var mixind_applied_rd 	"Average"
twoway (scatter mixind_applied_rd nber1ALT_1SIC if sample==1, xlabel(1[1]8) mcolor(blue) ///
		msize(large) msymbol(circle) yaxis(1))(lfit yhat nber1ALT_1SIC if sample==1 & outbg_totrdsales==0,  xlabel(1[1]8) ///
		lcolor(red) lwidth(medthick)), ytitle(BASIC RESEARCH INTENSITY) xtitle(NUMBER OF 1 DIGIT SIC INDUSTRIES)  ///
		xlabel(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6" 7 "7" 8 "8+", labels labcolor(black)) ///
		legend(on) graphregion(fcolor(white) lcolor(white))
gr export partialledoutwhitefinal.pdf, replace


***************************
*****  MICRO MOMENTS  *****
***************************

*  INTENSIVE MARGIN FOR EACH # OF INDUSTRIES
bysort nber1ALT_1SIC: sum mixbg_applied_rd if sample==1 & outbg_totrdsales==0

*  EXTENSIVE MARGIN FOR EACH # OF INDUSTRIES
bysort nber1ALT_1SIC: sum dbasic if sample==1 & outbg_totrdsales==0

*  RETURN ON SALES
sum bgros_wemp if sample==1 & outbg_totrdsales==0

*  INDUSTRY DISTRIBUTION
sum nber1ALT_1SIC if sample==1 & outbg_totrdsales==0
gen SQnber1ALT_1SIC=nber1ALT_1SIC^2
sum SQnber1ALT_1SIC if sample==1 & outbg_totrdsales==0

*  R&D EXPENDITURES OVER TOTAL LABOR PLUS INTERMEDIATE INPUT COSTS
gen rdlaborplus_bg=totbg_total_rd/totalbg_laborcosts_plus
sum rdlaborplus_bg if sample==1 & outbg_totrdsales==0

*  AGE CUTOFF AT MEDIAN EMPLOYMENT
gen amedian_bgemployment=0 if sample==1 & totalbg_employment!=.
_pctile totalbg_employment if sample==1 & totalbg_employment!=., percentiles(1 5 25 50 75 95 99)
replace amedian_bgemployment=1 if totalbg_employment>r(r4) & sample==1 & totalbg_employment!=.
bysort amedian_bgemployment: sum bg_avage_wemp if sample==1 & outbg_totrdsales==0

*  ORGANIZE DATA: EMPLOYMENT GROWTH
use tempfirm, replace
sort siren year
merge siren year using moments1_panel100102_20200630.dta, keep(fwd_* bg_exit)
drop _merge
sort siren year
merge siren year using moments1_panel102104_20200630.dta, keep(fwd_* bg_exit) update
drop _merge
sort siren year
merge siren year using moments1_panel104106_20200630.dta, keep(fwd_* bg_exit) update
*  EMPLOYMENT GROWTH
sum fwd_bggrowth_wemp if sample==1 & outbg_totrdsales==0 & year<106


*  ORGANIZE DATA: EXIT
use moments_economypanel100102_20200630.dta, replace
drop if year==102
append using moments_economypanel102104_20200630.dta
drop if year==104
append using moments_economypanel104106_20200630.dta
drop if year==106
*  FREQUENCIES
tab  bg_exit 
