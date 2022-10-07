
*  INITIALIZE DATA
use allcitingcitations_onerow_updatedappyear_restud.dta, replace
*  NOTE: COMPUTES FOR ALL PATENTS THE NUMBER OF CITATIONS WITHIN X YEARS

*  GENERATE SPLIT
gen public=.
replace public=0 if asscode==2 | asscode==3
replace public=1 if asscode==6 | asscode==7
drop if public==.

*  VARIABLE
gen totcites_15y=TotN_citations_15yrs if year<1991

*  KEEP
keep if country=="FR"

*  DESCIPTIVE STATISTICS		
latabstat TotN_citations_15yrs if public==0 & year<1991, stat(mean p25 p50 p75 sd min max count) col(stat) f(%10.2f) clabel(summary-statistics) 
latabstat TotN_citations_15yrs if public==1 & year<1991, stat(mean p25 p50 p75 sd min max count) col(stat) f(%10.2f) clabel(summary-statistics) 

*  RMS
gen x2=(TotN_citations_15yrs)^2
bysort public: egen mx2=mean(x2)
gen rms=sqrt(mx2)
bysort public: sum rms


