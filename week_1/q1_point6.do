cd "/Users/veronica/Dropbox/PhD/2023_2/EC_702_Macro_Theory/problem_set/week_1"

global overleaf "/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_Ding_Liu_Perez"

use pwt1001, clear
 
 
set scheme gta
* making growth rate gdp 
sort countrycode year


by countrycode: gen lag_gdp = rgdpna[_n-1]
br countrycode year rgdpna lag_gdp

bysort countrycode: g growth_rate=(rgdpna[_n]-lag_gdp)/lag_gdp

br countrycode year rgdpna lag_gdp growth_rate

corr growth_rate csh_i

gen ln_growth_rate = ln(growth_rate)
gen ln_csh_i = ln(csh_i)

scatter (csh_i growth_rate)
graph export "corr_growth_csh.png", replace


encode countrycode, gen(countrycode_codd)

** baseline, yearly growth rate
reghdfe csh_i growth_rate, absorb(country) cluster(countrycode_codd)

// outreg2 using "$overleaf/tables/PS1_regression_growth.tex", replace tex(fragment) /// 
// ctitle("Dependent Variable: Share gross K formation", "Cross-Section all years") ///
// addtext(Country FE, YES)

** run models by aggregated longer periods
			
** aggregated data, yearly growth rate



local values "10 15 20 30 50"


foreach x of local values {
	preserve
		
		// Sort the dataset by country and year
		sort countrycode year

		// Create a new variable for the decade
		gen time_frame = floor((year - 1950) / `x') * `x' + 1950

		collapse (mean) growth_rate csh_i, by(countrycode_codd time_frame)
		gen growth_rate_2 = growth_rate^2 
		
		reghdfe csh_i growth_rate, absorb(countrycode) cluster(countrycode_codd)
		
// 		outreg2 using "$overleaf/tables/PS1_regression_growth.tex", append tex(fragment) /// 
// 		ctitle("t = `x'") addtext(Country FE, YES)
		
		reghdfe csh_i growth_rate growth_rate_2, absorb(countrycode) cluster(countrycode_codd)
		
		twoway (scatter csh_i growth_rate) (qfit csh_i growth_rate), ///
		xtitle("avg. growth rate") ///
		title("t = `x'") legend(order(1 "avg. investement rate" 2 "" ))
		graph export "corr_growth_csh_`x'.png", replace
		
// 		outreg2 using "$overleaf/tables/PS1_regression_growth.tex", append tex(fragment) /// 
// 		ctitle("t = `x'") addtext(Country FE, YES)
		
	restore 
}



** calculations in the long run

