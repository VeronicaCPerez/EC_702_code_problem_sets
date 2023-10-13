cd "/Users/veronica/Dropbox/PhD/2023_2/EC_702_Macro_Theory/problem_set/week_2"

global overleaf "/Users/veronica/Dropbox/Apps/Overleaf/Macro_EC_702_PS_cd0814_jayinliu_vcperez"



************	question 1.5
****************************************************************
import delimited "ps2_q5_simulated_data_2", clear encoding("utf-8")

*each row is a country, each column is a period, reshape 

gen id_country = _n
gen country_str = "country_"
egen country = concat(country_str id_country)

drop id_country country_str

forvalues i = 1/100 {
	rename v`i' time`i' 
}

order country, first

reshape long time, i(country) j(year)
rename time gdp_pc

* gen ln gdp 
gen ln_gdp_pc = ln(gdp_pc)
label var ln_gdp_pc "$ ln(y_{t}) $"

sort country year

by country: gen ln_gdp_pc_plus1 = ln_gdp_pc[_n+1]

encode country, gen(countrycode_num)

* gen delta ln gdp
sort country year

by country: gen delta_ln_gdp_pc = ln_gdp_pc_plus1-ln_gdp_pc

* regression
********************
* eq 3 no fixed effects
reg delta_ln_gdp_pc ln_gdp_pc year, r
outreg2 using "$overleaf/tables/PS2_q5.tex", replace tex(fragment) /// 
ctitle("Eq. (3)") ///
addtext(Country FE, No, Year FE, No) label

* eq 3 drop year 

reg delta_ln_gdp_pc ln_gdp_pc, r
outreg2 using "$overleaf/tables/PS2_q5.tex", append tex(fragment) /// 
ctitle("Eq. (3) no Year") ///
addtext(Country FE, No, Year FE, No) label

* eq 3 year fe 

reghdfe delta_ln_gdp_pc ln_gdp_pc, absorb(year) 
outreg2 using "$overleaf/tables/PS2_q5.tex", append tex(fragment) /// 
ctitle("Year FE") ///
addtext(Country FE, No, Year FE, Yes) label

* eq 3 country + year fe 

reghdfe delta_ln_gdp_pc ln_gdp_pc, absorb(year country) cluster(country)

outreg2 using "$overleaf/tables/PS2_q5.tex", append tex(fragment) /// 
ctitle("Country Year FE") ///
addtext(Country FE, Yes, Year FE, Yes) label




************	question 1.6
****************************************************************
import delimited "ps2_q6_simulated_data_2", clear encoding("utf-8")

*each row is a country, each column is a period, reshape 

gen id_country = _n
gen country_str = "country_"
egen country = concat(country_str id_country)

drop id_country country_str


forvalues i = 1/100 {
	rename v`i' time`i' 
}

order country, first

reshape long time, i(country) j(year)
rename time gdp_pc

* gen ln gdp 
gen ln_gdp_pc = ln(gdp_pc)
label var ln_gdp_pc "$ ln(y_{t}) $"

sort country year

by country: gen ln_gdp_pc_plus1 = ln_gdp_pc[_n+1]

encode country, gen(countrycode_num)

* gen delta ln gdp
sort country year 
by country: gen delta_ln_gdp_pc = ln_gdp_pc_plus1-ln_gdp_pc

* regression
********************

* eq 3 no fixed effects
reg delta_ln_gdp_pc ln_gdp_pc year
outreg2 using "$overleaf/tables/PS2_q6.tex", replace tex(fragment) /// 
ctitle("Eq. (3)") ///
addtext(Country FE, No, Year FE, No) label

* eq 3 drop year 

reg delta_ln_gdp_pc ln_gdp_pc
outreg2 using "$overleaf/tables/PS2_q6.tex", append tex(fragment) /// 
ctitle("Eq. (3) no Year") ///
addtext(Country FE, No, Year FE, No) label



* eq 3 year fe 

reghdfe delta_ln_gdp_pc ln_gdp_pc, absorb(year)
outreg2 using "$overleaf/tables/PS2_q6.tex", append tex(fragment) /// 
ctitle("Year FE") ///
addtext(Country FE, No, Year FE, Yes) label

* eq 3 country + year fe 

reghdfe delta_ln_gdp_pc ln_gdp_pc, absorb(year country)
outreg2 using "$overleaf/tables/PS2_q6.tex", append tex(fragment) /// 
ctitle("Country Year FE") ///
addtext(Country FE, Yes, Year FE, Yes) label


************	question 1.8 - 1.9
****************************************************************


use pwt1001, clear
 
 
set scheme gta

****************************************************************
************			Preparing data				************
****************************************************************
sort countrycode year

gen gdp_pc = rgdpna/pop

gen ln_gdp_pc = ln(gdp_pc)
label var ln_gdp_pc "$ ln(y_{t}) $"

by countrycode: gen ln_gdp_pc_plus1 = ln_gdp_pc[_n+1]

encode countrycode, gen(countrycode_num)
label var countrycode_num "Number of countries"

xtset countrycode_num year

sort countrycode year

gen delta_ln_gdp_pc = ln_gdp_pc_plus1-ln_gdp_pc

************	question 1.8
****************************************************************


* eq 3 no fixed effects
reg delta_ln_gdp_pc ln_gdp_pc year
outreg2 using "$overleaf/tables/PS2_q8.tex", replace tex(fragment) /// 
ctitle("Eq. (3)") ///
addtext(Country FE, No, Year FE, No) label


* eq 3 country fixed effects
reghdfe delta_ln_gdp_pc ln_gdp_pc year, absorb(countrycode)
outreg2 using "$overleaf/tables/PS2_q8.tex", append tex(fragment) /// 
ctitle("Trend + Country FE") ///
addtext(Country FE, Yes, Year FE, No) label


* eq 3 random effects 
xtreg delta_ln_gdp_pc ln_gdp_pc year
outreg2 using "$overleaf/tables/PS2_q8.tex", append tex(fragment) /// 
ctitle("Random Effects") ///
addtext(Country FE, No, Year FE, No) label


* year fixed effects
reghdfe delta_ln_gdp_pc ln_gdp_pc, absorb(year)
outreg2 using "$overleaf/tables/PS2_q8.tex", append tex(fragment) /// 
ctitle("Year FE") ///
addtext(Country FE, No, Year FE, Yes) label

* eq 3 country & year fixed effects
reghdfe delta_ln_gdp_pc ln_gdp_pc, absorb(countrycode year)

outreg2 using "$overleaf/tables/PS2_q8.tex", append tex(fragment) /// 
ctitle("Both FE") ///
addtext(Country FE, Yes, Year FE, Yes) label


************	question 1.9
**************************************************************

// Plot the path for GDP per capita for a subset of countries you find interesting. 
// Are the patterns in line with those obtained in numeral 4 or 5?


* Colombia
* Argentina
* South Korea
* US
* UK



tw (line gdp_pc year if country == "Colombia") ///
(line gdp_pc year if country == "Argentina") ///
(line gdp_pc year if country == "Republic of Korea") ///
(line gdp_pc year if country == "United States") ///
(line gdp_pc year if country == "United Kingdom"), ///
legend(order(1 "Colombia" 2 "Argentina" 3 "Republic of Korea" ///
4 "United States" 5 "United Kingdom") cols(2)) ///
ytitle("GDP per capita{sub:t}") xscale(range(1950 2020)) ///
 xlabel(1950(10)2020)
 
graph export "$overleaf/figures/PS2_q9.png", replace



************	question 1.8 every 5 years
**************************************************************

// Sort the dataset by country and year
sort countrycode year

// Create a new variable for the decade
gen time_frame = floor((year - 1950) / 5) * 5 + 1950


collapse (mean) gdp_pc, by(countrycode_num time_frame)

sort countrycode time_frame

gen ln_gdp_pc = ln(gdp_pc)
label var ln_gdp_pc "$ ln(y_{t}) $"

by countrycode: gen ln_gdp_pc_plus1 = ln_gdp_pc[_n+1]

xtset countrycode_num time_frame

sort countrycode time_frame

gen delta_ln_gdp_pc = ln_gdp_pc_plus1-ln_gdp_pc


xtset countrycode_num time_frame

label var time_frame "time"

* eq 3 no fixed effects
reg delta_ln_gdp_pc ln_gdp_pc time_frame
outreg2 using "$overleaf/tables/PS2_q8_2.tex", replace tex(fragment) /// 
ctitle("Eq. (3)") ///
addtext(Country FE, No, Year FE, No) label


* eq 3 country fixed effects
reghdfe delta_ln_gdp_pc ln_gdp_pc time_frame, absorb(countrycode_num)
outreg2 using "$overleaf/tables/PS2_q8_2.tex", append tex(fragment) /// 
ctitle("Trend + Country FE") ///
addtext(Country FE, Yes, Year FE, No) label


* eq 3 random effects 
xtreg delta_ln_gdp_pc ln_gdp_pc time_frame
outreg2 using "$overleaf/tables/PS2_q8_2.tex", append tex(fragment) /// 
ctitle("Random Effects") ///
addtext(Country FE, No, Year FE, No) label


* year fixed effects
reghdfe delta_ln_gdp_pc ln_gdp_pc, absorb(time_frame)
outreg2 using "$overleaf/tables/PS2_q8_2.tex", append tex(fragment) /// 
ctitle("Year FE") ///
addtext(Country FE, No, Year FE, Yes) label

* eq 3 country & year fixed effects
reghdfe delta_ln_gdp_pc ln_gdp_pc, absorb(countrycode_num time_frame)

outreg2 using "$overleaf/tables/PS2_q8_2.tex", append tex(fragment) /// 
ctitle("Both FE") ///
addtext(Country FE, Yes, Year FE, Yes) label
		


























