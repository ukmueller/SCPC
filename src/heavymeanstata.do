capture program drop robttest
program robttest, eclass
//	capture noisily setscores
//	if _rc!=0 exit
	rtt_setBM
	rtt_setscores
	matrix b=e(b)
	local na : colnames e(V)
	matrix robCis =e(b)'*(1,1)
	matrix robCI =(1,2)
	matrix robpval=(1)
	matrix robpvals=e(b)
	tokenize `na'
	qui sum rtt_sel
	if r(sum)>=50{
		local k 8
		disp "proceeding with default of k=8"
		}
	else {
		if r(sum)>=25 local k 4
		else {
			disp "number of independent obs/clusters smaller than 25; aborting"
			exit(999)
		}
	}
	
	matrix WR = J(colsof(b),`k',.) 
	matrix rownames WR=`na'
	matrix WL = WR
//matrix myv=e(b)
	local i = 1
	while "``i''" != "" {
		capture drop w
		matrix coef=r(sum)*rtt_Bread[`i',1...]
		quietly matrix score w=coef if rtt_sel==1
		quietly replace w=w+b[1,`i'] if rtt_sel==1		
//qui sum w if rtt_sel==1
//matrix myv[1,`i']=r(Var)/r(N)
		mata setCIfromS(`k')
		matrix robCis[`i',1]=robCI[1,1..2]
		mata setpvalfromS(`k')
		matrix robpvals[1,`i']=robpval[1,1]
		matrix Wextremes=Wextremes'
		matrix WR[`i',1]=Wextremes[1,1...]
		matrix WL[`i',1]=Wextremes[2,1...]
		local ++i
	}
	matrix A=(b \ robpvals \ robCis')'
	matrix colnames A = Coef p-val "95% Conf" "Interval "
	local n = rowsof(A)
	local ands = `n'*"&"
	local rs &-`ands'
	matlist A, border(rows) title("Results using robust t-tst, k=`k'") cspec(o2& %12s | %9.0g o2 & o1 %5.3f & o2 %9.0g o1 &  o1 %9.0g o2&) rspec(`rs')
	matlist WR, title("Normalized largest k=`k' terms in right tail") names(rows)
	matlist WL, title("Normalized largest k=`k' terms in left tail") names(rows)
	ereturn matrix robCIs = robCis
	ereturn matrix robpvals = robpvals
end

capture program drop rtt_setBM
program rtt_setBM
	matrix rtt_Bread=e(V_modelbased)
	if( rtt_Bread[1,1]==.) {
		disp "e(V_modelbased) missing; aborting"
		exit(999)
	}
end program

capture program drop rtt_setscores
program rtt_setscores, sortpreserve
	tempvar e
	local slist ""
//	quietly{
	if(e(cmd)=="regress" | e(cmd)=="areg" | e(cmd)=="logit" | e(cmd)=="probit"){
		predict `e',sc
		local na : colnames e(V)
		tokenize `na'
		local i = 1
		while "``i''" != "" {
			capture drop rtt_score`i'
			if("``i''"=="_cons") gen rtt_score`i'=`e' if e(sample)
			else gen rtt_score`i'=`e'*``i'' if e(sample)
			local slist "`slist' rtt_score`i'"
			local ++i
		}
	}
	else if(e(cmd)=="ivregress"){
		tempvar stata_scr
		predict `stata_scr'*,sc
		
		local j=`: word count `e(instd)''+1
		gen `e'=`stata_scr'`j'
		di `j'
		local na : colnames e(V)
		tokenize `na'
		local i = 1
		while "``i''" != "" {
			capture drop rtt_score`i'
			if (`i'<`j') gen rtt_score`i'=`stata_scr'`i' if e(sample)
			else{
				if("``i''"=="_cons") gen rtt_score`i'=`e' if e(sample)
				else gen rtt_score`i'=`e'*``i'' if e(sample)
			}
			local slist "`slist' rtt_score`i'"
			local ++i
		}
	}
	else {
		di "rtt_setscores does not support " e(cmd)
		exit(999)
	}
	
	matrix colnames rtt_Bread = `slist'
	capture drop rtt_sel
	gen rtt_sel=0
	if( "`e(clustvar)'"=="") replace rtt_sel=1 if e(sample)
	else{
		tempvar TotScore
		sort `e(clustvar)'
		local i = 1
		while "``i''" != ""{
			capture drop `TotScore'
			by `e(clustvar)': egen `TotScore'=total(rtt_score`i') if e(sample)
			replace rtt_score`i'=`TotScore'
			by `e(clustvar)': replace rtt_sel = 1 if _n==1 & e(sample)
			local ++i
		}
	}
//	}
end

capture program drop rtt_analyticV
program rtt_analyticV
	qui corr rtt_score* if rtt_sel==1, cov 
	tempname V
	mat `V'=r(N)*rtt_Bread*r(C)*rtt_Bread'
	mat list `V'
	mat list e(V)
	mata Vx=st_matrix("`V'")
	mata V=st_matrix("e(V)")
	mata Vx:/V
//	end
end program




//rtt_setBM
//rtt_setscores
//rtt_analyticV
