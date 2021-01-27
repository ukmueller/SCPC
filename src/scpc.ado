

set matastrict on
mata mata clear
mata 

real scalar getavc(real scalar c, real vector dist)
// computes average correlation given vector of pairwise distances
{
	return(sum(exp(-c*dist))/length(dist))
}	

real matrix getdistmat(real matrix s)
// computes matrix of distances from locations
{
	real matrix mat
	real scalar n,i
	n=rows(s)
	mat=J(n,n,0)
	for(i=1;i<=n;i++){
		mat[.,i]=sqrt(rowsum((s:-s[i,.]):^2))
	}
	return(mat)
}	

real vector lvech(real matrix mat)
// takes lower triangular part of matrix and puts it into vector
{
	real scalar i,c
	real vector v
	v=J(rows(mat)*(rows(mat)-1)/2,1,0)
	c=1
	for(i=1;i<=rows(mat)-1;i++){
		v[c::c+rows(mat)-i-1]=mat[i+1::rows(mat),i]	
		c=c+rows(mat)-i
	}	
	return(v)
}	

real matrix demeanmat(real matrix mat)
// demeans matrix
{
	mat=mat:-(rowsum(mat)/cols(mat))
	mat=mat:-(colsum(mat)/rows(mat))
	return(mat)
}

real scalar getc0fromavc(real vector dist,real scalar avc0)
// solves for c0 from vector of pairwise distances given avc0
{
	real scalar c0,c1,c
	
	c0=10
	do{
		c0=.5*c0
	} while(getavc(c0,dist)<avc0)
	c1=20
	do{
		c1=2*c1		
	} while(getavc(c1,dist)>avc0)	
	do{
		c=sqrt(c1*c0)
		if(getavc(c,dist)>avc0){
			c0=c
			}
		else{
			c1=c
		}
	} while(c1-c0>0.001)
	return(c)
}

real matrix getW(real matrix distmat, real scalar c0)
// computes first qmax eigenvectors of demeaned Sigma(c0) and stores them in W, along with column vector of ones
// all columns of W are normalized to length 1
{
	external real scalar qmax
	real vector	evals
	real matrix W

	W=J(rows(distmat),qmax,1)
	evals=(1::qmax)
	symeigensystemselecti(demeanmat(exp(-c0*distmat)),(1,qmax),W,evals)
	W=(J(rows(distmat),1,1/sqrt(rows(distmat))),W)
	return(W)
}

struct mats{
// helper struct for vector of matrices Oms
		real matrix mat
}

real scalar getrp(real matrix Om, real scalar cv)
{
// computes rejection probability given Om and cv; cv here is cv/sqrt(q) in notation of paper
	complex vector cevals
	real vector evals
	real matrix Omx
	external real matrix GQxw
	real scalar i,u,val
	
	Omx=-cv^2*Om
	Omx[1,.]=Om[1,.]
	cevals=eigenvalues(Omx)
	evals=Re(cevals)	
	evals=-select(evals,evals:<0)/max(evals)		
	val=0
	
	for(i=1;i<=rows(GQxw);i++){
		u=GQxw[i,1]
		val=val+GQxw[i,2]/sqrt((1-u^2)*exp(sum(log(1:+evals/(1-u^2)))))
	}
	return(val*2.0/3.14159265359)
}

real scalar maxrp(struct mats vector Oms, real scalar q, real scalar cv)
{
// computes largest rejection probability given vector of Om matrices stored in Oms given q and cv
	real scalar i, val
	val=0
	for(i=1;i<=length(Oms);i++){
		val=max((val,getrp(Oms[i].mat[1::q+1,1::q+1],cv)))
	}
	return(val)
}

void setfinalW(struct mats vector Oms, real matrix W, real scalar cv)
// solves for optimal q and cv from Oms, stores results in W and cv
{
	external real scalar qmax
	real scalar i, cv0,cv1,q
	real vector cvs,lengths,qvec,rps
	real matrix kmat
	
	lengths=(1::qmax)
	cvs=lengths
	rps=(1::length(Oms))
	for(q=1;q<=qmax;q++){
		i=1			// compute cv1 first for Oms[1]; if that doesn't lead to size control, iterate
		cv0=invttail(q,0.025)/sqrt(q)   
		do{
			cv1=cv0
			while(1){
				if(getrp(Oms[i].mat[1::q+1,1::q+1],cv1)>0.05){
					cv0=cv1
					cv1=cv1+1/sqrt(q)
				}
				else{
					break
				}
			}
			
			while(cv1-cv0>0.01/sqrt(q)){
				cv=0.5*(cv1+cv0)
				if(getrp(Oms[i].mat[1::q+1,1::q+1],cv)>0.05){
					cv0=cv
					}
				else{
					cv1=cv
				}
			}
			
			for(i=1;i<=length(Oms);i++){
				rps[i]=getrp(Oms[i].mat[1::q+1,1::q+1],cv1)		
			}
			maxindex(rps,1,qvec,kmat)			
			i=qvec[1]				// set potential new i to grid index with largest rejection prob
			cv0=cv1
		}while(rps[i]>0.05)		
		cvs[q]=cv1
		lengths[q]=cv1*gamma(.5*(q+1))/(gamma(.5*q))
	} 
	minindex(lengths,1,qvec,kmat)
	q=qvec[1]
	cv=cvs[q]*sqrt(q)	
	W=W[.,1::q+1]
}

void setfinalWold(struct mats vector Oms, real matrix W, real scalar cv)
// solves for optimal q and cv from Oms, stores results in W and cv
{
	external real scalar qmax
	real scalar i, cv0,cv1,q
	real vector cvs,lengths,qvec
	real matrix kmat
	
	lengths=(1::qmax)
	cvs=lengths
	for(q=1;q<=qmax;q++){
		cv0=invttail(q,0.025)/sqrt(q)   
		cv1=cv0
		do{
			cv1=cv1+1/sqrt(q)
		} while(maxrp(Oms,q,cv1)>0.05)
		cv0=cv1-1/sqrt(q)
		do{
			cv=0.5*(cv1+cv0)
			if(maxrp(Oms,q,cv)>0.05){
				cv0=cv
				}
			else{
				cv1=cv
			}
		}while(cv1-cv0>0.01/sqrt(q))		
		cvs[q]=cv1
		lengths[q]=cv1*gamma(.5*(q+1))/(gamma(.5*q))
	} 
	minindex(lengths,1,qvec,kmat)
	q=qvec[1]
	cv=cvs[q]*sqrt(q)	
	W=W[.,1::q+1]
}

real scalar getnc(real scalar c0, real scalar cmax)
// computes number of c-values in grid so that largest c is at least cmax
{
	external real scalar cgridfac
	return(max((2,ceil(log(cmax/c0)/log(cgridfac)))))
}

struct mats vector getOms(real matrix distmat, real scalar c0, real scalar cmax, real matrix W)
// computes vector of qmax x qmax Om(c) matrices from distmat, c0 and W
{
	external real scalar cgridfac
	struct mats vector Oms
	real scalar i, c
	
	
	Oms=mats(getnc(c0,cmax))
	c=c0
	for(i=1;i<=length(Oms);i++){	
		Oms[i].mat=cross(W,exp(-c*distmat))*W 
		c=c*cgridfac			
	}
	return(Oms)
}

real scalar gettau(real vector y, real matrix W)
// computes SCPC t-stat
{
	return(sqrt(cols(W)-1)*cross(W[.,1],y)/norm(cross(W[.,2::cols(W)],y)))
}


void setGQxw()
{
// sets gaussian quadrature weights
	external real matrix GQxw
	GQxw=(  8.811451447204854404E-004,  4.636880650271735238E-003,  1.137002500811307160E-002,  2.104159039310410373E-002,  3.359359586066185122E-002,  4.895059651556293856E-002,  6.702024839387038524E-002,  8.769388458334448355E-002,  1.108471742867402909E-001,  1.363408724050364507E-001,  1.640216576929103831E-001,  1.937230551660097388E-001,  2.252664374524356861E-001,  2.584620991569107629E-001,  2.931103978141975097E-001,  3.290029545871209216E-001,  3.659239074963732685E-001,  4.036512096493146129E-001,  4.419579646623724711E-001,  4.806137912469746754E-001,  5.193862087530253246E-001,  5.580420353376275289E-001,  5.963487903506855536E-001,  6.340760925036268425E-001,  6.709970454128790784E-001,  7.068896021858023238E-001,  7.415379008430891261E-001,  7.747335625475638698E-001,  8.062769448339899281E-001,  8.359783423070896724E-001,  8.636591275949634383E-001,  8.891528257132595980E-001,  9.123061154166557385E-001,  9.329797516061295592E-001,  9.510494034844372280E-001,  9.664064041393380933E-001,  9.789584096068959518E-001,  9.886299749918872060E-001,  9.953631193497285423E-001,  9.991188548552796256E-001 \ 2.260638549266569437E-003,  5.249142265576826200E-003,  8.210529190953967313E-003,  1.112292459708323286E-002,  1.396850349001168591E-002,  1.673009764127359952E-002,  1.939108398723625448E-002,  2.193545409283656489E-002,  2.434790381753590763E-002,  2.661392349196848098E-002,  2.871988454969627222E-002,  3.065312124646369166E-002,  3.240200672830125667E-002,  3.395602290761685904E-002,  3.530582369564324446E-002,  3.644329119790197524E-002,  3.736158452898413751E-002,  3.805518095031324571E-002,  3.851990908212438863E-002,  3.875297398921226377E-002,  3.875297398921258990E-002,  3.851990908212385434E-002,  3.805518095031288489E-002,  3.736158452898429017E-002,  3.644329119790179483E-002,  3.530582369564296691E-002,  3.395602290761668557E-002,  3.240200672830127054E-002,  3.065312124646427452E-002,  2.871988454969589405E-002,  2.661392349196792934E-002,  2.434790381753596661E-002,  2.193545409283669673E-002,  1.939108398723569590E-002,  1.673009764127403667E-002,  1.396850349001163734E-002,  1.112292459708396665E-002,  8.210529190953890985E-003,  5.249142265576283231E-003,  2.260638549266727297E-003)'
}

/* 
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXX now routines for n>2000 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  */


void normalize_s(real matrix s)
// normalizes locations to ensure invariance given probabilistic algorithm
{
	real matrix evecs
	real vector evals
	external real vector permfin
	
	s=s:-mean(s)		
	symeigensystem(cross(s,s),evecs,evals)
	s=s*evecs		
	if(max(s[.,1])!=max(abs(s[.,1]))) s=-s
	permfin=order(s,1)	
	s=s[permfin,.]
	s=s:-colmin(s)
	s=s/max(s)	
}

real scalar nextU()
// simple pseudo-random number generator, state (=seed) in global variable random_t 
{
	external real scalar random_t
	random_t=mod(64389*random_t + 1,2^32)
	return(random_t/2^32)
}

void jumble_s(real matrix s, real scalar m)
// randomly jumbles first m rows of location matrix s
{
	real scalar i,j,n
	real vector sj
	n=rows(s)	
	for(i=1;i<=m;i++){
		j=floor(nextU()*n)+1		
		sj=s[j,.]
		s[j,.]=s[i,.]
		s[i,.]=sj
	}
}	

struct sms{
// helper struct for capN random sets of m locations: mat stores distmat, s stores locations
		real matrix mat
		real matrix s
}

void LNsetWc0(real matrix s,real scalar avc,real matrix W, real scalar c0, real scalar cmax)
{
// computes W and c0 when n is large
	external real scalar m,capN
	external real scalar qmax, minavc
	struct sms vector ms
	real vector distvec
	real rowvector evals
	real colvector v
	real scalar i,j,n
	real matrix r,W0,Wx,Wall,evecs
	
	n=rows(s)
	ms=sms(capN)
	distvec=J(capN*m*(m-1)/2,1,0) // vector of all pairwise distances across capN random locations
	
	r=s	// ensure s does not get jumbled
	for(i=1;i<=capN;i++){
		jumble_s(r,m)		
		ms[i].s=r[1::m,.]
		ms[i].mat=getdistmat(ms[i].s)
		distvec[(i-1)*m*(m-1)/2+1::i*m*(m-1)/2]=lvech(ms[i].mat)
	}
	c0=getc0fromavc(distvec,avc)
	cmax=getc0fromavc(distvec,minavc)
	
	Wall=J(n,capN*qmax,0)		// matrix of all eigenvectors 1,...,qmax across all capN random locations
	for(i=1;i<=capN;i++){
		symeigensystemselecti(demeanmat(exp(-c0*ms[i].mat)),(1,qmax),W0,evals)
		Wx=J(n,qmax,0)		// m x qmax eigenvectors from ith random set of m locations
		for(j=1;j<=m;j++){
			v=exp(-c0*sqrt(rowsum((s:-ms[i].s[j,.]):^2)))
			Wx=Wx+v*W0[j,.]		// Nystrom extension 
		}
		Wx=Wx:-mean(Wx)
		Wx=Wx:/sqrt(colsum(Wx:^2))				
		for(j=1;j<=qmax;j++){
			Wall[.,(j-1)*capN+i]=Wx[.,j] // store such that first capN vectors correspond to eigenvectors corresponding to largest eigenvalue, etc
		}
	}	
	W=J(n,qmax,0)
	for(i=1;i<=qmax;i++){
		Wx=Wall[.,1::capN*i]
		symeigensystemselecti(cross(Wx,Wx),(1,1),evecs,evals)	
		W[.,i]=Wx*evecs		// first (sample) principle component of first 1,...,q*capN columns of Wall
		W[.,i]=W[.,i]/sqrt(colsum(W[.,i]:^2))
		Wall=Wall-W[.,i]*cross(W[.,i],Wall) // project off new principle component
	}
	W=(J(n,1,1/sqrt(n)),W)
}

real vector raninds(real scalar n,real scalar capM)
{
// yields capM+1 vector of random indices in {1,...,n} such that each pair is distinct
 	real scalar i,j
	real vector v

	v=J(capM+1,1,0)
	j=floor(n*nextU())
	for(i=1;i<=capM+1;i++){
		v[i]=j+1
		j=mod(j+1+floor(nextU()*(n-1)),n)
	}
	return(v)
}

struct mats vector LNgetOms(real matrix s, real scalar c0, real scalar cmax,real matrix W)
{
// version of getOms that is faster if n is very large using a random selection of capM pairwise distances
	external real scalar capM, cgridfac
	struct mats vector Oms
	external real scalar qmax
	real scalar i, c, n
	real vector inds,dist,cd
	real matrix W1,W2
	
	n=rows(s)
	Oms=mats(getnc(c0,cmax))
	inds=raninds(n,capM)
	dist=sqrt(rowsum((s[inds[1::capM],.]-s[inds[2::capM+1],.]):^2)) // capM distance vector of pairwise distances
	W1=W[inds[1::capM],.]		
	W2=W[inds[2::capM+1],.]
	c=c0
	for(i=1;i<=length(Oms);i++){
		cd=exp(-c*dist)
		Oms[i].mat=I(qmax+1)+.5*(n*(n-1)/capM)*(cross(W1,cd,W2)+cross(W2,cd,W1))
		c=c*cgridfac			
	}
	return(Oms)
}


/* 
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXX driver routines XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  */

void setOmsWfin(real scalar avc0, string scalar sel)
{
	real matrix s
	real scalar c0,cmax,cv,n
	real matrix W,distmat
	external real scalar random_t	
	struct mats vector Oms
	external real scalar qmax, capM, m, cgridfac, capN, minavc

	external struct mats vector Omsfin
	external real matrix Wfin
	external real vector permfin
	external real scalar cvfin
	
	s=st_data(.,"s_*",sel)		// expects locations in s_1, s_2, s_3... etc. Dimension d is equal to the number of s_ variables
	n=rows(s)
	if(n<5){
		stata("disp as text"+char(34)+"no locations found in variables s_1, s_2 etc; aborting"+char(34))
		exit(999)
	}		
	setGQxw()
	minavc=0.00001			// minimal avc value for which size control is checked
	cgridfac=1.2			// factor in c-grid for size control
	
	if(avc0>=0.05){
		qmax=10
	}
	else{
		if(avc0>=0.01){
		qmax=20
		}
		else{
			if(avc0>=0.005){
				qmax=60
			}
			else{
				qmax=120
			}
		}
	}
	while(1){
		if(n<2000){
	// code for small n
			permfin=(1::n)
			distmat=getdistmat(s)
			c0=getc0fromavc(lvech(distmat),avc0)		
			cmax=getc0fromavc(lvech(distmat),minavc)	
			W=getW(distmat,c0)
			Oms=getOms(distmat,c0,cmax,W)		
			}
		else{
			capN=20			// number of random sets of locations (of size m) used to approximate eigenvectors
			capM=1000000		// number of location pairs used to approximate Om(c) (in production: capM=1000000 or so)
			m=1000			// size of set of locations for which eigenvector is computed, and then Nystrom extended from m to n (in production: m=1000 or so)
							// note: m cannot be larger than n

	//		stata("display c(current_time)")
			random_t=1
			normalize_s(s)	
			LNsetWc0(s,avc0,W,c0,cmax)
	//		"done with W and c0"
	//		stata("display c(current_time)")

			Oms=LNgetOms(s,c0,cmax,W)	
	//		"done with Oms"
	//		stata("display c(current_time)")
		}
		setfinalW(Oms,W,cv)
		if(cols(W)-1<qmax) break
		qmax=round(qmax+qmax/2,1)
	}
	stata("disp as text"+char(34)+"SCPC using "+strofreal(rows(s), "%6.0f")+" observations / clusters and "+strofreal(cols(s), "%3.0f")+"-dimensional locations in s_*"+char(34))
	stata("disp as text"+char(34)+"resulting optimal q="+strofreal(cols(W)-1, "%3.0f")+" and 5% two-sided critical value = "+strofreal(cv, "%6.3f")+char(34))
	
	Wfin=W
	Omsfin=Oms
	cvfin=cv
}
	
void set_scpcstats(string scalar w, string scalar sel)
{	
	real vector y
	real scalar tau, pval, SE
	external struct mats vector Omsfin
	external real matrix Wfin
	external real vector permfin
	external real scalar cvfin

	y=st_data(.,w,sel)

	y=y[permfin]	
	tau=gettau(y,Wfin)
	pval=maxrp(Omsfin,cols(Wfin)-1,abs(tau)/sqrt(cols(Wfin)-1)) 
	SE=norm(cross(Wfin[.,2::cols(Wfin)],y))/(sqrt(cols(Wfin)-1)*sqrt(rows(Wfin)))
	st_matrix("scpcstats", (mean(y), SE, tau, pval, mean(y)-cvfin*SE,mean(y)+cvfin*SE ))	// coef, SE, tstat, p-value, 95% CI
}


end

capture program drop scpc_setscores
program scpc_setscores, sortpreserve
	syntax, scpc_sel(name)
	tempvar e
	local slist ""
//	quietly{
	if(e(cmd)=="regress" | e(cmd)=="areg" | e(cmd)=="logit" | e(cmd)=="probit"){
		predict `e',sc
		local na : colnames e(V)
		tokenize `na'
		local i = 1
		while "``i''" != "" {
			capture drop scpc_score`i'
			if("``i''"=="_cons") gen scpc_score`i'=`e' if e(sample)
			else gen scpc_score`i'=`e'*``i'' if e(sample)
			local slist "`slist' scpc_score`i'"
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
			capture drop scpc_score`i'
			if (`i'<`j') gen scpc_score`i'=`stata_scr'`i' if e(sample)
			else{
				if("``i''"=="_cons") gen scpc_score`i'=`e' if e(sample)
				else gen scpc_score`i'=`e'*``i'' if e(sample)
			}
			local slist "`slist' scpc_score`i'"
			local ++i
		}
	}
	else {
		di "scpc_setscores does not support " e(cmd)
		exit(999)
	}
	
	matrix colnames scpc_Bread = `slist'
	gen `scpc_sel'=0
	if( "`e(clustvar)'"=="") qui replace `scpc_sel'=1 if e(sample)
	else{
		tempvar TotScore
		sort `e(clustvar)'
		local i = 1
		while "``i''" != ""{
			capture drop `TotScore'
			by `e(clustvar)': egen `TotScore'=total(scpc_score`i') if e(sample)
			replace scpc_score`i'=`TotScore'
			by `e(clustvar)': replace `scpc_sel' = 1 if _n==1 & e(sample)
			local ++i
		}
	}
//	}
end


capture program drop scpc
program scpc, eclass sortpreserve
	syntax [if] [in], ///
		[ ///
			avc(real -1) ///
		]
		//	syntax [if] [in], [
//		avc(real 1.0)
//		]
	
	matrix scpc_Bread=e(V_modelbased)
	if( scpc_Bread[1,1]==.) {
		disp as text "e(V_modelbased) missing; aborting"
		exit(999)
	}	
	if `avc'==-1 {
		local avc 0.05	// set default
	}
	if `avc'<0.001 | `avc'>0.99 {
		disp as text "average correlation must be between 0.001 and 0.99; aborting"
		exit(999)
	}
	tempvar scpc_sel
	scpc_setscores, scpc_sel(`scpc_sel')
	qui sum `scpc_sel', detail
	local neff r(sum)
	matrix b=e(b)
	local na : colnames e(V)
	matrix scpctab =e(b)'*(1,1,1,1,1,1)
	tokenize `na'
	mata setOmsWfin(`avc',"`scpc_sel'")
	
	// Loop through covariates
	local i = 1
	while "``i''" != "" {
		tempvar w
		matrix coef = r(sum)*scpc_Bread[`i', 1...]
		qui matrix score `w' = coef if `scpc_sel' == 1
		qui replace `w' = `w' + b[1,`i'] if `scpc_sel' == 1		
		mata set_scpcstats("`w'", "`scpc_sel'")
//matlist scpcstats 		
		matrix scpctab[`i', 1] = scpcstats[1, 1..6]
//		mata setpvalfromS(`k', "`w'", "`scpc_sel'", "`all_sd'")
		local ++i
	}
//matlist scpctab	
	matrix colnames scpctab = "Coef" "Std_Err" "  t  " "P>|t|" "95% Conf" "Interval"
	local n = rowsof(scpctab)
	local ands = `n'*"&"
	local rs &-`ands'
	local pavc : di %5.3f `avc'
	matlist scpctab, border(all) title("Results from SCPC, maximal average pairwise correlation = `pavc'") cspec(o2& %12s | %9.0g o2 & %9.0g o2 &o1 %5.2f o1& o2 %6.3f o1 & o2 %9.0g & o1 %9.0g o2&) rspec(`rs')
 	// Return results
	ereturn matrix scpcstats = scpctab
	cap drop scpc_score*
end

