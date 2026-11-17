* Run from src so Stata loads scpc.ado from this checkout.
clear all
set more off
which scpc

set obs 32
generate cluster = ceil(_n/4)
generate byte in_sample = mod(_n,7) != 0
generate double s_1 = _n
generate double z1 = sin(.9*_n) + .3*cos(_n)
generate double z2 = cos(.6*_n)
generate double w = sin(.3*_n)
generate double x = z1 + .4*z2 + .2*w + .5*cos(1.3*_n)
generate double y = .5*x + w + sin(.7*_n + .5)

quietly ivregress 2sls y w (x=z1) if in_sample, vce(robust)
quietly scpc, avc(.2) k(1) cvs

mata:
void check_iv_weights(string scalar instruments, real scalar clustered)
{
	external real matrix Wfin, Wx
	real colvector sites, locations, loading
	real matrix observed, Z, P_Z, bread, M, membership, Xi, Xi_s, R, expected, wrong
	real scalar i, n

	observed=st_data(.,("x","w"),"in_sample")
	observed=(observed,J(rows(observed),1,1))
	Z=st_data(.,tokens(instruments),"in_sample")
	Z=(Z,st_data(.,"w","in_sample"),J(rows(Z),1,1))
	P_Z=Z*luinv(Z'*Z)*Z'
	bread=luinv(observed'*P_Z*observed)
	if(cols(Z)==cols(observed)){
		// Section 5.3's just-identified residual maker.
		M=I(rows(observed))-observed*luinv(Z'*observed)*Z'
	}
	else{
		M=I(rows(observed))-observed*bread*observed'*P_Z
	}
	sites=st_data(.,"s_1","in_sample")
	locations=uniqrows(sites)
	n=rows(locations)
	membership=I(n)
	if(clustered){
		membership=J(rows(sites),n,0)
		for(i=1;i<=n;i++) membership[.,i]=(sites:==locations[i])
	}
	loading=n*(P_Z*observed*bread)[.,1]
	Xi=diag(loading)*membership
	Xi_s=Xi*diag(1:/sqrt(colsum(Xi:^2)))
	R=Wfin[.,2..cols(Wfin)]
	expected=(Xi_s'*Xi*Wfin[.,1], Xi_s'*M'*Xi*R)
	wrong=Xi_s'*M*Xi*R
	assert(max(abs(expected[.,2..cols(expected)]-wrong)) > .01)
	assert(max(abs(expected-Wx))/max(abs(expected)) < 1e-6)
}
end

mata: check_iv_weights("z1",0)

quietly replace s_1 = cluster
quietly ivregress 2sls y w (x=z1 z2) if in_sample, vce(cluster cluster)
quietly scpc, avc(.2) k(1) cvs
mata: check_iv_weights("z1 z2",1)

display "Conditional IV order tests passed."
