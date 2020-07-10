

set matastrict on
mata mata clear
set matalnum on

// Locate "all.sd" file
qui findfile "all.sd"
local all_sd `r(fn)'

mata:

struct tst{
	real vector lam_stage1, lam
	real matrix	thlist0_stage1, thlist0
	real scalar level, levelu, k
	real vector cv_zt
}

real scalar getdens_base(real vector msx, real vector Y)
{
	real scalar	m, sig, xsi, zk, k

	k=length(Y)
	m=msx[1]
	sig=msx[2]
	xsi=msx[3]
	zk=1+xsi*(Y[k]/sig-m)
	if(zk<0) return(0)
	if(1+xsi*(Y[1]/sig-m)<0) return(0)
	return(exp(-k*log(sig)-zk^(-1/xsi)-(1+1/xsi)*sum(log(1.0:+xsi*(Y/sig:-m)))))	
}	

real scalar getmeanadj(real vector msx, real scalar yk) 
{
	real scalar	m, sig, xsi, zk
	m=msx[1]
	sig=msx[2]
	xsi=msx[3]
	zk=1+xsi*(yk/sig-m)
	if(zk<0) return(0)
    return((sig/xsi)*zk^(-1/xsi)*(zk/(1-xsi)+xsi*m-1))
}

real scalar getdens0_stage1(real vector msx, real vector Y, real scalar Zsum, real scalar var)
{
	real scalar ma
	ma=getmeanadj(msx,Y[length(Y)])
	return(getdens_base(msx,Y)*exp(-0.5*(Zsum+ma)^2/var)/sqrt(var))
}

real scalar getdens0_stage2(real vector th, real vector Y1, real vector Y2, real scalar Zsum)
{
	real scalar ma, val,k
	k=length(Y1)
	ma=getmeanadj(th[1::3],Y1[k])-getmeanadj(th[4::6],Y2[k])
	val=getdens_base(th[1::3],Y1)*getdens_base(th[4::6],Y2)*exp(-0.5*(Zsum+ma)^2)
	ma=getmeanadj(th[4::6],Y1[k])-getmeanadj(th[1::3],Y2[k])
	val=val+getdens_base(th[4::6],Y1)*getdens_base(th[1::3],Y2)*exp(-0.5*(Zsum+ma)^2)
	return(0.5*val)
}

real scalar getf1(real vector Y)
{
	real scalar val,s,xsi,k
	real vector xs
	real scalar i,j
	real scalar nxsigrid
	
	external real matrix GQxw
	
	k=length(Y)
	if(Y[1]/Y[k]<1.05 & Y[k]>0) return(1E100)	
	nxsigrid=10
	xs=(Y:-Y[k])/(Y[1]-Y[k])
	val=0
	for(j=1;j<=nxsigrid;j++){
		xsi=-0.49+(j-1)/(nxsigrid-1)
		for(i=1;i<=40;i++){
			s=(GQxw[i,1]^(-xsi)-1)/xsi
			val=val+GQxw[i,2]*gamma(k-xsi)*GQxw[i,1]^(-xsi-1)*exp(-(1.0+1.0/xsi)*sum(log(1.0:+xsi*s*xs))+(k-1)*log(s))
		}
	}
	val=val/nxsigrid
	val=val/(Y[1]-Y[k])^k
	return(val)
} 

real scalar switchind(real vector Y)
{
	real scalar k,val
	k=length(Y)
	val=Y[1]-0.5
	if(Y[k]>0){
		if(k==8) val=min((val,Y[1]/Y[k]-1.8))
		if(k==4) val=min((val,Y[1]/Y[k]-1.5))
	}
	return(max((0.0,val)))
}

real scalar getlLR1(real vector Y1, real vector Y2, real scalar Y0, real scalar fa, struct tst scalar tst)
{
	real scalar i,h,ma,var,Zsum
	Zsum=Y0-sum(Y2)
	var=1+sum(Y2:^2)
	h=0
	for(i=1;i<=length(tst.lam_stage1);i++){
		ma=getmeanadj(tst.thlist0_stage1[.,i],Y1[length(Y1)])
		h=h+tst.lam_stage1[i]*getdens_base(tst.thlist0_stage1[.,i],Y1)*exp(-0.5*(Zsum+ma)^2/var)/sqrt(var)
	}
	if(h>1E-100){
		h=-5*switchind(Y1)+log(h/fa)
		}
	else{
		h=-100
	}
	return(h)	
}

real scalar getlLR2(real vector Y1, real vector Y2, real scalar Y0, real scalar fa, struct tst scalar tst)
{
	real scalar i,h
	h=0
	for(i=1;i<=length(tst.lam);i++){
		h=h+tst.lam[i]*getdens0_stage2(tst.thlist0[.,i],Y1,Y2,Y0)
	}
	return(log(h/fa))
}

real scalar gett2(real vector Y1, real vector Y2, real scalar Y0, struct tst scalar tst)
{
	real scalar var,cv,cvw

	var=1.0+sum(Y1:^2)+sum(Y2:^2)
	cvw=1.0/var
	cv=cvw*tst.cv_zt[1]+(1.0-cvw)*tst.cv_zt[2]
	return(cv^2-(Y0+sum(Y1)-sum(Y2))^2/var)
}
		
real colvector getLRs(real vector Y1, real vector Y2, real scalar Y0, struct tst scalar tst)
{
	real scalar fa1,fa2
	fa1=getf1(Y1); fa2=getf1(Y2)
	return((getlLR2(Y1,Y2,Y0,fa1*fa2,tst)\getlLR1(Y1,Y2,Y0,fa1,tst)\getlLR1(Y2,Y1,-Y0,fa2,tst)\gett2(Y1,Y2,Y0,tst)))
}
	
real scalar getcross(real vector x, real vector y)
{
	real scalar	 a,b,c,d,val
	d=((x[1] - x[2])*(x[1] - x[3])*(x[2] - x[3]))
	a=(x[3]*(-y[1] + y[2]) + x[2]*(y[1] - y[3]) + x[1]*(-y[2] + y[3]))/d
	b=(x[3]^2*(y[1] - y[2]) + x[1]^2*(y[2] - y[3]) + x[2]^2*(-y[1] + y[3]))/d
	c=(x[3]*(x[2]*(x[2] - x[3])*y[1] + x[1]*(-x[1] + x[3])*y[2]) + x[1]*(x[1] - x[2])*x[2]*y[3])/d
	d=b^2-4*a*c
	if(d<0) return(1000)
	val=(-b+sqrt(d))/(2*a)
	if(val>=x[1]-1E-8 & val<=x[3]+1E-8) return(val)
	val=(-b-sqrt(d))/(2*a)
	if(val>=x[1]-1E-8 & val<=x[3]+1E-8) return(val)
	return(1000)
}

real scalar getCIlb(real scalar n_sample, real vector Y1, real vector Y2, real scalar Y0, struct tst scalar tst)
{
	real scalar i,mu,s,maxLR,k,mustep,xp,j
	real vector cmus
	real matrix cLRs
	
	k=tst.k
	mu=(sum(Y1)-sum(Y2)+Y0)/n_sample
	s=sqrt(1+sum((Y1:-mu):^2)+sum((Y2:+mu):^2))/n_sample
	mustep=max((1.0/(n_sample-2*k),.1*s))

	cmus=J(3,1,.)
	cLRs=J(4,3,0)
	mu=mu-5*s
	do{
		mu=mu-s
		cLRs[.,2]=getLRs(Y1:-mu,Y2:+mu,Y0-(n_sample-2*k)*mu,tst)
	} while(max(cLRs[.,2])>-4)

	cmus[2]=mu
	mu=mu+mustep
	cmus[3]=mu
	cLRs[.,3]=getLRs(Y1:-mu,Y2:+mu,Y0-(n_sample-2*k)*mu,tst)
	do{
		maxLR=max(cLRs[.,3])
		if(maxLR<-4.0){
			mu=mu+2*mustep
		}
		else{
			if(maxLR<-1.5){
			mu=mu+mustep
			}
			else{
			mu=mu+.2*mustep
			}
		}
		cmus[1..2]=cmus[2..3]
		cmus[3]=mu
		cLRs[.,1..2]=cLRs[.,2..3]
		cLRs[.,3]=getLRs(Y1:-mu,Y2:+mu,Y0-(n_sample-2*k)*mu,tst)
		xp=1000
		for(i=1;i<=4;i++){
			xp=min((xp,getcross(cmus,cLRs[i,.])))
		}
		if(xp<1000) return(xp)
	} while(1)
}

real vector getCI(real scalar n_sample, real vector Y1, real vector Y2, real scalar Y0, struct tst scalar tst)
{
	return((getCIlb(n_sample,Y1,Y2,Y0,tst),-getCIlb(n_sample,Y2,Y1,-Y0,tst)))
}

real scalar reject(real vector Y1, real vector Y2, real scalar Zsum, struct tst scalar tst)
{
	return(max(getLRs(Y1,Y2,Zsum,tst))<0)
}

void setYsfromW(real vector W, real scalar k, real vector Y1, real vector Y2, real scalar Zsum, real scalar s)
{
	real matrix Wo
	real scalar n
	real vector ind

	n=length(W)
	Wo=J(n,1,0)
	Wo[.,1]=W
	Wo=sort(Wo,1)
	ind=(1::k)
	ind=(n+1):-ind
	Y1=Wo[ind,1]
	Y2=-Wo[(1::k),1]
	Zsum=sum(Wo[(k+1::n-k),1])
	s=sqrt(sum(Wo[(k+1::n-k),1]:^2)-Zsum^2/(n-2*k))
	Y1=Y1/s
	Y2=Y2/s
	Zsum=Zsum/s
}

void mytest(real scalar ind, string scalar all_sd)
{
	real matrix W
	real vector Y1,Y2,v
	real scalar Zsum,s,k
	external struct tst vector tsts
	
	st_view(W=.,.,ind,0)
	v=W[.,1]
	k=4
	loadall(all_sd)
	Y1
	Y2
	Zsum
	getCI(v,tsts[3])	
}

real matrix getCIfromW(real vector W, struct tst scalar tst)
{
	real vector Y1, Y2
	real scalar scale,n, k, Zsum
	real matrix CI
	
	setYsfromW(W,tst.k,Y1,Y2,Zsum,scale)
	n=length(W)
	CI=scale*getCI(n,Y1,Y2,Zsum,tst)
	return(CI)
}

void setCIfromS(real scalar k, string scalar w, string scalar rtt_sel, string scalar all_sd)
{
	external struct tst matrix tsts
	real matrix W
	loadall(all_sd)
	W=st_data(., w, rtt_sel)
	st_replacematrix("robCI",getCIfromW(W[.,1],tsts[k/4,1]))
}

void setpvalfromS(real scalar k, string scalar w, string scalar rtt_sel, string scalar all_sd)
{
	external struct tst matrix tsts
	real matrix W
	real matrix p
	real vector Y1, Y2,CI
	real scalar n,scale, Zsum,j,pu,pl,pc,level,tj

	loadall(all_sd)
	W=st_data(., w, rtt_sel)
	n=length(W)
	tj=k/4
	setYsfromW(W,k,Y1,Y2,Zsum,scale)	
	CI=getCI(n,Y1,Y2,Zsum,tsts[tj,1])
	
	if(CI[1]<0 & CI[2]>0){
		CI=getCI(n,Y1,Y2,Zsum,tsts[tj,3])		
		if(CI[1]<0 & CI[2]>0){
			pu=1000; pl=100	; pc=tsts[tj,3].levelu		
			}
		else{
			pu=100; pl=50 ; pc=tsts[tj,1].levelu
			}
		}
	else{
		CI=getCI(n,Y1,Y2,Zsum,tsts[tj,2])
		if(CI[1]<0 & CI[2]>0){
			pu=50; pl=10; pc=tsts[tj,2].levelu
			}
		else{
			pu=10; pl=2; pc=2
			}
		}

	for(j=4;j<=cols(tsts);j++){
		level=tsts[tj,j].level		
		if(level>pu) continue
		if(level<pl) continue	
		if(!reject(Y1,Y2,Zsum,tsts[tj,j])) pc=max((tsts[tj,j].levelu,pc))
	}	
	p=pc/1000
	st_replacematrix("robpval",p)
	st_matrix("Wextremes", (Y1 , Y2))

}


	
void addmat(string scalar all_sd)
{
	real scalar fh
	real matrix X,Xs
	st_view(X=.,.,.)
	fh=fopen(all_sd, "a")
	Xs=X[.,(2::cols(X))]
	fputmatrix(fh,Xs)
	fclose(fh)
}

void saveall(string scalar all_sd)
{
	real vector levellist
	string scalar fname, kname
	real scalar i,level,k
	
	fname="c:/dropbox/mystuff/heavymean/2020/swind/GQxw.txt"
	fname="import delimited "+fname+", delim("+char(34)+" "+char(34)+", collapse)"
	fname
	stata("clear")
	stata(fname)
	addmat(all_sd)

	levellist=(50,10,100,2,4,6,8,20,30,40,60,70,80,90,120,140,160,180,200,250,300,400,500)
	for(k=4;k<=8;k=k+4){
		if(k==4) kname="_4_25_"
		else kname="_8_50_"
		for(i=1;i<=length(levellist);i++){
			level=levellist[i]
			fname="c:/dropbox/mystuff/heavymean/2020/swind/lamth_stage1"+kname+strofreal(level)+".txt"
			fname="import delimited "+fname+", delim("+char(34)+" "+char(34)+", collapse)"
			fname
			stata("clear")
			stata(fname)
			addmat(all_sd)
			fname="c:/dropbox/mystuff/heavymean/2020/swind/lamth"+kname+strofreal(level)+".txt"
			fname="import delimited "+fname+", delim("+char(34)+" "+char(34)+", collapse)"
			fname
			stata("clear")
			stata(fname)
			addmat(all_sd)
		}
	}
}

void loadall(string scalar all_sd)
{
	real scalar i,level,fh,lx,j
	real vector levellist
	external struct tst matrix tsts
	external real matrix GQxw
	real matrix mdat
	levellist=(50,10,100,2,4,6,8,20,30,40,60,70,80,90,120,140,160,180,200,250,300,400,500,1000)
	tsts=tst(2,length(levellist))
	fh=fopen(all_sd, "r")
	GQxw=fgetmatrix(fh)
	for(j=1;j<=2;j++){
		for(i=1;i<=length(levellist)-1;i++){
			level=levellist[i]
			mdat=fgetmatrix(fh)
			tsts[j,i].lam_stage1=mdat[1,.]
			tsts[j,i].thlist0_stage1=mdat[2::4,.]
			mdat=fgetmatrix(fh)
			tsts[j,i].lam=mdat[1,.]
			tsts[j,i].thlist0=mdat[2::7,.]
			tsts[j,i].level=level
			tsts[j,i].levelu=levellist[i+1]
			tsts[j,i].k=4*j
			lx=level/1000
			tsts[j,i].cv_zt=(invnormal(1-lx/2),invttail(80+10*log(lx), lx/2))
		}
	}
	fclose(fh)
}
	

end



mata loadall("`all_sd'")

