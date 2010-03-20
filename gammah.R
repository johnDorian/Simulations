gammah<-function(h,modtyp,co,c1,a1,c2,a2){
	alpha=1.75
	sill=co+c1+c2

	switch(modtyp,
		iso.lin=return(co*h),
		iso.cir=
		if(h>a1)
			return(sill)
		else{
			hovera=h/a1
			angle=atan(sqrt(1-hovera^2)/(hovera+0.000001))
			h4=twoopi*hovera*sqrt(1-hovera^2)
			return(co+c1*(twoopi*angle+h4))
		},
		
		iso.sph=
		if(h>a1)
			return(sill)
		else{
			hovera=h/a1
			return(co+c1*(1.5*hovera-0.5*hovera^3))
		},
		iso.exp=return(co+c1*(1-exp(-h/a1))),
		
		iso.pen=if(h>a1)
			return(sill)
		else
			return(co+c1*(1.875*(h/a1)-1.25*(h/a1)^3+0.375*(h/a1)^5)),
		
		iso.dou=if(h>a2)
			return(sill)
		else
			if(h>a1)
				return(co+c1+c2*(1.5*(h/a2)-0.5*(h/a2)^3))
			else
				return(co+c1*(1.5*(h/a1)-0.5*(h/a1)^3)+c2*(1.5*(h/a2)-0.5*(h/a2)^3)),

		dou.exp=return(co+c1(1-exp((-h/a1)))+c2*(1-exp((-h/a2)))),
		
		stable=return(co+c1*(1-exp(-h^alpha/a1^alpha))),
		
		dou.sta=return(co+c1*(1-exp(-h^alpha/a1^alpha))+c2*(1-exp(-h/a2))),
		
		bro=return(co*h^a1)
	)
}
		
		