c	Subroutine that fits the linear model of coregionalisation to a set of 
c	autovariograms and a cross-variogram.  A weighted least-squares criterion 
c	is minimised by simulated annealing.
c
c	Caps text are my own amendments from RML's original code 'FITLMCR3.FOR'. 
c
	SUBROUTINE LMCR(SEMVAR,NOLAGS,WGT,ICVP,CPARF,MODTYP,COVAR,MAXDIST,
     $GUESSA,LOCK,C_OUT,A_OUT,vAIC)

	use numerical_libraries

c	nvar - no. variables, nstr - number of structures (excluding the nugget). 
C	amax - max. value of distance parameter.
c	icvp - set to 1 to constrain covariances to be positive, else 0.
c	a - holds distance parameters.
c	c - holds variances c(0:2 [structures],nvar,nvar).
c	am - holds maximum change of distance parameters.
c	dm - holds maximum change for each parameter (isomorphic to c).
c  
c	DATAF is (nvar,nvar, 3 [h,gamma,npairs,st. dev of gamma at h], nlag [lines]).

	INTEGER,PARAMETER::NVAR=2
	INTEGER::WGT,NOLAGS,MODTYP,I,J,IR,IC,IX,IPA,IO
	INTEGER::ISTR,NSTR,ICONTO,ICPO,IS,NUNCH,IMC,KOUNT,ICVP
	INTEGER::NLAGS(3),LOCK
	REAL,DIMENSION(500,5)::SEMVAR
	REAL::R1,R2,R3
	REAL::COVAR(2,2),SD(2)
	REAL::dataf(NVAR,NVAR,3,500),MAXDIST,AMAX
	REAL::F,FO,FIC,ACCREJ,REJ,ACC,RC,PACC,PACCIN
	REAL::am(2),GUESSA(2),A(2),AO,c(0:2,nvar,nvar),dm(0:2,nvar,nvar)
	real::cpar,cparf,alp,COLD,AMIN,vAIC
	integer iwopt,nmarkov,istop,ics
	DOUBLEPRECISION::R,C_OUT(0:2,nvar,nvar),A_OUT(2)
	

	DATAF=0	
	INDX=0
	KOUNT=0
	do 4 ir=1,nvar
		do 5 ic=ir,nvar	

			INDX=INDX+1

			do 6 ix=1,NOLAGS

				KOUNT=KOUNT+1

				dataf(ir,ic,1,ix)=SEMVAR(IX,1)

				IF(INDX==1)THEN
					dataf(ir,ic,2,ix)=SEMVAR(IX,2)
				ENDIF
				IF(INDX==2)THEN
					dataf(ir,ic,2,ix)=SEMVAR(IX,3)
				ENDIF
				IF(INDX==3)THEN
					dataf(ir,ic,2,ix)=SEMVAR(IX,4)
				ENDIF

				dataf(ir,ic,3,ix)=SEMVAR(IX,5)
c	print*,ir,ic,dataf(ir,ic,1,ix),dataf(ir,ic,2,ix),dataf(ir,ic,3,ix)
6			continue

5		continue
4	continue
	NLAGS=NOLAGS
	

C	Standard deviations of each variable.
	SD(1)=SQRT(COVAR(1,1))
	SD(2)=SQRT(COVAR(2,2))

c	Adjust 'AMAX' if necessary.
	AMAX=MAXDIST
	IF(MODTYP==1.OR.MODTYP==3.OR.MODTYP==4)THEN
		AMAX=AMAX/3
	ENDIF
	AMAX=AMAX*1.5


	CALL RNSET(0)
		
c	Enter initial guesses of the variogram parameters.

C	Distance parameters.
	A=GUESSA

C	Nugget of the first variable's auto-variogram.
2003	C(0,1,1)=COVAR(1,1)*0.25

C	Nugget of the second variable's auto-variogram.
	C(0,2,2)=COVAR(2,2)*0.25

C	Nugget of the cross-variogram.
	C(0,1,2)=COVAR(1,2)*0.05

	IF(C(0,1,1)==0.OR.C(0,2,2)==0)THEN
		C(0,1,2)=0
	ENDIF
	C(0,2,1)=C(0,1,2)

	IF(MODTYP<4)THEN

		IF(MODTYP==3)THEN
			NSTR=2
		ELSE
			NSTR=1
		ENDIF

C		C1 of the first variable's auto-variogram.
		C(1,1,1)=COVAR(1,1)*0.75

C		C1 of the second variable's auto-variogram.
		C(1,2,2)=COVAR(2,2)*0.75

C		C1 of the cross-variogram
		C(1,1,2)=COVAR(1,2)*0.25

		IF(C(1,1,1)==0.OR.C(1,2,2)==0)THEN
			C(1,1,2)=0
		ENDIF
		C(1,2,1)=C(1,1,2)

	ENDIF

	IF(MODTYP>=4)THEN

		NSTR=2

C		C1 of the first variable's auto-variogram.
		C(1,1,1)=COVAR(1,1)*0.375

C		C1 of the second variable's auto-variogram.
		C(1,2,2)=COVAR(2,2)*0.375

C		C1 of the cross-variogram
		C(1,1,2)=COVAR(1,2)*0.1

		IF(C(1,1,1)==0.OR.C(1,2,2)==0)THEN
			C(1,1,2)=0
		ENDIF
		C(1,2,1)=C(1,1,2)

C		C2 of the first variable's auto-variogram.
		C(2,1,1)=COVAR(1,1)*0.375

C		C2 of the second variable's auto-variogram.
		C(2,2,2)=COVAR(2,2)*0.375

C		C2 of the cross-variogram.
		C(2,1,2)=COVAR(1,2)*0.25

		IF(C(2,1,1)==0.OR.C(2,2,2)==0)THEN
			C(2,1,2)=0
		ENDIF
		C(2,2,1)=C(2,1,2)

	ENDIF

c	Check guesses for positive-definiteness.
	do 13 istr=0,nstr
		call check (c,nstr,nvar,istr,icpo)

		if (icpo.lt.0) then
			print*,'Values for semivariogram are not' 
			print*,'positive-definite.'
			print*,'Try remaking semivariogram with' 
			print*,'different specifications.'
			stop
		endif

13	continue


c	Enter the annealing parameters.
c
c	These constrain how much the variogram parameters can change with each iteration.
	do 14 ix=1,nstr

		IF(MODTYP==1.OR.MODTYP==4)THEN

			AMIN=DATAF(1,1,1,2)/3

			IF(A(IX)<33)THEN
				AM(IX)=1
			ELSE
				AM(IX)=2
			ENDIF

		ENDIF
		IF(MODTYP==2.OR.MODTYP==5)THEN

			AMIN=DATAF(1,1,1,2)

			IF(A(IX)<100)THEN
				AM(IX)=1
			ELSE
				AM(IX)=2
			ENDIF

		ENDIF

14	continue
	IF(MODTYP==3)THEN

		AMIN=DATAF(1,1,1,2)/(3-((A(2)-1)*1.26795))

		IF(A(1)<(100/(3-((A(2)-1)*1.26795))))THEN
			AM(1)=1
		ELSE
			AM(1)=2
		ENDIF

		AM(2)=0.05

	ENDIF

	do 18 is=0,nstr
		do 19 ir=1,nvar
			do 20 ic=ir,nvar

				IF(IR==IC)THEN
					DM(IS,IR,IC)=COVAR(IR,IC)*0.25
				ELSE
				DM(IS,IR,IC)=0.25*(SQRT(COVAR(IR,IR))*SQRT(COVAR(IC,IC)))
				ENDIF

20			continue
19		continue
18	continue

c	Parameter 'cparf' is the initial temperature in the input file.  If no 
c	previous attempt at a solution has been made then 'cpar' is set to 'cparf'.
	if(ics.eq.0)then
		cpar=cparf
	endif
	alp=0.975
	nmarkov=60
	istop=50
	iwopt=WGT

c	Initiate the fitting criterion.
	call fcn(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt,f)
	print*, 'Weighted sum-of-squares of guessed parameters is:',f
	PRINT*,' '
	WRITE(13,*),'WSS of guessed parameters is:',f
	write(13,*) ' '
	write(13,*) 'Fitting parameters'
	write(13,*) ' '
	write(13,*) 'Initial temperature:',cpar
	write(13,*) 'Cooling parameter:',alp
	write(13,*) 'Number of trials in a Markov chain:',nmarkov
	write(13,*) 'Number of M. chains returning no change to stop:'	
     :,istop
55	continue
	write(13,*) ' '
	write(13,*) 'Weighting option:',iwopt
	write(13,*) ' '
	write(13,*) ' '
	

c	**************************************************************************************************
c
c		Start simulated annealing
c
c	**************************************************************************************************

C	'fic' and 'nunch' are used to count the number of chains for which 'f'
c	has not been changed.

	fic=-10.0
	nunch=0

c	Iterate cooling step.

	do 310 ics=1,2000

		if (nunch.gt.istop) then
			goto 2001
		endif

		rej=0
		acc=0

c		Iterate within the step.
		do 320 imc=1,nmarkov

c			Adjust each parameter in turn.

c			Distance parameter(s) first.
			IF(LOCK==0)THEN

				do 330 istr=1,nstr

					ao=a(istr)

1965					r=drnunf()
 					rc=(r-0.5)*2.0*am(istr)

c					Reject inappropriate values.
					IF(MODTYP<4)THEN

						IF(MODTYP==3.AND.ISTR==2)THEN
							IF((ao+rc)<=1.OR.(ao+rc)>=2)THEN
								goto 1965
							ENDIF
						ELSE
							IF((ao+rc)<=AMIN.OR.(ao+rc)>AMAX)THEN
								goto 1965
							ENDIF
						ENDIF

					ENDIF
					IF(MODTYP==4.OR.MODTYP==5)THEN
						IF(istr==1)THEN
							IF((ao+rc)<=AMIN.OR.(ao+rc)>=A(2))THEN
								goto 1965
							ENDIF
						ELSE
							IF((ao+rc)<=A(1).OR.(ao+rc)>AMAX)THEN
								goto 1965
							ENDIF
						ENDIF

					ENDIF

					a(istr)=ao+rc
					fo=f

					call fcn(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt,f)

					if(f.le.fo) then
						goto 1969
					endif

					call metrop(f,fo,cpar,accrej)

					if(accrej.lt.0) then
						rej=rej+1
						a(istr)=ao
						f=fo
						goto 330
					endif

1969					acc=acc+1
330				continue

			ENDIF


c			Now adjust variances.
			do 335 istr=0,nstr

				IF(MODTYP==3.AND.ISTR==2)GOTO 335

				do 340 ir=1,nvar
					do 350 ic=ir,nvar

						cold=c(istr,ir,ic)

2065						r=drnunf()
						rc=(R-0.5)*2.0*dm(istr,ir,ic)

c						Reject inappropriate values.
						IF(IR==IC)THEN
							if((cold+rc)<0)then
								goto 2065
							endif
						ELSE

							IF(ICVP==1.AND.(COLD+RC)<0)THEN
								GOTO 2065
							ENDIF
							IF(ICVP==-1.AND.(COLD+RC)>0)THEN
								GOTO 2065
							ENDIF

						ENDIF

c						Zero value if less than tolerance.
						IF(ABS(cold+rc)>=ABS(COVAR(IR,IC)/1000))THEN
							c(istr,ir,ic)=cold+rc
						ELSE
							c(istr,ir,ic)=0
						ENDIF
						c(istr,ic,ir)=c(istr,ir,ic)

c						Check for positive-definiteness.
						call check(c,nstr,nvar,istr,icpo)

						If (icpo.lt.0) then
							c(istr,ir,ic)=cold
							c(istr,ic,ir)=c(istr,ir,ic)
							goto 2065
						endif

						IF(LOCK==1)THEN
							fo=f
						ENDIF

				call fcn(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt,f)

						if (f.le.fo) then
							goto 2069
						endif

						call metrop(f,fo,cpar,accrej)

						if (accrej.lt.0) then
							rej=rej+1
							c(istr,ir,ic)=cold
							c(istr,ic,ir)=c(istr,ir,ic)
							f=fo
							goto 350
						endif
						
2069						acc=acc+1

350					continue
340				continue
335			continue			

			pacc=acc/(acc+rej)

			if (ics.eq.1) then
				paccin=pacc
			endif

320		continue

		if (fic.eq.f) then
			nunch=nunch+1
		else
			nunch=0
			fic=f
		endif

		if (ics.eq.1) then
		print*,'Proportion of changes accepted in 1st chain:',paccin
			PRINT*,' '
		PRINT*,'(Between 0.90-0.99 is optimal. If equal to one, reject'
		PRINT*,'and decrease the initial temperature. If less than 0.9,'
		PRINT*,'reject and increase the initial temperature).'
			print*,' '
			print*,'To accept press 0, to exit press 1:'
			read*,iconto

			IF(ICONTO==0)THEN
				PRINT*,' '
				PRINT*,'Minimising weighted sum-of-squares'
				PRINT*,'(please be patient)...'
				PRINT*,' '
			ENDIF
			IF(ICONTO==1)THEN
				STOP
			ENDIF

		endif
		
C		Cool things.
		cpar=cpar*alp

c		Parameter 'ics' is the Markov chain number, 'f' is the criterion minimised
c		and 'pacc' is how much it has changed since the last time.
		write(14,*)ics,f,pacc

310	continue

	
c	Write out results.
2001	print*,'Solution at:'
	write(13,*) 'Solution after ',ics-1,' Markov chains at:'
	
	Do 100 istr=0,nstr
		print*, 'Structure ',istr,' (0 denotes the nugget)'
		write(13,*) 'Structure ',istr,' (0 denotes the nugget)'

		if(istr.gt.0)then
			print*, 'Distance parameter:',a(istr)
			write(13,*) 'Distance parameter:',a(istr)
			A_OUT(ISTR)=DBLE(A(ISTR))
		endif

		do 110 ir=1,nvar
			do 120 ic=ir,nvar
				print*,ir,ic,c(istr,ir,ic)
				write(13,*)ir,ic,c(istr,ir,ic)
				C_OUT(istr,ir,ic)=DBLE(c(istr,ir,ic))
120			continue
110		continue
	
		IF(ISTR==2)EXIT

		print*, ' '
		print*, ' '
		write(13,*) ' '
		write(13,*) ' '

100	continue

	print*, ' '
	write(13,*) ' '

	print*,'Weighted sum-of-squares of final solution is:',f
	print*,'************************************************'
	print*,' '

	write(13,*)'Final WSS:',f


C	Overall goodness-of-fit according to Akaike Information Criterion
C	(Webster & McBratney, 1989).
	IF(MODTYP==1.OR.MODTYP==2)THEN
		vAIC=SUM(NLAGS)*LOG(F)+14
	ENDIF
	IF(MODTYP==3)THEN
		vAIC=SUM(NLAGS)*LOG(F)+16
	ENDIF
	IF(MODTYP==4.OR.MODTYP==5)THEN
		vAIC=SUM(NLAGS)*LOG(F)+22
	ENDIF
	
	WRITE(13,*)'Variable portion of Akaike Information Criterion:',vAIC
	WRITE(13,*)' '


C	Effective ranges (if applicable).
	IF(MODTYP==1.OR.MODTYP==3.OR.MODTYP==4)THEN

		IF(MODTYP==1)THEN
			WRITE(13,*)'Effective range (m) =',3*A(1)
		ENDIF
		IF(MODTYP==3)THEN
		WRITE(13,*)'Effective range (m) =',(3-((A(2)-1)*1.26795))*A(1)
		ENDIF
		IF(MODTYP==4)THEN
		WRITE(13,*)'Effective range of first structure (m) =',3*A(1)
		WRITE(13,*)'Effective range of second structure (m) =',3*A(2)
		ENDIF

		WRITE(13,*)' '
		
	ENDIF

	WRITE(13,*)'************************************************'	

	DATAF=0
	A=0
	C=0
	ICS=0

	RETURN

	END SUBROUTINE
c
c
c	**************************************************************************************************
c	**************************************************************************************************
c
c	Subroutines and Functions
c
c	**************************************************************************************************
c	**************************************************************************************************
	SUBROUTINE metrop(f,fo,c,accrej)

c	Metropolis criterion.
	REAL::R

	pracc=exp((fo-f)/c)

	CALL RNUN(1,R)
	ran=R
	if (ran.le.pracc) then
		accrej=1.0
	else
		accrej=-1.0
	endif
	return
	end subroutine

c	**************************************************************************************************
c	**************************************************************************************************
	SUBROUTINE check(c,nstr,nvar,istr,icpo)

c	Checks that the coregionalization matrix for structure 'istr' is positive-definite.

	REAL::c(0:2,nvar,nvar),cm(nvar,nvar),eval(nvar)

	do 10 ir=1,nvar
		do 20 ic=ir,nvar
			cm(ir,ic)=c(istr,ir,ic)
			cm(ic,ir)=c(istr,ic,ir)
20		continue
10	continue

c	Eigenvalues.
	call evlsf(nvar,cm,nvar,eval)

	rmin=10000.0
	do 30 ix=1,nvar
		if (rmin.gt.eval(ix)) then
			rmin=eval(ix)
		endif
30	continue

	if (rmin.le.-0.00001) then
		icpo=-1
	else
		icpo=1
	endif

	return
	end subroutine

c	**************************************************************************************************
c	**************************************************************************************************
	SUBROUTINE FCN(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt,f)

c	Returns 'f', which is the weighted sum of squares given LMCR parameters in 'a' and 'c'.

	integer::MODTYP,NLAGS(3),N,NSTR,NVAR,IWOPT,NLG
	real::dataf(NVAR,NVAR,3,500),sd(2),c(0:2,nvar,nvar),a(2),of,f,h
	DOUBLEPRECISION::C0,C1,A1,C2,A2,A3


	of=0.0
	rnpt=0.0

	do 10 ir=1,nvar
		do 20 ic=ir,nvar

			IF(IR==1.AND.IC==1)THEN
				NLG=NLAGS(1)
			ENDIF
			IF(IR==1.AND.IC==2)THEN
				NLG=NLAGS(2)
			ENDIF
			IF(IR==2.AND.IC==2)THEN
				NLG=NLAGS(3)
			ENDIF

			C0=c(0,ir,ic)
			C1=c(1,ir,ic)
			a1=a(1)

			if (nstr.gt.1) then
				c2=c(2,ir,ic)
				a2=a(2)
			else
				c2=0.0
				a2=0.0
			endif

			do 30 lag=1,NLG

				h=dataf(ir,ic,1,lag)
				gam=dataf(ir,ic,2,lag)
				rnp=dataf(ir,ic,3,lag)

c				Call 'GAMMAH' function.
				prgam=gammah(h,modtyp,C0,C1,a1,C2,a2)

c				Weighted sums-of-squares.
				if (ir.eq.ic) then
					wt=(sd(ir))**(-4.0)
				else
					wt=((sd(ir))**(-4.0))+((sd(ic))**(-4.0))
				endif
				wdevs=((gam-prgam)**2.0)*wt
			
				if (iwopt.eq.1) then
					wdevs=wdevs
				else if (iwopt.eq.2) then
					wdevs=wdevs*rnp
				else if (iwopt.eq.3) then
					wdevs=wdevs*(rnp/(prgam**2.0))
				else if (iwopt.eq.4) then
					wdevs=wdevs*(rnp/(h**2.0))
				endif

				of=of+wdevs
	
30			continue
			
20		continue
10	continue

	f=of

	return
	end subroutine