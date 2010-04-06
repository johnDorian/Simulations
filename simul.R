##First part takes ~8 secs
start<-Sys.time()
x=170000
y1_0=rnorm(x,0,1)
y1_1=grf(x,cov.pars=c(1,13.59084),mean=0,nug=0)$data

y2_0=rnorm(x,0,1)
y2_1=grf(x,cov.pars=c(1,13.59084),mean=0,nug=0)$data
Sys.time()-start



#Nugget
a11_0
a12_0
#sill
a11_1
a12_1
#Nugget
a21_0
a22_0
#sill
a21_1
a22_1
#means
mu1
mu2

Z1=a11_0*y1_0+a12_0*y2_0+a11_1*y1_1+a121*y1_2+mu1
Z2=a21_0*y1_0+a22_0*y2_0+a21_1*y1_1+a22_1*y1_2+mu2
