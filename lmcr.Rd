\name{lmcr}
\alias{lmcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a linear model of coregionalisation to a gstat object.
}
\description{
Given a gstat model object with two variables and a gstat variogram with two variograms and a cross co variogram, this function will fit a lmcr.
}
\usage{
lmcr(g, v, covar, guessa, modtyp, cpar)

lmcr(g, v, covar, guessa, modtyp, cpar, 
	wgt=1, icvp=1, lock=0, istop = 50, 
	plot.wss.change = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
A gstat object with no model parameters and two variables.
}
  \item{v}{
A gstat cross variogram with with the same two variables and the cross variogram
}
  \item{covar}{
The nugget of each covariogram. In the order, auto.1,cross,auto.2. Must be three in length
}
  \item{guessa}{
Initial guesses of the distance of the variogram model. If modtyp is <6 only one number is required, else two variables must be specified.
}
  \item{modtyp}{
An integer between 1 and 8 indicating which type of model is to be used. See details for more information.
}
  \item{cpar}{
Cooling paparameter for the simulating annealing process. See Lark's 2003 paper for more details.
}
  \item{wgt}{
The weighted sum of squares method to use. Defualt is set to 1. Equal wieghting. See details below for more information.
}
  \item{icvp}{
Should covariance matrix be positive definitive. 1 indicates yes. defualt is 1.
}
  \item{lock}{
If the distance parameter should be fitted. 0 = fit range, 1 = accept guessa for distance. Default = 0.
}
  \item{istop}{
The length of each markov chain. See Lark's 2003 paper for more details. Defualt is 50.
}
  \item{plot.wss.change}{
Wether a plot of each iterative weighted sum of squares should be given. Default is TRUE
}
}
\details{
Need to talk about modtyp, cpar and wgt, istop, and mention plotting function.
This function is an adaptation of fortran code, based on Lark and Papritz (2003). This function has been modified to use a cross co-variogram of two variables from the gstat package, and a gstat model.
The function can use 8 different models.
\arguments{
	\item{1}
An isotropic linear model.
	\item{2}
iso.cir
iso.sph
iso.exp
iso.pen
iso.dou
dou.exp
stable
dou.sta
bro
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
