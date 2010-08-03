\name{lmcr}
\alias{lmcr}
\title{
Fit a linear model of coregionalisation to a gstat object.
}
\description{
Given a gstat model object with two variables and a gstat variogram with two variograms and a cross co variogram, this function will fit a lmcr, using the method of Lark and Papritz (2003).
}
\usage{
lmcr(g, v, covar, guessa, modtyp, cpar,...)

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

This function is an adaptation of fortran code, based on Lark and Papritz (2003). The orginal fortran code has been modified to fit a LMCR to a gstat cross-covariogram object (Please see the example for further details).

The function can use 10 different semi-variogram models;
\describe{
	\item{1 }{An isotropic linear model.}
	\item{2 }{An istopric circular model.}
	\item{3 }{An istropic spherical model.}
	\item{4 }{An istropic exponential model.}
	\item{5 }{An iso.pen}
	\item{6 }{An istropic double model.}
	\item{7 }{A double exponential model.}
	\item{8 }{A stable model.}
	\item{9 }{A double stable model.}
	\item{10 }{bro}
}
The optimisation of the model is based on a weighted sums of squares approach and can accomidate four different weighting methods (TODO: complete the wss weighting methods);
\describe{
	\item{1 }{Equal weighting}
	\item{2 }{wdevs=wdevs*rnp}
	\item{3 }{wdevs=wdevs*(rnp/(prgam^2.0))}
	\item{4 }{wdevs=wdevs*(rnp/(h^2.0))}
}
}
\value{

	\item{WSSGuessedParameter }{The intial WSS of the input parameters.}
	\item{initialTemperature }{The initial Cooling parameter used.}
	\item{coolingParameter }{}
	\item{numberTrialMarkovChain }{}
	\item{numberMarkovChainReturningNoChange }{}
	\item{weightingOption }{The weighting option set for the WSS.}
	\item{other.wss }{The WSS of each fitted model from the annealing process.}
	\item{other.pacc }{The amount of change from previous model.}
	\item{solutionAfter }{The amount of models fitted.}
	\item{structure }{}
	\item{distance }{The distance of the final LMCR.}
	\item{ir }{}
	\item{ic }{}
	\item{c }{The LMCR fitted nugget and psill (auto1,cross,auto2)}
	\item{finalWSS }{The WSS of the returned LMCR}
	\item{variogram }{A gstat cross-covariogram object of the initial input}
	\item{gstat }{A gstat formatted LMCR of the final model.}
	\item{variableAIC }{}
	\item{effectiveRange1stStructure }{The practical range of the LMCR}
	\item{effectiveRange2ndStructure }{Depending on the model specified The 2nd practical range of the LMCR.}

}
\references{
Lark, RM and Papritz, A (2003) Fitting a linear model of coregionalization for soil properties using simulated annealing, \bold{Geoderma},115, 245-260.
}
\author{
JS LESSELS
}
\note{
The function is not the most stabile, and care should be taken to insure the order of the parameters are provided in the correct manor.
}


