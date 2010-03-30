\name{lmcr}
\alias{lmcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a linear model of coregionalisation to a gstat object using Lark's simulated annealing method.
}
\description{
Given a gstat model object with two variables and a gstat variogram with two variograms and a cross co variogram, this function will fit a lmcr.
}
\usage{
lmcr(g, v, wgt=1, icvp=1, cpar, modtyp, covar, maxdist, guessa, lock, istop = 50, plot.wss.change = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
A gstat object with no model parameters and two variables.
}
  \item{v}{
A gstat cross variogram with with the same two variables and the cross variogram
}
  \item{wgt}{
The weighted sum of squares method to use. Defualt is set to 1. Equal wieghting. See details below for more information.
}
  \item{icvp}{
Should covariance matrix be positive definitive. 1 indicates yes. defualt is 1.
}
  \item{cpar}{
Cooling paparameter for the simulating annealing process. See Lark's 2003 paper for more details.
}
  \item{modtyp}{
An integer between 1 and 8 indicating which type of model is to be used. See details for more information.
}
  \item{covar}{
A 2*2 matrix with the initial nugget values. See details for more information.
}
  \item{maxdist}{
The maximum distance of the variogram model. This is the same for each model.
}
  \item{guessa}{
Initial guesses of the distance of the variogram model. If modtyp is <6 only one number is required, else two variables must be specified.
}
  \item{lock}{
If the distance parameter should be fitted. 0 = fit range.
}
  \item{istop}{
The length of each markov chain. See Lark's 2003 paper for more details. Defualt is 50.
}
  \item{plot.wss.change}{
Wether a plot of each iterative weighted sum of squares should be given. Default is TRUE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
