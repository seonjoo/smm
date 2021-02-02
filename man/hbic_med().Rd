\name{hbic_med()}
\alias{hbic_med()}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to do HBIC computation for multiple mediations model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hbic_med(fit, fit.n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{x} here~~
Fitted model uses sparse.mediation.sgrplasso.largep_omega() from the smm package.
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
hbic_med = hbic_med(fit, fit.n)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
