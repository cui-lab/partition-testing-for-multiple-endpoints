\name{clrt}
\alias{clrt}
\title{Consonant Likelihood Ratio Test for Primary and Secondary Endpoints}
\description{
Calculate critical constants and perform consonant likelihood ratio test for primary and secondary endpoints problem.
}
\usage{
clrt(alpha=NULL,cv=NULL,n=NULL,grpmean=NULL,grpstdev=NULL,distr.null=NULL,
     distr.null.size=NULL,distr.alter=NULL,distr.alter.size=NULL,rho=0,
     delta=c(0,0),scenario=NULL,savenull=FALSE,savealter=FALSE,
     showpvalues=FALSE,showLRT=FALSE)
}
\arguments{
\item{alpha}{Significance level for hypothesis testing and critical constants calculation, between 0 and 1. Alpha needs to be specified unless the user only interested in generating Monte-Carlo samples for the null or alternative distirbution.}
\item{cv}{A list of critical constants saved from previous computation. If missing, new values will be calculated based on provided significance level.}
\item{n}{Vector of sample sizes of each treatment group. The last element should be the sample size of the control group.}
\item{grpmean}{A matrix with two rows of effect size for each treatment group. The first row are effect sizes for primary endpoint and the second row are effect sizes for the secondary endpoint. The effect size of the control group is in the last column.}
\item{grpstdev}{A vector of length two, the common standard deviation of each endpoint.}
\item{distr.null}{Samples of null distribution previously generated and saved by clrt.}
\item{distr.null.size}{Number of Monte-Carlo samples generated for the null distribution.}
\item{distr.alter}{Samples of alternative distribution previously generated and saved by clrt.}
\item{distr.alter.size}{Number of Monte-Carlo samples generated for the alternative distribution.}
\item{rho}{The known correlation coefficient between primary and secondary endpoints that is used for the alternative distribution. Default value is 0.}
\item{delta}{A vector of length two, provide the effecacy margin specified in the hypotheses for primary and secondary endpoints. Default value is 0.}
\item{scenario}{A matrix with two rows and number of columns equals the number of dose levels compared with the control. The matrix provide true effect size for each dose-endpoint combination for the alternative distribution.}
\item{savenull}{A logical variable indicating if the null distribution sample will be saved. Saving the null distribution samples let the user to compute critical constants of different significance level in the future. In order to return p-values when using 'clrt' for testing, the null distribution samples are required to safe time. Users need to make sure enough space is available to save the file.}
\item{savealter}{A logical variable indicating if the alternative distribution sample will be saved. Saving the alternative distribution samples let the user to run the testing procedure and calculate power under different significance level in the future. Users need to make sure enough space is available to save the file.}
\item{showpvalues}{A logical variable indicating if the pvalues should be returned when perform the hypotheses testing. Returning the pvalues requires to have Monte-Carlo samples form the null distribution. If no such samples are provided from previously saved files, new samples will be generated which will require extra time to run the program.}
\item{showLRT}{A logical variable indicating if the likelihood ratio statistics should be returned for each partition hypothesis.}
}
\details{
This function can generate Monte-Carlo samples for null distribution or alternative distribution for the hypotheses testing problems involving multiple doses and primary, secondary endpoints. To generate null distribution samples, the parameters of n, distr.null.size are required and savenull must be set to TRUE. To generate alternative distribution, the parameters of n, distr.alter.size, delta and rho are required and savealter must be TRUE.

This function can also compute the critical constants for consonant likelihood ratio test. Alpha is required to be specified in this case and either distr.null or distr.null.size should be provided. If distr.null is available from previously saved calculations, the computation of critical constants will be very simple. Otherwise, if distr.null.size is provided, new Monte-Carlo samples of the null distribution will be generated and the critical constants will be calculated afterwards. The user can also specify at this time if the newly generated null distribution samples should be saved by giving the savenull parameter.

This function can perform hypotheses testing for multiple dose and primary-secondary endpoints problems using the consonant likelihood ratio test. The user can either provide the mean and variance from a real problem or provide a group of means and variances generated from a simulation study. In the first case, one testing result will be returned. In the latter case, the testing results for each sample is obtained and the power of all the results are provided. To provide the mean and variance of one sample, the user need to specify grpmean and grpstdev. To provide samples from an alternative distribution, the user need to either borrow result from a previously saved samples of alternative distribution with distr.alter, or give the number of samples need to be generated by specifying distr.alter.size. To perform testing, alpha is required to be specified in advance. The cv (critical constants) saved from previous computation must match the currently specified alpha. Otherwise, the null distribution samples will be generated again. In order to return pvalues from the consonant likelihood ratio test, the null ditribution samples must be provided. Otherwise, it will be generated before the testing procedure can start. It is recommended to always save the null distribution samples for future computations if enough space is available.
}

\value{
\item{alpha}{Significance level for hypothesis testing and critical constants calculation.}
\item{cv}{A list of critical constants for the consonant likelihood ratio test.}
\item{distr.null}{Monte-Carlo samples from the null distribution.}
\item{distr.alter}{Monte-Carlo samples from the alternative distribution.}
\item{Power.Decision}{The testing result (pvalue) or the power of consonant likelihood ratio test.}
\item{FWER}{The family-wise error rate from simulation study using the provided alternative distribution samples.}
\item{Partition.power}{The consonant likelihood ratio test power for each partition hypothesis.}
\item{Partition.pvalue}{The pvalue for each partition hypothesis, when returning the pvalues is requested.}
\item{Partition.rejection}{The testing result for each partition hypothesis}
\item{Hypotheses}{The indices of the partition hypotheses, refer to Wang and Cui (2011).}
\item{ComponentLRT}{The likelihood ratio statistic for each partition hypothesis.}
}
\references{
   Wang, B. and Cui, X. (2011)
}
\examples{
## Generate samples for null distribution
null.sample<-clrt(n=c(30,30,30,30),distr.null.size=10^4,savenull=TRUE)
## Generate samples for alternative distribution
alter.sample<-clrt(n=c(30,30,30,30),distr.alter.size=10^3,delta=c(0,0),
rho=0.3,scenario=matrix(c(0.5,0.7,0.9,0.6,0.6,0.6),nrow=2,ncol=3,
byrow=TRUE),savealter=TRUE)
## Compute critical constants for alpha=0.05
cv005<-clrt(alpha=0.05,n=c(30,30,30,30),distr.null=null.sample$distr.null)
## Perform hypotheses testing using saved null distribution samples.
clrt(alpha=0.025,n=c(30,30,30,30),distr.null=null.sample$distr.null,
distr.alter=alter.sample$distr.alter,showLRT=FALSE)
}
\keyword{file}
