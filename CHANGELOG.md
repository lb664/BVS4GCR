
# News/Changelog

Changes in **BVS4GCR** package version **0.2.1** (25/10/2026):
<ul>
    <li> Fixed the assignment of the beta-binomial prior hyperparameter in the Bayesian Variable Selection step ('bPi <- hyperPar$aPi' is now 'bPi <- hyperPar$bPi'). The prior on the predictor inclusion probability now has the mean implied by 'EVgamma', instead of being symmetric
    <li> Sim_GCRdata() now uses its 'seed' argument ('set.seed(seed)'): previously the seed was only recorded in the output and the simulated data were different at every call
    <li> The argument 'tau' of BVS4GCR() is now used as the prior variance of the regression coefficients. Previously it was ignored and the prior variance was always set at 1. When 'std = TRUE' (the predictors are standardised) 'tau' is still set at 1, whatever value is supplied
    <li> Fixed the posterior mean of the standard deviations of the responses, 'postMean$D', which was computed from the last Markov chain Monte Carlo iteration only (and divided by the number of responses)
    <li> Fixed the size of the stored output when '(niter - burnin) / thin' is not an integer, which stopped the analysis with a "subscript out of bounds" error (for example 'niter = 20, burnin = 5, thin = 2'). The output is now sized from the exact set of saved iterations
    <li> The initial set of predictors included in the model now uses the p-values of the one at-a-time regressions, as intended. Previously the wrong field name was used and no predictor was included initially
    <li> The p-value used to build the proposal distribution of an ordinal response is now based on the t-value of the predictor (previously on the regression coefficient). Each binomial response now uses its own number of trials 'nTrial' (previously the first value was used for all binomial responses)
    <li> For binomial and negative binomial responses, the inclusion indicators are now reset when the proposed model is empty, as for the other responses. Missing values of Gaussian responses are now initialised with independent random draws
    <li> BVS4GCR() stops with an informative message when an ordinal response has fewer than three categories (a binary response must be specified as 'binary')
    <li> The internal functions are now named with a leading dot ('.BVS', '.Sample_G', '.Sample_gamma', '.Sim_HIW', '.mod_BDgraph', '.mod_BDgraphInit') and all their arguments are documented. A package help page ('BVS4GCR-package') has been added
    <li> Corrected the documentation: default values of 'EVgamma', 'probAdd' and 'std', the elements of the output of BVS4GCR() and Sim_GCRdata() (including the missing 'Sigma' and the 'NBPar' entry), a Greek letter used in place of 'B' that prevented the manual from being built, and the definition of the marginal posterior inclusion probability (mPIP) and posterior edge inclusion probability (PEIP)
    <li> The README has been updated (installation with the 'pak' package, references, definitions of mPIP and PEIP). Added the author's ORCID, 'URL', 'BugReports' and 'Language' to DESCRIPTION, declared the import of 'stats' and removed the obsolete C++11 request from 'src/Makevars'
</ul>
