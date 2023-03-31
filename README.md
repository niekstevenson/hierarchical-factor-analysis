# hierarchical-factor-analysis
A bayesian implementation and sampler to estimate hierarchical factor analysis developped by Reilly Innes and Niek Stevenson.

In 'toy-example.R' we showcase how to implement a joint hierarchical factor model and how to perform posterior-inference.
In 'test-identifiability.R' we showcase how a researcher can specify a design equal in principle to their own design to test whether their design can be recovered. 
'plots.R' contains some function that can be used to plot the factor loadings.

Furthermore the other files show how to implement the other group-level distributions, while sampling for the individual level using the likelihood defined by the user.
1. Diagonal: No relations estimated between parameters at the group-level.
2. Blocked: Only relations estimated between blocks of the parameter space.
3. Single Subject: No group level parameters estimated, just a sampler for bayesian inference
4. Standard: Covariances estimated between all parameters included in the model

I'm happy to help no strings attached with anyone interested in applying this approach to their own models, just shoot me a question at:
niek.stevenson@gmail.com
