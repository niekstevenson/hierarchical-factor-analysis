# hierarchical-factor-analysis
A bayesian implementation and sampler to estimate hierarchical factor analysis developped by Reilly Innes and Niek Stevenson

The example files implement the following group-level distributions, whilst also sampling for the individual level using the likelihood defined by the user.
1. Diagonal: No relations estimated between parameters at the group-level.
2. Blocked: Only relations estimated between blocks of the parameter space.
3. Single Subject: No group level parameters estimated, just a sampler for bayesian inference
4. Standard: Covariances estimated between all parameters included in the model
5. Factor analysis: Using a factor analysis decomposition of the covariance matrix
6. Multitasks: Similar to standard, but then with an example how to include multiple models in the group level.

I'm happy to help no strings attached with anyone interested in applying this approach to their own models, just shoot me a question at:
niek.stevenson@gmail.com
