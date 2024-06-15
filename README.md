# SCR-Classify-Pairwise
MCMC samplers for latent identity spatial capture-recapture models that use pairwise scores between samples

Minimal information for now, but there are two models here that are latent ID SCR models where we observe pairwise match scores between all samples. 

1. Scores are similar to what you get from the SURF or SIFT algorithms as implemented in Hotspotter, etc. A toy example of a model for this type of data is "pairwise Poisson", where the correct and incorrect match score distributions are assumed to have Poisson distributions with different parameters. The correct match lambda should be larger than the incorrect match lambda for SURF/SIFT data that count the matching key points between photos. Currently, I assume only 1 pairwise score between samples, but multiple classifiers could be used.

2. Scores are categorical with multiple classifiers. Here, the categories could be human-assigned correct-incorrect, or correct-uncertain-incorrect. Or they could be something else entirely. This model is set up for multiple classifiers, e.g., multiple humans assigning pairwise scores between photos.

These models will be sensitive to correctly specifying the match score distributions, which will not be something nice like a Poisson for SURF/SIFT-type data or most data sources, probably. Those could be converted to categorical scores, possibly. Also, known-ID validation data with pairwise scores could be useful in determining the most appropriate score distribution. This is uncharted territory.

Finally, I will direct everyone to Amanda Ellis's PhD dissertation on this topic in a nonspatial context. See here:

https://uknowledge.uky.edu/statistics_etds/31/