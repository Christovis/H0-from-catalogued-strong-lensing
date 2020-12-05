# Inferring H0 with Catalogued Strong Lensing Time-delays

Code repository for the paper **The impact of line-of-sight structures on measuring the $$H_0$$ with Strong Lensing Time-delays**
by Nan Li, Christoph Becker, and Simon Dye


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-2006.08540%20-green.svg)](https://arxiv.org/abs/2006.08540)

## Abstract

Measurements of The Hubble-Lemaitre constant from early- and local-universe observations show a significant discrepancy.  In an attempt to understand the origin of this mismatch, independent techniques to measure $$H_0$$ are required. One such technique, strong lensing time delays, is set to become a leading contender amongst the myriad methods due to forthcoming large strong lens samples. It is therefore critical to understand the systematic effects inherent in this method. In this paper, we quantify the influence of additional structures along the line-of-sight by adopting realistic lightcones derived from the _CosmoDC2_ semi-analytical extra-galactic catalogue. Using multiple lens plane ray-tracing to create a set of simulated strong lensing systems, we have investigated the impact of line-of-sight structures on time-delay measurements and in turn, on the inferred value of $$H_0$$. We have also tested the reliability of existing procedures for correcting for line-of-sight effects. We find that if the integrated contribution of the of line-of-sight structures is close to a uniform mass sheet, the bias in $$H_0$$ can be adequately corrected by including a constant external convergence $$\kappa_{\rm ext}$$ in the lens model. However, for realistic line-of-sight structures comprising many galaxies at different redshifts, this simple correction over-estimates the bias by a factor of approximately three. We therefore conclude that lens modelling must incorporate multiple lens planes to account for line-of-sight structures for accurate and precise inference of $$H_0$$.

## Code

For the inferrence of $$H_0$$ three different software were evaluated: [_glafic_](https://ascl.net/1010.012), [_Gravlens_](http://ascl.net/1102.003), [_lenstronomy_](https://github.com/sibirrer/lenstronomy)

As _lenstronomy_ is the only open-source software and allowed for the most robust optimization procedure, it was the code on which the published results are based on.
