# VAF fdr

As statistic we define the $VAF_{LR}$ or VAF log-ratio. Formally defined as

$$
VAF_{LR} = | log2 \bigg( \frac{VAF_{i}}{VAF_{j}} \bigg) |
$$

where $VAF_{j}$ and $VAF_{i}$ are 2 VAF values corresponing to a pair of mutations.

## True Positive Set

The true positive set consist in the $VAF_{LR}$ of the clustered positions if $VAF_{i} = VAF_{j}$.
We define a true positive set as the $VAF_{LR}$ values of the actual vaf value in each clustered mutation and the corresponding simulated vaf ($VAF'$).
This artifial VAF value tries to simulate the noise expected in the detection of a mutation that will decrease as the covrage increases.
We define the simulated VAF as

$$
VAF' = \frac{sr \cdot ar + (1-sr) ar'}{tr}
$$

where $sr$ is the number of shared reads between 2 positions.
This correction avoids the overestimation of noise when the pair of mutations are very close and therefore they share a proportion of reads where the difference of VAF values have to be 0.
The proportion is estimated using

$$
sr = \frac{50 - x}{50 + x}
$$

for $sr \subset [0,1]$  and 50 representing the read length.

Then, $AR'$ is the estimated number of reads using the original number of alternate reads and the vaf as parameters for a Binomial distribution.
For each mutation only one vaf value is sampled.
This is because we expect deviations in a single value will be corrected by other cases.
In this case, the distribution is the overall meaningfull element.

$$
ar' \sim \operatorname{B} \left({tr, vaf}\right)
$$