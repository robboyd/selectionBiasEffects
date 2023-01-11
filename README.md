# selectionBiasEffects
## Introduction

The sample mean is the most common estimator of the population mean. Where there is a selection bias, however, it is likely to be a biased one. In this document, I demonstrate this unfortunate property of the sample mean using an example from biodiversity science. 

## Calluna vulgaris example 

In the example, the estimand is mean occupancy {1,0} of the plant Calluna vulgaris in Britain. The population comprises 229,772 1 km grid cells, each of which being occupied by C. vulgaris (1) or not (0). True occupancy was estimated from the union of the UKCEH Land Cover Map [1], 1 km Atlas data (2000-2019; [2]) and the UK Countryside Survey 2007 plot data [3], constrained using the plant's 10 km distribution 2000-2019 [2]. 19,419 1 km grid cells were sampled by colunteers, who submitted records of vascular plants to iRecord (www.irecord.org.uk) via the smartphone app, the website, or indirectly via the iSpot initiative (https://www.ispotnature.org/) from 2000-2019. 

The data, in a slimmed down form, can be accessed via the R package selectionBiasEffects. In this form, the data consists of two vectors: R, which is the sample inclusion indicator (1 if the grid cell is in the sample and 0 otherwise); and Y, which is the value of true occupancy. Each row corresponds to one 1 km grid cell.

```{r}
# install selectionBiasEffects
#devtools::install_github("https://github.com/robboyd/selectionBiasEffects")

library(ggplot2)
library(reshape2)
library(selectionBiasEffects)

# load calluna data 
data("Calluna")

head(Calluna)

```

We can use the function howBadIsIt in selectionBiasEffects to evaluate the performance of the sample mean as an estimator of the population mean. Here, d1 is the correlation between sample membership and the true occupancy [4]. N is the population size and n is the sample size. ybar_sample and ybar_pop are the sample and population means, respectively. nEff is the effective sample size, defined as the size of the simple random sample (SRS) whose estimate the population mean would have the same Mean Squared Error (MSE) on average [4]. The other elements of calStats needn't be defined here; type ?howBadIsIt into the R console for details.

```{r}

calStats <- howBadIsIt(Calluna)

calStats

```

Whilst a SRS of size nEff and biased sample of size n might have the same MSE on average, they are not equally good for inferential purposes. To see this, we can simulate 100 realisations of each, using the function simulateSamples in selectionBiasEffects. This biased samples are simulated to have the same correlation between R and Y as the empirical sample. simulateSamples returns various statistics from each sample, including an estimate of the population mean and its error. 

```{r}

stats <- simulateSamples(dat = Calluna,
                         N = calStats$N,
                         n = CalStats$n,
                         n_eff = calStats$nEff,
                         ddc = calStats$d1,
                         truth = calStats$ybar_pop,
                         nSim = 100)

head(stats)

```

Plotting the distributions of these statistics provides insight into the differences between the SRSs of size nEff and the biased samples of size n. The small SRSs estimate the population mean accurately on average, whereas the larger biased samples always underestimate it. Nevertheless, the samples have similar MSEs, as shown in panel B, and these would converge on the same value if we were to increase the number of simulated samples. 

```{r}
# format the data for plotting
stats <- stats[,-c(3,4,5,6)]

stats2 <- melt(stats, id.vars = c("sim"))

stats2$samp <- ifelse(startsWith(as.character(stats2$variable), "small"), "small SRS", "large biased")

stats2$stat <- NA

stats2$stat[endsWith(as.character(stats2$variable), "Estimate")] <- "A) estimate"

stats2$stat[endsWith(as.character(stats2$variable), "SquaredError")] <- "B) squared error"

head(stats2)

stats2$int <- ifelse(stats2$stat == "A) estimate", 0.3, NA)

## add MSEs 
m1 <- mean(stats$smallSRSSquaredError)

m2 <- mean(stats$largeBiasedSquaredError)

stats2$int2 <- ifelse(stats2$stat == "B) squared error", m1, NA)

stats2$int3 <- ifelse(stats2$stat == "B) squared error", m2, NA)

# plot
ggplot(data = stats2, aes(x = value, fill = samp), colour = c("red", "blue")) +
  geom_density(colour = NA, alpha = 0.5) + 
  labs(x = "",
       fill = "") +
  theme_minimal() +
  facet_wrap(~stat, 
             scales = "free",
             nrow = 2) +
  geom_vline(aes(xintercept = int), linetype = "dotdash") +
  geom_vline(aes(xintercept = int2), linetype = "dotdash", colour = "blue") +
  geom_vline(aes(xintercept = int3), linetype = "dotdash", colour = "red") +
  scale_fill_manual(values = c("red", "blue"))


```

## References

1. Morton, D. et al. (2014) Land Cover Map 2007 (1km percentage target class, GB) v1.2.11112
2. Stroh, P.A. et al. (2023) Plant Atlas 2020: Mapping Changes in the Distribution of the British and Irish Flora, Princeton Univ. Press
3. Carey, P. et al. (2008) Countryside Survey: UK Results from 2007 CHAPTER 1 Countryside Survey Methodology Secretary of State for Environment, Food and Rural Aff airs1-119
4. Meng, X.L. (2018) Statistical paradises and paradoxes in big data (I): Law of large populations, big data paradox, and the 2016 us presidential election. Ann. Appl. Stat. 12, 685-726
