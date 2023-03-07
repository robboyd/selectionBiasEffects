library(raster)
library(ggplot2)
library(survey)
library(rstanarm)
library(PracTools)
library(reshape2)

## load population data
pop <- read.csv("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/NERC_exploring_frontiers/Data/all_data.csv")

pop <- pop[complete.cases(pop), ]

## transform auxiliary data

#pop$postcode_density_299_neighbours[pop$postcode_density_299_neighbours == 0] <- 0.1

#pop$postcode_density_299_neighbours <- log(pop$postcode_density_299_neighbours)

#pop$UKelv[pop$UKelv == 0 | pop$UKelv < 0] <- 0.1

#pop$UKelv <- log(pop$UKelv)

## pull out the auxiliary data for the whole population
pop_aux_cont <- pop[,5:9]

## for reasons that will become clear later, we need to discretize the auxiliary data in pop
for (i in c(5,8,9)) {

  q <- as.numeric(quantile(pop[,i], probs = c(0, 0.33, 0.66, 1)))

  print(q)

  pop[,i] <- cut(pop[,i],
                 breaks = q,
                 labels = FALSE,
                 include.lowest = TRUE,
                 right = TRUE)

  pop[,i] <- as.numeric(pop[,i])

}

## we'll make PA and open access land coverage binary

pop$openAccessGB <- ifelse(pop$openAccessGB > 0, 1, 0)

pop$allPACoverage <- ifelse(pop$allPACoverage > 0, 1, 0)

## pull out columns relevant to period 1
pop_p1 <- pop[,c(1,3,5,6,7,8,9,10)]

## and period 2
pop_p2 <- pop[,c(2,4,5,6,7,8,9,11)]

## pull out sampled rows for periods 1 and 2
samp_p1 <- pop_p1[pop_p1$sampled_units_1987.1999 == 1, ]

samp_p2 <- pop_p2[pop_p2$sampled_units_2010.2019 == 1, ]

## pull out the auxiliary data for the whole population
pop_aux <- pop[,5:9]

## calculate the population means
pop_mean_p1 <- mean(pop_p1$heather_true_dist_1987.1999)

pop_mean_p2 <- mean(pop_p2$heather_true_dist_2010.2019)

## now set up survey designs using the survey package to construct estimators
design_p1 <- svydesign(ids=~0,
                       data = samp_p1)

design_p2 <- svydesign(ids=~0,
                       data = samp_p2)

## calculate sample means
samp_mean_p1 <- svymean(design = design_p1,
                        x=~heather_true_dist_1987.1999)

samp_mean_p2 <- svymean(design = design_p2,
                        x=~heather_true_dist_2010.2019)

## and their confidence intervals
samp_mean_p1_conf <- confint(object = samp_mean_p1,
                             level = 0.95)

samp_mean_p2_conf <- confint(object = samp_mean_p2,
                             level = 0.95)

## now calculate a weighted mean
## first create a new survey design with the estimated inclusion probabilities
weighted_design_p1 <- svydesign(ids=~0,
                                data = samp_p1,
                                probs=~inclusionProbs_1987.1999)

weighted_design_p2 <- svydesign(ids=~0,
                                data = samp_p2,
                                probs=~inclusionProbs_2010.2019)

## then get the weighted sample means
weighted_samp_mean_p1 <- svymean(design = weighted_design_p1,
                                 x=~heather_true_dist_1987.1999)

weighted_samp_mean_p2 <- svymean(design = weighted_design_p2,
                                 x=~heather_true_dist_2010.2019)

## and their confidence intervals
weighted_samp_mean_p1_conf <- confint(object = weighted_samp_mean_p1,
                                      level = 0.95)

weighted_samp_mean_p2_conf <- confint(object = weighted_samp_mean_p2,
                                      level = 0.95)

## next we want to postratify. We use the auxiliary data from earlier
## first, cross the covariates to get the poststrata
cells <- data.frame(table(pop_aux))

## now poststratify
ps_design_p1 <- postStratify(design = design_p1,
                             strata = samp_p1[,3:7],
                             population = data.frame(table(pop_aux)),
                             partial = T)

ps_design_p2 <- postStratify(design = design_p2,
                             strata = samp_p2[,3:7],
                             population = data.frame(table(pop_aux)),
                             partial = T)

## and get the weighted mean across poststrata
ps_samp_mean_p1 <- svymean(design = ps_design_p1,
                           x=~heather_true_dist_1987.1999,
                           na.rm = T)

ps_samp_mean_p2 <- svymean(design = ps_design_p2,
                                 x=~heather_true_dist_2010.2019)


## and their confidence intervals
ps_samp_mean_p1_conf <- confint(object = ps_samp_mean_p1,
                                level = 0.95)

ps_samp_mean_p2_conf <- confint(object = ps_samp_mean_p2,
                                level = 0.95)

## another estimator is regression-based
## first, calculate sums of the auxiliary variables

aux_tots <- c(nrow(pop_aux_cont), colSums(pop_aux_cont))

names(aux_tots)[1] <- "(Intercept)"

## now calibrate
calib_design_p1 <- calibrate(design = design_p1,
                             formula = ~ postcode_density_299_neighbours +
                               openAccessGB + allPACoverage +
                               road_length_299_neighbours + UKelv,
                             population = aux_tots,
                             calfun="linear")

sum(weights(calib_design_p1))

summary(weights(calib_design_p1))

sp_samp_mean_p1 <- svymean(~heather_true_dist_1987.1999, design=calib_design_p1)

calib_design_p2 <- calibrate(design = design_p2,
                             formula = ~ postcode_density_299_neighbours +
                               openAccessGB + allPACoverage +
                               road_length_299_neighbours + UKelv,
                             population = aux_tots,
                             calfun="linear")

sum(weights(calib_design_p2))

summary(weights(calib_design_p2))

sp_samp_mean_p2 <- svymean(~heather_true_dist_2010.2019, design=calib_design_p2)

## and the confidence intervals
sp_samp_mean_p1_conf <- confint(object = sp_samp_mean_p1,
                                level = 0.95)

sp_samp_mean_p2_conf <- confint(object = sp_samp_mean_p2,
                                level = 0.95)

## next is multiple regression and poststratification
#fit <- stan_glmer(
# heather_true_dist_1987.1999 ~ 1 + (1 | postcode_density_299_neighbours) +
#    (1 | road_length_299_neighbours) + (1 | allPACoverage) +
#    (1 | openAccessGB) + (1 | UKelv),
#  family = binomial(link = "logit"),
#  data = samp_p1,
#  chains = 2
#)

#print(fit)

#posterior_prob <- posterior_linpred(fit, transform = TRUE, newdata = data.frame(table(pop_aux)))

#poststrat_prob <- posterior_prob %*% cells$Freq / sum(cells$Freq)

poststrat_prob <- read.csv("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/NERC_exploring_frontiers/Data/poststrat_prob.csv")[,1]

#write.csv(poststrat_prob,
#          "W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/NERC_exploring_frontiers/Data/poststrat_prob.csv",
#         row.names = F)

model_popn_pref <- c(mean = mean(poststrat_prob),
                     lower = quantile(poststrat_prob, probs = 0.025),
                     upper = quantile(poststrat_prob, probs = 0.975))

#round(model_popn_pref, 3)

#fit2 <- stan_glmer(
#  heather_true_dist_2010.2019 ~ 1 + (1 | postcode_density_299_neighbours) +
#    (1 | road_length_299_neighbours) + (1 | allPACoverage) +
#    (1 | openAccessGB) + (1 | UKelv),
#  family = binomial(link = "logit"),
#  data = samp_p2[sample(1:nrow(samp_p2), 30000), ],
#  chains = 2
#)

#print(fit2)

#posterior_prob_p2 <- posterior_linpred(fit2, transform = TRUE, newdata = cells)

#poststrat_prob_p2 <- posterior_prob_p2 %*% cells$Freq / sum(cells$Freq)

poststrat_prob_p2 <- read.csv("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/NERC_exploring_frontiers/Data/poststrat_prob_p2.csv")[,1]

#write.csv(poststrat_prob_p2,
#          "W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/NERC_exploring_frontiers/Data/poststrat_prob_p2.csv",
#          row.names = F)

model_popn_pref_p2 <- c(mean = mean(poststrat_prob_p2),
                        lower = quantile(poststrat_prob_p2, probs = 0.025),
                        upper = quantile(poststrat_prob_p2, probs = 0.975))

#round(model_popn_pref_p2, 3)

#plotDat <- read.csv("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/NERC_exploring_frontiers/Data/plotDat.csv")

## let's plot some of the estimates so far
plotDat <- data.frame(p = c(1,2,1,2,1,2,1,2,1,2,1,2),
                      est = c(pop_mean_p1[1], pop_mean_p2[1], samp_mean_p1[1], samp_mean_p2[1], ps_samp_mean_p1[1], ps_samp_mean_p2[1], weighted_samp_mean_p1[1], weighted_samp_mean_p2[1], sp_samp_mean_p1[1], sp_samp_mean_p2[1], model_popn_pref[1],model_popn_pref_p2[1]),
                      type = c("Population", "Population", "Sample", "Sample", "Poststratified", "Poststratified", "Weighted", "Weighted", "Superpopulation", "Superpopulation", "MRP", "MRP"),
                      lower = c(pop_mean_p1[1], pop_mean_p2[1], samp_mean_p1_conf[1], samp_mean_p2_conf[1], ps_samp_mean_p1_conf[1], ps_samp_mean_p2_conf[1], weighted_samp_mean_p1_conf[1], weighted_samp_mean_p2_conf[1], sp_samp_mean_p1_conf[1], sp_samp_mean_p2_conf[1], model_popn_pref[2], model_popn_pref_p2[2]),
                      upper = c(pop_mean_p1[2], pop_mean_p2[2], samp_mean_p1_conf[2], samp_mean_p2_conf[2], ps_samp_mean_p1_conf[2], ps_samp_mean_p2_conf[2], weighted_samp_mean_p1_conf[2], weighted_samp_mean_p2_conf[2], sp_samp_mean_p1_conf[2], sp_samp_mean_p2_conf[2], model_popn_pref[3], model_popn_pref_p2[3]))

#write.csv(plotDat,
#          "W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/NERC_exploring_frontiers/Data/plotDat.csv",
#          row.names = F)

png("Fig1.png", width = 5, height = 5, units = "in", res = 500)
ggplot(data = plotDat, aes(x = p, y = est, colour = type, fill = type)) +
  geom_point() +
  geom_line() +
  theme_linedraw() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  labs(x = "",
       y = "Estimate",
       fill = "",
       colour = "") +
  scale_x_continuous(breaks = c(1,2), labels = c("1987-1999", "2010-2019"))
dev.off()

## now let's look at the trends
trends <- lapply(unique(plotDat$type),
                 function(x) {
                   data.frame(difference = plotDat$est[plotDat$p == 2 & plotDat$type == x] - plotDat$est[plotDat$p == 1 & plotDat$type == x],
                              estimator = x)
                 })

trends <- do.call("rbind", trends)

## and the confidence intervals for the trends

## sample mean
sample_mean_se <- sqrt(SE(ps_samp_mean_p1)^2 + SE(ps_samp_mean_p2)^2)

sample_mean_upper <- (samp_mean_p2 - samp_mean_p1) + 1.96 * sample_mean_se

sample_mean_lower <- (samp_mean_p2 - samp_mean_p1) - 1.96 * sample_mean_se

## quasi-randomisation
weighted_mean_se <- sqrt(SE(weighted_samp_mean_p1)^2 + SE(weighted_samp_mean_p2)^2)

weighted_mean_upper <- (weighted_samp_mean_p2 - weighted_samp_mean_p1) + 1.96 * weighted_mean_se

weighted_mean_lower <- (weighted_samp_mean_p2 - weighted_samp_mean_p1) - 1.96 * weighted_mean_se

## poststratification
ps_mean_se <- sqrt(SE(ps_samp_mean_p1)^2 + SE(ps_samp_mean_p2)^2)

ps_mean_upper <- (ps_samp_mean_p2 - ps_samp_mean_p1) + 1.96 * ps_mean_se

ps_mean_lower <- (ps_samp_mean_p2 - ps_samp_mean_p1) - 1.96 * ps_mean_se

## superpopulation
sp_mean_se <- sqrt(SE(sp_samp_mean_p1)^2 + SE(sp_samp_mean_p2)^2)

sp_mean_upper <- (sp_samp_mean_p2 - sp_samp_mean_p1) + 1.96 * sp_mean_se

sp_mean_lower <- (sp_samp_mean_p2 - sp_samp_mean_p1) - 1.96 * sp_mean_se

## MRP

diffs <- poststrat_prob_p2 - poststrat_prob

#diffs <- diffs[,1]

MRP_up <- quantile(diffs, probs = 0.975)

MRP_low <- quantile(diffs, probs = 0.025)

head(trends)

trends$lower <- c(NA, sample_mean_lower, ps_mean_lower, weighted_mean_lower, sp_mean_lower, MRP_low)

trends$upper <- c(NA, sample_mean_upper, ps_mean_upper, weighted_mean_upper, sp_mean_upper, MRP_up)

png("Fig1.png", width = 5, height = 4, units = "in", res = 500)
ggplot(data = trends, aes(x = difference, y = estimator)) +
         geom_point() +
         theme_linedraw() +
         geom_vline(xintercept = 0) +
  labs(x = "Trend",
       y = "") +
  geom_errorbar(aes(xmin = lower, xmax = upper, width = .2))
dev.off()

## demonstrating the effects of weighting

pop <- read.csv("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/NERC_exploring_frontiers/Data/all_data.csv")

pop <- pop[complete.cases(pop), ]

colnames(pop) <- c("truth_p1", "truth_p2", "sampled_p1", "sampled_p2",
                   "postcode_density", "open_access", "protected_area",
                   "road_length", "elevation", "inclusion_prob_p1",
                   "inclusion_prob_p2")

relFreqPlot <- function(pop,
                        R,
                        x,
                        weights,
                        breaks,
                        RNames) {

  dat <- lapply(1:length(R),
                function(y) {

                  stats <- lapply(x,
                                  function(z) {

                                    pop$bin <- cut(pop[,z], breaks = breaks, labels = FALSE)

                                    samp <- pop[pop[R[y]]==1,]

                                    samp$weights = weights[[y]]

                                    weightedFreq <- lapply(unique(pop$bin),
                                                           function(x) {
                                                             data.frame(weighted_sample = sum(samp$weights[samp$bin==x]) / sum(samp$weights),
                                                                        sample = nrow(samp[samp$bin== x,]) / nrow(samp),
                                                                        population = nrow(pop[pop$bin == x,]) / nrow(pop),
                                                                        bin = x,
                                                                        var = z,
                                                                        id = RNames[y])
                                                           })

                                    weightedFreq <- do.call("rbind", weightedFreq)

                                    melt(weightedFreq, id = c("bin", "var", "id"))

                                  })

                  if (length(x) > 1 | length(R) > 1) stats <- do.call("rbind", stats)

                })


    if (length(x) > 1 | length(R) > 1) dat <- do.call("rbind", dat)

    p <- ggplot(data=dat,aes(y = value, x = bin, colour = variable)) +
                geom_line() +
                theme_linedraw() +
                labs(colour = "",
                x = "",
                y = "Relative frequency")

    if (length(x) > 1 | length(R) > 1) p <- p + facet_grid(id~var,
                                                           scales = "free_y")

}

p <- relFreqPlot(pop = pop,
            x = c("road_length", "elevation"),
            R = c("sampled_p1", "sampled_p2"),
            RNames = c("Period_1", "Period_2"),
            weights = list(1/weighted_design_p1$prob, 1/weighted_design_p2$prob),
            breaks = 50)

png("Fig2.png", width = 7, height = 4, units = "in", res = 500)
p + theme_dark() + theme(legend.position = "bottom")
dev.off()

## token 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
