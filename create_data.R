library(cowfit)


# LMM ---------------------------------------------------------------------

set.seed(1)
myherd <- mastitis[c("herd", "id", "sire")]
df <- expand.grid(protein = 1:5,
            animal = unique(mastitis$id))
df <- merge(df, myherd, by.x = "animal", by.y = "id")
myf <- y ~ protein + (1|herd) + (protein|sire)
(myvar <- cowfit_var_comp(formula = myf, data = df))
myvar$vcov <- c(500, 400, 300, -200, 50)
dat <- sim_lmer(formula = myf, data = df
         , pedigree = list(sire = pedSires), var_comp = myvar$vcov
         , beta = c(6000, 200), return_ranef = FALSE)
sim_milk <- dat[c("animal", "sire", "herd", "protein", "y")]
save(sim_milk, file="data/sim_milk.RData")


# GLMM --------------------------------------------------------------------

set.seed(1)
myherd <- mastitis[c("herd", "id", "sire")]
df <- expand.grid(lact = 1:5,
                  animal = unique(mastitis$id))
df <- merge(df, myherd, by.x = "animal", by.y = "id")
myf <- y ~ lact + (lact|sire)
(myvar <- cowfit_var_comp(formula = myf, data = df))
myvar$vcov <- c(4, 3, -0.5, 1)
dat <- sim_glmer(formula = myf, data = df,
                 family = "binomial", pedigree = list(sire = pedSires),
                 var_comp = myvar$vcov, beta = c(-3, 0.5), return_ranef = FALSE)
sim_twin <- dat[c("animal", "sire", "herd", "lact", "y")]
save(sim_twin, file="data/sim_twin.RData")
