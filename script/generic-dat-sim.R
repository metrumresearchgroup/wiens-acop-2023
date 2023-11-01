# This script generates data for a DCM without reference to specific PK related covariates
# The covariate relationships include non-linear and interaction effect

library(tidyverse)
library(mrgsolve)
library(mrggsave)
library(here)
library(data.table)
library(dmutate)

dataDir <- here("data", "derived")
thisScript <- "generic-dat-sim.R"

### data simulation model ----------------------------
runno= '108'

simmod <- mread(here(glue::glue("model/simmod/{runno}b"))) 

n = 10000

set.seed(1901)
cov_dat <- tibble(ID = 1:n,
                  DOSE = sample(c(5, 20, 50, 100), size = n, replace = T),
                  SEX = sample(c(0, 1), size = n, replace = T),
                  WT = round(runif(n, 60, 70)),
                  X1 = runif(n, -2, 2),
                  X2 = runif(n, -2, 2),
                  X3 = runif(n, -2, 2),
                  X4 = runif(n, -2, 2),
                  X5 = runif(n, -2, 2),
                  X6 = runif(n, -2, 2))

params <- list(THETA1 = log(1), ## ka
               THETA2 = log(70), 
               THETA3 = log(4), 
               THETA4 = log(70), 
               THETA5 = log(3), 
               THETA6 = 0.05, 
               THETA7 = -1, 
               THETA8 = -0.2, 
               THETA9 = 0.5, 
               THETA10 = 0.1)


simmod <- update(simmod,
                 param = params)

simdat <-   mutate(cov_dat,
                   TIME=0,
                   EVID=1,
                   AMT=DOSE,
                   CMT=2)


set.seed(1)
simout <- mrgsim(simmod, 
                data = simdat,
                carry_out = "EVID,AMT,CMT",
                recover = "CL",
                output = "df", 
                # obsonly=T, ## retaining for NONMEM data set
                Req="Y,IPRED,TVLCL,CL,V2,V3,KA,Q",
                tgrid=c(0.5,seq(1,12,1),seq(24,72,24), 120, 168),
                quiet  = TRUE) 

simout2 <- cov_dat %>%
  left_join(simout %>%
              filter(EVID == 0) %>% ## manually filter to just obs
              select(ID,CL,TVLCL,IPRED,V2,V3,KA,Q,TIME,DV=Y))

###  check covariate relationships ----------------------
simout2 %>%
  distinct(ID, TVLCL, CL) %>%
  slice_sample(n = 600) %>%
  ggplot(aes(x = TVLCL, y = log(CL))) +
  geom_point(size = 0.75) +
  geom_abline(slope = 1, linetype = 2)

simout2 %>%
  distinct(ID, X2, CL, TVLCL) %>%
  slice_sample(n = 600) %>%
  ggplot(aes(x = X2, y = TVLCL)) +
  geom_point(size = 0.75) 

pset <- simout2 %>% 
  distinct(ID, TVLCL, SEX, X1, X2, X3, X4, X5, X6) %>% 
  slice_sample(n = 600) %>%
  pivot_longer(cols=X1:X6) %>% 
  mutate(
    Sex = if_else(SEX == 1, "Male", "Female"),
    CL = exp(TVLCL)
  ) 

p1 <- pset %>% 
  ggplot(aes(x=value, y=TVLCL)) + 
    geom_point(size = 0.75) + 
    facet_grid(Sex~name, scales="free_x") + 
    theme_bw()

p2 <- pset %>% 
  ggplot(aes(x=value, y=CL)) + 
    geom_point(size = 0.75) + 
    facet_grid(Sex~name, scales="free_x") + 
    theme_bw()

mrggsave(
  list(p1,p2),
  dir = here("deliv/figure"),
  stem = "cl-v-covs", 
  width=7, height=5, 
  script="generic-dat-sim.R"
)

map(seq(1, 45, 9), ~{
  simout2 %>% filter(between(ID, .x, .x+8)) %>% 
    ggplot() + 
    geom_point(aes(x=TIME, y=DV, group=ID)) +
    geom_line(aes(x=TIME, y=IPRED, group=ID)) +
    geom_smooth(aes(x=TIME, y=DV, group=ID)) +
    scale_y_log10() +
    facet_wrap(~ID, scales="free_y")
})

## huge tail on half-life, need to pull in
simout2 %>% 
  mutate(t12 = mrgmisc::half_life(CL,V2,V3,Q)) %>% 
  distinct(ID, t12) %>% 
  ggplot() + 
    geom_boxplot(aes(y=t12)) 

simout2 %>% 
  mutate(t12 = mrgmisc::half_life(CL,V2,V3,Q)) %>% 
  distinct(ID, t12, CL,V2) %>% 
  ggplot() + 
  geom_point(aes(x=CL,y=t12,colour=V2)) + 
  scale_x_log10() + 
  theme_bw() + theme(legend.position = "bottom") +
  labs(y="Half-life (hours)", x="CL (L/hour)", colour="V2 (L)")

### Output ------------------------------------------

data.table::fwrite(
  x = simout2,
  file = here::here("data", "derived", "pk-sim-N10000-generic.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  na = "."
)

data.table::fwrite(
  x = simout %>% 
    transmute(C=".", ID, TIME, MDV=EVID, EVID, AMT, CMT=1,
              DV=if_else(EVID == 1, 0, Y), 
              NUM=row_number()) %>% 
    left_join(cov_dat, by="ID"),
  file = here::here("data", "derived", "nm-N10000-generic.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  na = "."
)
