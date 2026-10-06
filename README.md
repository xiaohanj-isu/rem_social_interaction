# A Genetic Relational Event Model for Dynamic Analysis of Aggressive Interactions

This study uses a relational event model to analyse directed aggression among
pigs, accounting for interaction history, genetic effects, permanent
environmental effects and group effects. 

## Files

- `rem_table_all.rds`: prepared analysis table with 3,151,438 event-by-dyad rows,
  representing 18,230 directed events in 57 groups and 772 animals.
- `G.rds`: 772 by 772 genomic relationship matrix, with animal identifiers as
  row and column names.
- `rem_poisson_gr_cor.R`: reads the data, compiles the model, fits it and saves
  the results.
- `rem_poisson_gr_cor.cpp`: TMB implementation of the Poisson REM.

## Variables

Each row in `rem_table_all.rds` represents a directed dyad in the risk set. 
Read both data files with `readRDS()`.

- `group` and `event_number`: group identifier and event index within that group.
- `giver` and `receiver`: animal identifiers for the aggression giver and
  receiver. 
- `y`: 1 for the observed dyad, 0 for other dyads in the risk set.
- `offset`: logarithm of the inter-event exposure time in minutes.
- `nursery` and `litter`: whether the two pigs shared a nursery pen or litter,
  respectively; 1 = shared, 2 = not share.
- `giver_wt` and `receiver_wt`: body weight in kg.
- `outdegreeSender`, `indegreeSender`, `outdegreeReceiver` and
  `indegreeReceiver`: prior outgoing or incoming aggression for the giver
  or receiver.
- `inertia` and `reciprocity`: prior aggression in the same or reverse dyadic
  direction, respectively.
- `otp`, `itp`, `osp` and `isp`: outgoing two-paths, incoming two-paths,
  outgoing shared partners and incoming shared partners, respectively.
- `RF`, `PP` and `IP`: prior dyadic frequencies of reciprocal fighting,
  parallel pressing and inverse pressing. 
- `baseline`: constant 1 for intercept.


## Software

R version 4.5.0
TMB 1.9.21

## Run

The script reads the two RDS files, compiles the C++ model and fits it.
