NETMED Tutorial
================
Ritoban Kundu
2024-12-13

## Introduction

# NETMED: Network-Based Mediation Analysis

## Overview

Social network or spatial interference induces spillover effects from
neighbors’ exposures, increasing the complexity of statistical analysis,
especially when mediators are involved.

The `NETMED` package addresses these challenges by: - Providing a
theoretical framework employing a structural graphical modeling
approach. - Investigating both mediation and interference effects within
network or spatial data. - Enabling the estimation of mechanistic
pathways through which neighboring units’ exposures and mediators exert
direct and indirect influences on individual outcomes. - Extending the
exposure mapping paradigm using Random-Effects Network Structural
Equation Models (REN-SEM) for spillover effects estimation. - Offering
maximum likelihood estimation and inference procedures with theoretical
guarantees.

## Features

The package estimates six distinct effects and provides variance
estimates for each:

1.  **Individual Direct Effect ($\tau_1$)**: Influence of self-exposure
    ($A$) on own outcome ($Y$).
2.  **Mediation Effect ($\tau_2$)**: Self-exposure ($A$) impacts own
    outcome ($Y$) through own mediator ($M$).
3.  **Neighbor-Mediation Effect ($\tau_3$)**: Self-exposure ($A$)
    affects own outcome ($Y$) via neighbors’ mediators.
4.  **Spillover Effect ($\tau_4$)**: Neighbors’ exposures
    ($\boldsymbol{A}_{N^\dagger}$) directly influence own outcome ($Y$).
5.  **Spillover-Mediation Effect ($\tau_5$)**: Neighbors’ exposures
    ($\boldsymbol{A}_{N^\dagger}$) impact own outcome ($Y$) through own
    mediator ($M$).
6.  **Complex Spillover-Mediation Effect ($\tau_6$)**: Neighbors’
    exposures ($\boldsymbol{A}_{N^\dagger}$) influence own outcome ($Y$)
    through neighbors’ mediators.

## Theoretical Guarantees

The framework provides: - Consistent maximum likelihood estimation for
REN-SEM. - Asymptotic variance estimators under a non-i.i.d. asymptotic
theory.

## Installation

Install the package directly from GitHub:

``` r
#devtools::install_github("Ritoban1/NETMED")
```

## Data Generation for Network Data

``` r
set.seed(100)
## Sample Size
n=400


# Exposure
prob_a=0.5

# Confounder
mean_c=0
sd_c=1



###
## Adjacency Matrix
E=matrix(0,n,n)

prob_p=0.05

for(i in 1:n){
  for(j in 1:n){
    if(j < i){
      E[i,j]=E[j,i]=rbinom(1,1,prob_p)
    }
    E[i,i]=1
  }
}


## Random effects
mean_b_y=0
mean_b_m=0
sd_b_m=0.1
sd_b_y=0.1

## Mediator
sigma_m=1
gamma=c(-1,2,0.9,1.8,0.7)


## Response
beta=c(-2,1.5,0.8,1.2,0.4,2.1,1.3)
sigma_y=1
E_rowsum=apply(E,2,sum)
## tilde E
E_t=E
for(i in 1:nrow(E)){
  E_t[i,]=E[i,]/(E_rowsum[i]-E[i,i])
  E_t[i,i]=0
}

# Exposure
A <- rbinom(n, size = 1, prob = prob_a)

# Confounder
C <- rnorm(n, mean_c, sd_c)

# Random Effects
b_m <- rnorm(n, mean_b_m, sd_b_m)
b_y <- rnorm(n, mean_b_y, sd_b_y)

# Mediator
M <- as.numeric(
  gamma[1] + gamma[2] * A + gamma[3] * (E_t %*% A) +
  gamma[4] * C + gamma[5] * (E_t %*% C) +
  E %*% b_m + rnorm(n, 0, sigma_m)
)

# Outcome
Y <- as.numeric(
  beta[1] + beta[2] * A + beta[3] * (E_t %*% A) +
  beta[4] * M + beta[5] * (E_t %*% M) +
  beta[6] * C + beta[7] * (E_t %*% C) +
  E %*% b_y + rnorm(n, 0, sigma_y)
)
```

## Estimation using NETMED

``` r
res=NETMED(A=A,M=M,Y=Y,E=E,C=as.matrix(C),C_names=as.vector("C"))


## All six effects
res$Effects_res
```
