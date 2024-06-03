# *Drosophila melanogaster* Bottleneck Simulations

## Overview
This project involves the simulation of genetic differentiation in a neutrally evolving population of *Drosophila* after consecutive winter bottlenecks. The goal is to study the effects of these bottlenecks on the genetic diversity of a wild-caught population in Charlotteville, Virginia, USA. The scripts and data in this repository is part of the earlier work published in this article: (https://academic.oup.com/genetics/article/226/2/iyad207/7459204).

## Author
Started by Connor S. Murray  
Date: October 31, 2022

## Workflow

### 1. Generate the Burn-in VCF
This step generates the VCF for a metapopulation to achieve a theta pi ~ 0.01, representing the genetic diversity of a wild *Drosophila* metapopulation.

**Tools Used:**  
- `msprime` (Python)

### 2. Run the SLiM Scripts and Output Statistics
This step involves running SLiM scripts to independently replicate the simulated universe over 50 generations, including 2 overwintering bottleneck events.

**Tools Used:**  
- `SLiM`

### 3. Merge Parsed Data
This step collects data from each simulation run and concatenates it for the analysis step.

**Tools Used:**  
- `bash`
- `python`
- `R`

### 4. Analyze, Perform ABC, and Plot Data
This step involves performing Approximate Bayesian Computation (ABC) and plotting the results.

**Tools Used:**  
- `R`
