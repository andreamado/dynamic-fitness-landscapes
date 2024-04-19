# Simulations
This folder contains:
- the source code to the Rust programs for the simulations. Contains three programs:
  - create_landscapes, that creates landscapes to be used by the other two programs
  - convergence, that runs simulations and records detailed information and graphs of the fitness landscapes
  - ecoevo_landscapes, that runs simulations and records statistical information on population and fitness landscapes
- run_simulations.py - An example script to run a batch of simulations in parallel

To compile the Rust programs, install Rust following the instructions in the [Rust webpage](https://www.rust-lang.org/tools/install). Then, open a terminal in the `simulations` folder and run the command `cargo build --release`. This will create the three executables described above in the folder `target/release/`. For instructions on how to run them use the `--help` option, e.g., `target/release/ecoevo_landscapes --help`.

## Examples
- The following command generates 5 Rough Mount Fuji landscapes with no additive effects and epistatic effects with a variance of 0.1 and a covariance of 0.05 [0.05 = 0.1 (variance) * 0.5 (correlation)] of the effect between resources `target/release/create_landscape --landscapes 5 --rmf 0 0 0 0.1 0.05`

- After, actual simulations could be run on these landscapes. The following command runs simulations on the landscapes with indexed 0 to 4, 10 replicates per landscape, with a mutation probability of 0.01 per locus per generation, population size of 20 and an equal amount of resources for each of the two resource types `target/release/ecoevo_landscapes --landscapes 0 5 --replicates 10 --mutation_rate 0.01 --size 20 --resources 1 1 --rmf 0 0 0 0.1 0.09`

## Parameters
The command `--help` lists the parameters.

The parameter `mu` corresponds to the mean value of the additive contribution.

The `c` parameters correspond to the entries of the covariance matrix used for the multivariate normal distribution. `c_diag` represents the diagonal entries of the covariance matrix, i.e., the variance, and the `c_offdiag` the offdiagonal entries, i.e., the covariance between different resources. The covariance between resources is given by the product of their correlation and the variance. `a` is used for the additive contribution and `b` for the epistatic contribution.