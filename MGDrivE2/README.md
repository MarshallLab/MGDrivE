# **MGDrivE 2.0**: Mosquito Gene Drive Explorer 2

## Installation

To install **MGDrivE 2.0**, please use `devtools::install_github(repo = "MarshallLab/MGDrivE",subdir = "MGDrivE2")`.

## Brief Description

**MGDrivE 2.0** is a new simulation platform which extends capabilities from the **MGDrivE** simulation package in a new mathematical and computational framework. For more information about **MGDrivE**, see our [publication](https://doi.org/10.1111/2041-210X.13318). Some of the notable capabilities of **MGDrivE 2.0** include incorporation of human populations, epidemiological dynamics, time-varying parameters, and a continuous-time simulation framework with various sampling algorithms for both deterministic and stochastic interpretations.

## Computational Features

The design feature that has made it possible to include a rich set of functionality into **MGDrivE 2.0** is using a Petri net (PN) formalism to structure the software. By separating how a model is specified (state space, and events that are allowed to change state), and numerical methods which draw trajectories from a model, numerical methods in **MGDrivE 2.0** are independent of any particular model, as long as it is expressed as a PN. Additionally, the well understood PN formalism, based on a bipartite graph describing how states and events are allowed to affect one another, means that quickly adding new features (such as human populations) can be done without needing to rewrite large parts of the existing codebase.

When assigning to each event a function that gives the current rate that event will occur (fire) at, given the present state and time, and assuming the Markov property, the simulation becomes a continuous-time Markov chain, and the PN is referred to as a stochastic Petri net (SPN). If we assume that rather than describing exponentially distributed random variables, those rate functions describe deterministic (continuous) rates of events firing, the SPN mathematically becomes a set of ordinary differential equations (ODE). In either case we can rely on the wealth of algorithms developed for solving ODEs or simulating CTMCs.

**MGDrivE 2.0** uses the cube data structure from **MGDrivE** [(CRAN link)](https://CRAN.R-project.org/package=MGDrivE) to parameterize genetic inheritance.

## Vignettes

There are a large number of vignettes included in **MGDrivE 2.0** to describe its functionality. Aa a suggestion, we recommend new users to read them in the following rough order to acquaint themselves with the software.

  1. <a href="https://marshalllab.github.io/MGDrivE/docs_v2/articles/lifecycle-node.html">MGDrivE2: One Node Lifecycle Dynamics</a>: this vignette describes the basic use of **MGDrivE 2.0** in a single node, with only mosquitoes present (no epidemiology). It introduces how to build the Petri net objects, specify parameters, make a vector of hazard functions, and simulate from the SPN using several of the provided simulation algorithms.
  2. <a href="https://marshalllab.github.io/MGDrivE/docs_v2/articles/lifecycle-network.html">MGDrivE2: Metapopulation Network Lifecycle Dynamics</a>: this describes how to set up and simulate metapopulation networks of mosquitoes, without introducing epidemiological dynamics yet.
  3. <a href="https://marshalllab.github.io/MGDrivE/docs_v2/articles/epi-node.html">MGDrivE2: One Node Epidemiological Dynamics</a>: introduces the SIS-SEI (Susceptible Infected Susceptible - Susceptible Exposed Infectious) human-mosquito epidemiological dynamics in a single node.
  4. <a href="https://marshalllab.github.io/MGDrivE/docs_v2/articles/epi-network.html">MGDrivE2: Metapopulation Network Epidemiological Dynamics</a>: shows the user how to simulate SIS-SEI dynamics on a metapopulation of nodes.
  5. <a href="https://marshalllab.github.io/MGDrivE/docs_v2/articles/epi-SEIR.html">MGDrivE2: SEIR Epidemiological Dynamics</a>: demonstrates the SEIR-SEI (Susceptible Exposed Infectious Recovered) option for the human population, to model strongly immunizing arboviruses.
  6. <a href="https://marshalllab.github.io/MGDrivE/docs_v2/articles/output_storage.html">MGDrivE2: Data Storage and Analysis</a>: after a user is familiar with the basic simulation workflow, this vignette introduces the simulation functions that write output to external .CSV files on the user's system, and other functions to analyze those output; writing to external .CSV files is useful for large-scale Monte Carlo simulations.
  7. <a href="https://marshalllab.github.io/MGDrivE/docs_v2/articles/inhomogeneous.html">MGDrivE2: Simulation of Time-inhomogeneous Stochastic Processes (Seasonality)</a>: shows the user how to write custom time-varying hazard functions to simulate inhomogeneous (where rates depend on time) models, using temperature-driven adult mosquito mortality as an example.
  8. <a href="https://marshalllab.github.io/MGDrivE/docs_v2/articles/advanced_topics.html">MGDrivE2: Advanced Topics</a>: shows advanced users how to write their own numerical methods to sample from the SPN; uses a simple explicit Euler method for ODEs as an example.
