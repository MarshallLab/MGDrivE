# MGDrivE: Mosquito Gene Drive Explorer

## Brief Description

MGDrivE is a framework designed to serve as a testbed in which gene-drive releases for mosquito-borne diseases control can be tested. It is being developed to accommodate various mosquito-specific gene drive systems within a population dynamics model that allows migration of individuals between nodes in a spatial landscape.

## Demonstration

In this demo, we are releasing a total of 100 mosquitoes homozygous for the CRISPR/CAS9 and one with a mutation that makes the mosquito resistant to the construct. Each node in the network represents a mosquito population laid down in a spatial scenario (this could be though of as a household, house block or even city if needed). We simulate how the genetic construct would propagate across the nodes of the network if mosquitoes were slowly migrating between populations with a probability based on proximity. To watch more videos take look at our [youtube](https://www.youtube.com/watch?v=sZXuUtToszw&list=PLRzY6w7pvIWqFJi94ZfhPkSVnazlUylpN) playlist.

## How does it work?

The main idea behind this model is to consider the inheritance matrix of genotypes a three-dimensional structure in which each intersection point determines the ratio/probability of a specific offspring genotype (z axis) provided that a certain combination of male-female genotypes (x and y axis). This allows us to use tensors as the basis for our calculations which has many advantages, some of them being: computational speed, model's transparency and extendability.

The second novel idea in our framework is to consider the spatial layout as a network of inter-connected breeding habitats. By performing this abstraction we are able to transform these landscapes into distances matrices, and then into transition probabilities matrices (through the use of movement kernels). This allows our framework to be able to model arbitrary topologies in which we can simulate mosquito populations mating and migrating in realistic geographical settings.

## Project's Status: Writing the paper

Most of the software development has been finished. Testing and debugging are almost done too. Writing is our last step towards the release of the package.
