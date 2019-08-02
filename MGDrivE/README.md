# **MGDrivE**: Mosquito Gene Drive Explorer

## Brief Description

**MGDrivE** is a framework designed to serve as a testbed in which gene-drive releases for mosquito-borne diseases control can be tested. It is being developed to accommodate various mosquito-specific gene drive systems within a population dynamics model that allows migration of individuals between nodes in a spatial landscape.

## How does it work?

The main idea behind this model is to consider genetic inheritance a three-dimensional cube in which each element determines the probability of a specific offspring genotype (z axis) given a certain combination of male-female parent genotypes (x and y axis). This allows us to use tensors as the basis for computation which has many advantages, some of them being: computational speed, model transparency and extendability.

The second novel idea in our framework is to consider the spatial layout as a network of inter-connected breeding habitats. By performing this abstraction we are able to transform these landscapes into distances matrices, and then into transition probabilities matrices. This allows our framework to be able to model arbitrary topologies in which we can simulate mosquito populations mating and migrating in realistic geographical settings.

For more details, please read the vignettes accompanying this package!

## Demonstration

In this demo, we release 100 mosquitoes homozygous for the CRISPR/Cas-9 homing drive system, and one with a mutation that makes the mosquito resistant to the construct. Each node in the network represents a mosquito population laid down in a spatial scenario (this could be though of as a household, house block or even city if needed). We simulate how the genetic construct would propagate across the nodes of the network if mosquitoes were slowly migrating between populations with a probability based on proximity. 

<br>
<div align="center">
<a href="https://youtu.be/sZXuUtToszw" target="_blank"><img src="http://img.youtube.com/vi/sZXuUtToszw/0.jpg" width="640" height="480" border="10"/></a>
</div>
<br>

To watch more videos please take a look at our [youtube playlist!](https://www.youtube.com/playlist?list=PLRzY6w7pvIWqFJi94ZfhPkSVnazlUylpN)
