# Structural Controllability of Large-Scale Hypergraphs

This repository contains the code and experimental results accompanying the paper: *Structural Controllability of Large-Scale Hypergraphs*.

---

## Overview

Many real-world networked systems exhibit **higher-order interactions**, where the evolution of a state depends on multiple other states simultaneously. Examples arise in ecological systems, biological regulation, and engineered infrastructure networks. While structural controllability is well understood for graph-based systems with pairwise interactions, analogous theory for **hypergraph systems** remains largely undeveloped.

This work develops a **structural controllability framework for hypergraphs** by modeling hypergraph dynamics using polynomial dynamical systems. The resulting framework extends classical graph-based structural control theory to systems with higher-order interactions. The key idea is to represent the sparsity structure of polynomial dynamics using **directed hypergraphs**, allowing controllability properties to be analyzed using topology rather than precise parameter values. This work extends graph based definitions of vertex accessibility, dilations, and structural controllability to hypergraphs as shown here:

<br>

<p align="center">
  <img src="https://github.com/Jpickard1/hypergraph-structural-control/blob/main/figures/fig1.png" width="450">
</p>

<p align="left">
<em>
Figure 1: Hypergraph walk, accessibility, and dilation. This hypergraph is not structurally controllable.
Hyperedges are uniquely colored, and nodes are shaded to match the color of their incident hyperedge head.
Hyperedges are numbered according to the order they may appear in a walk originating at the control node.
The black and dark purple nodes are inaccessible, and the light purple nodes form a dilation.
</em>
</p>

## Main Contributions

The paper makes the following contributions:

- **Structural controllability framework for hypergraphs**  
  Extends classical structural control theory to nonlinear systems with higher-order interactions represented as hypergraphs. This work establishes conditions under which a polynomial system is structurally controllable using properties of its associated directed hypergraph.

- **Topology-based driver node selection**  
  Derives a structural lower bound on the number of driver nodes required to achieve controllability based on vertex accessibility and dilations found in the hypergraph structure.

- **Large-scale experiments**  
  Demonstrates scalability on hypergraphs with tens of thousands of nodes and higher-order interactions.

## Repository Structure
```
├── environment.yml        # Conda environment for reproducing experiments  
├── LICENSE  
├── README.md  
├── codes/                 # Python scripts and Jupyter notebooks for experiments and figure generation  
├── data/                  # Experimental results and intermediate outputs
└── figures/               # Figures included in the paper
```
---

## Contributors
- [Joshua Pickard](https://github.com/Jpickard1), Broad Institute of MIT and Harvard
- [Xin Mao](https://scholar.google.com/citations?user=uMLBc7gAAAAJ&hl=en), UNC School of Data Science and Society
- [Can Chen](https://tarheels.live/canc/), UNC School of Data Science and Society
