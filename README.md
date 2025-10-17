# Seam Carver

A **content-aware image resizer** that removes or inserts low-energy pixel paths (called seams) to resize images without distorting key content.

---

## Overview

Each pixel’s **energy** is computed using the gradient magnitude of color changes in both directions.  
Low-energy seams (least visually important paths) are found and removed to reduce the image’s width or height.

This implementation follows Avidan & Shamir’s 2007 paper *“Seam Carving for Content-Aware Image Resizing.”*

---

## Features
- Computes pixel energy using RGB gradients  
- Finds vertical and horizontal seams with the least total energy  
- Removes seams to resize images dynamically  
- Uses Princeton’s `algs4` library for data structures


---
---

# WordNet

A Java implementation of a simplified **WordNet** semantic network, inspired by the Princeton Algorithms II assignment.

## Overview

WordNet represents relationships between English nouns using two main concepts:
- **Synsets** – groups of synonymous words.
- **Hypernyms** – “is-a” relationships between synsets (e.g., `dog → animal`).

This project builds a directed acyclic graph (DAG) from given input files, where:
- each vertex represents a synset
- each directed edge represents a hypernym relationship

It supports finding the **shortest ancestral path (SAP)** and **semantic distance** between two nouns.

## Features

- Parse `synsets.txt` and `hypernyms.txt` input files  
- Validate graph structure as a rooted DAG  
- List all nouns in the dataset  
- Check if a word is a valid WordNet noun  
- Compute:
  - `distance(nounA, nounB)` → shortest path length between two nouns  
  - `sap(nounA, nounB)` → lowest common ancestor synset of two nouns

## Example

```java
WordNet wn = new WordNet("synsets.txt", "hypernyms.txt");
System.out.println(wn.sap("group_action", "action"));
System.out.println(wn.distance("group_action", "action"));
