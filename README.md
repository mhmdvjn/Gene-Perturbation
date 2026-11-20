# Gene-Perturbation
A computational research project focused on predicting transcriptional responses to genetic perturbations using deep learning, inspired by the GEARS framework. The repository explores gene perturbation modeling, gene-gene interactions, and the prediction of cellular responses to single- and multi-gene perturbations using single-cell RNA-seq data.


## Projects

### 1. PerturbPair 
[README: Perturbation analysis](./PerturbPair/)

A computational framework for analyzing
gene perturbation effects using machine learning
and downstream biological interpretation.

**Contents:**
- Machine learning-based modeling
- Bootstrap analysis
- Gene Ontology enrichment analysis

### 2. Prior Biology
[README: Does Prior Biology Help?](./PriorBiology/)

Test whether prior biology knowledge improves GEARS generalization
to held-out combinatorial gene pairs by comparing three model variants:
full GEARS (GO + coexpression), GEARS without GO, and GEARS without coexpression.

**Contents:**
- Three trained GEARS variants
- Held-out prediction error comparison
- Conclusion on prior biology contribution

### 3. Druggable Hypotheses
[README: Find Druggable Hypotheses](./DruggableHypotheses/)

Nominate gene targets whose perturbation suggests a druggable hypothesis by
combining GEARS-predicted transcriptional responses with gene dependency and
drug-target annotation evidence.

**Contents:**
- GEARS-predicted response per candidate target
- Dependency + drug-target annotation overlay
- Prioritized testable-target list with caveat
