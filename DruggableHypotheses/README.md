## Find Druggable Hypotheses

> #### Research question: Which gene targets, when perturbed, suggest a druggable hypothesis?

### Deliverables
- Predicted transcriptional response per candidate target
- Dependency + drug-target annotation overlay per candidate
- Prioritized testable-target list annotated with evidence and caveat

### Problem definition
Use GEARS-predicted transcriptional responses to nominate gene targets whose perturbation suggests a druggable hypothesis, cross-referenced with gene dependency and known drug-annotation evidence. Prioritize testable targets while noting that a hypothesis is not a validated drug.

---

A GEARS-based analysis of <mark style="background: #ADCCFFA6;">single-gene perturbations</mark> in the Norman et al. (2019) dataset. The notebook <mark style="background: #FFF3A3A6;">loads the Norman dataset, prepares the simulation split, extracts candidate single-gene targets, loads or trains GEARS, predicts each target's response, overlays dependency and drug-target annotations, and ranks candidates into a testable-target list.</mark>


---

## Pipeline

- **Candidate targets** — single-gene perturbations kept only if present in the GEARS perturbation graph.
- **GEARS model** — loads `models/drug/` pretrained if present; otherwise trains from scratch (`hidden_size=16`, 20 epochs) and saves.
- **Drug annotation** — flag each target with `drug_annotations.py`: approved drugs, clinical/indirect coverage, or no coverage.
- **Rank** — sort by score descending; save `results/ranked_targets.csv`.


## References

- Norman et al. (2019), *Exploring genetic interaction manifolds constructed from rich single-cell phenotypes*. The notebook uses the Norman perturbation dataset.
- GEARS: the Gene Expression Analysis and Regulation using Systems biology framework used for perturbation response prediction.