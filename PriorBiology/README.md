# Does Prior Biology Help?

### Problem definition
Test whether prior biology knowledge improves GEARS generalization to held-out combinatorial pairs by comparing three model variants: full GEARS (GO + coexpression graphs), GEARS without GO (no functional prior), and GEARS without coexpression (no expression graph). Report held-out prediction error for each.

### Problem info
The baseline GEARS conditions a graph-based gene expression model on two biological graphs: a functional prior derived from Gene Ontology (GO) and an expression graph derived from coexpression. Ablating each graph isolates its contribution to predicting novel multi-gene perturbations. Held-out pairs (never perturbed together in training) measure generalization.

### Solution
Modify the baseline GEARS pipeline as follows:
1. Train three GEARS variants on the same single- and multi-gene perturbation data: 
   - (a) both graphs (full GEARS)
   - (b) coexpression graph only (remove GO) 
   - (c) GO graph only (remove coexpression)

2. For each held-out combinatorial gene pair, predict the response profile with each variant.
3. Compute held-out prediction error per variant (e.g., mean squared error against observed response).
4. Compare error across the three models to quantify the value of each prior.
5. Report whether full GEARS, GO-only, or coexpression-only achieves the lowest held-out error.
