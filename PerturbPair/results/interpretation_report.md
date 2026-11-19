# Prioritize Gene Pairs — Results Summary

Dataset: norman (Norman et al. 2019)
Number of held-out pairs evaluated: 71
Pairs with 95% CI entirely above zero: 71

## Top 10 ranked held-out pairs

| Rank | Gene pair    | Bootstrap mean score | 95% CI           | Significant |
| ---- | ------------ | -------------------- | ---------------- | ----------- |
| 1    | CEBPB+CEBPA  | 0.0991               | [0.0991, 0.0991] | Yes         |
| 2    | JUN+CEBPA    | 0.0820               | [0.0820, 0.0820] | Yes         |
| 3    | JUN+CEBPB    | 0.0594               | [0.0594, 0.0594] | Yes         |
| 4    | CEBPE+CEBPB  | 0.0589               | [0.0589, 0.0589] | Yes         |
| 5    | FOXA3+HOXB9  | 0.0486               | [0.0486, 0.0486] | Yes         |
| 6    | FOSB+CEBPB   | 0.0483               | [0.0483, 0.0483] | Yes         |
| 7    | FEV+ISL2     | 0.0462               | [0.0462, 0.0462] | Yes         |
| 8    | AHR+KLF1     | 0.0434               | [0.0434, 0.0434] | Yes         |
| 9    | LYL1+CEBPB   | 0.0431               | [0.0431, 0.0431] | Yes         |
| 10   | MAP2K3+IKZF3 | 0.0415               | [0.0415, 0.0415] | Yes         |

## Pathway interpretation

Significant GO-term associations (FDR < 0.05) for top 15 pairs: 884

### Top enriched GO terms

| Pair     | GO term    | p-value (adj) | N in pair / N in background |
| -------- | ---------- | ------------- | --------------------------- |
| AHR+KLF1 | GO:0000987 | 0.0003124     | 2 / 11                      |
| AHR+KLF1 | GO:0000976 | 0.003224      | 2 / 58                      |
| AHR+KLF1 | GO:0006355 | 0.003224      | 2 / 68                      |
| AHR+KLF1 | GO:0003700 | 0.004759      | 2 / 96                      |
| AHR+KLF1 | GO:0045893 | 0.005052      | 2 / 134                     |
| AHR+KLF1 | GO:0001094 | 0.005052      | 1 / 1                       |
| AHR+KLF1 | GO:0002841 | 0.005052      | 1 / 1                       |
| AHR+KLF1 | GO:0034752 | 0.005052      | 1 / 1                       |
| AHR+KLF1 | GO:1904613 | 0.005052      | 1 / 2                       |
| AHR+KLF1 | GO:0034751 | 0.005052      | 1 / 2                       |
| AHR+KLF1 | GO:0030888 | 0.005052      | 1 / 2                       |
| AHR+KLF1 | GO:0017025 | 0.005052      | 1 / 2                       |
| AHR+KLF1 | GO:0000785 | 0.005388      | 2 / 192                     |
| AHR+KLF1 | GO:0000981 | 0.005388      | 2 / 196                     |
| AHR+KLF1 | GO:1904322 | 0.005388      | 1 / 3                       |

## Files produced

- ranked_pairs.csv : full ranking with bootstrap mean, SE, 95% CI
- pair_scores_bootstrap.csv : per-pair bootstrap replicates summary
- pair_scores_raw.csv : single-model point-estimate scores
- go_enrichment_top_pairs.csv : GO enrichment for top pairs
- assets/top_pairs.png : bar plot of top pairs with CI error bars
- assets/distribution.png : histogram of all pair scores
