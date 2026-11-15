import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    _HAS_PLOTLY = True
except Exception:
    _HAS_PLOTLY = False

try:
    import scanpy as sc
    _HAS_SCANPY = True
except Exception:
    _HAS_SCANPY = False

sns.set_theme(style="whitegrid", font_scale=1.1)


def _gene_idx(model, gene_names):
    gene_list = model.gene_list
    gene2idx = {g: i for i, g in enumerate(gene_list)}
    return np.array([gene2idx[g] for g in gene_names if g in gene2idx], dtype=int)


def _pert_key(pert):
    if isinstance(pert, (list, tuple)):
        return "_".join(pert)
    return pert


def _ctrl_mean(model):
    return model.ctrl_expression.detach().cpu().numpy()


def _top_de_genes(model, pert, n=20):
    adata = model.adata
    cond2name = dict(adata.obs[["condition", "condition_name"]].values)
    name2cond = {v: k for k, v in cond2name.items()}
    full = name2cond.get(pert, pert)
    key = cond2name.get(full, full)
    de = adata.uns.get("top_non_dropout_de_20", {}).get(key, None)
    if de is None:
        de = adata.uns.get("rank_genes_groups_cov_all", {}).get(key, [])[:n]
    return list(de)[:n]


def plot_pert_pred_vs_truth(model, pert, top_n=20, save_file=None, show=True):
    """Boxplot of predicted vs true delta-expression over control for top DE genes."""
    adata = model.adata
    gene_names = _top_de_genes(model, pert, n=top_n)
    idx = _gene_idx(model, gene_names)
    if len(idx) == 0:
        raise ValueError(f"No DE genes found for perturbation {pert}")

    cond = pert if pert in adata.obs["condition"].values else None
    if cond is None:
        cond2name = dict(adata.obs[["condition", "condition_name"]].values)
        cond = [k for k, v in cond2name.items() if v == pert][0]

    truth = adata[adata.obs.condition == cond].X.toarray()[:, idx]
    ctrl = adata[adata.obs["condition"] == "ctrl"].to_df().mean()[gene_names].values
    truth = truth - ctrl

    pred = model.predict([list(pert.split("+")) if "+" in pert else [pert]])[_pert_key(pert)][idx]
    ctrl_mean = _ctrl_mean(model)[idx]
    pred = pred - ctrl_mean

    df = pd.DataFrame({
        "gene": np.repeat(gene_names, truth.shape[0]),
        "truth": truth.ravel(),
        "pred": np.repeat(pred, truth.shape[0]),
    })

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    sns.boxplot(x="gene", y="truth", data=df, ax=axes[0], showfliers=False,
                medianprops=dict(linewidth=0), color="#4c72b0")
    axes[0].set_title(f"{pert} - Truth (delta over ctrl)")
    axes[0].axhline(0, linestyle="--", color="green")
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=90)

    df2 = pd.DataFrame({
        "gene": np.repeat(gene_names, truth.shape[0]),
        "pred": np.repeat(pred, truth.shape[0]),
    })
    sns.boxplot(x="gene", y="pred", data=df2, ax=axes[1], showfliers=False,
                medianprops=dict(linewidth=0), color="#dd8452")
    axes[1].set_title(f"{pert} - Prediction (delta over ctrl)")
    axes[1].axhline(0, linestyle="--", color="green")
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=90)

    fig.suptitle(f"Top {top_n} DE genes: {pert}", fontsize=14)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return fig


def plot_pred_vs_truth_scatter(model, pert, gene_subset=None, n_genes=200,
                               save_file=None, show=True, interactive=False):
    """Scatter of predicted vs true expression with correlation metrics."""
    adata = model.adata
    if gene_subset is None:
        gene_subset = _top_de_genes(model, pert, n=n_genes)
    idx = _gene_idx(model, gene_subset)
    if len(idx) == 0:
        raise ValueError(f"No genes found for perturbation {pert}")

    cond = pert if pert in adata.obs["condition"].values else None
    if cond is None:
        cond2name = dict(adata.obs[["condition", "condition_name"]].values)
        cond = [k for k, v in cond2name.items() if v == pert][0]

    truth = adata[adata.obs.condition == cond].X.toarray()[:, idx].mean(axis=0)
    pred = model.predict([list(pert.split("+")) if "+" in pert else [pert]])[_pert_key(pert)][idx]
    ctrl = _ctrl_mean(model)[idx]
    truth_d = truth - ctrl
    pred_d = pred - ctrl

    r, _ = pearsonr(pred_d, truth_d)
    s, _ = spearmanr(pred_d, truth_d)
    mse = float(np.mean((pred_d - truth_d) ** 2))

    if interactive and _HAS_PLOTLY:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=truth_d, y=pred_d, mode="markers",
                                 marker=dict(size=6, color="#4c72b0", opacity=0.7),
                                 name="genes"))
        lims = [min(truth_d.min(), pred_d.min()), max(truth_d.max(), pred_d.max())]
        fig.add_trace(go.Scatter(x=lims, y=lims, mode="lines", line=dict(dash="dash", color="green"),
                                 name="y=x"))
        fig.update_layout(
            title=f"{pert}: Pred vs Truth (r={r:.3f}, rho={s:.3f}, MSE={mse:.4f})",
            xaxis_title="True delta over ctrl", yaxis_title="Predicted delta over ctrl",
            width=700, height=600)
        if save_file:
            fig.write_html(save_file)
        if show:
            fig.show()
        return fig

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(truth_d, pred_d, s=20, alpha=0.7, color="#4c72b0", edgecolor="k", linewidth=0.3)
    lims = [min(truth_d.min(), pred_d.min()), max(truth_d.max(), pred_d.max())]
    ax.plot(lims, lims, ls="--", color="green", label="y=x")
    ax.set_xlabel("True delta over ctrl")
    ax.set_ylabel("Predicted delta over ctrl")
    ax.set_title(f"{pert}: r={r:.3f}, rho={s:.3f}, MSE={mse:.4f}")
    ax.legend()
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return {"fig": fig, "pearson": r, "spearman": s, "mse": mse}


def plot_de_heatmap(model, perts, top_n=20, save_file=None, show=True, method="mean"):
    """Heatmap of predicted vs true delta-expression for top DE genes across perts."""
    adata = model.adata
    all_genes = []
    for p in perts:
        all_genes.extend(_top_de_genes(model, p, n=top_n))
    genes = list(dict.fromkeys(all_genes))
    idx = _gene_idx(model, genes)
    if len(idx) == 0:
        raise ValueError("No DE genes found")

    ctrl = _ctrl_mean(model)[idx]
    rows, labels = [], []
    for p in perts:
        key = list(p.split("+")) if "+" in p else [p]
        pred = model.predict([key])[_pert_key(p)][idx] - ctrl
        cond2name = dict(adata.obs[["condition", "condition_name"]].values)
        cond = [k for k, v in cond2name.items() if v == p]
        if cond:
            truth = adata[adata.obs.condition == cond[0]].X.toarray()[:, idx].mean(axis=0) - ctrl
        else:
            truth = np.zeros_like(pred)
        if method == "mean":
            rows.append(np.stack([pred, truth], axis=0))
        else:
            rows.append(pred)
            rows.append(truth)
        labels.extend([f"{p} (pred)", f"{p} (truth)"])

    data = np.concatenate(rows, axis=0) if method == "mean" else np.stack(rows, axis=0)

    fig, ax = plt.subplots(figsize=(max(8, len(genes) * 0.35), max(4, len(labels) * 0.4)))
    im = ax.imshow(data, aspect="auto", cmap="vlag", vmin=-np.abs(data).max(), vmax=np.abs(data).max())
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_title("Top DE genes: predicted vs true (delta over ctrl)")
    fig.colorbar(im, ax=ax, label="Delta expression")
    fig.tight_layout()
    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return fig


def plot_dotplot_top_de(model, pert, top_n=20, save_file=None, show=True):
    """Dotplot comparing predicted vs true mean delta-expression per gene."""
    adata = model.adata
    genes = _top_de_genes(model, pert, n=top_n)
    idx = _gene_idx(model, genes)
    if len(idx) == 0:
        raise ValueError(f"No DE genes found for {pert}")

    cond2name = dict(adata.obs[["condition", "condition_name"]].values)
    cond = [k for k, v in cond2name.items() if v == pert][0]
    truth = adata[adata.obs.condition == cond].X.toarray()[:, idx].mean(axis=0)
    pred = model.predict([list(pert.split("+")) if "+" in pert else [pert]])[_pert_key(pert)][idx]
    ctrl = _ctrl_mean(model)[idx]
    truth_d = truth - ctrl
    pred_d = pred - ctrl

    df = pd.DataFrame({
        "gene": genes + genes,
        "value": np.concatenate([truth_d, pred_d]),
        "type": ["Truth"] * len(genes) + ["Prediction"] * len(genes),
    })

    fig, ax = plt.subplots(figsize=(max(8, len(genes) * 0.4), 5))
    sns.stripplot(x="gene", y="value", hue="type", data=df, ax=ax,
                  dodge=True, jitter=0.15, size=5, alpha=0.8,
                  palette={"Truth": "#4c72b0", "Prediction": "#dd8452"})
    ax.axhline(0, linestyle="--", color="green")
    ax.set_title(f"{pert}: Top {top_n} DE genes - Truth vs Prediction")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_ylabel("Delta expression over ctrl")
    ax.legend(title="", loc="best")
    fig.tight_layout()
    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return fig


def plot_fold_change(model, pert, save_file=None, show=True, threshold=0.5):
    """Volcano-style fold-change comparison: predicted vs true."""
    adata = model.adata
    genes = model.gene_list
    idx = np.arange(len(genes))
    cond2name = dict(adata.obs[["condition", "condition_name"]].values)
    cond = [k for k, v in cond2name.items() if v == pert][0]
    truth = adata[adata.obs.condition == cond].X.toarray()[:, idx].mean(axis=0)
    pred = model.predict([list(pert.split("+")) if "+" in pert else [pert]])[_pert_key(pert)]
    ctrl = _ctrl_mean(model)

    truth_fc = truth / np.maximum(ctrl, 1e-8)
    pred_fc = pred / np.maximum(ctrl, 1e-8)
    truth_fc[np.isnan(truth_fc)] = 0
    pred_fc[np.isnan(pred_fc)] = 0
    truth_fc[np.isinf(truth_fc)] = 0
    pred_fc[np.isinf(pred_fc)] = 0

    mask = np.maximum(truth, ctrl) > threshold
    x = np.log2(np.maximum(pred_fc, 1e-3))[mask]
    y = np.log2(np.maximum(truth_fc, 1e-3))[mask]
    g = np.array(genes)[mask]

    r, _ = pearsonr(x, y)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=10, alpha=0.6, color="#4c72b0")
    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lims, lims, ls="--", color="green", label="y=x")
    ax.axvline(0, color="k", lw=0.5)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("log2 predicted fold change over ctrl")
    ax.set_ylabel("log2 true fold change over ctrl")
    ax.set_title(f"{pert}: Fold change (r={r:.3f})")
    ax.legend()
    fig.tight_layout()
    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return {"fig": fig, "pearson": r}


def plot_gi_synergy(model, combo, save_file=None, show=True):
    """Bar chart of GI coefficients for a combinatorial perturbation."""
    from .utils import get_GI_params
    gi = get_GI_params(model.GI_predict(combo), combo)
    keys = ["c1", "c2", "mag", "dcor", "dcor_singles", "dcor_first", "dcor_second",
            "corr_fit", "dominance", "eq_contr"]
    vals = [gi.get(k, np.nan) for k in keys]

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["#4c72b0", "#4c72b0", "#dd8452", "#55a868", "#55a868", "#55a868",
              "#55a868", "#8172b3", "#8172b3", "#8172b3"]
    ax.bar(keys, vals, color=colors)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticklabels(keys, rotation=45, ha="right")
    ax.set_ylabel("Value")
    ax.set_title(f"Genetic Interaction: {' + '.join(combo)}")
    fig.tight_layout()
    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return {"fig": fig, "gi": gi}


def plot_uncertainty(model, pert, top_n=30, save_file=None, show=True):
    """Plot per-gene uncertainty (logvar) vs absolute prediction error."""
    if not model.config.get("uncertainty"):
        raise ValueError("Uncertainty mode is not enabled for this model.")
    preds, logvar_sum = model.predict([list(pert.split("+")) if "+" in pert else [pert]])
    key = _pert_key(list(pert.split("+")) if "+" in pert else [pert])
    unc = model.saved_logvar_sum.get(key, None)
    if unc is None:
        raise ValueError("No uncertainty saved for this perturbation.")

    adata = model.adata
    genes = model.gene_list
    cond2name = dict(adata.obs[["condition", "condition_name"]].values)
    cond = [k for k, v in cond2name.items() if v == pert][0]
    truth = adata[adata.obs.condition == cond].X.toarray().mean(axis=0)
    pred = preds[key]
    ctrl = _ctrl_mean(model)
    err = np.abs((pred - ctrl) - (truth - ctrl))

    order = np.argsort(err)[::-1][:top_n]
    x = unc[order]
    y = err[order]
    g = np.array(genes)[order]

    fig, ax = plt.subplots(figsize=(max(8, top_n * 0.3), 5))
    ax.scatter(x, y, s=40, c=y, cmap="viridis", edgecolor="k", linewidth=0.3)
    for i in range(len(g)):
        ax.annotate(g[i], (x[i], y[i]), fontsize=6, ha="left", va="bottom")
    ax.set_xlabel("Uncertainty (exp(-mean logvar))")
    ax.set_ylabel("Absolute prediction error (delta over ctrl)")
    ax.set_title(f"{pert}: Uncertainty vs Error (top {top_n} genes by error)")
    fig.tight_layout()
    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return fig


def plot_umap_by_condition(model, color_by="condition", save_file=None, show=True):
    """UMAP of cells colored by condition or prediction error (requires scanpy)."""
    if not _HAS_SCANPY:
        raise ImportError("scanpy is required for UMAP plotting.")
    adata = model.adata.copy()
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    sc.pl.umap(adata, color=color_by, show=False, save=save_file)
    if show:
        plt.show()
    return adata


def plot_training_history(train_losses, val_losses=None, save_file=None, show=True):
    """Plot training and validation loss curves."""
    fig, ax = plt.subplots(figsize=(7, 5))
    epochs = np.arange(1, len(train_losses) + 1)
    ax.plot(epochs, train_losses, marker="o", label="Train loss", color="#4c72b0")
    if val_losses is not None:
        ax.plot(epochs, val_losses, marker="o", label="Val loss", color="#dd8452")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.set_title("Training history")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return fig


def plot_pert_network(model, save_file=None, show=True, layout="spring"):
    """Plot the gene similarity / perturbation network (GO or co-express)."""
    import networkx as nx
    G = None
    if model.config.get("G_go") is not None:
        edge_index = model.config["G_go"].cpu().numpy()
        weight = model.config["G_go_weight"].cpu().numpy() if model.config.get("G_go_weight") is not None else None
        pert_list = model.pert_list
        G = nx.DiGraph()
        for i in range(edge_index.shape[1]):
            src, dst = int(edge_index[0, i]), int(edge_index[1, i])
            G.add_edge(pert_list[src], pert_list[dst], weight=float(weight[i]) if weight is not None else 1.0)
    elif model.config.get("G_coexpress") is not None:
        edge_index = model.config["G_coexpress"].cpu().numpy()
        weight = model.config["G_coexpress_weight"].cpu().numpy() if model.config.get("G_coexpress_weight") is not None else None
        gene_list = model.gene_list
        G = nx.Graph()
        for i in range(edge_index.shape[1]):
            src, dst = int(edge_index[0, i]), int(edge_index[1, i])
            G.add_edge(gene_list[src], gene_list[dst], weight=float(weight[i]) if weight is not None else 1.0)

    if G is None or G.number_of_nodes() == 0:
        raise ValueError("No network graph available in model config.")

    fig, ax = plt.subplots(figsize=(10, 8))
    if layout == "spring":
        pos = nx.spring_layout(G, seed=0)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    else:
        pos = nx.kamada_kawai_layout(G)
    weights = [d.get("weight", 1.0) for _, _, d in G.edges(data=True)]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=100, node_color="#4c72b0", alpha=0.8)
    nx.draw_networkx_edges(G, pos, ax=ax, width=[max(0.5, w * 3) for w in weights], alpha=0.4, edge_color="#888")
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=6)
    ax.set_title("Gene / Perturbation similarity network")
    ax.axis("off")
    fig.tight_layout()
    if save_file:
        fig.savefig(save_file, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return fig


def plot_all(model, perts, save_dir="gears_viz_output", top_n=20, interactive=False):
    """Generate all standard visualizations for a list of perturbations."""
    os.makedirs(save_dir, exist_ok=True)
    for p in perts:
        safe = p.replace("/", "_").replace("+", "_plus_")
        plot_pert_pred_vs_truth(model, p, top_n=top_n,
                                save_file=os.path.join(save_dir, f"{safe}_boxplot.png"), show=False)
        plot_pred_vs_truth_scatter(model, p, n_genes=top_n,
                                   save_file=os.path.join(save_dir, f"{safe}_scatter.png"),
                                   show=False, interactive=interactive)
        plot_dotplot_top_de(model, p, top_n=top_n,
                            save_file=os.path.join(save_dir, f"{safe}_dotplot.png"), show=False)
        plot_fold_change(model, p,
                         save_file=os.path.join(save_dir, f"{safe}_foldchange.png"), show=False)
    combos = [p for p in perts if "+" in p]
    for c in combos[:5]:
        safe = c.replace("/", "_").replace("+", "_plus_")
        plot_gi_synergy(model, c.split("+"),
                        save_file=os.path.join(save_dir, f"{safe}_gi.png"), show=False)
    plot_pert_network(model, save_file=os.path.join(save_dir, "network.png"), show=False)
    print(f"All visualizations saved to {save_dir}/")
