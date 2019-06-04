import pandas as pd
import numpy as np
import os
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import sklearn.decomposition
import traceback


def setup_figure(aw=4.5, ah=3, xspace=[0.75,0.25], yspace=[0.75,0.25]):
    fw = aw + np.sum(xspace)
    fh = ah + np.sum(yspace)
    fig = plt.figure(facecolor=(1,1,1), figsize=(fw,fh))
    ax = fig.add_axes([xspace[0]/fw, yspace[0]/fh, aw/fw, ah/fh])
    return ax

def metrics_plot(v, cohort_s, ylim, ms=12, alpha=1, title='', cohort_colors=None,
    date_s=None, threshold=None, threshold_dir=None, show_legend=False, show_xticklabels=False):

    dl = 0.75
    aw = 4.5
    ds = 0.1
    daw = 0.5
    dr = 0.2
    if show_xticklabels:
        db = 1.5
    else:
        db = 0.75
    ah = 3
    dt = 0.25
    fw = dl + aw + ds + daw + dr
    fh = db + ah + dt

    fig = plt.figure(facecolor=(1,1,1), figsize=(fw,fh))
    ax = fig.add_axes([dl/fw, db/fh, aw/fw, ah/fh])
    dax = fig.add_axes([(dl+aw+ds)/fw, db/fh, daw/fw, ah/fh])

    ns = len(v)

    if date_s is not None:
        date_ix = date_s.sort_values(na_position='first').index
        xlabel = 'Samples, ordered by sequencing date'
    else:
        date_ix = v.index
        xlabel = 'Samples'
    xpos = pd.Series(np.arange(1,ns+1), index=date_ix)

    for t in np.unique(cohort_s):
        ix = cohort_s[cohort_s==t].index
        ax.scatter(xpos[ix], v[ix], s=ms, edgecolor='none', label=t,
            c=[cohort_colors[t]], alpha=alpha, clip_on=False, rasterized=True)

    if threshold is not None:
        ax.plot([-0.02*ns, 1.02*ns], 2*[threshold], '--', color=[0.6,0.6,0.6], lw=1, alpha=0.8)
        if threshold_dir=='gt':
            ix = v>threshold
        elif threshold_dir=='lt':
            ix = v<threshold
        ax.scatter(xpos[np.where(ix)[0]], v[ix], c='none', edgecolor='k', s=ms, lw=1, label=None)

    try:
        sns.kdeplot(v, ax=dax, vertical=True, legend=False, shade=True)
    except np.linalg.LinAlgError:
        print(traceback.format_exc())
        print("Unable to render KDE for", title)

    ax.set_xlim([1-0.02*ns, 1.02*ns])
    if ylim is None:
        ax.set_ylim([0, ax.get_ylim()[1]])
    else:
        ax.set_ylim(ylim)

    dax.set_ylim(ax.get_ylim())
    dax.set_yticks([])

    if show_xticklabels:
        ax.set_xticks(xpos)
        ax.set_xticklabels(date_ix, rotation=45, ha='right', va='top')
    else:
        # hc.locator = ticker.MaxNLocator(nbins=5)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    ax.set_ylabel(title, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.tick_params(labelsize=12)
    dax.set_xticks([])
    dax.set_xlabel('Frequency', fontsize=12, labelpad=7)

    if show_legend:
        ax.legend(fontsize=8)
    return (fig, ax)

def plot_qc_figures(metrics_df, cohort_s=None, cohort_colors=None, date_s=None, ms=12, output_dir=None,
        intergenic_rate=None, exonic_rate=None, million_mapped_reads=None,
        million_mapped_reads_qc=None,
        mapping_rate=None, end1_mismatch_rate=None, end2_mismatch_rate=None,
        rrna_rate=None, unique_rate=None, dpi=300, name='', alpha=0.8, show_legend=False,
        show_xticklabels=False):
    """
    metrics_df: output from RNA-SeQC
    cohort_s: mapping of sample ID to cohort/cluster/etc.
    """
    figure_outputs = []
    if cohort_s is None:
        cohort_s = pd.Series('all_samples', index=metrics_df.index)

    if 'Mapped' in metrics_df and 'Mapped Reads' not in metrics_df:  # v1 data
        metrics_df = metrics_df.rename(columns={
            'Mapped':'Mapped Reads',
            'Mapped Unique':'Mapped Unique Reads',
            'rRNA rate':'rRNA Rate',
            'Duplication Rate of Mapped':'Duplicate Rate of Mapped'
        })

    cohorts = np.unique(cohort_s)
    nc = len(cohorts)
    cohort_index_dict = {s:i for i,s in enumerate(cohorts)}

    if date_s is not None:
        date_ix = date_s.sort_values(na_position='first').index
        xlabel = 'Samples, ordered by sequencing date'
    else:
        date_ix = metrics_df.index
        xlabel = 'Samples'
    xpos = pd.Series(np.arange(metrics_df.shape[0]), index=date_ix)

    ns = metrics_df.shape[0]
    c = [cohort_index_dict[t] for t in cohort_s.loc[metrics_df.loc[date_ix].index]]

    if cohort_colors is None:
        cohort_colors = {i:j for i,j in zip(cohorts, cm.get_cmap('tab10', 10)(np.arange(nc)))}

    metrics_args = {'show_xticklabels':show_xticklabels, 'show_legend':show_legend, 'cohort_colors':cohort_colors, 'ms':ms}

    v = metrics_df.loc[date_ix, 'Intergenic Rate']
    figure_outputs.append(
        metrics_plot(v, cohort_s, ylim=[0, 0.55], title='Intergenic rate', threshold=intergenic_rate, threshold_dir='gt', **metrics_args)
    )
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, 'intergenic_rate_{}.pdf'.format(name)), dpi=dpi)

    v = metrics_df.loc[date_ix, 'Exonic Rate']
    figure_outputs.append(
        metrics_plot(v, cohort_s, ylim=[0, 1], title='Exonic rate', threshold=exonic_rate, threshold_dir='lt', **metrics_args)
    )
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, 'exonic_rate_{}.pdf'.format(name)), dpi=dpi)

    v = metrics_df.loc[date_ix, 'Mapped Reads'].copy()/1e6
    v[v>250] = 250
    figure_outputs.append(
        metrics_plot(v, cohort_s, ylim=[0, 250], title='Mapped reads (millions)', threshold=million_mapped_reads, threshold_dir='lt', **metrics_args)
    )
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, 'mapped_reads_{}.pdf'.format(name)), dpi=dpi)

    v = metrics_df.loc[date_ix, 'Unique Mapping, Vendor QC Passed Reads'].copy()/1e6
    v[v>250] = 250
    figure_outputs.append(
        metrics_plot(v, cohort_s, ylim=[0, 250], title='QC-passed reads (millions)', threshold=million_mapped_reads_qc, threshold_dir='lt', **metrics_args)
    )
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, 'qc_passed_reads_{}.pdf'.format(name)), dpi=dpi)

    v = metrics_df.loc[date_ix, 'Mapping Rate'].copy()
    if np.sum(v<0.5)/len(v)<0.25:
        v[v<0.5] = 0.5
        ylim = [0.5,1.03]
    else:
        ylim = [0,1]
    figure_outputs.append(
        metrics_plot(v, cohort_s, ylim=ylim, title='Mapping Rate', threshold=mapping_rate, threshold_dir='lt', **metrics_args)
    )
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, 'mapping_rate_{}.pdf'.format(name)), dpi=dpi)


    # 4) Mismatch rate thresholds

    ax = setup_figure(3,3, xspace=[1,0.25], yspace=[0.75,0.25])
    x = metrics_df['End 1 Mismatch Rate'].copy()
    y = metrics_df['End 2 Mismatch Rate'].copy()
    x[x>0.01] = 0.01
    y[y>0.025] = 0.025
    for t in cohorts:
        ix = cohort_s[cohort_s==t].index
        ax.scatter(x[ix], y[ix], s=ms, edgecolor='none', label=t,
            c=[cohort_colors[t]], alpha=alpha, clip_on=False, rasterized=True)
    if end1_mismatch_rate is not None:
        ax.plot(2*[end1_mismatch_rate], [0,0.2], '--', color=[0.6]*3, zorder=0, lw=1, alpha=0.8)
    if end2_mismatch_rate is not None:
        ax.plot([0,0.02], 2*[end2_mismatch_rate], '--', color=[0.6]*3, zorder=0, lw=1, alpha=0.8)
    if end1_mismatch_rate is not None or end2_mismatch_rate is not None:
        ix = (x>end1_mismatch_rate) | (y>end2_mismatch_rate)
        ax.scatter(x[ix], y[ix], c='none', edgecolor='k', s=ms, lw=1, label=None)
    ax.set_xlim([0, 0.01])
    ax.set_ylim([0, 0.025])
    ax.set_xlabel('End 1 mismatch rate', fontsize=16)
    ax.set_ylabel('End 2 mismatch rate', fontsize=16)
    ax.tick_params(labelsize=12)
    if show_legend:
        ax.legend(fontsize=8)
    figure_outputs.append((ax.figure, ax))
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, 'end_mismatch_rate_{}.pdf'.format(name)), dpi=dpi)


    v = metrics_df.loc[date_ix, 'rRNA Rate'].copy()
    figure_outputs.append(
        metrics_plot(v, cohort_s, ylim=[0, 0.4], title='rRNA Rate', threshold=rrna_rate, threshold_dir='gt', **metrics_args)
    )
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, 'rRNA_rate_{}.pdf'.format(name)), dpi=dpi)

    v = metrics_df.loc[date_ix, 'Duplicate Rate of Mapped']
    if not np.all(v==0):
        figure_outputs.append(
            metrics_plot(v, cohort_s, ylim=[0, 1], title='Duplicate rate (of mapped)', **metrics_args)
        )
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, 'duplication_rate_{}.pdf'.format(name)), dpi=dpi)

        v = metrics_df.loc[date_ix, 'Unique Rate of Mapped']
        figure_outputs.append(
            metrics_plot(v, cohort_s, ylim=[0, 1], title='Unique rate of mapped', **metrics_args)
        )
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, 'unique_rate_{}.pdf'.format(name)), dpi=dpi)

    if "Median 3' bias" in metrics_df.columns:
        v = metrics_df.loc[date_ix, "Median 3' bias"]
        figure_outputs.append(
            metrics_plot(v, cohort_s, ylim=[0, 1], title="Median 3' bias", **metrics_args)
        )
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, 'bias3p_{}.pdf'.format(name)), dpi=dpi)
    return figure_outputs


def normalize_counts(gct_df, C=None):

    gct_norm_df = gct_df / deseq2_size_factors(gct_df)
    gct_norm_df = np.log10(1+gct_norm_df)

    # threshold low expressed genes
    mask = np.mean(gct_norm_df > 1, axis=1) > 0.1  # >=10 counts in >10% of samples
    gct_norm_df = gct_norm_df[mask]

    if C is not None:
        gct_norm_df = remove_covariates(gct_norm_df, C, center=False)

    # gct_norm_std_df = center_normalize(gct_norm_df)
    gct_norm_std_df = gct_norm_df - gct_norm_df.mean(axis=0)
    gct_norm_std_df = gct_norm_std_df / np.sqrt(gct_norm_std_df.pow(2).sum(axis=0))

    return gct_norm_std_df

def get_pcs(gct_df, normalize=True, C=None, n_components=5):
    """
    Scale input GCT, threshold, normalize and calculate PCs
    """
    if normalize:
        gct_norm_std_df = normalize_counts(gct_df, C=C)
    else:
        gct_norm_std_df = gct_df

    pca = sklearn.decomposition.PCA(n_components=n_components)
    pca.fit(gct_norm_std_df.T)
    P = pca.transform(gct_norm_std_df.T)
    P_df = pd.DataFrame(P, index=gct_norm_std_df.columns)

    return P_df, pca

def plot_pca(P_df, pca, c=None, cohort_s=None, cohort_colors=None, cohort_args=None, order=[1,2,3], outliers=None, title='',
    vmin=None, vmax=None, alpha=1, lw=0, s=30, cmap=plt.cm.Spectral_r, cticks=None, cticklabels=None, clabel='',
    show_legend=True, show_ax2=True):
    """
    cohort_s: Series encoding cohorts
    cohort_colors: dict

    Modes:
    """
    if cohort_s is not None:
        cohorts = cohort_s.unique()
        nc = len(cohorts)
        if cohort_colors is None and cohort_args is None:
            # cohort_colors = {i:j for i,j in zip(cohorts, cm.get_cmap(cmap, nc)(np.arange(nc)))}
            cohort_colors = {i:j for i,j in zip(cohorts, sns.husl_palette(nc, s=1, l=0.6))}
        if cohort_args is None:
            cohort_args = {}
            for k in np.unique(cohort_s):
                cohort_args[k] = {'color': cohort_colors[k], 'marker':'o', 'edgecolor':'none', 's':s}

    if show_ax2:
        fig = plt.figure(facecolor=(1,1,1), figsize=(10.5,5.5))
        ax1 = fig.add_axes(np.array([1/10.5, 0.75/5.5, 4/10.5, 4/5.5]))
    else:
        fig = plt.figure(facecolor=(1,1,1), figsize=(5.5,5.5))
        ax1 = fig.add_axes(np.array([1/5.5, 0.75/5.5, 4/5.5, 4/5.5]))
    if cohort_s is None:  # c[P_df.index]
        sa = ax1.scatter(P_df[order[1]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
    else:
        for k in np.unique(cohort_s):
        # for k in cohort_s.unique():
            i = cohort_s[cohort_s==k].index
            ax1.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])
    format_plot(ax1, fontsize=10)
    ax1.set_xlabel('PC {0} ({1:.2f}%)'.format(order[1], pca.explained_variance_ratio_[order[1]-1]*100), fontsize=12)
    ax1.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if show_ax2:
        ax2 = fig.add_axes(np.array([6/10.5, 0.75/5.5, 4/10.5, 4/5.5]))
        if cohort_s is None:
            ax2.scatter(P_df[order[2]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
        else:
            for k in np.unique(cohort_s):
                i = cohort_s[cohort_s==k].index
                ax2.scatter(P_df.loc[i,order[2]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])
            # ax2.legend(loc=3, fontsize=10, scatterpoints=1, handletextpad=0.1, framealpha=0.5, bbox_to_anchor=(-0.5,-0.1))

        format_plot(ax2, fontsize=10)
        ax2.set_xlabel('PC {0} ({1:.2f}%)'.format(order[2], pca.explained_variance_ratio_[order[2]-1]*100), fontsize=12)
        ax2.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if outliers is not None:
        ax1.scatter(P_df.loc[outliers, order[1]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)
        if show_ax2:
            ax2.scatter(P_df.loc[outliers, order[2]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)

    fig.suptitle(title, fontsize=12)

    if cohort_s is not None and show_legend:
        # ax2.legend(loc=0, fontsize=10, scatterpoints=1, handletextpad=0.1, framealpha=0.5, bbox_to_anchor=(-0.5,-0.1))
        leg = ax1.legend(loc=0, fontsize=9, scatterpoints=1, handletextpad=0.1, framealpha=1, labelspacing=0.35)
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    # if cohort_s is None and c is not None and not isinstance(c, list) and not isinstance(c, str):
    if cohort_s is None and c is not None and len(c)==P_df.shape[0]:
        if show_ax2:
            cax = fig.add_axes(np.array([3.5/10.5, 5/5.5, 1.5/10.5, 0.15/5.5]))
        else:
            cax = fig.add_axes(np.array([3.5/5.5, 5/5.5, 1.5/5.5, 0.15/5.5]))
        # cax = fig.add_axes(np.array([3.5/10.5, 4.85/5.5, 1.5/10.5, 0.15/5.5]))
        hc = plt.colorbar(sa, cax=cax, orientation='horizontal')
        if cticks is not None:
            hc.set_ticks(cticks)
        if cticklabels is not None:
            # hc.set_ticks([0,0.5,1])
            hc.ax.tick_params(labelsize=9)
            # cax.invert_xaxis()
            cax.set_xticklabels(cticklabels, fontsize=10)

        hc.locator = ticker.MaxNLocator(integer=True, min_n_ticks=2, nbins=5)
        hc.update_ticks()

        cax.set_ylabel(clabel, rotation=0, ha='right', va='center', fontsize=12)
    return fig

def deseq2_size_factors(counts_df):
    """
    Calculate DESeq size factors
    median of ratio to reference sample (geometric mean of all samples)

    References:
     [1] Anders & Huber, 2010
     [2] R functions:
          DESeq::estimateSizeFactorsForMatrix
    """
    idx = np.all(counts_df>0, axis=1)
    tmp_df = np.log(counts_df[idx])
    s = np.exp(np.median(tmp_df.T - np.mean(tmp_df, axis=1), axis=1))
    return s

def remove_covariates(df, C, center=False, fail_colinear=False):
    """
    Residualizes rows of M relative to columns of C
    """

    # transform input
    if isinstance(df, pd.DataFrame) or isinstance(df, pd.Series):
        M = df.values
    else:
        M = df

    isvector = False
    if isinstance(M, list) or (hasattr(M, 'shape') and len(M.shape)==1):
        M = np.array(M).reshape(1,-1)
        isvector = True

    if isinstance(C, list) or (hasattr(C, 'shape') and len(C.shape)==1):
        C = np.array(C).reshape(-1,1)

    Q = orthogonalize_covariates(C, fail_colinear=fail_colinear)

    # residualize M relative to C
    M0 = (M.T - np.mean(M,axis=1)).T
    if center:
        M0 = M0 - np.dot(np.dot(M0, Q), Q.T)
    else:
        M0 = M - np.dot(np.dot(M0, Q), Q.T)  # retain original mean

    if isvector:
        M0 = M0[0]

    if isinstance(df, pd.DataFrame):
        M0 = pd.DataFrame(M0, index=df.index, columns=df.columns)
    elif isinstance(df, pd.Series):
        M0 = pd.Series(M0, index=df.index, name=df.name)

    return M0

def orthogonalize_covariates(C, fail_colinear=True):
    """
    C: covariates (columns)
    """
    # center and orthogonalize
    Q,R = np.linalg.qr(C-np.mean(C,axis=0))

    # check for colinearity
    colinear_ix = np.abs(np.diag(R)) < np.finfo(np.float64).eps * C.shape[1]
    if np.any(colinear_ix):
        if fail_colinear:
        # if np.min(np.abs(np.diag(R))) < np.finfo(np.float64).eps * C.shape[1]:
            raise ValueError("Colinear or zero covariates detected")
        else:  # drop colinear covariates
            print('  * Colinear covariates detected. {} covariates dropped.'.format(np.sum(colinear_ix)))
            Q = Q[:, ~colinear_ix]

    return Q

def format_plot(ax, tick_direction='out', tick_length=4, hide=['top', 'right'], hide_spines=True, lw=1, fontsize=9):

    for i in ['left', 'bottom', 'right', 'top']:
        ax.spines[i].set_linewidth(lw)

    # ax.axis["left"].major_ticklabels.set_ha("left")
    ax.tick_params(axis='both', which='both', direction=tick_direction, labelsize=fontsize)

    # set tick positions
    if 'top' in hide and 'bottom' in hide:
        ax.get_xaxis().set_ticks_position('none')
    elif 'top' in hide:
        ax.get_xaxis().set_ticks_position('bottom')
    elif 'bottom' in hide:
        ax.get_xaxis().set_ticks_position('top')
    else:
        ax.get_xaxis().set_ticks_position('both')

    if 'left' in hide and 'right' in hide:
        ax.get_yaxis().set_ticks_position('none')
    elif 'left' in hide:
        ax.get_yaxis().set_ticks_position('right')
    elif 'right' in hide:
        ax.get_yaxis().set_ticks_position('left')
    else:
        ax.get_yaxis().set_ticks_position('both')

    if hide_spines:
        for i in hide:
            ax.spines[i].set_visible(False)


    # adjust tick size
    for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
    #for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(tick_length) # tick length
        line.set_markeredgewidth(lw) # tick line width

    for line in (ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True)):
        line.set_markersize(tick_length/2) # tick length
        line.set_markeredgewidth(lw/2) # tick line width
