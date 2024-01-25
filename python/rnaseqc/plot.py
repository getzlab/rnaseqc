import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm, ListedColormap, hsv_to_rgb
import seaborn as sns
import os
import qtl.plot


def get_cohort_colors(cohorts):
    nc = len(cohorts)
    if nc > 5:
        cohort_colors = {i:j for i,j in zip(cohorts, plt.cm.get_cmap('Spectral', nc)(np.random.permutation(np.arange(nc)))[:,:-1])}
    else:
        cohort_colors = {i:j for i,j in zip(cohorts, plt.cm.get_cmap('tab10', 10)(np.arange(nc))[:,:-1])}
    return cohort_colors


def sort_samples(sample_ix, cohort_s=None, cohort_order=None, date_s=None):
    """Sort samples by date and cohort label"""
    if cohort_s is None and date_s is None:
        return sample_ix

    if cohort_s is not None:
        assert sample_ix.isin(cohort_s.index).all()
        cohort_s = cohort_s[sample_ix]
    if date_s is not None:
        assert sample_ix.isin(cohort_s.index).all()
        date_s = date_s[sample_ix]

    if date_s is not None:
        if cohort_s is not None:  # sort samples by date and cohort
            sorted_ix = pd.concat([
                pd.to_datetime(date_s).rename('date'), cohort_s.rename('cohort')], axis=1
            ).sort_values(['date', 'cohort'], na_position='first').index
        else:
            sorted_ix = pd.to_datetime(date_s).sort_values(na_position='first').index
    else:  # sort by cohort only
        if cohort_order is None:
            sorted_ix = cohort_s.sort_values(na_position='first').index
        else:
            sorted_ix = cohort_s.map({j:i for i,j in enumerate(cohort_order)}).sort_values(na_position='first').index

    return sorted_ix


def mismatch_rates(metrics_df, cohort_s=None, cohort_order=None, cohort_colors=None, ms=12, alpha=1, aw=2,
                   end1_threshold=None, end2_threshold=None,
                   end1_limit=0.01, end2_limit=0.025):
    """Plot base mismatch rates ('NM' tag) for read mate 1 vs read mate 2."""

    if cohort_s is not None:
        assert metrics_df.index.isin(cohort_s.index).all()
        cohort_s = cohort_s.loc[metrics_df.index]
    else:
        cohort_s = pd.Series('NA', index=metrics_df.index)

    ax = qtl.plot.setup_figure(aw, aw)

    x = metrics_df['End 1 Mismatch Rate'].copy()
    y = metrics_df['End 2 Mismatch Rate'].copy()
    x[x>end1_limit] = end1_limit
    y[y>end2_limit] = end2_limit

    sorted_ix = sort_samples(metrics_df.index, cohort_s=cohort_s, cohort_order=cohort_order)
    cohorts = cohort_s.loc[sorted_ix].unique()
    if cohort_colors is None:
        cohort_colors = get_cohort_colors(cohorts)

    for t in cohorts:
        ix = cohort_s[cohort_s==t].index
        ax.scatter(x[ix], y[ix], s=ms, edgecolor='none', label=t,
            c=[cohort_colors[t]], alpha=alpha, clip_on=False, rasterized=True)

    if end1_threshold is not None:
        ax.plot(2*[end1_threshold], [0,0.2], '--',  color=[0.6]*3, zorder=0, lw=1, alpha=0.8)
    if end2_threshold is not None:
        ax.plot([0,0.02], 2*[end2_threshold], '--', color=[0.6]*3, zorder=0, lw=1, alpha=0.8)
    if end1_threshold is not None or end2_threshold is not None:
        ix = (x > end1_threshold) | (y > end2_threshold)
        ax.scatter(x[ix], y[ix], c='none', edgecolor='k', s=ms, lw=1, label=None, clip_on=False, rasterized=True)

    qtl.plot.format_plot(ax, fontsize=10)
    ax.set_xlim([0, end1_limit])
    ax.set_ylim([0, end2_limit])

    ax.spines['left'].set_position(('outward', 6))
    ax.spines['bottom'].set_position(('outward', 6))
    ax.plot([0, end1_limit], [0, end1_limit], '--', c=[0.6]*3, lw=1, zorder=0)

    ax.set_xlabel('End 1 mismatch rate', fontsize=12)
    ax.set_ylabel('End 2 mismatch rate', fontsize=12)


def metrics(metric_s, cohort_s=None, cohort_order=None, cohort_colors=None, date_s=None,
            threshold=None, threshold_dir=None, outlier_method='threshold', plot_density=True, show_legend=False,
            ms=12, alpha=1, ylim=None, ylabel=None,
            show_xticklabels=False, highlight_ids=None,
            dl=0.85, aw=6, ds=0.2, daw=0.5, dr=0.25,
            db=0.75, ah=2, dt=0.25, fontsize=10, rasterized=True):
    """Plot a single QC metric sorted by cohort and/or date"""

    if ylabel is None:
        ylabel = metric_s.name

    if metric_s.median() > 1e5:
        metric_s = metric_s.copy() / 1e6
        if threshold is not None:
            threshold = threshold / 1e6
        ylabel += ' (millions)'

    if cohort_s is not None:
        assert metric_s.index.isin(cohort_s.index).all()
        cohort_s = cohort_s.loc[metric_s.index]
    else:
        cohort_s = pd.Series('NA', index=metric_s.index)

    if show_xticklabels:
        db += 0.75

    if plot_density:
        fw = dl + aw + ds + daw + dr
    else:
        fw = dl + aw + dr
    fh = db + ah + dt
    fig = plt.figure(facecolor=(1,1,1), figsize=(fw,fh))
    ax = fig.add_axes([dl/fw, db/fh, aw/fw, ah/fh])
    if plot_density:
        dax = fig.add_axes([(dl+aw+ds)/fw, db/fh, daw/fw, ah/fh], sharey=ax)

    if date_s is not None:
        xlabel = 'Samples, ordered by date'
    else:
        xlabel = 'Samples'

    sorted_ix = sort_samples(metric_s.index, cohort_s=cohort_s, cohort_order=cohort_order, date_s=date_s)
    cohorts = cohort_s.loc[sorted_ix].unique()
    if cohort_colors is None:
        cohort_colors = get_cohort_colors(cohorts)

    ns = len(metric_s)
    xpos = pd.Series(np.arange(1,ns+1), index=sorted_ix)

    # plot
    for t in cohorts:
        ix = cohort_s[cohort_s==t].index
        ax.scatter(xpos[ix], metric_s[ix], s=ms, edgecolor='none', label=t,
                   c=[cohort_colors[t]], alpha=alpha, clip_on=False, rasterized=rasterized)

    if highlight_ids is not None:
        ax.scatter(xpos[highlight_ids], metric_s[highlight_ids], marker='s',
                   edgecolor='k', facecolor='none', clip_on=False, rasterized=rasterized)

    if threshold is not None:
        ax.plot([-0.02*ns, 1.02*ns], 2*[threshold], '--', color=[0.6,0.6,0.6], lw=1, alpha=0.8)

    if outlier_method.lower() == 'iqr':
        p = np.percentile(metric_s, [25, 75])
        if threshold_dir == 'gt':
            ix = metric_s[metric_s > p[1] + 1.5*(p[1]-p[0])].index
        elif threshold_dir == 'lt':
            ix = metric_s[metric_s < p[0] - 1.5*(p[1]-p[0])].index
        if any(ix):
            ax.scatter(xpos[ix], metric_s[ix], c='none', edgecolor='k', s=ms, lw=1, label=None, clip_on=False, rasterized=rasterized)
    elif outlier_method.lower() == 'threshold' and threshold is not None:
        if threshold_dir == 'gt':
            ix = metric_s[metric_s > threshold].index
        elif threshold_dir == 'lt':
            ix = metric_s[metric_s < threshold].index
        if any(ix):
            ax.scatter(xpos[ix], metric_s[ix], c='none', edgecolor='k', s=ms, lw=1, label=None, clip_on=False, rasterized=rasterized)

    # plot density
    if plot_density:
        sns.kdeplot(y=metric_s, ax=dax, legend=False, fill=True, lw=1.5)
        dax.set_ylabel(None)
        qtl.plot.format_plot(dax, fontsize=fontsize, hide=['top', 'right', 'bottom'])
        plt.setp(dax.get_yticklabels(), visible=False)
        dax.set_xticks([])
        dax.set_xlabel('Freq.', ha='left', x=0, fontsize=fontsize, labelpad=7)

    qtl.plot.format_plot(ax, fontsize=fontsize)
    ax.spines['left'].set_position(('outward', 8))

    ax.set_xlim([1, ns])
    if ylim is None:
        ax.set_ylim([0, ax.get_ylim()[1]])
    else:
        ax.set_ylim(ylim)

    if show_xticklabels:
        ax.set_xticks(xpos)
        ax.set_xticklabels(sorted_ix, rotation=45, ha='right', va='top')
    else:
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    ax.set_ylabel(ylabel, fontsize=fontsize+2)
    ax.set_xlabel(xlabel, fontsize=fontsize+2)

    if show_legend:
        ax.legend(fontsize=9, handlelength=1, labelspacing=0.5, title=cohort_s.name)

    if plot_density:
        return ax, dax
    else:
        return ax


def detection_bias(metrics_df, bias_metric="Median 3' bias", c='Duplicate Rate of Mapped',
                   ah=2, aw=2, ct=0, rasterized=False):
    """Plot genes detected vs a bias metric (e.g., Median Exon CV)"""

    ax, cax = qtl.plot.setup_figure(ah, aw, xspace=[0.75, 0.75],
                                    colorbar=True, ds=0.05, cw=0.1, ct=ct)

    ix = metrics_df[c].sort_values().index
    h = ax.scatter(metrics_df.loc[ix, 'Genes Detected'], metrics_df.loc[ix, bias_metric],
                   c=metrics_df.loc[ix, c], cmap=plt.cm.GnBu,
                   clip_on=False, s=36, edgecolor='k',lw=0.5,
                   vmin=0, vmax=1, rasterized=rasterized)

    ax.set_xlabel('Genes detected', fontsize=12)
    ax.set_ylabel(bias_metric, fontsize=12)
    qtl.plot.format_plot(ax, fontsize=10)
    ax.autoscale(True)
    ax.spines['left'].set_position(('outward', 6))
    ax.spines['bottom'].set_position(('outward', 6))
    hc = plt.colorbar(h, cax=cax)
    hc.set_label('Duplicate Rate', fontsize=12, labelpad=6)
    return ax, cax


def mapping_sense(metrics_df, cohort_s=None, cohort_order=None, cohort_colors=None, date_s=None, width=0.8,
                  dl=0.75, aw=4, dr=1.5, db=0.5, ah=2, dt=0.25, ds=0.066, dc=0.1):
    """Summary of sense/antisense alignments.

    For stranded protocols, most reads should be 'End 1 Antisense' and 'End 2 Sense',
    or vice versa, depending on protocol.
    For unstranded protocols, the 4 categories are expected to be of equal proportion (~0.25).
    """
    sorted_ix = sort_samples(metrics_df.index, cohort_s=cohort_s,
                             cohort_order=cohort_order, date_s=date_s)
    df = metrics_df.loc[sorted_ix, ['End 1 Sense', 'End 1 Antisense', 'End 2 Sense', 'End 2 Antisense']]
    df = df / np.sum(df.values, axis=1, keepdims=True)

    fw = dl + aw + dr
    fh = db + ah + dt
    fig = plt.figure(facecolor=(1,1,1), figsize=(fw,fh))
    ax = fig.add_axes([dl/fw, db/fh, aw/fw, ah/fh])
    df.reset_index(drop=True).plot(kind='bar', width=width,
                                   stacked=True, xticks=[], ax=ax,
                                   color=hsv_to_rgb([
                                       [0.1, 0.6, 1],
                                       [0.4, 0.7, 0.75],
                                       [0.25, 0.4, 0.85],
                                       [0.15, 0.55, 1],
                                   ])
                                   )
    ax.set_ylim([0,1])
    ax.legend(loc='upper left', handlelength=0.66, bbox_to_anchor=(1,1))
    ax.set_ylabel('Proportion of mapped reads', fontsize=12)
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_xlim([-width/2, metrics_df.shape[0]-width/2])

    if cohort_s is not None:
        cax = fig.add_axes([dl/fw, (db+ah+ds)/fh, aw/fw, dc/fh], sharex=ax)
        cax.set_yticks([])
        _plot_cohort_labels(cax, cohort_s.loc[df.index], cohort_colors=cohort_colors,
                            lax=ax, legend=False, orientation='horizontal')


def calculate_expression_cdfs(tpm_df):
    """Sort and compute CDF for each sample independently"""
    cdf_df = tpm_df.reset_index(drop=True).copy()
    if 'Description' in cdf_df:
        cdf_df.drop('Description', axis=1, inplace=True)
    for c in cdf_df:
        cdf_df[c] = np.cumsum(cdf_df[c].sort_values(ascending=False).values) / 1e6
    return cdf_df


def cumulative_expression(cdf_df, cohort_s=None, cohort_colors=None, ax=None, cmap=plt.cm.Spectral_r, c=[0.6,0.6,0.6],
                          reference_df=None, reference_name=None, alpha=0.5, mode='lines', lw=1, legend=False, rasterized=False):
    """
    Plot cumulative gene expression for each sample.
    This enables identification of samples with dominant expression of few genes.

    With mode='ci', median and confidence intervals are shown instead of individual samples.
    """
    if cohort_s is None:
        cohort_s = pd.Series('_NA', index=cdf_df.columns)

    if cohort_colors is None:
        cohorts = cohort_s.unique()
        nc = len(cohorts)
        if nc==1:
            cohort_colors = {cohorts[0]: c}
        else:
            cohort_colors = {i:j for i,j in zip(cohorts, plt.cm.get_cmap(cmap.name, nc)(np.arange(nc)))}

    if ax is None:
        ax = qtl.plot.setup_figure(4, 2.5)
    ax.set_xscale('log')

    if reference_df is not None:  # plot reference distribution
        # mu = reference_df.mean(axis=1)
        s = reference_df.std(axis=1)
        mu = reference_df.median(axis=1)
        # s = np.median(np.abs(gtex_cdf-mu), axis=0) / 0.6745
        if mode=='ci':
            ax.fill_between(np.arange(reference_df.shape[0])+1, mu-1.96*s, mu+1.96*s, facecolor='k', edgecolor='k', alpha=0.2, label=None, zorder=20)
            ax.plot(mu, 'k', lw=2, alpha=0.8, rasterized=rasterized, label=reference_name, zorder=30)
        else:
            ax.fill_between(np.arange(reference_df.shape[0])+1, mu-1.96*s, mu+1.96*s, facecolor='k', edgecolor='none', alpha=0.2, label=f'{reference_name} 95% CI', zorder=20)
            ax.plot(mu, 'k', lw=1.5, alpha=0.8, rasterized=rasterized, label=f'{reference_name} mean', zorder=30)

    for c in cohort_s.unique():
        x = np.arange(1, cdf_df.shape[0]+1)
        ix = cohort_s[cohort_s==c].index
        if mode == 'ci':  # plot confidence intervals
            mu = cdf_df[ix].median(axis=1)
            s = cdf_df[ix].std(axis=1)  # replace with MAD?
            fc = cohort_colors[c]
            pc = ax.fill_between(x, mu-1.96*s, mu+1.96*s, facecolor=fc, edgecolor=fc, alpha=0.2, label=None, zorder=20, lw=1)
            ax.plot(mu, '-', color=cohort_colors[c], lw=2, alpha=0.8, rasterized=rasterized, label=c, zorder=30)
        else:
            ax.plot(x, cdf_df[ix[0]], color=cohort_colors[c], alpha=alpha, lw=lw, rasterized=rasterized, label=c)  # plot first one w/ label
            if len(ix)>1:
                ax.plot(x, cdf_df[ix[1:]], color=cohort_colors[c], alpha=alpha, lw=lw, rasterized=rasterized)

    ax.set_ylim([0,1])
    ax.set_xlim([1,10000])
    qtl.plot.format_plot(ax, fontsize=10)
    ax.set_xlabel('Number of genes', fontsize=12)
    ax.set_ylabel('Cumulative transcriptional output', fontsize=12)
    ax.spines['left'].set_position(('outward', 6))

    if legend and not (cohort_s == '_NA').all():
        leg = ax.legend(loc=4, handlelength=1, fontsize=10)
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    return ax


def _plot_cohort_labels(ax, cohort_s, cohort_colors=None, lax=None, legend=True, orientation='vertical'):
    """Internal function for adding a cohort color legend to a figure (in a separate axis)"""

    cohort_index_dict = {i:k for k,i in enumerate(np.unique(cohort_s))}
    if cohort_colors is None:
        n = len(cohort_index_dict)
        cmap = ListedColormap(plt.cm.get_cmap('Spectral', n)(np.arange(n)), 'indexed')
    else:
        cmap = ListedColormap(pd.Series(cohort_index_dict).sort_values().index.map(cohort_colors))

    if orientation == 'vertical':
        ax.imshow(cohort_s.map(cohort_index_dict).values.reshape(-1,1), aspect='auto', cmap=cmap)
    else:
        ax.imshow(cohort_s.map(cohort_index_dict).values.reshape(1,-1), aspect='auto', cmap=cmap)

    if lax is None:
        lax = ax
    for k,i in cohort_index_dict.items():
        lax.scatter([], [], marker='s', c=[cmap(i)], label=f'{k}')
    if legend:
        lax.legend(loc='upper left', borderaxespad=None, bbox_to_anchor=(1,1), handlelength=1, title='Cohort')


def insert_sizes(insertsize_df, cohort_s=None, cohort_colors=None,
                 cohort_order=None, sort_order='mean', max_size=1000,
                 legend=False, dl=0.75, aw=3, dr=0.5, db=0.5, ah=2, dt=0.25):
    """Plot heat map of insert size distributions"""

    # expand to 'max_size' bp
    df = insertsize_df.reindex(np.arange(1,max_size+1)).fillna(0).astype(np.int32).T

    # sort by mean if > 100000 reads
    mu = df.mul(df.columns.values, axis=1).sum(1)
    n = df.sum(1).sort_values()
    si = n[n<100000].index.tolist() + mu.loc[n[n>=100000].index].sort_values().index.tolist()

    if cohort_s is not None and sort_order == 'cohort':  # sort within each cohort
        if cohort_order is None:
            cohort_order = cohort_s.value_counts().index
        sort_s = pd.Series(cohort_s[si], index=si)
        si = []
        for c in cohort_order:
            si.extend(sort_s[sort_s==c].index)

    # set up figure
    if cohort_s is not None:
        cw = 0.1
        ds = 0.05
    else:
        cw = 0
        ds = 0
    fw = dl + cw + ds + aw + dr
    fh = db + ah + dt
    fig = plt.figure(facecolor=(1,1,1), figsize=(fw,fh))
    ax = fig.add_axes([(dl+cw+ds)/fw, db/fh, aw/fw, ah/fh])

    # add cohort information and legend
    if cohort_s is not None:
        # set up axes
        cax = fig.add_axes([dl/fw, db/fh, cw/fw, ah/fh], sharey=ax)
        plt.setp(ax.get_yticklabels(), visible=False);
        for line in ax.yaxis.get_ticklines():
            line.set_markersize(0)
            line.set_markeredgewidth(0)
        cax.set_xticks([])
        cax.set_ylabel('Sample', fontsize=12)

        # plot labels
        _plot_cohort_labels(cax, cohort_s[si], cohort_colors=cohort_colors, lax=ax, legend=legend)

    ax.imshow(df.loc[si], interpolation='none', aspect='auto', norm=LogNorm())
    ax.set_xlabel('Insert size (bp)', fontsize=12)
    ax.set_xlim([1, max_size])
    return ax


def gc_content(gc_content_df, cohort_s=None, cohort_colors=None,
               cohort_order=None, sort_order='mean', legend=False,
               dl=0.75, aw=3, dr=0.5, db=0.5, ah=2, dt=0.25):
    """Plot heat map of GC content distributions"""

    # sort by mean
    x = gc_content_df.index.values
    mu = (gc_content_df * x.reshape(-1,1)).sum()
    si = mu.sort_values().index

    if cohort_s is not None and sort_order == 'cohort':  # sort within each cohort
        if cohort_order is None:
            cohort_order = cohort_s.value_counts().index
        sort_s = pd.Series(cohort_s[si], index=si)
        si = []
        for c in cohort_order:
            si.extend(sort_s[sort_s==c].index)

    # set up figure
    if cohort_s is not None:
        ch = 0.1
        ds = 0.05
    else:
        ch = 0
        ds = 0
    fw = dl + aw + dr
    fh = db + ah + ch + ds + dt
    fig = plt.figure(facecolor=(1,1,1), figsize=(fw,fh))
    ax = fig.add_axes([dl/fw, db/fh, aw/fw, ah/fh])

    # add cohort information and legend
    if cohort_s is not None:
        # set up axes
        cax = fig.add_axes([dl/fw, (db+ah+ds)/fh, aw/fw, ch/fh], sharex=ax)
        plt.setp(cax.get_xticklabels(), visible=False);
        for line in cax.xaxis.get_ticklines():
            line.set_markersize(0)
            line.set_markeredgewidth(0)
        cax.set_yticks([])

        # plot labels
        _plot_cohort_labels(cax, cohort_s[si], orientation='horizontal',
                            cohort_colors=cohort_colors, lax=ax, legend=legend)

    ax.imshow(gc_content_df[si], origin='lower', interpolation='none', aspect='auto', norm=LogNorm())
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Fragment GC Content', fontsize=12)
    y = np.arange(0, 120, 20)
    ax.set_yticks(y)
    ax.set_yticklabels(y/100)
    return ax


def xy_expression(tpm_df, sex_s=None, flag_klinefelter=True, highlight_ids=None,
                  x_threshold=5, y_threshold=30, s=24, verbose=True, rasterized=False, **kwargs):
    """Expression of sex-specific genes (XIST and RPS4Y1) to identify sample swaps.

    sex_s: pd.Series annotating the sex of each sample, as Male/Female.
    """

    x_id = tpm_df.index[tpm_df.index.str.startswith('ENSG00000229807')][0]  # XIST
    y_id = tpm_df.index[tpm_df.index.str.startswith('ENSG00000129824')][0]  # RPS4Y1
    x_s = tpm_df.loc[x_id].rename('XIST')
    y_s = tpm_df.loc[y_id].rename('RPS4Y1')

    ax = qtl.plot.setup_figure(3, 3, xspace=[0.75, 1.75])
    ax.set_xscale('symlog')
    ax.set_yscale('symlog')

    if sex_s is not None:  # flag potential swaps based on thresholds
        assert tpm_df.columns.isin(sex_s.index).all()
        res_s = pd.Series('NA', index=sex_s.index[sex_s.index.isin(tpm_df.columns)], name='inferred_sex')

        args =  {'ec':'none', 'lw':0, 'rasterized':rasterized, 'clip_on':False, 's':s, 'alpha':0.33}
        args.update(kwargs)
        args2 = {**args, 'ec':'k', 'lw':1, 's':s+6, 'alpha':1}

        # infer missing labels based on thresholds
        ix = sex_s[sex_s.isnull() & (x_s <= x_threshold) & (y_s > y_threshold)].index
        if len(ix) > 0:
            ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0.6,0.8,0.7]).reshape(1,-1), **args, label=f"Male* ({len(ix)})")
            res_s[ix] = 'Male'
        ix = sex_s[sex_s.isnull() & (x_s > x_threshold) & (y_s <= y_threshold)].index
        if len(ix) > 0:
            ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0,0.8,0.7]).reshape(1,-1), **args, label=f"Female* ({len(ix)})")
            res_s[ix] = 'Female'
        ix = sex_s[sex_s.isnull() & (x_s > x_threshold) & (y_s > y_threshold)].index
        if len(ix) > 0:
            ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0.75,0.8,0.7]).reshape(1,-1), **args, label=f"XXY* ({len(ix)})")
            res_s[ix] = 'Klinefelter (XXY)'
        ix = sex_s[sex_s.isnull() & (x_s <= x_threshold) & (y_s <= y_threshold)].index
        if len(ix) > 0:
            ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0.0,0,0.7]).reshape(1,-1), **args, label=f"? ({len(ix)})")
            res_s[ix] = np.NaN

        # matching samples
        ix = sex_s[(sex_s == 'Male') & (x_s <= x_threshold)].index
        if len(ix) > 0:
            res_s[ix] = 'Male'
            ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0.6,0.8,0.7]).reshape(1,-1), label=f"Male ({len(ix)})", **args)
        ix = sex_s[(sex_s == 'Female') & (y_s <= y_threshold)].index
        if len(ix) > 0:
            res_s[ix] = 'Female'
            ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0,0.8,0.7]).reshape(1,-1), label=f"Female ({len(ix)})", **args)

        # mismatches
        if flag_klinefelter:
            ix = sex_s[(sex_s == 'Male') & (x_s > x_threshold) & (y_s <= y_threshold)].index
            if len(ix) > 0:
                ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0.6,1,0.9]).reshape(1,-1), label=f'M > F swap ({len(ix)})', **args2)
                if verbose:
                    print(f'F mislabeled as M:\n{ix.tolist()}')
                res_s[ix] = 'Female'
            ix = sex_s[(sex_s == 'Female') & (y_s > y_threshold) & (x_s <= x_threshold)].index
            if len(ix) > 0:
                ax.scatter(x_s[ix], y_s[ix], c=[[0.9, 0, 0, 1]], label=f'F > M swap ({len(ix)})', **args2)
                if verbose:
                    print(f'M mislabeled as F:\n{ix.tolist()}')
                res_s[ix] = 'Male'

            # Klinefelter
            ix = sex_s[(sex_s == 'Male') & (x_s > x_threshold) & (y_s > y_threshold)].index
            if len(ix) > 0:
                ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0.75,1,0.9]).reshape(1,-1), label=f'XXY? ({len(ix)})', **args2)
                if verbose:
                    print(f'Possible Klinefelter (XXY): {ix.tolist()}')
                res_s[ix] = 'Possible Klinefelter (XXY)'
            ix = sex_s[(sex_s == 'Female') & (y_s > y_threshold) & (x_s > x_threshold)].index
            if len(ix) > 0:
                ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0,1,0.9]).reshape(1,-1), label=f'XXY? ({len(ix)})', **args2)
                if verbose:
                    print(f'Possible Klinefelter (XXY): {ix.tolist()}')
                res_s[ix] = 'Possible Klinefelter (XXY)'

        else:
            ix = sex_s[(sex_s == 'Male') & (x_s > x_threshold)].index
            if len(ix) > 0:
                ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0.6,1,0.9]).reshape(1,-1), label=f'M > F swap ({len(ix)})', **args2)
                if verbose:
                    print(f'F mislabeled as M:\n{ix.tolist()}')
                res_s[ix] = 'Female'
            ix = sex_s[(sex_s == 'Female') & (y_s > y_threshold)].index
            if len(ix) > 0:
                ax.scatter(x_s[ix], y_s[ix], c=hsv_to_rgb([0,1,0.9]).reshape(1,-1), label=f'F > M swap ({len(ix)})', **args2)
                if verbose:
                    print(f'M mislabeled as F:\n{ix.tolist()}')
                res_s[ix] = 'Male'
    else:
        ax.scatter(x_s, y_s, s=s, alpha=0.5, edgecolors='none', lw=0.5, rasterized=True, clip_on=False)

    if highlight_ids is not None:  # highlight selected samples
        ax.scatter(x_s[highlight_ids], y_s[highlight_ids], c=[hsv_to_rgb([0.075,1,1])], s=s+12, alpha=1, edgecolors='k', lw=1, zorder=50, rasterized=False, clip_on=False, label=None)

    qtl.plot.format_plot(ax, fontsize=12)
    ax.spines['left'].set_position(('outward', 6))
    ax.spines['bottom'].set_position(('outward', 6))
    ax.set_xlabel('XIST expression (TPM)', fontsize=14)
    ax.set_ylabel('RPS4Y1 expression (TPM)', fontsize=14)

    xlim = list(ax.get_xlim())
    ylim = list(ax.get_ylim())
    xlim[0] = 0
    ylim[0] = 0
    ax.plot(2*[x_threshold], ylim, '--', c=[0.75]*3)
    ax.plot(xlim, 2*[y_threshold], '--', c=[0.75]*3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if sex_s is not None:
        leg = ax.legend(loc='upper left', fontsize=12, handlelength=0.5, labelspacing=0.2, bbox_to_anchor=(1,1))
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        return res_s
