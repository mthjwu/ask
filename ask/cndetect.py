################################################################################
# infer segments from breakpoint pairs
################################################################################
#------------------------------------------------------------------------------#
import numpy as np
import pandas as pd
import math
import scipy.stats as stats
import statsmodels.api as sm
from itertools import groupby
from itertools import combinations
from operator import itemgetter
from collections import defaultdict
from statsmodels.formula.api import glm
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------#
import cbs
import misc
from grange import GRange
from bpjoint import join_breakpoint


#------------------------------------------------------------------------------#
# get bias in binsize
#------------------------------------------------------------------------------#
def get_bin_bias(bin_count, biasfile, binsize = 10000):
    """calculate bias factor in binsize
    """
    df = pd.read_csv(biasfile, sep = '\t')
    df.columns = ['Chrom', 'Start', 'End', 'GCcontent', 'Mappability']
    df.Start = (np.floor(df.Start/binsize)*binsize).astype('int')
    df_agg = df.groupby(['Chrom', 'Start'])[['GCcontent', 'Mappability']]\
        .mean().reset_index().rename(columns = {"Start":"Coord"})
    dfm = pd.merge(bin_count, df_agg, on = ['Chrom', 'Coord'])
    return dfm


#------------------------------------------------------------------------------#
# glm remove GC and mappability biasa
#------------------------------------------------------------------------------#
def norm_glm(dfm, cutoff = [0, 0], method = 'nb', fold_bound = 5, plot = False):
    """fit glm to remove biases from bin counts

    Parameters
    ----------
    cutoff : select data points used in the model, default use all
    method : 'nb' - negative binomial
             'poisson' - poisson
    fold_bound : bound of bias factor to consider
        to avoid getting super small or large normalized counts
    """
    keep = (dfm.GCcontent >= cutoff[0]) & (dfm.Mappability >= cutoff[1])
    if (method == 'nb'):
        fit = glm('Count ~ GCcontent + Mappability', data = dfm,
            subset = keep, family = sm.families.NegativeBinomial()).fit()
    else:
        fit = glm('Count ~ GCcontent + Mappability', data = dfm,
            subset = keep, family = sm.families.Poisson()).fit()

    coef = fit.params
    fitted_y = np.exp(coef[0] + coef[1] * dfm.GCcontent +
                      coef[2] * dfm.Mappability)
    fac = fitted_y / np.median(fitted_y)

    fac[fac > fold_bound] = fold_bound
    fac[fac < 1/fold_bound] = 1/fold_bound

    norm_y = round(dfm.Count / fac, 4)

    bin_norm = dfm.iloc[:,0:4].copy()
    bin_norm['Count'] = norm_y + 1
    bin_norm['CN'] = 2 * bin_norm['Count'] / get_mode(bin_norm['Count'])

    if (plot):
        plt.plot(dfm.GCcontent, dfm.CN, 'o', alpha = 0.01)
        plt.yscale('log')
        plt.show()

        plt.plot(dfm.GCcontent, norm_y, 'o', alpha = 0.01)
        plt.yscale('log')
        plt.show()

        plt.plot(dfm.Mappability, dfm.CN, 'o', alpha = 0.01)
        plt.yscale('log')
        plt.show()

        plt.plot(dfm.Mappability, norm_y, 'o', alpha = 0.01)
        plt.yscale('log')
        plt.show()

    return bin_norm


#------------------------------------------------------------------------------#
# 2D bias plot
#------------------------------------------------------------------------------#
def plot_bias_2D(x, y, z, npts = 100, cmap = 'gist_rainbow',
                 xlab = 'GCcontent', ylab = 'Mappability'):
    """ 2D bias plot ('viridis')
    """

    # bin data and take average
    zi = stats.binned_statistic_2d(y, x, z, 'mean', bins=npts).statistic
    print(np.nansum(abs(zi-np.nanmedian(zi))))

    # make heatmap
    plt.imshow(zi, cmap=cmap, origin = 'lower', extent=(0,1,0,1))
    plt.colorbar()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.show()


#------------------------------------------------------------------------------#
# apply triangular smoothing on the normalized bin counts
#------------------------------------------------------------------------------#
def smooth_cn(bin_norm, k = 10):
    """smooth normalized bin counts
    """
    df = bin_norm.copy()
    # calculate log2 ratio
    df['Log2Ratio'] = np.log2(np.array(df.CN / 2))

    df_smooth = []
    for chr in misc.unique(df['Chrom']):
        dfsub = df[df['Chrom'] == chr]
        x = np.array(dfsub['Log2Ratio'])
        x_smooth = trismooth_1d(x, k = k)
        dfsub = dfsub.assign(Log2Ratio = x_smooth)
        df_smooth.append(dfsub)
    return pd.concat(df_smooth)


#------------------------------------------------------------------------------#
def trismooth_1d(x, k = 10, outlier_sd = 5,
                 edge_repeat = False):
    """triangular smoothing on a vector, only smooth outliers
    k : # of points to consider on the left and the right of a point
    outlier_sd : mean +- outlier_sd * sd is the outlier in the region
    edge_repeat : False - boundaries extended by filling nan
        True - boundaries extended by repeating edge points
    """

    x = np.array(x)
    y = np.zeros((len(x), 2*k + 1), dtype = x.dtype)
    y[:, k] = x
    for i in range(k):
        j = k - i
        y[j:,i] = x[:-j]
        y[:-j,-(i+1)] = x[j:]
        if (edge_repeat):
            y[:j,i] = x[0]
            y[-j:,-(i+1)] = x[-1]
        else:
            y[:j,i] = np.nan
            y[-j:,-(i+1)] = np.nan

    # identify outlier
    ym = y[:,k]
    ys = y[:,[i for i in range(2*k+1) if i != k]]
    ys_mean = np.nanmean(ys, axis = 1)
    ys_sd = np.nanstd(ys, ddof=1, axis = 1)
    ys_upper = ys_mean + outlier_sd * ys_sd
    ys_lower = ys_mean - outlier_sd * ys_sd

    # calculate on outlier
    outlier = (ym > ys_upper) | (ym < ys_lower)
    y[outlier, k] = np.nan
    y_outlier = y[outlier,]

    # weighted average on outliers
    tt = list(range(1, k+1))
    w = tt +[k+1] + tt[::-1]
    masked_y = np.ma.masked_array(y_outlier, np.isnan(y_outlier))
    avg = np.ma.average(masked_y, axis = 1, weights = w)
    avg = avg.filled(np.nan)

    # assign back
    smooth_x = x
    smooth_x[outlier] = avg

    return smooth_x


#------------------------------------------------------------------------------#
def trismooth_1d_plus(x, k = 10, outlier_sd = 5,
                      trim_quantile = 0.1, edge_repeat = False):
    """triangular smoothing on a vector, only smooth outliers
    k : # of points to consider on the left and the right of a point
    outlier_sd : mean +- outlier_sd * sd is the outlier in the region
    trim_quantile : remove trim_quantile extreme data points
                    to calcuate mean and sd
                    for both upper and lower quantile
    edge_repeat : False - boundaries extended by filling nan
        True - boundaries extended by repeating edge points
    """

    x = np.array(x)
    y = np.zeros ((len (x), 2*k + 1), dtype = x.dtype)
    y[:, k] = x
    for i in range(k):
        j = k - i
        y[j:,i] = x[:-j]
        y[:-j,-(i+1)] = x[j:]
        if (edge_repeat):
            y[:j,i] = x[0]
            y[-j:,-(i+1)] = x[-1]
        else:
            y[:j,i] = np.nan
            y[-j:,-(i+1)] = np.nan

    # weights
    tt = list(range(1, k+1))
    w = np.array(tt +[k+1] + tt[::-1])

    # smooth
    smooth_x = []
    for i in range(len(x)):
        z = y[i,:]
        zm = z[k]
        z[k] = np.nan

        with np.errstate(invalid='ignore'):
            tf = (z < np.nanquantile(z, trim_quantile)) |\
                 (z > np.nanquantile(z, 1 - trim_quantile))
        z[tf] = np.nan

        z_mean = np.nanmean(z)
        z_sd = np.nanstd(z, ddof=1)
        z_upper = z_mean + outlier_sd * z_sd
        z_lower = z_mean - outlier_sd * z_sd

        if (zm > z_upper or zm < z_lower):
            idx = np.where(~np.isnan(z))[0]
            smooth_x.append(np.average(z[idx], weights = w[idx]))
        else:
            smooth_x.append(zm)

    return smooth_x


#------------------------------------------------------------------------------#
def get_mode(x):
    """get mode from density
    """
    x = np.array(x)
    x = x[np.isfinite(x)]
    density = stats.kde.gaussian_kde(x)
    v = np.linspace(0, np.quantile(x, 0.99), 500)
    density = density(v)
    return v[density.argmax()]


#------------------------------------------------------------------------------#
# output function
#------------------------------------------------------------------------------#
def get_cn_segment(bin_count_clean, binsize = 10000, nperm = 1000, p = 0.01,
        method = 'cutoff', min_cn = None, std_scale = 8):
    """ get copy number segmentations bin count data
    output a dataframe of amplicons

    Parameters
    ----------
    df : counts in genomic bins
    method : method used to perform segmentation
        'cutoff' : use hard cutoff of min_cn to choose amplified regions
        'cbs' : use cbs method to choose amplified regions

    std_scale : mean + std_scale * std as min_cn
    binsize : target genomic bin size to calculate read counts

    output
    ----------
    segmented regions in dataframe
    """
    df = bin_count_clean.copy()

    # set min_cn by mean + n * std
    if (min_cn is None):
        xx = np.array(df['CN'])
        min_cn = np.mean(xx) + std_scale * np.std(xx, ddof = 1)
        print(min_cn)

    # calculate log2 ratio if use cbs
    if (method == 'cbs'):
        df['L2R'] = np.log2(np.array(df.CN / 2))

    # call copy number segmentations in each chromosome
    cnsegs = []
    for chr in set(df['Chrom']):
        # print(chr)
        dfsub = df[df['Chrom'] == chr]
        # seg = []
        if (method == 'cutoff'):
            x = np.array(dfsub['CN'])
            seg = simple_amplicon_segment(x, min_cn)
            seg = [list(dfsub.iloc[i[0],0:2]) +
                   list(dfsub.iloc[i[1]-1,1:2]+binsize) +
                   [i[2]] for i in seg]
        else:
            x = np.array(dfsub['L2R'])
            seg = cbs.cbs_segment(x, nperm=nperm, p=p)
            seg = cbs.cbs_recheck(x, seg, nperm=nperm, p=p)
            # cbs.cbs_plot_segment(x, seg) # plot function
            seg = [list(dfsub.iloc[i[0],0:2]) +
                   list(dfsub.iloc[i[1]-1,1:2]+binsize) +
                   [2**i[2]*2] for i in seg]
        cnsegs.append(pd.DataFrame(seg,
            columns=['Chrom', 'Start', 'End', 'CN']))
        seg_df = pd.concat(cnsegs)

    return seg_df[seg_df.CN >= min_cn]



#------------------------------------------------------------------------------#
# simple segmentation
#------------------------------------------------------------------------------#
def simple_seg(bin_count, binsize = 10000, min_cn = 5, std_scale = 8):
    """ get copy number segmentations bin count data
    output a dataframe of amplicons
    use hard cutoff of min_cn to choose amplified regions

    Parameters
    ----------
    bin_count : counts in genomic bins
    std_scale : mean + std_scale * std as min_cn
    binsize : target genomic bin size to calculate read counts
    min_cn : minimal copy number to consider

    output
    ----------
    segmented regions in dataframe
    """
    df = bin_count.copy()

    # set min_cn by mean + n * std
    if (min_cn is None):
        xx = np.array(df['CN'])
        min_cn = np.mean(xx) + std_scale * np.std(xx, ddof = 1)
        print(min_cn)

    # call copy number segmentations in each chromosome
    cnsegs = []
    for chr in set(df['Chrom']):
        # print(chr)
        dfsub = df[df['Chrom'] == chr]

        # segmentation
        x = np.array(dfsub['CN'])
        seg = simple_amplicon_segment(x, min_cn)
        seg = [list(dfsub.iloc[i[0],0:2]) +
               list(dfsub.iloc[i[1]-1,1:2]+binsize) +
               [i[2]] for i in seg]

        cnsegs.append(pd.DataFrame(seg,
            columns=['Chrom', 'Start', 'End', 'CN']))
        seg_df = pd.concat(cnsegs)

    return seg_df[seg_df.CN >= min_cn]


#------------------------------------------------------------------------------#
# simple calling of amplicons
#------------------------------------------------------------------------------#
def simple_amplicon_segment(x, min_cn = 5, min_dist_to_merge = 2):
    """
    Simple way to identify amplified segments
    score of genomic bins higher than a certain score
    let's say bins with more than 5 copies

    min_dist_to_merge: merge closeby intervals
    """

    xx = np.where(x > min_cn)[0].tolist()
    segments =[]

    # get the continuous numbers in a list with a tolerance
    grp = groupby(zip(xx, xx[1:] + xx[:0]), \
        lambda x: x[1] - x[0] <= min_dist_to_merge)
    for k,g in grp:
        if (k == True):
            gg = [i for i in g]
            segments.append([gg[0][0], gg[-1][1]+1])

    # get the segmentation mean table
    seg = []
    for j in segments:
        seg_start = j[0]
        seg_end = j[1]
        seg_mean = np.mean(x[seg_start:seg_end])
        seg.append([seg_start, seg_end, seg_mean])

    return seg


#------------------------------------------------------------------------------#
def segment_filter_by_size(cn_amp, binsize = 10000, fold = 5):
    """
    filter copy number segments by size (remove small spikes)

    """
    # old version
    # return cn_amp[cn_amp['End'] - cn_amp['Start'] >= fold * binsize]

    cn_amp_merged = misc.merge_bed(cn_amp, gap = 100000)
    cn_amp_drop = pd.DataFrame([
        row for row in cn_amp_merged if (row[2] - row[1] < fold * binsize
        )], columns = cn_amp.columns[0:3])
    df = pd.merge(cn_amp, cn_amp_drop, indicator = True, how = 'left'
        ).query('_merge == "left_only"').drop('_merge', axis = 1)
    return df


#------------------------------------------------------------------------------#
def bin_filter_by_blacklist(bin_count, blacklistfile, \
    binsize = 10000, ext = 0, set_value = 2):
    """
    set the bin CN to 2 (diploid), if overlap with blacklist

    ext - extend both end on ext * binsize
    """
    bin_count_ = bin_count.copy()

    # convert bed file to gr
    gr = GRange(blacklistfile, 'bedfile')

    # init blacklist bins
    bl_bin = defaultdict(lambda: False)

    # save blacklist bins to dict
    for _gr in gr.gr:
        _start = math.floor(_gr[1].start/binsize) - ext
        _end = math.floor(_gr[1].stop/binsize) + 1 + ext
        for i in range(_start, _end):
            bl_bin[(_gr[0], i*binsize)] = True

    # get the bool vector of blacklist bins
    tf = [bl_bin[row[1:3]] for row in bin_count_.itertuples()]

    # set blacklist bins to 0
    # bin_count_.loc[tf, 'Count'] = set_value
    # bin_count_['CN'] = 2*bin_count_['Count']/np.mean(bin_count_['Count'])

    # set blacklist bins CN to 2, left count unchanged
    bin_count_.loc[tf, 'CN'] = set_value

    return bin_count_

#------------------------------------------------------------------------------#
# identify cn segment boundaries by using clip reads
#------------------------------------------------------------------------------#
def finemap_cn_segment_boundary(cn_amp, bp_all, binsize = 10000):
    """
    search break points for single base resolution
    boundaries of amplified segments
    cn_amp: copy number high amplicons
    bp_all: bp_all
    """
    cn_amp = cn_amp[['Chrom', 'Start', 'End', 'CN']]
    seg_bdry = []
    for idx, seg in cn_amp.iterrows():
        bdry = list(seg) + [0, 0]

        # left boundary
        left = bp_all['left']
        left = left[left['Chrom'] == seg[0]]
        left = left.assign(dist = abs(seg[1]-left['Coord']))
        cand = left[left['dist'] <= binsize]
        if (not cand.empty):
            cand = cand.sort_values(['dist'])
            cand = cand.sort_values('Count', ascending=False)
            bdry[1] = list(cand['Coord'])[0]
            bdry[4] = list(cand['Count'])[0]

        # right boundary
        right = bp_all['right']
        right = right[right['Chrom'] == seg[0]]
        right = right.assign(dist = abs(seg[2]-right['Coord']))
        cand = right[right['dist'] <= binsize]
        if (not cand.empty):
            cand = cand.sort_values(['dist'])
            cand = cand.sort_values('Count', ascending=False)
            bdry[2] = list(cand['Coord'])[0]
            bdry[5] = list(cand['Count'])[0]

        # output
        seg_bdry.append(bdry)

    colnames = ['Chrom', 'Start', 'End', 'CN', 'ClipLeft', 'ClipRight']
    return pd.DataFrame(seg_bdry, columns=colnames)


#------------------------------------------------------------------------------#
def get_bppair_from_cnamp(cn_amp, bamfile, min_nt = 4):
    """
    get the bppairs from the refined cnamp
    as addition to the bppairs identifed by bp_cand
    """
    op = []
    for row in cn_amp.itertuples():
        op.append([row[1], row[2], 'L'])
        op.append([row[1], row[3], 'R'])
    op

    tt = [pair[0] + pair[1] + list(
            join_breakpoint(bamfile, *pair, min_nt = min_nt)) \
        for pair in combinations(op, 2)]

    t2 = [row[0:6] + [row[6][1] + row[7][1]] \
          for row in tt \
              if row[6][1] > 0 and row[7][1] > 0 \
              and row[6][0] == row[7][0]]

    colnames = ["Chrom1", "Coord1", "Clip1", "Chrom2", "Coord2", "Clip2", 'Count']
    bp_pair_extra = pd.DataFrame(t2, columns = colnames)

    return bp_pair_extra


#------------------------------------------------------------------------------#
def segment_annotation(cn_amp, refgenef):
    """
    annotate the copy number segments by gene
    """
    return GRange.map_gene(cn_amp, refgenef)






#------------------------------------------------------------------------------#
# filter copy number segments by blacklist intervals (currently not in use)
#------------------------------------------------------------------------------#
def segment_filter_by_blacklist(cn_amp, blacklistfile, f_sort=True):
    """
    filter copy number segments by blacklist intervals

    """
    def filter_blacklist(df, blacklistfile):
        gr1 = GRange(df, 'dataframe_hasend')
        gr2 = GRange(blacklistfile, 'bedfile')
        # extend 100 bp on both end of the blacklist
        gr = gr1.intersect(gr2, a_extend = 0, b_extend = 100, invert = True)

        return pd.DataFrame(\
            [[row[0]] + [row[1].start, row[1].stop] + list(row[2]) for row in gr.gr]\
            , columns=df.columns)

    df = filter_blacklist(cn_amp, blacklistfile)

    if (f_sort):
        df = df.sort_values(['Chrom', 'Start', 'End'])

    return df
