################################################################################
# circular binary segmentation implementation
# (ugly max t-statistic searching method,
# might only work for high focal amplification)
# currently input should be absolute copy number
# todo:
#   include 1 bin spikes
################################################################################
#------------------------------------------------------------------------------#
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import groupby, combinations
from operator import itemgetter

#------------------------------------------------------------------------------#
from grange import GRange
import misc



#------------------------------------------------------------------------------#
# select amplified segments on cbs results
#------------------------------------------------------------------------------#
def cbs_amplicon(seg_df, min_cn = 5):
    """get amplicon from cbs segmentation results
    """

    return seg_df[seg_df.CN >= min_cn]



#------------------------------------------------------------------------------#
# cbs segmentation
#------------------------------------------------------------------------------#
def cbs(bin_count, binsize = 10000, nperm = 1000, p = 0.01):
    """ get copy number segmentations from bin count data
    output a dataframe of amplicons

    Parameters
    ----------
    df : counts in genomic bins
    binsize : target genomic bin size to calculate read counts
    min_cn : minimal copy number

    output
    ----------
    segmentation results in dataframe
    """
    df = bin_count.copy()

    # calculate log2 ratio if use cbs
    if ('Log2Ratio' not in df.columns):
        df['Log2Ratio'] = np.log2(np.array(df.CN / 2))

    # call copy number segmentations in each chromosome
    cnsegs = []
    for chr in misc.unique(df['Chrom']):
        # print(chr)
        dfsub = df[df['Chrom'] == chr]

        # perform cbs on each chromosome
        x = np.array(dfsub['Log2Ratio'])
        seg = cbs_segment(x, nperm = nperm, p = p)
        seg = cbs_recheck(x, seg, nperm = nperm, p = p)
        # cbs_plot_segment(x, seg, (3000, 4000)) # plot function
        seg = [list(dfsub.iloc[i[0],0:2]) +
               list(dfsub.iloc[i[1]-1,1:2]+binsize) +
               [i[2]] for i in seg]

        cnsegs.append(pd.DataFrame(seg,
            columns=['Chrom', 'Start', 'End', 'Log2Ratio']))
        seg_df = pd.concat(cnsegs)
        seg_df['CN'] = 2**seg_df.Log2Ratio*2

    return seg_df[['Chrom', 'Start', 'End', 'CN', 'Log2Ratio']]


#------------------------------------------------------------------------------#
# filter copy number amplified segments by blacklist intervals
#------------------------------------------------------------------------------#
def amplicon_filter_by_blacklist(cn_amp, blacklistfile, f_sort=True):
    """
    filter amplicon by blacklist intervals
    for left and right clip, respectively
    (use pandas dataframe to process and store)

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


#------------------------------------------------------------------------------#
# simple calling of amplicons
#------------------------------------------------------------------------------#
def simple_amplicon_segment(x, min_cn = 5):
    """
    Simple way to identify amplified segments
    score of genomic bins higher than a certain score
    let's say bins with more than 5 copies
    """

    xx = np.where(x > min_cn)[0].tolist()
    segments =[]

    # get the combine the continuous numbers in a list
    for k,g in groupby(enumerate(xx),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        segments.append((group[0],group[-1]+1))

    # get the segmentation mean table
    seg = []
    for j in segments:
        seg_start = j[0]
        seg_end = j[1]
        seg_mean = np.mean(x[seg_start:seg_end])
        seg.append([seg_start, seg_end, seg_mean])

    return seg


#------------------------------------------------------------------------------#
# segmentation function
#------------------------------------------------------------------------------#
def cbs_segment(x, nperm = 1000, p = 0.01, k = 1000):
    """
    segmentation function

    k : split to k data points per block
    """
    start = 0
    end = len(x)
    segments = []

    # get segments
    for chunk in make_window(end):
        segments = cbs_segment_helper(
            x, chunk[0], chunk[1], segments, nperm = nperm, p = p)

    # get the segmentation mean table
    seg = []
    for j in segments:
        seg_start = j[0]
        seg_end = j[1]
        seg_mean = np.nanmean(x[seg_start:seg_end])
        seg.append([seg_start, seg_end, seg_mean])

    return seg

#------------------------------------------------------------------------------#
def cbs_segment_helper(x, start, end, segments, nperm = 1000, p = 0.01):
    """
    Recursive segmentation helper function
    """
    # print(start, end)

    # if (end - start <= 3): # no split need for <= 3 length interval
    #     return segments

    # print(start, end, segments)
    ij = cbs_determine_ij(x, start, end, nperm = nperm, p = p)
    # print(start, end, ij, segments)

    # if (ij[1] - ij[0] == 1 and end - start > 3):
    #     ij[2] = True

    if (ij[2] ==  False):
        segments.append((start, end))

    else:
        if (ij[0] - start >= 2):
            cbs_segment_helper(x, start, ij[0]+1, segments, nperm=nperm, p=p)
        elif (ij[0] - start < 2): # left edge spike
            segments.append((start, ij[0]+1))

        if (ij[1]-ij[0] >= 3):
            cbs_segment_helper(x, ij[0]+1, ij[1]+1, segments, nperm=nperm, p=p)
        elif (ij[1]-ij[0] < 3): # middle spike
            segments.append((ij[0]+1, ij[1]+1))

        if (end - ij[1] >= 4):
            # print(x, end, ij, segments)
            cbs_segment_helper(x, ij[1]+1, end, segments, nperm=nperm, p=p)
        elif (end - ij[1] < 4 and end - ij[1] > 1): # right edge spike
            segments.append((ij[1]+1, end))

    return segments

#------------------------------------------------------------------------------#
def make_window(n, k = 1000, l = 100):
    """bin n numbers into windows,
    k points per window with l overlap
    """
    op = [[i-l, i+k+l] for i in range(0, n, k)]
    op[0][0] = 0
    op[-1][1] = n
    if (len(op) > 1 and op[-1][1]-op[-1][0] < k/2 + l):
        op.pop()
        op[-1][1] = n
    return op


#------------------------------------------------------------------------------#
# recheck function
#------------------------------------------------------------------------------#
def cbs_recheck(x, segments, nperm = 1000, p = 0.01, tmin = 1.5):
    """
    recheck whether every three splits combination are significant

    tmin : minimal t stats to run permutation, too low t stats means
        no difference of copy number between segments
    min_fc : minimal fold change between segments
    """
    sp_cand = [i[0] for i in segments]+[len(x)] # split points
    sp_cand = misc.unique(sp_cand)
    sp_cand.sort()
    # print(sp_cand)
    sp = [0]
    while (len(sp_cand) >= 3):
        start = sp_cand[0]
        mid = sp_cand[1]
        end = sp_cand[2]
        i = mid - start - 1

        xs = x[start:end]
        S = np.cumsum(xs) # cumulative sum of sliced x

        tmax = cbs_tstats_ij(S, i) # t-statistic

        if (tmax >= tmin):
            tf = cbs_permute(xs, tmax, i, nperm = nperm, p = p) # permutation
        else:
            tf = False
        # print(start, mid, end, tmax, tf)

        if (tf == True):
            sp.append(mid)
            # sp.append(end)
            sp_cand.remove(start)
        else:
            # if mid in sp: sp.remove(mid)
            sp_cand.remove(mid)

    if (sp[-1] != sp_cand[-1]):
        sp.append(sp_cand[-1]) # add chrom end

    sp = misc.unique(sp) # make unique breaks

    seg = []
    for j in range(1, len(sp)):
        seg_start = sp[j - 1]
        seg_end = sp[j]
        seg_mean = np.mean(x[seg_start:seg_end])
        seg.append([seg_start, seg_end, seg_mean])

    return seg


#------------------------------------------------------------------------------#
# sub functions
#------------------------------------------------------------------------------#
def cbs_tstats_ij(S, i, j = None):
    """
    calculate the t-statistic for i, j breaks or one i break

    S : np.array of cumulative sum of x
    i : np.array of segment 1 end
    j : np.array of segment 2 end (optional)
    """
    if (j is not None):
        Si = S[j] - S[i]
        Sj = S[-1] - Si
        k = j - i
        n = len(S)
        Tn = (Si/k - Sj/(n-k))
        Td = (1/k+1/(n-k))**(1/2) # assume equal variance, cancel np.std(x, ddof=1)
        T = abs(Tn/Td)
    else:
        Si = S[i]
        Sj = S[-1] - Si
        k = i + 1
        n = len(S)
        Tn = (Si/k - Sj/(n-k))
        Td = (1/k+1/(n-k))**(1/2) # assume equal variance, cancel np.std(x, ddof=1)
        T = abs(Tn/Td)

    return T # return t-statistic and fold change

#------------------------------------------------------------------------------#
def cbs_permute(x, tmax, i, j = None, nperm = 1000, p = 0.01):
    """
    permutation test for t-statistic

    x: copy number data per genomic bin
    i: segment 1 end
    j: segment 2 end
    tmax: max t-statistic between two segments
    nperm: # of permutations
    p: p-value cutoff
    """
    h0_count = 0
    alpha = nperm * p
    xp = x.copy()
    for p in range(nperm):
        seed = p
        np.random.seed(seed)
        np.random.shuffle(xp)
        S = np.cumsum(xp)
        if (j is not None): # two split points
            h0 = cbs_tstats_ij(S, i, j)
        else:  # one split point
            h0 = cbs_tstats_ij(S, i)
        if h0 >= tmax:
            h0_count += 1
        if h0_count > alpha:
            return False
    return True

#------------------------------------------------------------------------------#
def cbs_determine_ij(x, start, end, nperm = 1000, p = 0.01, tmin = 1.5):
    """
    Determine i and j at max t-statistic

    x: copy number data per genomic bin
    start: start index in x
    end: end index in x
    tmin : minimal t stats to run permutation, too low t stats means
        no difference of copy number between segments

    slice x by start and end to perform the analysis
    and output the i, j split achieve max t-statistic
    """
    if (end - start <= 3):
        return [start, end,  False]
    xs = x[start:end]
    S = np.cumsum(xs) # cumulative sum of sliced x

    ## refine the i j choice by t-statistic
    ii = []
    jj = []
    xs_len = len(xs)
    for k1 in range(1, xs_len - 2):
        for k2 in range(k1 + 2, xs_len):
            ii.append(k1)
            jj.append(k2)
    ii = np.array(ii)
    jj = np.array(jj)
    # comb = list(combinations(range(len(xs)), 2))
    # ii = np.array([xx[0] for xx in comb])
    # jj = np.array([xx[1] for xx in comb])

    # determine max T
    T = cbs_tstats_ij(S, ii, jj)

    ## get the max T case
    maxT = np.argmax(T)
    i, j, tmax = ii[maxT], jj[maxT], T[maxT]

    ## Test significance by permutation
    # tf = cbs_permute(xs, tmax, i, j, nperm = nperm, p = p)
    if (tmax >= tmin):
        tf = cbs_permute(xs, tmax, i, j, nperm = nperm, p = p)
    else:
        tf = False

    # print(i + start, j + start, tmax, tf)

    return [i + start, j + start, tf]


#------------------------------------------------------------------------------#
# plot functions
#------------------------------------------------------------------------------#
def cbs_plot_segment(x, segments, xlim = (0, None)):
    '''
    segmentation plot
    '''
    split_points = [s[0] for s in segments]+[len(x)]
    plt.figure()
    p = sns.scatterplot(range(len(x)), x, color='grey', size = 0.1, legend=None)
    for i in split_points:
        p.axvline(i-0.5, lw = 0.3, alpha = 0.5)
    for j in segments:
        seg_start = j[0]
        seg_end = j[1]
        seg_mean = j[2]
        p.hlines(seg_mean, seg_start-0.5, seg_end-0.5, color='red')
    p.set_xlim(xlim[0], xlim[1])
    p.get_figure().set_size_inches(16, 4)
    return p
