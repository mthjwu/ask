################################################################################
# breakpoint detection by analyzing soft and hard clip reads
################################################################################
#------------------------------------------------------------------------------#
import pysam
import numpy as np
import pandas as pd
import math

#------------------------------------------------------------------------------#
from grange import GRange
import misc

#------------------------------------------------------------------------------#
# function to check read
#------------------------------------------------------------------------------#
def check_read(read, mapq = 20, nmmax = 1):
    """
    true if the read meets certian conditions
    """
    return not read.is_unmapped \
        and read.mapping_quality >= mapq \
        and read.get_tag('NM') <= nmmax \
        and not read.is_duplicate

#------------------------------------------------------------------------------#
# get all clip depth
#------------------------------------------------------------------------------#
def clip_from_bam(bamfile, mapq = 20, nmmax = 1, n = 0):
    """
    get all clip positions from bam files
    -- count softclip and hardclip reads per chromosome coordinates
    for left and right clip, respectively

    usage:
        clip_from_bam(bamfile)

    input: bamfile
    output: pandas dataframe of split read counts in a dict
    nmmax: max number of mismatches
    mapq: minimal mapping quality
    readq: minimal read mean quality (readq = 20, not necessary)
    """

    with pysam.AlignmentFile(bamfile, "rb") as bamf:

        left = list() # init left clip
        right = list() # init right clip

        for read in bamf.fetch():
            if check_read(read, mapq = mapq, nmmax = nmmax):

                # left clip read
                if (read.cigartuples[0][0] == 4): # soft clip reads
                    left.append((read.reference_name, \
                                 read.reference_start, "L"))
                elif (read.cigartuples[0][0] == 5): # hard clip reads
                    left.append((read.reference_name, \
                                 read.reference_start, "L"))

                # right clip reads
                if (read.cigartuples[-1][0] == 4): # soft clip reads
                    right.append((read.reference_name, \
                                  read.reference_end-1, "R"))
                elif (read.cigartuples[-1][0] == 5): # hard clip reads
                    right.append((read.reference_name, \
                                  read.reference_end-1, "R"))

    # convert all clip read counts to df
    left = list2df_count(left, n = n, f_sort = False)
    right = list2df_count(right, n = n, f_sort = False)

    # return left and right clip points by chrom.position
    return dict(left=left, right=right)

#------------------------------------------------------------------------------#
def list2df_count(lst, colnames=["Chrom", "Coord", "Clip"], n=10, f_sort=True):
    """
    convert list to dataframe, count and sort

    usage:
        list2df_count(lst)

    input: list of softclip reads with "Chrom", "Coord" and "Clip"
    output: pandas dataframe of split read counts

    """
    df = pd.DataFrame(lst, columns=colnames)
    df = df.groupby(colnames).size().reset_index(name='Count')
    if (n > 0):
        df = df[df['Count']>=n]
    if (f_sort):
        df = df.sort_values('Count', ascending=False)
    return df

#------------------------------------------------------------------------------#
def clip2bedgraph(clip, n = 2):
    """
    clip to bedgraph for visualization
    """
    # init
    left = clip['left'].copy()
    right = clip['right'].copy()

    # filter
    left = left[left['Count'] >= n]
    right = right[right['Count'] >= n]

    # output bedgraph format to a pandas dataframe
    left['Start'] = left['Coord']
    left['End'] = left['Start'] + 1
    right['Start'] = right['Coord']
    right['End'] = right['Start'] + 1
    left['Count'] = -left['Count']
    df = pd.concat([left[['Chrom', 'Start', 'End', 'Count']],
                    right[['Chrom', 'Start', 'End', 'Count']]])
    return df.sort_values(['Chrom', 'Start', 'End'])



#------------------------------------------------------------------------------#
# apply heavy smoothing on sub bin counts to remove peaks in ATAC-seq
#------------------------------------------------------------------------------#
def smooth_count(bin_count, k = 20, outlier_sd = 2, npass = 10):
    """smooth bin counts
    """
    df = bin_count.copy()

    df_smooth = []
    for chr in misc.unique(df['Chrom']):
        dfsub = df[df['Chrom'] == chr]
        x_smooth = smooth_count_1d(dfsub['Count'], k, outlier_sd, npass)
        dfsub = dfsub.assign(Count = x_smooth)
        df_smooth.append(dfsub)
    return pd.concat(df_smooth)


#------------------------------------------------------------------------------#
def smooth_count_1d(x, k = 20, outlier_sd = 2, npass = 10, edge_repeat = False):
    """smoothing on a vector, only smooth outliers
    k : # of points to consider on the left and the right of a point
    outlier_sd : mean +- outlier_sd * sd is the outlier in the region
    edge_repeat : False - boundaries extended by filling nan
        True - boundaries extended by repeating edge points
    """

    x = np.array(x, dtype = 'float64')
    n = 0

    while n < npass:
        # contruct working matrix
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
        outlier = (ym > ys_upper)

        # assign outlier to na
        x[outlier] = np.nan

        # iteration counter
        n = n + 1

    return x


#------------------------------------------------------------------------------#
# count reads in genomic bins
#------------------------------------------------------------------------------#
def region_count(bamfile, genomesizefile = None,
                 binsize = 10000, mapq = 20,
                 nmmax = 1, sort_ = True):
    """
    count reads in genome bins
    """
    cnt = list()

    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        for read in bamf.fetch():
            if check_read(read, mapq = mapq, nmmax = nmmax):
                pos = math.floor(read.reference_start/binsize)*binsize
                cnt.append([read.reference_name, pos])

    # convert all read counts to df
    colnames = ['Chrom', 'Coord']
    df = pd.DataFrame(cnt, columns=colnames)
    df = df.groupby(colnames).size().reset_index(name='Count')
    df['CN'] = 2*df['Count']/np.mean(df['Count'])

    if (sort_): df = df.sort_values(colnames)
    return df


#------------------------------------------------------------------------------#
def region_count_2pass(bamfile, genomesizefile = None,
                   binsize = 10000, mapq = 20,
                   nmmax = 1, sort_ = True,
                   sub_binsize = 100, q = 0.5,
                   colnames = ['Chrom', 'Coord']):
    """wrapper of calculating robust read counts
    in target genomic bins

    Parameters
    ----------
    binsize : target genomic bin size to calculate read counts
    sub_binsize : smaller bin size within target bin size to
        calculate robust score
    q : quantile as robust score, if None auto determine
    """
    # calculate read counts in smaller bin size
    df_sub = region_count(bamfile, genomesizefile, sub_binsize,
                          mapq, nmmax, False)

    # smooth the sub bin counts to remove spikes
    df_sub_smooth = smooth_count(df_sub, k = 2*int(binsize/sub_binsize))

    # calculate read counts in target bin size
    df = region_count_bias(df_sub_smooth, binsize, sub_binsize, q, colnames)

    # fill NA with linear average of nearby counts
    df = df.assign(
        Count = df.groupby(['Chrom'])['Count'].apply(
            lambda group: group.interpolate()
        )
    )

    # calculate the absolute copy number in each bin
    df['CN'] = 2*df['Count']/np.mean(df['Count'])

    if (sort_): df = df.sort_values(colnames)
    return df

#------------------------------------------------------------------------------#
def region_count_bias(df_sub, binsize = 10000, sub_binsize = 100,
                  q = 0.5, colnames = ['Chrom', 'Coord']):
    """calculate robust read counts in target genomic bins
    by using quantiles of read counts in smaller bins
    """
    df = df_sub.copy()
    df = df.assign(Coord = (
        np.floor(np.array(df['Coord'])/binsize)*binsize).astype(int))
    df = df.groupby(colnames)['Count'].quantile(q = q).reset_index(name='Count')
    return df
