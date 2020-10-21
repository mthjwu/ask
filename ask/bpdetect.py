################################################################################
# breakpoint detection by analyzing soft and hard clip reads
################################################################################
#------------------------------------------------------------------------------#
import pysam
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------#
from grange import GRange
from covbam import check_read


#------------------------------------------------------------------------------#
# filter breakpoint by "n" supporting reads
#------------------------------------------------------------------------------#
def bp_filter_by_n(bpall, n=5, f_sort=True):
    """
    filter breakpoint by # of softclip and hardclip reads
    for left and right clip, respectively
    (use pandas dataframe to process and store)

    usage:
        breakpoint_filter_by_n(bpall)

    input: output of breakpoint_from_bam
    output: pandas dataframe of split read counts in a dict

    """
    df_left = bpall['left']
    df_right = bpall['right']

    if (n > 0):
        df_left = df_left[df_left['Count']>=n]
        df_right = df_right[df_right['Count']>=n]
    if (f_sort):
        df_left = df_left.sort_values('Count', ascending=False)
        df_right = df_right.sort_values('Count', ascending=False)

    ## return a dict of chrom.position with "left" and "right" keys
    return dict(left=df_left, right=df_right)


#------------------------------------------------------------------------------#
# at most include 200 bp candidates on each side
#------------------------------------------------------------------------------#
def bp_top_n(bp_cand, topn = 200):
    """
    at most include 200 bp candidates on each side
    """
    left = bp_cand['left']
    if (len(left.index) > topn):
        cutoff = left['Count'][topn-1]
        left = left[left['Count'] > cutoff]

    right = bp_cand['right']
    if (len(right.index) > topn):
        cutoff = right['Count'][topn-1]
        right = right[right['Count'] > cutoff]

    return dict(left=left, right=right)

#------------------------------------------------------------------------------#
# filter breakpoints by blacklist intervals
#------------------------------------------------------------------------------#
def bp_filter_by_blacklist(bpall, blacklistfile, f_sort=True):
    """
    filter breakpoint by blacklist intervals
    for left and right clip, respectively
    (use pandas dataframe to process and store)

    usage:
        breakpoint_filter_by_blacklist(bpall)

    input: output of breakpoint_from_bam
    output: pandas dataframe of split read counts in a dict

    """
    def filter_blacklist(df, blacklistfile):
        gr1 = GRange(df, 'dataframe')
        gr2 = GRange(blacklistfile, 'bedfile')
        # extend 100 bp on both end of the blacklist
        gr = gr1.intersect(gr2, a_extend = 0, b_extend = 100, invert = True)
        return pd.DataFrame(\
            [[row[0]] + [row[1].start] + list(row[2]) for row in gr.gr]\
            , columns=df.columns)

    df_left = bpall['left']
    df_right = bpall['right']

    df_left = filter_blacklist(df_left, blacklistfile)
    df_right = filter_blacklist(df_right, blacklistfile)

    if (f_sort):
        df_left = df_left.sort_values('Count', ascending=False)
        df_right = df_right.sort_values('Count', ascending=False)

    ## return a dict of chrom.position with "left" and "right" keys
    return dict(left=df_left, right=df_right)


#------------------------------------------------------------------------------#
# filter breakpoints by amplified segments
#------------------------------------------------------------------------------#
def bp_filter_by_amplicon(bpall, cn_amplicon, binsize = 10000, f_sort=True):
    """
    filter breakpoint by amplified segments
    for left and right clip, respectively
    (use pandas dataframe to process and store)

    bin_extend: extend cn_amplicon by <N> bp on both sides
    """
    def filter_amplicon(df, cn_amplicon):
        gr1 = GRange(df, 'dataframe')
        gr2 = GRange(cn_amplicon, 'dataframe_hasend')
        # extend binsize bp on both end of the cn_amplicon
        gr = gr1.intersect(gr2, a_extend = 0, b_extend = 2*binsize, invert = False)
        gr.gr
        return pd.DataFrame(\
            [[row[0]] + [row[1].start] + list(row[2]) for row in gr.gr]\
            , columns=df.columns)

    df_left = bpall['left']
    df_right = bpall['right']

    df_left = filter_amplicon(df_left, cn_amplicon)
    df_right = filter_amplicon(df_right, cn_amplicon)

    if (f_sort):
        df_left = df_left.sort_values('Count', ascending=False)
        df_right = df_right.sort_values('Count', ascending=False)

    ## return a dict of chrom.position with "left" and "right" keys
    return dict(left=df_left, right=df_right)


#------------------------------------------------------------------------------#
# Breakpoint evaluation by sequence depth
#------------------------------------------------------------------------------#
def bp_seq_depth(bp_cand, bamfile, perc = 0.2, mapq = 20, nmmax = 1):
    """
    calculate the sequence depth on the breakpoint (InDepth)
    and 1-bp out of breakpoint (OutDepth)

    Clean breakpoint (CleanBP) is defined as the
    ClipDepth/InDepth >= perc,
    and means it can't be set as alternative segment boundary

    """
    bp_cand_df = pd.concat([bp_cand['left'], bp_cand['right']])\
        .sort_values('Count', ascending=False).reset_index(drop=True)

    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        op = []
        for index, row in bp_cand_df.iterrows():
            # depth_in = bamf.count(row[0], row[1], row[1]+1, read_callback = check_read)

            depth_in = len([
                read for read in bamf.fetch(row[0], row[1], row[1]+1) \
                if check_read(read, mapq = mapq, nmmax = nmmax)])

            depth_out = depth_in - row[3]
            # if (depth_out < 0):
            #     print([[row], depth_in])
            depth_bp = [depth_in, depth_out, row[3]/depth_in >= perc]
            op.append(depth_bp)

    colnames=['InDepth', 'OutDepth', 'CleanBP']
    df = pd.concat([bp_cand_df, pd.DataFrame(op, columns=colnames)], axis=1)
    df = df[['Chrom', 'Coord', 'Clip', 'CleanBP', 'Count', 'InDepth', 'OutDepth']]
    df.rename(columns={'Count':'ClipDepth'}, inplace=True)
    return df
