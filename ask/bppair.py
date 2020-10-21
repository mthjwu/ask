################################################################################
# breakpoint pair detection and filtering
################################################################################
#------------------------------------------------------------------------------#
import pysam
import re
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------#
from grange import GRange
import misc


#------------------------------------------------------------------------------#
# filter breakpoint pairs
#------------------------------------------------------------------------------#
def bp_pair_filter(bp_duo, short_distance = 200, min_junc_cnt = 5,
                   max_offset_overhang = 1000, min_match_len = 5,
                   max_offset_insertion = 20, max_srp = 0.7):
    """
    filter breakpoint pairs

    max_srp : max simple repeat proportion
    """
    bp_pair = bp_duo[bp_duo['Count'] >= min_junc_cnt]

    ind = (bp_pair['Chrom1'] == bp_pair['Chrom2']) \
        & (abs(bp_pair['Coord2'] - bp_pair['Coord1']) < short_distance)
    # also control for max insertion size (+) and max overhang size (-)
    ind = ind | (bp_pair['offset'] > max_offset_insertion) \
              | (bp_pair['offset'] < -max_offset_overhang)

    # get the matching sequence which in upper case
    match_seq = [re.search('([A-Z]+)', i).group(0) for i in bp_pair['Seq']]
    # length of the matched sequence
    match_len = [len(s) for s in match_seq]
    # get the whether is a simple repeat or not
    srp = [simple_repeat_proportion(s) for s in match_seq]
    # remove simple repeat bp pairs which is artificial
    ind2 = (np.array(srp) > max_srp) | \
           (np.array(match_len) < min_match_len)

    # include the PE-Support pairs
    ind = (~ind & (bp_pair['Seq'] == 'PE_Support')) | \
          (~ind & ~ind2)

    return bp_pair.loc[ind].reset_index(drop=True)

#------------------------------------------------------------------------------#
def simple_repeat_proportion(s, n_iter = 4):
    """search simple repeat in a string
    return the proportion
    """
    match_ = [i for i in re.compile(r"(.+?)\1+").finditer(s)]
    return sum([len(i.group()) for i in match_
                if len(i.group()) >= n_iter])/len(s)


#------------------------------------------------------------------------------#
# get breakpoint from bp pairs and cn segments
#------------------------------------------------------------------------------#
def bp_refine(bp_pair, bp_cand_stats, cn_amp):
    """
    get breakpoints from breakpoint pairs and amplified segments
    """

    # add sequence depth to stats
    col1 = ['Chrom1', 'Coord1', 'Clip1']
    col2 = ['Chrom2', 'Coord2', 'Clip2']
    newcol = ['Chrom', 'Coord', 'Clip']
    bp_df = pd.concat([bp_pair[col1].set_axis(newcol, axis=1, inplace=False), \
                    bp_pair[col2].set_axis(newcol, axis=1, inplace=False)])
    bp_df = bp_df.drop_duplicates()

    # add clip depth to stats
    bp_stats = pd.merge(bp_df, bp_cand_stats, on=['Chrom', 'Coord', 'Clip'])

    # breakpoints from amplified segments
    ## merge adjacent segments into one
    cn_amp_merged = misc.merge_bed(cn_amp)
    cn_amp_clip = dict()
    for i in cn_amp.itertuples():
        cn_amp_clip[(i[1], i[2], 'L')] = i[5]
        cn_amp_clip[(i[1], i[3], 'R')] = i[6]
    ## make breakpoints dataframe
    bplist = [i[1:4] for i in bp_stats.itertuples()]
    op = []
    for row in cn_amp_merged:
        L = (row[0], row[1], 'L')
        R = (row[0], row[2], 'R')
        if (L not in bplist):
            op.append([row[0], row[1], 'L', True, cn_amp_clip[L], 0, 0])
        if (R not in bplist):
            op.append([row[0], row[2], 'R', True, cn_amp_clip[R], 0, 0])
    cn_seg_df = pd.DataFrame(op, columns=bp_stats.columns)
    cn_seg_df

    # merge and output
    df = pd.concat([bp_stats, cn_seg_df])
    df = df.drop_duplicates().sort_values(['Chrom', 'Coord', 'Clip'])
    return df.reset_index(drop=True)


#------------------------------------------------------------------------------#
# infer segments from breakpoint pairs
#------------------------------------------------------------------------------#
def get_segment(bp_fine, bp_pair, restrict = False, \
        min_segment_size = 100, max_segment_size = 100000000):
    """
    infer segments from breakpoint table

    usage:
        get_segment(bp_fine, bp_pair)

    input:
    bp_fine: dataframe of finally refined breakpoints
    output: pandas dataframe of segments
    restrict: only use breakpoints with True cleanBP to infer segments if True
    remove segments with size < min_segment_size and > max_segment_size
    """

    df = bp_fine.sort_values(['Chrom', 'Coord'])

    if (restrict): # only use CleanBP True
        df = df[df['CleanBP']]

    # init
    seg = []
    colnames = ['Chrom', "Start", "End"]

    # infer potential segments from breakend patterns (<-...->)
    for chrom in set(df['Chrom']):
        ind = df['Chrom'] == chrom
        df_all = df.loc[ind].copy()
        # df_all.loc[df_all['InDepth'] == 0, 'CleanBP'] = False
        df_l = df_all.loc[df_all['Clip'] == "L"].reset_index(drop=True)
        df_r = df_all.loc[df_all['Clip'] == "R"].reset_index(drop=True)

        if (df_l.size > 0 and df_r.size > 0):
            for d,row in df_l.iterrows():
                dst = (row['Coord']-df_r['Coord'])
                dst = np.where(dst < 0, -dst, np.inf)
                sort_index = np.argsort(dst)
                try:
                    ind = list(df_r['CleanBP'][sort_index]).index(True)
                    end = [df_r['Coord'][sort_index[i]] for i in range(ind+1)]
                    df_seg = pd.DataFrame(
                        [list(row[0:2]) + [i] for i in end],
                        columns = colnames)
                    seg.append(df_seg)
                except:
                    print(row, '\n')

            for d,row in df_r.iterrows():
                dst = (df_l['Coord']-row['Coord'])
                dst = np.where(dst < 0, -dst, np.inf)
                sort_index = np.argsort(dst)
                try:
                    ind = list(df_l['CleanBP'][sort_index]).index(True)
                    end = [df_l['Coord'][sort_index[i]] for i in range(ind+1)]
                    df_seg = pd.DataFrame(
                        [[row[0]] + [i] + [row[1]] for i in end],
                        columns = colnames)
                    seg.append(df_seg)
                except:
                    print(row, '\n')

    # add segment from direct loop
    # in case detailed structure can't be found
    # return this simple circle
    op = []
    for row in bp_pair.itertuples():
        if (row[1] == row[4]):
            if (row[2] < row[5] and row[3] == 'L' and row[6] == 'R'):
                op.append([row[1], row[2], row[5]])
            elif (row[2] > row[5] and row[3] == 'R' and row[6] == 'L'):
                op.append([row[1], row[5], row[2]])
    seg.append(pd.DataFrame(op, columns = colnames))

    # merge all segments
    seg = pd.concat(seg).drop_duplicates()
    seg = seg.drop_duplicates()

    # remove oversized segments
    seg_size = seg['End'] - seg['Start']
    seg = seg[(seg_size >= min_segment_size) & (seg_size <= max_segment_size)]

    return seg.sort_values(['Chrom', 'Start', 'End']).reset_index(drop=True)


#------------------------------------------------------------------------------#
# get breakpoint from bp pairs and cn segments
#------------------------------------------------------------------------------#
def add_cn(df, cn_amp):
    """
    map copy number to breakpoint segments
    """

    # create GRange objects
    gr1 = GRange(df[['Chrom', 'Start', 'End']], 'dataframe_hasend')
    gr2 = GRange(cn_amp, 'dataframe_hasend')

    # only use the mid point for the segments
    op = []
    for i in gr1.gr:
        mid = int((i[1].start + i[1].stop)/2)
        op.append((i[0], range(mid, mid+1), i[2]))
    gr1.gr = op

    # extend binsize on both end of the cn_seg
    map_list = gr1.gmap(gr2, a_extend = 0, b_extend = 0)

    # output
    df['CN'] = [round(i[0]) if (i is not None) else None for i in map_list]
    return df
