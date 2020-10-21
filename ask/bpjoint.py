################################################################################
# breakpoint detection by analyzing soft and hard clip reads
################################################################################
#------------------------------------------------------------------------------#
import pysam
import numpy as np
import pandas as pd
import re
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from itertools import combinations

#------------------------------------------------------------------------------#
from covbam import check_read
import misc


#------------------------------------------------------------------------------#
def check_read_clip(read, clip_type = 'soft_left'):
    """
    true if the read is the clip read
    """
    tf = check_read(read)

    if (tf):
        clip_tf = False
        if (clip_type == 'soft_left' and read.cigartuples is not None):
            clip_tf = read.cigartuples[0][0] == 4
        elif (clip_type == 'soft_right' and read.cigartuples is not None):
            clip_tf = read.cigartuples[-1][0] == 4
        tf = tf and clip_tf

    return tf

#------------------------------------------------------------------------------#
def get_consensus_sequence(bamfile, contig, start, stop):
    """
    get the consensus sequence of a given region
    """
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        seq = bamf.count_coverage(contig, start, stop)
        nn = ['A', 'C', 'G', 'T']
        return ''.join([nn[i] for i in np.argmax(seq, axis=0)])

#------------------------------------------------------------------------------#
def rev_compl(seq):
    nn = defaultdict(
        lambda: 'N', {
            'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
            'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'})
    return ''.join([nn[a] for a in seq][::-1])


#------------------------------------------------------------------------------#
# main function to joint any two breakpoints
#------------------------------------------------------------------------------#
def join_breakpoint(bamfile, a, b, min_nt = 2, seq_len = 200, \
    offset = 0, match_method = 'fuzzy_match'):
    """
    detect whether the two breakpoints are joined in the genome

    a, b - breakpoints in a tuple, such as:
        a = ('chr7', 54830975, 'L')
        b = ('chr7', 56117062, 'R')

    pattern -  mean the clip patterns of a and b breakpoints
        'LR': a - left clip; b - right clip
        'LL': a - left clip; b - left clip
        'RR': a - right clip; b - right clip

    min_qlen - minimal length of clipped bases to compare
        remove read with <= min_qlen

    seq_len - # of bases within the breakpoints to compare

    match_method - perfect_match, partial_match, fuzzy_match

    usage:
        join_breakpoint(bamfile, ('chr7', 55194959, 'L'), ('chr7', 56117062, 'R'))
    """

    if (a[2] == 'L'):
        a_seq, a_clipped = left_clip(bamfile, a, seq_len)
    else:
        a_seq, a_clipped = right_clip(bamfile, a, seq_len)

    if (b[2] == 'L'):
        b_seq, b_clipped = left_clip(bamfile, b, seq_len)
    else:
        b_seq, b_clipped = right_clip(bamfile, b, seq_len)

    # remove too short reads
    if (min_nt > 0):
        a_clipped = [i for i in a_clipped if len(i) >= min_nt]
        b_clipped = [i for i in b_clipped if len(i) >= min_nt]

    # format clip reads and breakend sequence to the same format
    if (a[2]+b[2] == 'LR'):
        _b2a = R2L(b_clipped, a_seq)
        _a2b = L2R(a_clipped, b_seq)
    elif (a[2]+b[2] == 'RL'):
        _b2a = L2R(b_clipped, a_seq)
        _a2b = R2L(a_clipped, b_seq)
    elif (a[2]+b[2] == 'LL'):
        _b2a = L2L(b_clipped, a_seq)
        _a2b = L2L(a_clipped, b_seq)
    elif (a[2]+b[2] == 'RR'):
        _b2a = R2R(b_clipped, a_seq)
        _a2b = R2R(a_clipped, b_seq)
    else:
        raise(Exception("Invalid pattern"))

    # map clip reads to breakend sequence
    if (match_method == 'perfect_match'):
        b2a = perfect_match(_b2a[0], _b2a[1], offset = offset, min_nt = min_nt)
        a2b = perfect_match(_a2b[0], _a2b[1], offset = offset, min_nt = min_nt)
    elif (match_method == 'partial_match'):
        b2a = partial_match(_b2a[0], _b2a[1], min_nt = min_nt)
        a2b = partial_match(_a2b[0], _a2b[1], min_nt = min_nt)
    elif (match_method == 'fuzzy_match'):
        b2a = fuzzy_match(_b2a[0], _b2a[1], offset = offset, min_nt = min_nt)
        a2b = fuzzy_match(_a2b[0], _a2b[1], offset = offset, min_nt = min_nt)
    else:
        raise(Exception("Invalid method"))

    return [b2a, a2b]

#------------------------------------------------------------------------------#
def left_clip(bamfile, l, seq_len):
    # get the clip read for breakpoint
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        l_read = [
            read for read in bamf.fetch(l[0], l[1], l[1]+1)
            if read.reference_start == l[1]
            and check_read_clip(read, 'soft_left')]

    # get the sequence bases within the breakpoint
    l_seq = get_consensus_sequence(bamfile, l[0], l[1], l[1]+seq_len)

    # get the clipped reads of the breakpoint
    l_clipped = [read.query_sequence[0:read.query_alignment_start]
                 for read in l_read]

#     print([read.query_sequence for read in l_read])
    return [l_seq, l_clipped]

#------------------------------------------------------------------------------#
def right_clip(bamfile, r, seq_len):
    # get the clip read for the breakpoint
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        r_read = [
            read for read in bamf.fetch(r[0], r[1], r[1]+1)
            if read.reference_end == r[1]+1
            and check_read_clip(read, 'soft_right')]

    # get the sequence bases within the breakpoint
    r_seq = get_consensus_sequence(bamfile, r[0], r[1]-seq_len+1, r[1]+1)

    # get the clipped reads of the breakpoint
    r_clipped = [read.query_sequence[read.query_alignment_end:]
                 for read in r_read]

#     print([read.query_sequence for read in r_read])
    return [r_seq, r_clipped]

#------------------------------------------------------------------------------#
def L2R(l_clipped, r_seq):
    # compare b_clipped to a_seq (L2R)
    _clipped = [ax[::-1] for ax in l_clipped]
    _seq = r_seq[::-1]
    return (_clipped, _seq)

def R2L(r_clipped, l_seq):
    # compare b_clippled to a_seq (R2L)
    _clipped = [bx for bx in r_clipped]
    _seq = l_seq
    return (_clipped, _seq)

def L2L(l_clipped, l_seq):
    # compare b_clippled to a_seq (L2L)
    _clipped = [rev_compl(bx) for bx in l_clipped]
    _seq = l_seq
    return (_clipped, _seq)

def R2R(r_clipped, r_seq):
    # compare b_clippled to a_seq (R2R)
    _clipped = [rev_compl(bx)[::-1] for bx in r_clipped]
    _seq = r_seq[::-1]
    return (_clipped, _seq)

#------------------------------------------------------------------------------#
def perfect_match(_clipped, _seq, offset = 0, min_nt = 2):
    """
    perfect match between clipped reads and breakend sequence
    for whole length

    min_nt : minimal nt match is better to be larger than 2,
        1 nt match is 25% change by random,
        2 nt match is 0.0625 change by random,
        3 nt match is 0.0156 change by random,
    """
    if (offset > 0):
        m = [re.search(_x[offset:], _seq)
             for _x in _clipped if (len(_x) > offset) and (len(_x) >= min_nt)]
    else:
        m = [re.search(_x, _seq) for _x in _clipped if (len(_x) >= min_nt)]

    mm = [i.start() for i in m if i is not None]
    if (mm):
        # get the most common split point
        return Counter(mm).most_common()[0]
    else:
        return (0, 0)

#------------------------------------------------------------------------------#
def partial_match(_clipped, _seq, min_nt = 2):
    """
    perfect match between clipped reads and breakend sequence
    for last n bases, mean the leading bases can be absent for potential
    insertions between two breakend

    """
    m = [SequenceMatcher(None, _x, _seq, autojunk = False)\
        .find_longest_match(0, len(_x), 0, len(_seq))
        for _x in _clipped if (len(_x) >= min_nt)]
    _clipped_len = [len(_x) for _x in _clipped]
    mm = [m[i].b for i in range(len(m)) \
        if (m[i].size + m[i].a == _clipped_len[i])]
    if (mm):
        # get the most common split point
        return Counter(mm).most_common()[0]
    else:
        return (0, 0)

#------------------------------------------------------------------------------#
def fuzzy_match(_clipped, _seq, offset = 0, min_nt = 2, error_rate = 0.1):
    """fuzzy match between clipped reads and breakend sequence
    use levenshtein distance to search complete sequence in _clipped
    that match part of _seq from start with error rate of 10%

    """
    match_start = 0
    if (offset > 0):
        m = [levenshtein_distance(_x[offset:], _seq)
             for _x in _clipped if (len(_x) > offset) and (len(_x) >= min_nt)]
    else:
        match_start = -offset
        _sub_seq = _seq[match_start:]
        m = [levenshtein_distance(_x, _sub_seq)
             for _x in _clipped if (len(_x) >= min_nt)]
    # print(m)
    return (match_start, sum(np.array(m) <= error_rate))


#------------------------------------------------------------------------------#
def levenshtein_distance(str1, str2, rate = True, _sub = 1, _del = 1, _ins = 1):
    """levenshtein distance
    calculate the distance between two strings
    """
    # init zero
    rows = len(str1) + 1
    cols = len(str2) + 1
    dist = np.zeros((rows, cols), dtype = int)

    # init first row and column
    for i in range(1, rows):
        for j in range(1, cols):
            dist[i][0] = i
            dist[0][j] = j

    # dynamic programming
    for col in range(1, cols):
        for row in range(1, rows):
            if str1[row - 1] == str2[col - 1]:
                mm = 0 # match
            else:
                mm = _sub # substitution
            dist[row][col] = min( # minimal penalty
                dist[row - 1][col - 1] + mm, # substitution
                dist[row - 1][col] + _del, # deletion
                dist[row][col - 1] + _ins, # insertion
            )
    # print(dist)
    if (rate):
        return dist[row].min()/len(str1)
    else:
        return dist[row].min()



#------------------------------------------------------------------------------#
# main function to search for breakpoint duo in a large list of breakpoints
#------------------------------------------------------------------------------#
def get_breakpoint_duo(bamfile, bplist, \
    seq_len = 50, seed_len = None, min_nt = 5):
    """
    get breakpoint duo from a list of breakpoints
    use pairwise assembly

    parameters:
    seq_len - # bases within breakend
    seed_len - # of bases up and down stream of the breakend
        in assembled sequence. up: seed_len; down: seed_len respectively
    """

    bpass = [breakpoint_assemble(bamfile, a, seq_len) for a in bplist]

    op = []
    for i,j in combinations(range(len(bplist)), 2):
        mtf, offset, ab = breakpoint_matcher(bpass[i], bpass[j], \
            bplist[i][2] + bplist[j][2], seq_len, seed_len, min_nt)
        if (mtf):
            op.append(bplist[i] + bplist[j] + (mtf, offset, ab))
    return op

#------------------------------------------------------------------------------#
def breakpoint_assemble(bamfile, a, seq_len = 50):
    """
    assemble sequence arround clipped breakpoint
    """
    if (a[2] == 'L'):
        a_seq, a_clipped = left_clip(bamfile, a, seq_len)
        if (a_clipped != []):
            # a_asmb = ''
            a_asmb = consensus_from_listofseq(a_clipped, a[2]) + a_seq
        else:
            a_asmb = ''
    else:
        a_seq, a_clipped = right_clip(bamfile, a, seq_len)
        if (a_clipped != []):
            # a_asmb = ''
            a_asmb = a_seq + consensus_from_listofseq(a_clipped, a[2])
        else:
            a_asmb = ''
    return a_asmb

#------------------------------------------------------------------------------#
def consensus_from_listofseq(rc, pattern = 'R'):
    """
    generate consensus sequence from a list of sequences
    """
    op = []
    length = max([len(x) for x in rc])

    if (pattern == 'R'):
        for j in range(length):
            op.append(Counter([s[j]
                      for s in rc if len(s) > j]).most_common(1)[0][0])
        return "".join(op)
    else:
        for j in range(length):
            op.append(Counter([s[::-1][j]
                      for s in rc if len(s) > j]).most_common(1)[0][0])
        return "".join(op)[::-1]


#------------------------------------------------------------------------------#
def breakpoint_matcher(a_ass, b_ass, pattern,
                       seq_len = 50, seed_len = None, min_nt = 5):
    """
    breakpoint match given two sequence

    output:
        offset:
            0  - perfect join ['ATTGC', 'GTCAT'] -> 'ATTGCGTCAT'
            1  - 1 bp gap (insertion) ['ATTGC', 'GTCAT'] -> 'ATTGCAGTCAT'
            -1 - 1 bp overlap ['ATTG', 'GTCAT'] -> 'ATTGTCAT'
    """
    # format assembled sequences to match
    # b->a = 5'->3'pattern
    if (pattern == 'RL'):
        a_ass, b_ass = b_ass, a_ass
    elif (pattern == 'RR'):
        a_ass = rev_compl(a_ass)
    elif (pattern == 'LL'):
        b_ass = rev_compl(b_ass)

    if (seed_len is not None): # match seed sequence first to increase speed
        seed_seq = b_ass[seq_len-seed_len:seq_len+seed_len]
        # seed_match = re.search(seed_seq, a_ass)
        # if (seed_match is not None):
        if (seed_seq in a_ass):
            m = SequenceMatcher(None, a_ass, b_ass, autojunk = False)
            m = m.find_longest_match(0, len(a_ass), 0, len(b_ass))
            ab = b_ass[0: m.b].lower() \
                + b_ass[m.b: m.b + m.size] \
                + a_ass[m.a + m.size:].lower() # sequence of the joint point
            is_paired = seq_len - min_nt >= m.b \
                and seq_len + min_nt <= m.b + m.size - 1
            offset = (len(a_ass) - seq_len - m.a) - (seq_len - 1 - m.b) - 1
        else:
            is_paired = False
            offset = None
            ab = None
    else: # don't use seed
        m = SequenceMatcher(None, a_ass, b_ass, autojunk = False)
        m = m.find_longest_match(0, len(a_ass), 0, len(b_ass))
        ab = b_ass[0: m.b].lower() \
            + b_ass[m.b: m.b + m.size] \
            + a_ass[m.a + m.size:].lower() # sequence of the joint point
        is_paired = seq_len - min_nt >= m.b \
            and seq_len + min_nt <= m.b + m.size - 1
        offset = (len(a_ass) - seq_len - m.a) - (seq_len - 1 -m.b) - 1

    return (is_paired, offset, ab)


#------------------------------------------------------------------------------#
def get_bppair(bamfile, bp_cand_df, \
    seq_len = 50, seed_len = 5, min_nt = 5,
    match_method = 'fuzzy_match'):
    """
    get the bppairs from bp_cand_stats (a list of bps)

    parameters:
    seq_len - # bases within breakend
    seed_len - # of bases up and down stream of the breakend
        in assembled b sequence. up: seed_len; down: seed_len respectively
    """

    bp_cand_df_sorted = bp_cand_df.sort_values(['Chrom', 'Coord', 'Clip'])
    bplist = [x[1:4] for x in bp_cand_df_sorted.itertuples()]

    # get the breakpoint pairs
    # note: searching use a different (shorter) seq_len parameter
    # to increase running efficiency
    bpduo = get_breakpoint_duo(bamfile, bplist, seq_len, seed_len, min_nt)

    # count the # of supporting reads (clip reads)
    # note: counting use a different (longer) seq_len parameter
    # to ensure the reads are fully coverred in the breakend
    tt = [list(row[0:3] + row[3:6] \
        + tuple(join_breakpoint(bamfile, row[0:3], row[3:6], \
                offset = row[7], match_method = match_method)) \
            + row[7:9])\
            for row in bpduo]

    # format output
    t2 = [row[0:6] + [row[6][1] + row[7][1]] + row[8:10] \
          for row in tt \
              if row[6][1] > 0 and row[7][1] > 0 \
              and row[6][0] == row[7][0]]

    colnames = ["Chrom1", "Coord1", "Clip1",
                "Chrom2", "Coord2", "Clip2", 'Count', 'offset', 'Seq']
    bp_pair_df = pd.DataFrame(t2, columns = colnames)

    return bp_pair_df.sort_values(
        'Count', ascending=False).reset_index(drop=True)



#------------------------------------------------------------------------------#
# extract break joint sequences
#------------------------------------------------------------------------------#
def left_clip_seq(bamfile, l, seq_len = 100):
    # get the clip read for breakpoint
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        l_read = [
            read for read in bamf.fetch(l[0], l[1], l[1]+1) \
            if check_read_clip(read, 'soft_left') \
            and read.reference_start == l[1]]

    # get the clipped reads of the breakpoint
    l_clipped_max = max([read.query_alignment_start for read in l_read])
    l_clipped = [' ' * (l_clipped_max-read.query_alignment_start) +
        read.query_sequence for read in l_read]
    l_clipped.sort()

    # get the sequence bases within the breakpoint
    l_seq = get_consensus_sequence(bamfile, l[0], l[1], l[1]+seq_len)
    l_seq = ' ' * (l_clipped_max) + l_seq

    return [l_seq] + l_clipped, l_clipped_max

#------------------------------------------------------------------------------#
def right_clip_seq(bamfile, r):
    # get the clip read for the breakpoint
    with pysam.AlignmentFile(bamfile, "rb") as bamf:
        r_read = [
            read for read in bamf.fetch(r[0], r[1], r[1]+1) \
            if check_read_clip(read, 'soft_right') \
            and read.reference_end == r[1]+1]

    # get the clipped reads of the breakpoint
    r_clipped_max = max([read.query_alignment_end for read in r_read])
    r_clipped = [' ' * (r_clipped_max-read.query_alignment_end) +
        read.query_sequence for read in r_read]
    r_clipped.sort(reverse = True)

    # get the sequence bases within the breakpoint
    r_seq = get_consensus_sequence(bamfile, r[0], r[1]-r_clipped_max+1, r[1]+1)

    return [r_seq] + r_clipped, r_clipped_max

#------------------------------------------------------------------------------#
def get_alignment(bamfile, a, b, offset = 0):
    """get alignment of bp pairs
    """
    if (a[2] == 'R'):
        a_seq, a_clip_max = right_clip_seq(bamfile, a)
    else:
        a_seq, a_clip_max = left_clip_seq(bamfile, a)

    if (b[2] == 'R'):
        b_seq, b_clip_max = right_clip_seq(bamfile, b)
    else:
        b_seq, b_clip_max = left_clip_seq(bamfile, b)

    pattern = a[2] + b[2]

    if (pattern == 'RL'):
        b_seq = [' ' * (a_clip_max - b_clip_max + offset) + read
                 for read in b_seq]
    elif (pattern == 'LR'):
        a_seq, b_seq = b_seq, a_seq
        a_clip_max, b_clip_max = b_clip_max, a_clip_max
        b_seq = [' ' * (a_clip_max - b_clip_max + offset) + read
                 for read in b_seq]
    elif (pattern == 'LL'):
        a_clipped_max = max([len(z) for z in a_seq])
        a_seq = [' ' * (a_clipped_max - len(z)) + rev_compl(z.strip())
                 for z in a_seq]
        b_seq = [' ' * (a_clipped_max - a_clip_max - b_clip_max + offset)
                 + read for read in b_seq]
    else:
        b_clipped_max = max([len(z) for z in b_seq])
        b_seq = [' ' * (b_clipped_max - len(z)) + rev_compl(z.strip())
                 for z in b_seq]
        b_seq = [' ' * (len(a_seq[0]) + offset - b_seq[0].count(' '))
                 + read for read in b_seq]

    return [a_seq[0]] + [b_seq[0]] + [] + a_seq[1:] + b_seq[1:]

#------------------------------------------------------------------------------#
def ouput_alignment(bp_pair, bamfile, output_align):
    """output alignment on pairs of breakpoints
    """
    f = open(output_align, 'w')
    for row in bp_pair.itertuples():
        if (row[9] != 'PE_Support'):
            a = row[1:4]
            b = row[4:7]
            offset = row[8]
            alignment = get_alignment(bamfile, a, b, offset)

            f.write('>' + '\t'.join(map(str, list(row[1:9]))) + "\n")
            for read in alignment:
                f.write("%s\n" % read)
            f.write("\n\n")
    f.close()




#------------------------------------------------------------------------------#
# search for breakpoint pairs candidates from improper paired reads
#------------------------------------------------------------------------------#
def get_bppair_peread(bamfile, bp_cand_df, \
    seq_len = 500, min_reads = 5):
    """
    get the bppairs from improper paired reads

    parameters:
    seq_len : # bases within breakend to search for
        improper paired reads
    min_reads : # paired reads supporting the breakpoint
    """

    # get the improper paired reads for each breakend
    bplist = []
    for row in bp_cand_df.itertuples():
        if (row[3] == 'L'):
            a = [row[1], row[2], row[2] + seq_len]
        elif (row[3] == 'R'):
            a = [row[1], row[2] - seq_len, row[2]]
        with pysam.AlignmentFile(bamfile, "rb") as bamf:
            a_read = [
                read.query_name for read in bamf.fetch(a[0], a[1], a[2])
                if (check_read(read) and (read.is_proper_pair == False))]
        if (len(a_read) >= min_reads):
            bplist.append([row[1], row[2], row[3], a_read])

    # search pairs
    op = []
    for i,j in combinations(range(len(bplist)), 2):
        if (bplist[i][0] == bplist[j][0] and
            abs(bplist[i][1] - bplist[j][1]) < seq_len):
            isc = 0
        else:
            isc = len(misc.intersect(bplist[i][3], bplist[j][3]))

        if (isc >= min_reads):
            op.append(bplist[i][0:3] + bplist[j][0:3] + [isc, 0, 'PE_Support'])
    colnames = ["Chrom1", "Coord1", "Clip1",
                "Chrom2", "Coord2", "Clip2", 'Count', 'offset', 'Seq']
    bp_pair_df = pd.DataFrame(op, columns = colnames)

    return bp_pair_df.sort_values(
        'Count', ascending=False).reset_index(drop=True)


#------------------------------------------------------------------------------#
def bp_cn_boundary(cn_amp, bp_cand_stats = None):
    """
    get breakpoints from breakpoint candidates and amplified segments
    """
    if (bp_cand_stats is not None):
        ## make breakpoints dataframe
        bplist = [i[1:4] for i in bp_cand_stats.itertuples()]
        op = []
        for row in cn_amp.itertuples():
            if ((row[1], row[2], 'L') not in bplist):
                op.append([row[1], row[2], 'L', True, row[5], 0, 0])
            if ((row[1], row[3], 'R') not in bplist):
                op.append([row[1], row[3], 'R', True, row[6], 0, 0])
        cn_seg_df = pd.DataFrame(op, columns = bp_cand_stats.columns)
        cn_seg_df

        # merge and output
        df = pd.concat([bp_cand_stats, cn_seg_df])
        df = df.drop_duplicates().sort_values(['Chrom', 'Coord', 'Clip'])
        return df.reset_index(drop = True)
    else:
        op = []
        for row in cn_amp.itertuples():
            op.append([row[1], row[2], 'L', True, row[5], 0, 0])
            op.append([row[1], row[3], 'R', True, row[6], 0, 0])
        columns = ['Chrom', 'Coord', 'Clip', 'CleanBP',
                   'ClipDepth', 'InDepth', 'OutDepth']
        cn_seg_df = pd.DataFrame(op, columns = columns)
        df = cn_seg_df.drop_duplicates().sort_values(['Chrom', 'Coord', 'Clip'])
        return df.reset_index(drop = True)
