################################################################################
# misc functions
################################################################################
#------------------------------------------------------------------------------#
import grange

#------------------------------------------------------------------------------#
def unique(seq):
    """
    unique a list by preserve the order
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


#------------------------------------------------------------------------------#
def all_max_index(a):
    """
    get all the index of the max
    """
    if (not a):
        return []
    idx = [0]
    max_ = a[0]
    for i in range(1, len(a)):
        if (a[i] > max_):
            idx = [i]
            max_ = a[i]
        elif (a[i] == max_):
            idx.append(i)
    return idx


#------------------------------------------------------------------------------#
def intersect(lst1, lst2):
    """intersection of two lists
    """
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


#------------------------------------------------------------------------------#
def map_cancergene(df, cgfile):
    # read cancer gene file
    with open(cgfile, 'rt') as f:
        next(f)
        cgc = []
        for line in f:
            word = line.strip().split()
            gene = word[0]
            cgc.append(gene)

    # intersect
    return df.assign(
        CancerGene = ['; '.join(intersect(row['Gene'].split("; "), cgc))
        for idx, row in df.iterrows()])



#--------------------------------------------------------------------------#
def merge_bed(amp_df, gap = 100000):
    """merge adjacent from a data frame of segments
    """
    # convert to list
    rois = [[row[0], row[1], row[2]] for idx, row in amp_df.iterrows()]
    rois_previous = []
    while rois != rois_previous:
        rois_previous = rois
        rois = merge_bed_helper(rois, gap)

    return rois


#--------------------------------------------------------------------------#
def merge_bed_helper(regions, gap = 100000):
    # search overlapping segments in a dataframe for 1 round
    rois = []
    for row in regions:
        a = (row[0], range(row[1], row[2]))
        a_ext = gap

        if (rois == []):
            rois.append([row[0], row[1], row[2]])
        else:
            is_hit = 0
            for i, roi in enumerate(rois):
                b = (roi[0], range(roi[1], roi[2]))
                b_ext = gap
                if grange.GRange.interval_overlap(a, b, a_ext, b_ext):
                    rois[i] = [roi[0], min(row[1], roi[1]), max(row[2], roi[2])]
                    is_hit = 1
                    break
            if (is_hit == 0):
                rois.append([row[0], row[1], row[2]])
    return rois
