################################################################################
# GRange class to analyze genomic intervels
################################################################################
import gzip
from itertools import compress
import pandas as pd

#------------------------------------------------------------------------------#
class GRange:
    """
    store and analyze genomic intervels
    """

    #--------------------------------------------------------------------------#
    def __init__(self, genomic_range, type = 'grange'):
        self.gr = self.create_gr(genomic_range, type)

    #--------------------------------------------------------------------------#
    def create_gr(self, genomic_range, type):
        if (type == 'bedfile'):
            return self.create_gr_from_bedf(genomic_range)
        elif (type == 'dataframe'):
            return self.create_gr_from_df(genomic_range)
        elif (type == 'dataframe_hasend'):
            return self.create_gr_from_df(genomic_range, has_end = True)
        elif (type == 'grange'):
            return genomic_range
        else:
            raise(Exception("Invalid intervel type"))

    #--------------------------------------------------------------------------#
    @staticmethod
    def create_gr_from_bedf(bedfile):
        op = []
        if (bedfile.endswith('gz')):
            bf = gzip.open(bedfile, 'rt')
        else:
            bf = open(bedfile, 'rt')
        for line in bf:
            word = line.strip().split()
            chrom = word[0]
            start = int(word[1])
            end = int(word[2])
            ext = tuple(word[3:])
            op.append((chrom, range(start, end), ext))
        bf.close()
        return op

    #--------------------------------------------------------------------------#
    @staticmethod
    def create_gr_from_df(df, has_end = False):
        op = []
        for ind, word in df.iterrows():
            chrom = word[0]
            start = int(word[1])
            if (has_end):
                end = int(word[2])
                ext = tuple(word[3:])
            else:
                end = start+1
                ext = tuple(word[2:])

            op.append((chrom, range(start, end), ext))
        return op

    #--------------------------------------------------------------------------#
    @staticmethod
    def interval_overlap(a, b, a_extend = 0, b_extend = 0):
        """
        whether a overlap with b
        extend: extend 'n' bp on both end (default: 0)

        usage:
            a = ('chr5', range(300, 500))
            b = ('chr5', range(100, 600))
            gr_overlap(a, b)
        """
        # get the range object
        x = a[1]
        y = b[1]

        # extend range on both end
        x = range(x.start-a_extend, x.stop+a_extend)
        y = range(y.start-b_extend, y.stop+b_extend)

        if (a[0] != b[0]):
            return False
        elif (x.start == x.stop or y.start == y.stop):
            return False
        else:
            return ((x.start < y.stop  and x.stop > y.start) or
                (x.stop  > y.start and y.stop > x.start))

    #--------------------------------------------------------------------------#
    def intersect(self, b, a_extend = 0, b_extend = 0, invert = False):
        """
        which intervals in a overlap with intervals in b
        extend: extend 'n' bp on both end (default: 0)
        """
        a = self.gr
        b = b.gr
        bool_list = []
        for x in a:
            tf = False
            for y in b:
                if (self.interval_overlap(x, y, a_extend, b_extend)):
                    tf = True
                    break
            bool_list.append(tf)
        return self.gr_subset(bool_list, invert)

    #--------------------------------------------------------------------------#
    def gmap(self, b, a_extend = 0, b_extend = 0, multi_hit = False):
        """
        map score in b to a if they overlap
        extend: extend 'n' bp on both end (default: 0)
        multi_hit: False if only expect one hit for each grange in a
                True if expect multiple hits for each grange in a

        todo: specify which column to map
        """
        a = self.gr
        b = b.gr
        map_list = []
        for x in a:
            map_item = None

            if (multi_hit == False):
                for y in b:
                    if (self.interval_overlap(x, y, a_extend, b_extend)):
                        map_item = y[2]
                        break
            else:
                map_item = [y[2] for y in b \
                    if (self.interval_overlap(x, y, a_extend, b_extend))]

            map_list.append(map_item)

        return map_list

    #--------------------------------------------------------------------------#
    def gr_subset(self, bool_list, invert = False):
        """
        subset a GRange by a bool list
        """

        if (invert == False):
            return GRange(list(compress(self.gr, bool_list)))
        else:
            return GRange(list(compress(self.gr, [not i for i in bool_list])))

    #--------------------------------------------------------------------------#
    @staticmethod
    def range_intersect_2list(a, b):
        """
        get intersections between two list of ranges
        range in format of [0, 10] refers to range(0, 10)
        """
        ranges = []
        i = j = 0
        while i < len(a) and j < len(b):
            a_left, a_right = a[i]
            b_left, b_right = b[j]

            if a_right < b_right:
                i += 1
            else:
                j += 1

            if a_right >= b_left and b_right >= a_left:
                end_pts = sorted([a_left, a_right, b_left, b_right])
                middle = [end_pts[1], end_pts[2]]
                ranges.append(middle)

        ri = 0
        while ri < len(ranges)-1:
            if ranges[ri][1] == ranges[ri+1][0]:
                ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]

            ri += 1

        return ranges

    #--------------------------------------------------------------------------#
    @staticmethod
    def map_gene(df, bed12):
        """
        search genes contained segments

        df: bedlike dataframe
        bed12: gene annotation in bed12 format
        """
        gr1 = GRange(df, 'dataframe_hasend')
        gr2 = GRange(bed12, 'bedfile')
        gr = gr1.gmap(gr2, multi_hit = True)
        return df.assign(Gene = ['; '.join(set([i[0] for i in row])) for row in gr])

    #--------------------------------------------------------------------------#
    def gr_to_df(self):
        return self.gr[0:n]

    #--------------------------------------------------------------------------#
    def gr_head(self, n=5):
        return self.gr[0:n]

    #--------------------------------------------------------------------------#
    def __repr__(self):
        return '# of intervals: ' + str(len(self.gr)) + '\n' \
            + str([i for i in self.gr_head()])
    def __str__(self):
        return '# of intervals: ' + str(len(self.gr)) + '\n' \
            + str([i for i in self.gr_head()])
