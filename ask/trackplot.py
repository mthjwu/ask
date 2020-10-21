################################################################################
# amplicon plot
################################################################################
#------------------------------------------------------------------------------#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Polygon, Rectangle
from matplotlib import cm, colors
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy
import numpy as np
import pandas as pd
import grange
import gzip

#------------------------------------------------------------------------------#
# processing functions
#------------------------------------------------------------------------------#
def prepare_plot(amp_df, gbeds, bincnt, binsize = 10000, ext = 0.3,
                 genre = 'circular'):
    """prepare data for track plot
    ext : surrounding regions (%) to be plotted
    genre : 'circular', 'linear' or 'cn'
    """

    # # get roi
    # print(amp_df)
    # roi = [amp_df.Chrom.iloc[0], min(amp_df.Start), max(amp_df.End)]
    # roi_len = roi[2] - roi[1]
    # offset = max([roi_len * ext, binsize * 10])
    # roi[1] = int(roi[1] - offset)
    # roi[2] = int(roi[2] + offset)
    # print(roi)

    # get rois
    rois = get_rois(amp_df, binsize * 10, ext)

    # get gene annotation
    genes = set([x for g in list(filter(None, amp_df['Gene']))
        for x in g.split('; ')])
    beds = [g for g in gbeds if (g.name in genes)]

    # get segments
    if (genre == 'cn'):
        segs = [[row[1], row[2], row[3], row[0], row[4], '.']
            for row in amp_df.itertuples()]
    else:
        segs = [[row[1], row[2], row[3], row[0], row[6], row[4]]
            for row in amp_df.itertuples()]
    segs = [Bed(x) for x in segs]

    # get links
    if (genre == 'cn'):
        links = [Link(['virtual', 0, 'virtual', 1, 0])]
    else:
        s = amp_df['Strand'] == '-'
        amp_df.loc[s,['Start', 'End']] = amp_df.loc[s,['End', 'Start']].values
        # idx_1 = list(amp_df.index)
        idx_1 = list(range(0, len(amp_df.index)))
        idx_2 = idx_1[1:] + [idx_1[0]]
        if (genre == 'linear'):
            idx_1 = idx_1[:-1]
            idx_2 = idx_2[:-1]
        links = [[amp_df.iloc[v, 0], amp_df.iloc[v, 2],
            amp_df.iloc[idx_2[i], 0], amp_df.iloc[idx_2[i], 1],
            amp_df.iloc[v, 4]] for i, v in enumerate(idx_1)]
        links = [Link(x) for x in links]
    if (links == []):
        links = [Link(['virtual', 0, 'virtual', 1, 0])]

    # get bin counts
    bg_list = []
    for roi in rois:
        bcnt_sub = bincnt[(bincnt['Chrom'] == roi[0]) &
                          (bincnt['Coord'] >= roi[1]) &
                          (bincnt['Coord'] <= roi[2])]
        op = bcnt_sub[['Chrom', 'Coord', 'CN']]
        bg_list.append(op.assign(Coord = op.Coord + binsize / 2))
    bgs = pd.concat(bg_list)

    return rois, segs, links, bgs, beds

#--------------------------------------------------------------------------#
def plot_amp(rois, segs, links, bgs, beds, fig_out = 'test.pdf',
             fig_width = 15, fontsize = 12):
    """plot amplicon structure
    """

    # get roi total length
    roi_len = sum([roi[2]-roi[1] for roi in rois])

    # get length per char
    char_len = get_font_length(fig_width, roi_len, fontsize)

    # get # of rows to plot
    segs_nrows = get_max_row(segs, char_len)
    beds_nrows = get_max_row(beds, char_len)

    # set the height for segs and beds
    segs_height = segs_nrows / 3
    beds_height = beds_nrows / 5

    # set figure height
    fig_height = 3.3 + segs_height + beds_height

    # layout
    gridspec = dict(hspace=0,
        height_ratios=[1, 0, 1, segs_height, 0, beds_height])
    fig, axs = plt.subplots(6, 1, figsize=(fig_width, fig_height),
                            gridspec_kw = gridspec)

    # plot tracks
    Bgtrack(axs[0], rois, bgs, char_len).plot() # bin count track
    axs[1].set_visible(False) # 10% space
    Arctrack(axs[2], rois, links, char_len, 'Reds').plot() # arc track
    Bedtrack(axs[3], rois, segs, char_len, 'Dark2').plot() # segment track
    axs[4].set_visible(False) # 10% space
    Genetrack(axs[5], rois, beds, char_len).plot() # gene track

    # aligh tracks
    poss = [ax.get_position() for ax in axs]
    x_pos = poss[0].x0
    x_width = poss[0].width
    for i in range(1, len(axs)):
        axs[i].set_position([x_pos, poss[i].y0, x_width, poss[i].height])

    # output
    plt.savefig(fig_out)

#--------------------------------------------------------------------------#
def read_bedf(bedfile):
    op = []
    if (bedfile.endswith('gz')):
        bf = gzip.open(bedfile, 'rt')
    else:
        bf = open(bedfile, 'rt')
    for line in bf:
        word = line.strip().split('\t')
        op.append(Bed(word))
    bf.close()
    return op

#--------------------------------------------------------------------------#
def get_font_length(fig_width, region_width, fontsize):
    """get the length of one char
    """
    inches_per_pt = 1.0 / 72.27
    font_in_inches = fontsize * inches_per_pt
    unit_per_inch = region_width / (fig_width * 0.8)
    font_in_unit = font_in_inches * unit_per_inch
    return font_in_unit

#--------------------------------------------------------------------------#
def get_max_row(beds, char_len, max_row = 30):
    """Get the # of rows to plot for bed or gene track
    """

    # arrange ypos
    ys = []
    ends = [0] * max_row # end position of each row
    for bed in beds:
        offset = char_len * (len(bed.name))
        for idx, end in enumerate(ends):
            if (end == 0 or bed.start > end + offset):
                # bed.ypos = idx
                ys.append(idx)
                ends[idx] = bed.end
                break

    if (ys == []):
        n_row = 1
    else:
        n_row = max(ys) + 1
    return n_row


#--------------------------------------------------------------------------#
def get_rois(amp_df, offset = 100000, ext = 0.3):
    """get rois from a data frame of segments
    """
    # convert to list
    rois = [[row[0], row[1], row[2]] for idx, row in amp_df.iterrows()]
    rois_previous = []
    while rois != rois_previous:
        rois_previous = rois
        rois = get_rois_helper(rois, offset, ext)

    # output extended rois
    rois_ext = []
    for roi in rois:
        ext_len = max([(roi[2] - roi[1]) * ext, offset])
        rois_ext.append([roi[0], round(roi[1] - ext_len),
                         round(roi[2] + ext_len)])

    return rois_ext


#--------------------------------------------------------------------------#
def get_rois_helper(regions, offset = 100000, ext = 0.3):
    # search overlapping segments in a dataframe for 1 round
    rois = []
    for row in regions:
        a = (row[0], range(row[1], row[2]))
        a_ext = round(max([(row[2]-row[1]) * ext, offset]))

        if (rois == []):
            rois.append([row[0], row[1], row[2]])
        else:
            is_hit = 0
            for i, roi in enumerate(rois):
                b = (roi[0], range(roi[1], roi[2]))
                b_ext = round(max([(roi[2]-roi[1]) * ext, offset]))
                if grange.GRange.interval_overlap(a, b, a_ext, b_ext):
                    rois[i] = [roi[0], min(row[1], roi[1]), max(row[2], roi[2])]
                    is_hit = 1
                    break
            if (is_hit == 0):
                rois.append([row[0], row[1], row[2]])
    return rois


#--------------------------------------------------------------------------#
def get_sub_df(amp_df, roi):
    """get sub dataframe from a roi
    """

    a = (roi[0], range(roi[1], roi[2]))
    tf = [grange.GRange.interval_overlap(
            a, (row[0], range(row[1], row[2]))
        ) for idx, row in amp_df.iterrows()]
    return amp_df[tf]


#--------------------------------------------------------------------------#
def get_rois_with_score(amp_df, offset = 100000, ext = 0.3):
    """get rois from a data frame of segments
    """
    # convert to list
    i = amp_df.columns.get_loc("CN")
    rois = [[row[0], row[1], row[2], row[i]] for idx, row in amp_df.iterrows()]
    rois_previous = []
    while rois != rois_previous:
        rois_previous = rois
        rois = get_rois_with_score_helper(rois, offset, ext)

    # output extended rois
    rois_ext = []
    for roi in rois:
        ext_len = max([(roi[2] - roi[1]) * ext, offset])
        rois_ext.append([roi[0], round(roi[1] - ext_len),
                         round(roi[2] + ext_len), roi[3]])

    return rois_ext


#--------------------------------------------------------------------------#
def get_rois_with_score_helper(regions, offset = 100000, ext = 0.3):
    # search overlapping segments in a dataframe for 1 round
    rois = []
    for row in regions:
        a = (row[0], range(row[1], row[2]))
        a_ext = round(max([(row[2]-row[1]) * ext, offset]))
        a_score = row[3]

        if (rois == []):
            rois.append([row[0], row[1], row[2], row[3]])
        else:
            is_hit = 0
            for i, roi in enumerate(rois):
                b = (roi[0], range(roi[1], roi[2]))
                b_ext = round(max([(roi[2]-roi[1]) * ext, offset]))
                b_score = row[3]
                if grange.GRange.interval_overlap(a, b, a_ext, b_ext):
                    rois[i] = [roi[0], min(row[1], roi[1]),
                               max(row[2], roi[2]), max(a_score, b_score)]
                    is_hit = 1
                    break
            if (is_hit == 0):
                rois.append([row[0], row[1], row[2], row[3]])
    return rois


# #--------------------------------------------------------------------------#
# def get_rois(amp_df, offset = 100000, ext = 0.3):
#     """get rois from a data frame of segments
#     """

#     rois = []
#     for df in group_segment(amp_df, offset, ext):
#         row = [df.iloc[0, 0], min(df.Start), max(df.End)]
#         row_len = row[2] - row[1]
#         ext_len = max([row_len * ext, offset])
#         rois.append([row[0], int(row[1] - ext_len), int(row[2] + ext_len)])

#     return rois

# #--------------------------------------------------------------------------#
# def group_segment(amp_df, offset = 100000, ext = 0.3):
#     """group segments by proximity
#     """

#     df_sort = amp_df.sort_values(['Chrom', 'Start', 'End']).copy()
#     grp = [[0]]
#     for i, idx in enumerate(list(df_sort.index)[:-1]):
#         a = (df_sort.Chrom.iloc[i],
#             range(df_sort.Start.iloc[i], df_sort.End.iloc[i]))
#         b = (df_sort.Chrom.iloc[i+1],
#             range(df_sort.Start.iloc[i+1], df_sort.End.iloc[i+1]))
#         if (df_sort.Chrom.iloc[i] == df_sort.Chrom.iloc[i+1] and
#             abs(df_sort.Start.iloc[i+1] - df_sort.End.iloc[i]) < offset) \
#             or grange.GRange.interval_overlap(a, b):
#             for g in grp:
#                 if (i in g):
#                     g.append(i+1)
#         else:
#             grp.append([i+1])

#     return [df_sort.iloc[g] for g in grp]

#--------------------------------------------------------------------------#
def map_region(igr, ori, cvt):
    """map regions between original and converted
    """
    return [cvt[0] + x - ori[0] for x in igr]

# def group_segment(amp_df, offset = 100000):
#     """group segments by proximity
#     """

#     grp = [[0]]
#     for i in list(amp_df.index)[:-1]:
#         if (amp_df.Chrom.iloc[i] == amp_df.Chrom.iloc[i+1]
#             and abs(amp_df.Start.iloc[i+1] - amp_df.End.iloc[i]) < offset):
#             for g in grp:
#                 if (i in g):
#                     g.append(i+1)
#         else:
#             grp.append([i+1])

#     return [amp_df.iloc[g] for g in grp]


#------------------------------------------------------------------------------#
# bed class
#------------------------------------------------------------------------------#
class Bed:

    #--------------------------------------------------------------------------#
    def __init__(self, bed):
        # bed3
        if (len(bed) >= 3):
            self.chrom = bed[0]
            self.start = int(bed[1])
            self.end = int(bed[2])
            self.type = 'bed3'
        else:
            raise(Exception("Invalid bed type"))

        # bed6
        if (len(bed) >= 6):
            self.name = str(bed[3])
            self.score = int(bed[4])
            self.strand = bed[5]
            self.type = 'bed6'

        # bed12
        if (len(bed) >= 12):
            self.thick_start = int(bed[6])
            self.thick_end = int(bed[7])
            self.item_rgb = bed[8]
            self.block_count = int(bed[9])
            self.block_sizes = [int(x) for x in bed[10].strip(',').split(',')]
            self.block_starts = [int(x) for x in bed[11].strip(',').split(',')]
            self.type = 'bed12'

    #--------------------------------------------------------------------------#
    def print(self):
        print(self.chrom, self.start, self.end, self.name, self.score,
            self.strand, self.thick_start, self.thick_end, self.item_rgb,
            self.block_count, self.block_sizes, self.block_starts)



#------------------------------------------------------------------------------#
# link class
#------------------------------------------------------------------------------#
class Link:

    #--------------------------------------------------------------------------#
    def __init__(self, link):
        self.chrom1 = link[0]
        self.start1 = int(link[1])
        self.chrom2 = link[2]
        self.start2 = int(link[3])
        self.score = int(link[4])



#------------------------------------------------------------------------------#
# plot bedgraph class
#------------------------------------------------------------------------------#
class Bgtrack:

    #--------------------------------------------------------------------------#
    def __init__(self, ax, roi, bgs, char_len):
        self.ax = ax
        self.roi = roi
        self.bgs = bgs
        self.ymax = int(max(bgs.iloc[:,2])) + 1

        # init converted axis
        self.bgs_cvt = self.bgs
        self.roi_cvt = self.roi

        # get the cutpoints of split segments
        self.cutpoints = []

        # get the gap width between split segments
        self.gap = char_len * 2

        # convert bgs to converted axis
        if (len(roi[0]) > 1):
            self.bgs_cvt, self.roi_cvt, self.cutpoints = \
                self.convert_bgs(bgs, roi, self.gap)

        # plot features
        self.edgecolor = 'black'
        self.facecolor = 'gray'
        self.alpha = 0.8

    #--------------------------------------------------------------------------#
    def plot(self):
        """plot bedgraph track
        """
        ax = self.ax
        bgs = self.bgs_cvt
        ymax = self.ymax
        roi = self.roi_cvt
        cutpoints = self.cutpoints
        facecolor = self.facecolor
        alpha = self.alpha

        # set axis
        ax.set_xlim(xmin = roi[1], xmax = roi[2])
        ax.set_ylim(ymin = 0, ymax = ymax)
        ax.set_xticks([])
        ax.set_yticks([0, ymax])
        ax.set_ylabel('Copy #', rotation = 'horizontal', ha = 'right')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        [ax.axvline(x, linestyle = '--', alpha = 0.5,
                    zorder = -1) for x in cutpoints]

        # plot track
        ax.fill_between(bgs.iloc[:,1], bgs.iloc[:,2], linewidth = 0,
                        color = facecolor, alpha = alpha)
        # ax.bar(bgs.iloc[:,1], bgs.iloc[:,2], linewidth = 5, width=0, ec="k")


    #--------------------------------------------------------------------------#
    @staticmethod
    def convert_bgs(bgs, rois, gap):
        """ convert bgs to new x-axis
        """
        # init
        bgs_cvt = []
        current_start = gap
        current_end = gap
        cutpoints = []

        # convert
        for roi in rois:
            current_end = current_start + roi[2] - roi[1]
            cutpoints.append(current_start)
            cutpoints.append(current_end)

            bgs_sub = bgs[(bgs['Chrom'] == roi[0]) &
                          (bgs['Coord'] >= roi[1]) &
                          (bgs['Coord'] <= roi[2])]
            cvt = np.interp(bgs_sub.Coord, [roi[1], roi[2]],
                            [current_start, current_end])
            bgs_cvt.append(bgs_sub.assign(Coord = cvt))
            bgs_cvt.append(pd.DataFrame([[roi[0], current_end, 0]],
                columns = ['Chrom', 'Coord', 'CN']))
            bgs_cvt.append(pd.DataFrame([[roi[0], current_start, 0]],
                columns = ['Chrom', 'Coord', 'CN']))
            current_start = current_end + gap
        bgs_cvt = pd.concat(bgs_cvt).sort_values(['Chrom', 'Coord'])
        return bgs_cvt, [rois[0][0], 0, current_end + gap], cutpoints



#------------------------------------------------------------------------------#
# plot arc class
#------------------------------------------------------------------------------#
class Arctrack:

    #--------------------------------------------------------------------------#
    def __init__(self, ax, roi, links, char_len, colormap = None):
        self.ax = ax
        self.roi = roi
        self.links = links
        self.roi_cvt = roi
        self.links_cvt = links
        self.char_len = char_len
        self.gap = char_len * 2
        self.cutpoints = []

        # convert links to new axis
        if (len(roi[0]) > 1):
            self.links_cvt, self.roi_cvt, self.cutpoints = \
                self.convert_links(links, roi, self.gap)

        # set color
        if (colormap is not None):
            self.sm = self.map_color([x.score for x in links], colormap)
        else:
            self.sm = None

        # plot features
        self.color = 'red'
        self.linewidth = 1
        self.fontsize = 12

    #--------------------------------------------------------------------------#
    def plot(self):
        """plot arc track
        """
        ax = self.ax
        links = self.links_cvt
        roi = self.roi_cvt
        color = self.color
        fontsize = self.fontsize
        linewidth = self.linewidth
        cutpoints = self.cutpoints
        char_len = self.char_len

        ax.set_xlim(xmin = roi[1], xmax = roi[2])
        ax.set_ylabel('Breakjoin', rotation = 'horizontal', ha = 'right')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set(frame_on=False)
        [ax.axvline(x, linestyle = '--', alpha = 0.5,
                    zorder = -1) for x in cutpoints]

        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes('top', size='10%', pad='5%')

        # # draw straight line
        # ax.plot(0, 0, zorder = 0)
        # plt.colorbar(self.sm, cax = cax, orientation='horizontal')
        # cax.xaxis.set_ticks_position("top")

        # draw straight line
        ax.plot(0, 0, zorder = 0)

        # draw arcs
        heights = []
        for lk in links:
            center = (lk.start1 + lk.start2)/2
            width = (lk.start2 - lk.start1)
            if (self.sm is not None):
                color = self.sm.to_rgba(lk.score)
            # height = np.log2(abs(width))*10
            height = np.sqrt(abs(width) + 10000)
            heights.append(height / 2)

            # plot arc
            ax.add_patch(Arc((center, 1), width, height, 0, 0, 180,
                              color = color, linewidth = linewidth))

            # plot arc score
            ax.text(center, height / 2, lk.score,
                    horizontalalignment = 'center',
                    verticalalignment = 'bottom',
                    fontsize = fontsize)

        # set the y-axis
        ax.set_ylim(ymin = 0, ymax = max(heights) * 1.3)

        # cbar = plt.colorbar(self.sm, ax = ax, orientation='vertical',
        #     pad = 0, fraction = 0.1, shrink = 0.9, aspect = 10,
        #     ticks = [0, max([x.score for x in self.links])])
        # cbar.ax.tick_params(labelsize = self.fontsize * 0.7)

    #--------------------------------------------------------------------------#
    @staticmethod
    def map_color(scores, colormap):
        """map color to given score
        """
        if (scores == []):
            scores = [1]
        norm = colors.Normalize(0, max(scores))
        cmap = cm.get_cmap(colormap)
        return cm.ScalarMappable(norm, cmap)

    #--------------------------------------------------------------------------#
    @staticmethod
    def convert_links(links, rois, gap):
        """ convert links to new x-axis
        """
        # init
        links_cvt = deepcopy(links)
        current_start = gap
        current_end = gap
        cutpoints = []

        # convert
        for roi in rois:
            current_end = current_start + roi[2] - roi[1]
            cutpoints.append(current_start)
            cutpoints.append(current_end)

            for i, lk in enumerate(links_cvt):
                lk_old = links[i]
                if (lk_old.chrom1 == roi[0] \
                        and lk_old.start1 >= roi[1] \
                        and lk_old.start1 <= roi[2]):
                    lk.start1 = np.interp(
                        lk_old.start1, [roi[1], roi[2]],
                        [current_start, current_end])
                if (lk_old.chrom2 == roi[0] \
                        and lk_old.start2 >= roi[1] \
                        and lk_old.start2 <= roi[2]):
                    lk.start2 = np.interp(
                        lk_old.start2, [roi[1], roi[2]],
                        [current_start, current_end])
            current_start = current_end + gap

        return links_cvt, [rois[0][0], 0, current_end + gap], cutpoints



#------------------------------------------------------------------------------#
# plot gene class
#------------------------------------------------------------------------------#
class Genetrack:

    #--------------------------------------------------------------------------#
    def __init__(self, ax, roi, beds, char_len):
        self.ax = ax
        self.roi = roi
        self.beds = beds
        self.char_len = char_len
        self.roi_cvt = roi
        self.beds_cvt = beds
        self.cutpoints = []
        self.gap = char_len * 2

        # convert beds to new axis
        if (len(roi[0]) > 1):
            self.beds_cvt, self.roi_cvt, self.cutpoints = \
                self.convert_beds(beds, roi, self.gap)
        else:
            self.roi_cvt = self.roi = [roi]

        # plot features
        self.height = 0.7
        self.edgecolor = 'black'
        self.facecolor = 'lightgrey'
        self.linewidth = 0.5
        self.fontsize = 12
        self.max_row = 30
        self.arrow_length = self.char_len

    #--------------------------------------------------------------------------#
    def plot(self):
        """plot gene track
        """
        ax = self.ax
        beds = self.beds_cvt
        roi_cvt = self.roi_cvt
        rois = self.roi
        cutpoints = self.cutpoints
        gap = self.gap

        # set axis
        ax.set_xlim(xmin = roi_cvt[1], xmax = roi_cvt[2])
        ax.set_ylabel('Gene', rotation = 'horizontal', ha = 'right')
        # ax.set_xlabel(roi_cvt[0])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        [ax.axvline(x, linestyle = '--', alpha = 0.5,
                    zorder = -1) for x in cutpoints]

        # set x-axis ticks to actual genomic region
        cvt_start = gap
        cvt_end = gap
        x_actual = []
        x_convert = []
        for rg in rois:
            cvt_end = cvt_start + rg[2] - rg[1] - 1
            x_convert.append(cvt_start)
            x_convert.append((cvt_start + cvt_end) / 2)
            x_convert.append(cvt_end)
            x_actual.append(f'{rg[1]:,}')
            x_actual.append('\n\n' + rg[0])
            x_actual.append('\n' + f'{rg[2]:,}')
            cvt_start = cvt_end + 1 + gap
            # ax.text((rg[1] + rg[2]) / 2, 1, rg[0])
        ax.set_xticks(x_convert)
        ax.set_xticklabels(x_actual)

        # get y position of each gene
        ys = self.get_gene_ypos()

        # plot track
        for idx, bed in enumerate(beds):
            # print(bed.start, bed.end, bed.name, ys[idx])
            self.draw_gene(bed, ys[idx])

    #--------------------------------------------------------------------------#
    def draw_gene(self, bed, ypos):
        """draw gene mode using one row of bed12
        """
        ax = self.ax
        height = self.height
        edgecolor = self.edgecolor
        facecolor = self.facecolor
        linewidth = self.linewidth
        fontsize = self.fontsize

        # draw straight line
        ax.plot([bed.start, bed.end], [ypos + height / 2, ypos + height / 2],
                color = edgecolor, linewidth = linewidth, zorder = 0)

        # get the first exon
        if (bed.strand == '-'):
            block_starts = bed.block_starts[::-1].copy()
            block_sizes = bed.block_sizes[::-1].copy()
            label_pos = bed.start - self.arrow_length * 1.3
        else:
            block_starts = bed.block_starts.copy()
            block_sizes = bed.block_sizes.copy()
            label_pos = bed.start - self.arrow_length * 0.3
        exon1_start = bed.start + block_starts.pop()
        exon1_end = exon1_start + block_sizes.pop()

        # get arrow xy
        arrow_xy = self.get_arrow_xy(
            exon1_start, exon1_end, bed.strand, ypos)

        # draw arrow head
        ax.add_patch(Polygon(arrow_xy,
                            edgecolor = edgecolor,
                            facecolor = facecolor,
                            linewidth = linewidth))

        # draw other parts of gene
        for idx, rg in enumerate(block_starts):
            start_pos = bed.start + rg
            end_pos = start_pos + block_sizes[idx]
            xy = [(start_pos, ypos), (start_pos, ypos + height),
                (end_pos, ypos + height), (end_pos, ypos)]
            ax.add_patch(Polygon(xy, edgecolor = edgecolor,
                facecolor = facecolor, linewidth = linewidth))

        # # plot gene symbol
        # print(label_pos, bed.start, bed.end, self.arrow_length, ypos)
        ax.text(label_pos, ypos + height / 2,
                bed.name, horizontalalignment = 'right',
                verticalalignment = 'center', fontsize = fontsize)

    #--------------------------------------------------------------------------#
    def get_arrow_xy(self, start, end, strand, ypos):
        """get arrow xy locations to indicate gene direction

        Parameters
        ----------
        ypos : y-axis position of the arrow
        length : length of the arrow head,
            relative to total genomic region to plot
        height : height of the arrow
        """
        height = self.height
        length = self.arrow_length

        # calculate y
        y0 = ypos
        y1 = ypos + height

        # calculate x
        if (strand == '+'):
            x0 = start
            x1 = end
            x2 = end + length
        elif (strand == '-'):
            x0 = end
            x1 = start
            x2 = start - length

        # get the xy
        xy = [(x0, y0), (x0, y1), (x1, y1),
            (x2, y0 + height / 2), (x1, y0)]

        return xy

    #--------------------------------------------------------------------------#
    def get_gene_ypos(self):
        """determine ypos for all genes to avoid overlap
        """
        beds = self.beds_cvt
        max_row = self.max_row
        char_len = self.char_len

        ys = []
        ends = [0] * max_row # end position of each row
        for bed in beds:
            offset = char_len * (len(bed.name))
            for idx, end in enumerate(ends):
                if (end == 0 or bed.start > end + offset):
                    # bed.ypos = idx
                    ys.append(idx)
                    ends[idx] = bed.end
                    break
                # any genes not fit in plot in the last row
                elif (idx == max_row - 1):
                    ys.append(max_row - 1)
                    ends[idx] = bed.end
        return [max(ys) - y for y in ys]

    #--------------------------------------------------------------------------#
    @staticmethod
    def convert_beds(beds, rois, gap):
        """ convert beds to new x-axis
        """
        # init
        # beds_cvt = deepcopy(beds)
        current_start = gap
        current_end = gap
        cutpoints = []

        # convert
        beds_cvt = []

        for roi in rois:
            current_end = current_start + roi[2] - roi[1]
            cutpoints.append(current_start)
            cutpoints.append(current_end)

            for i, bed in enumerate(beds):
                bed_cvt = deepcopy(bed)
                bed_old = beds[i]
                a = (bed_old.chrom, range(bed_old.start, bed_old.end))
                b = (roi[0], range(roi[1], roi[2]))
                if grange.GRange.interval_overlap(a, b):
                    bed_cvt.start, bed_cvt.end = map_region(
                        [bed_old.start, bed_old.end],
                        [roi[1], roi[2]],
                        [current_start, current_end])
                    # if (hasattr(bed_old, 'block_starts')):
                    #     bed_cvt.block_starts = map_region(
                    #     bed_old.block_starts,
                    #     [roi[1], roi[2]],
                    #     [current_start, current_end])
                    beds_cvt.append(bed_cvt)
                    # print(bed_cvt.start, bed_cvt.end,
                    #     bed_old.start, bed_old.end, roi[1], roi[2],
                    #     current_start, current_end, bed_old.name)

            current_start = current_end + gap

        return beds_cvt, [rois[0][0], 0, current_end + gap], cutpoints


#------------------------------------------------------------------------------#
# plot bed class
#------------------------------------------------------------------------------#
class Bedtrack(Genetrack):

    #--------------------------------------------------------------------------#
    def __init__(self, ax, roi, beds, char_len, colormap = 'Dark2'):
        super().__init__(ax, roi, beds, char_len)

        # # set specific object
        # self.beds = beds

        # set color
        if (colormap is not None):
            self.lc = cm.get_cmap(colormap)
        else:
            self.lc = None

    #--------------------------------------------------------------------------#
    def plot(self):
        """plot bed track
        """
        ax = self.ax
        lc = self.lc
        beds = self.beds_cvt
        roi = self.roi_cvt
        cutpoints = self.cutpoints

        # set axis
        ax.set_ylabel('Segment', rotation = 'horizontal', ha = 'right')
        ax.set_xlim(xmin = roi[1], xmax = roi[2])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set(frame_on=False)
        [ax.axvline(x, linestyle = '--', alpha = 0.5,
                    zorder = -1) for x in cutpoints]

        # get y positions
        ys = self.get_gene_ypos()

        # plot track
        for idx, bed in enumerate(beds):
            if (lc is None):
                self.draw_region(bed, ys[idx])
            else:
                self.draw_region(bed, ys[idx], lc(idx))

    #--------------------------------------------------------------------------#
    def draw_region(self, bed, ypos, color = None):
        """draw region track
        """
        ax = self.ax
        height = self.height
        edgecolor = self.edgecolor
        facecolor = self.facecolor
        linewidth = self.linewidth
        fontsize = self.fontsize

        # set ractangle color
        if (color is not None):
            facecolor = color

        # draw straight line
        ax.plot([bed.start, bed.end], [ypos + height / 2, ypos + height / 2],
                color = edgecolor, linewidth = linewidth, zorder = 0)

        # draw bed region
        if (bed.strand not in ['+', '-']):
            ax.add_patch(Rectangle((bed.start, ypos),
                        bed.end - bed.start, height,
                        edgecolor = edgecolor,
                        facecolor = facecolor,
                        linewidth = linewidth))
        else:
            arrow_xy = self.get_arrow_xy(bed.start, bed.end, bed.strand, ypos)
            ax.add_patch(Polygon(arrow_xy,
                                edgecolor = edgecolor,
                                facecolor = facecolor,
                                linewidth = linewidth))

        # plot region name
        ax.text((bed.start + bed.end) / 2, ypos + height / 2,
                bed.name, horizontalalignment = 'center',
                verticalalignment = 'center', fontsize = fontsize)

    #--------------------------------------------------------------------------#
    def get_gene_ypos(self):
        return super().get_gene_ypos()

    #--------------------------------------------------------------------------#
    def get_arrow_xy(self, start, end, strand, ypos):
        return super().get_arrow_xy(start, end, strand, ypos)
