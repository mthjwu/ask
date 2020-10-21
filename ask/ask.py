################################################################################
# APIs to AmpliconSeeK
################################################################################
import os
import pandas as pd
import covbam
import cndetect
import bpdetect
import bpjoint
import bppair
import ggraph
import grange
import trackplot
import misc
import cbs


#------------------------------------------------------------------------------#
def process_alignment(bamfile, gsfile = None, binsize = 10000,
                      mapq = 1, nmmax = 3,
                      mode = 'standard', sub_binsize = 1000,
                      seg_robost_quant = 0.5):
    """from bam file to amplified segments and breakpoint list

    Parameters
    ----------
    bamfile : input bam file with index in the same folder
        better bwa mem aligned, coordinate sorted and duplicates marked
    gsfile : genome size file (not yet in use)
    binsize : genomic bin size to count reads in
    mapq : minimal mapping quality (include)
    nmmax : maximal number of mismatches (include)
    mode : bin read count mode
        'standard' : calculate total read counts in target bins
                     used for data with low bias, such as Input
        'bias' : calculate robust statistics in target bins
                 used for data with bias coverage, such as ChIP
    sub_binsize : smaller bin size within target bin size to
        calculate robust statistics
    seg_robost_quant : quantile as robust statistics
        if None, if None auto determine for different library size
    """

    # extract all clip reads
    bp_all = covbam.clip_from_bam(bamfile, mapq = mapq, nmmax = nmmax)

    # output bedgraph
    clip_bg = covbam.clip2bedgraph(bp_all)

    # get read counts in genomic bins
    if (mode == 'standard'):
        bin_count = covbam.region_count(bamfile, gsfile,
            binsize = binsize, mapq = mapq, nmmax = nmmax)
    else:
        bin_count = covbam.region_count_2pass(bamfile, gsfile,
            binsize = binsize, mapq = mapq, nmmax = nmmax,
            sub_binsize = sub_binsize, q = seg_robost_quant)

    ## return
    return bp_all, clip_bg, bin_count


#------------------------------------------------------------------------------#
def detect_amplified_segment(bin_count, bp_all, blfile = None, genefile = None,
                             gsfile = None, cgfile = None, biasfile = None,
                             binsize = 10000,
                             min_cn = None, std_scale = 8,
                             min_segsize_scale = 5,
                             min_nt = 5):
    """detect amplified segments and breakpoint list

    Parameters
    ----------
    bin_count : bin counts
    blfile : blacklist file in bed format
        include regions having universal high read counts in any samples
    genefile : gene annotation file in bed12 format
    gsfile : genome size file (not yet in use)
    cgfile : cancer gene annotation file, gene symbol in 1st column
    binsize : genomic bin size to count reads in
    min_cn : minimal copy number to treat bin as amplification
        score based on assuming no whole genome duplication
        (None - use mean + std_scale * std as cutoff)
    std_scale : determin min_cn by mean + std_scale * std
    min_segsize_scale: minimal segment size factor
        min_segsize_scale * binsize
    min_nt : minimal nt required to be matched to bridge two breakend
    """

    # get the bias in binsize and remove bias from count data
    if (biasfile is not None):
        dfm = cndetect.get_bin_bias(bin_count, biasfile, binsize)
        bin_norm = cndetect.norm_glm(dfm)
    else:
        bin_norm = bin_count

    # set blacklist bins to diploid CN
    if (blfile is not None):
        bin_norm = cndetect.bin_filter_by_blacklist(
            bin_norm, blfile, binsize = binsize)

    # smoothing on CN column
    bin_norm = cndetect.smooth_cn(bin_norm)

    # segmentation
    cn_seg = cbs.cbs(bin_norm, binsize = binsize)

    # get amplicon
    cn_amp_raw = cbs.cbs_amplicon(cn_seg, min_cn = min_cn)

    # filter out small spikes (fold of binsize)
    cn_amp_raw = cndetect.segment_filter_by_size(cn_amp_raw, \
        binsize = binsize, fold = min_segsize_scale)

    # cn_seg boundary refine
    cn_amp = cndetect.finemap_cn_segment_boundary(cn_amp_raw, bp_all, \
        binsize = binsize)

    # annotation and output
    if (genefile is not None):
        cn_amp = cndetect.segment_annotation(cn_amp, genefile)
        if (cgfile is not None):
            cn_amp = misc.map_cancergene(cn_amp, cgfile)

    ## return
    return cn_amp, cn_seg, bin_norm


#------------------------------------------------------------------------------#
def detect_bp_pair(bamfile, bp_all, cn_amp,
                   blfile = None, bp_min_clip = 10,
                   binsize = 10000, max_n_bp = 2000, clean_bp_perc = 0.2,
                   only_keep_clean_bp = True, min_nt = 5,
                   mapq = 1, nmmax = 3):
    """ from breakpoint list to breakpoint pairs

    Parameters
    ----------
    bamfile : input bam file with index in the same folder
        better bwa mem aligned, coordinate sorted and duplicates marked
    bp_all : all breakpoints in a dict
        dict{'left' = pd_df, 'right' = pd_df}
    cn_amp : amplified segments in pd_df
    blfile : blacklist file in bed format
        include regions having universal high read counts in any samples
    bp_min_clip : minimal clipped reads of breakpoint to include
    binsize : genomic bin size
    max_n_bp : max # of breakpoints to included in the analysis
        if the sequence depth is too high, the clip points will be too many
        in this case, only the top N breakpoints will be analyzed
        if tie @ N all breakpoints <= counts @ N will be removed
    clean_bp_perc : percentage cutoff of the "CleanBP"
        CleanBP (breakpoint) is defined as the ClipDepth/InDepth >= perc
        - calculate the sequence depth on the breakpoint (InDepth)
        - # of clip reads (ClipDepth)
    only_keep_clean_bp: remove non CleanBP
    min_nt : minimal nt required to be matched to bridge two breakend
    mapq : minimal mapping quality (include)
    nmmax : maximal number of mismatches (include)

    """
    # get the breakpoint candidates >= "n" read counts
    bp_cand = bpdetect.bp_filter_by_n(bp_all, n = bp_min_clip)

    # filter breakpoints in blacklist regions
    if (blfile is not None):
        bp_cand = bpdetect.bp_filter_by_blacklist(bp_cand, blfile)

    # only include the breakpoints in amplified segments
    bp_cand = bpdetect.bp_filter_by_amplicon(bp_cand, cn_amp, binsize)

    # in this case, only the top N breakpoints will be analyzed
    bp_cand = bpdetect.bp_top_n(bp_cand, topn = max_n_bp)

    # put bp candidates into one table
    bp_cand_all = bpdetect.bp_seq_depth(bp_cand, bamfile, \
        perc = clean_bp_perc, mapq = mapq, nmmax = nmmax)

    # remove noise
    if (only_keep_clean_bp):
        bp_cand_stats = bp_cand_all[bp_cand_all['CleanBP'] == True]
    else:
        bp_cand_stats = bp_cand_all

    # get bp pairs from bp_cand_stats
    bp_duo = bpjoint.get_bppair(bamfile, bp_cand_stats, min_nt = min_nt)

    # get bp pairs from paired reads
    bp_cn = bpjoint.bp_cn_boundary(cn_amp)
    bp_duo_pe = bpjoint.get_bppair_peread(bamfile, bp_cn)

    # merge bp pairs
    bp_duo = pd.concat([bp_duo_pe, bp_duo])

    # return
    return bp_duo, bp_cand_stats, bp_cand_all


#------------------------------------------------------------------------------#
def output_bppair_alignment(bp_duo, bamfile, output_align):
    """output breakpoint pair alignments
    """
    if (bp_duo.empty == False):
        bpjoint.ouput_alignment(bp_duo, bamfile, output_align)


#------------------------------------------------------------------------------#
def construct_amplicon(bp_duo, bp_cand_stats, cn_amp,
                       genefile = None, cgfile = None,
                       min_junc_cnt = 10, bpp_min_dist = 200,
                       segment_restrict = True):
    """construct amplicon structure

    Parameters
    ----------
    bp_duo : final breakpoint pairs used to construct amplicon
    bp_cand_stats: breakpoint candidate stats table
    cn_amp : copy number amplified segments
    genefile : gene annotation file
    cgfile : cancer gene file with first column gene symbol
    min_junc_cnt : minimal supporting read counts of breakpoint pairs
    bpp_min_dist : remove breakpoint pairs with distance <= N
    segment_restrict : True - remove non CleanBP
        False - use non CleanBP as alternative segment stop
    """

    #--------------------------------------------------------------------------#
    # return empty construct if no bp pairs detected
    if (bp_duo.empty):
        colnames = ['Chrom', 'Start', 'End', 'Strand',
                    'SplitCount', 'CN', 'AmpliconID', 'Gene']
        circ_anno = pd.DataFrame(columns = colnames)
        line_anno = pd.DataFrame(columns = colnames)
        return circ_anno, line_anno, bp_duo

    #--------------------------------------------------------------------------#
    # remove breakpoint pairs in a short distance or have few split reads
    bp_pair = bppair.bp_pair_filter(bp_duo, bpp_min_dist, min_junc_cnt)

    # breakpoint stats table
    bp_fine = bppair.bp_refine(bp_pair, bp_cand_stats, cn_amp)

    # infer segments
    seg = bppair.get_segment(bp_fine, bp_pair, restrict = segment_restrict)

    # add cn to segments
    seg = bppair.add_cn(seg, cn_amp)

    #--------------------------------------------------------------------------#
    # create a new ggraph
    gg = ggraph.GGraph()

    # build ggraph from breakpoints
    gg.build_ggraph_from_bp(bp_pair, bp_fine, seg)

    # get all circles in the ggraph
    all_circ_path = gg.get_all_path_contain_circle()

    #--------------------------------------------------------------------------#
    # remove redundant
    circ_uniq = gg.path_unique(all_circ_path)

    # choose a representitive circle from each circle cluster
    circ = gg.get_representitive_path(circ_uniq)

    # make interpretable circular amplicon
    circ_df = gg.make_amplicon_df(circ)

    # annotate circular amplicon
    if (genefile is not None):
        circ_anno = circ_df.assign(
            Gene = grange.GRange.map_gene(circ_df, genefile)['Gene'])
        if (cgfile is not None):
            circ_anno = misc.map_cancergene(circ_anno, cgfile)
    else:
        circ_anno = circ_df

    #--------------------------------------------------------------------------#
    # remove redundant
    line_uniq = gg.path_unique(all_circ_path, type = 'linear')

    # choose a representitive one
    line = gg.get_representitive_path(line_uniq, type = 'linear')

    # make interpretable amplicon
    line_df = gg.make_amplicon_df(line, type = 'linear')

    # annotate circular amplicon
    if (genefile is not None):
        line_anno = line_df.assign(
            Gene = grange.GRange.map_gene(line_df, genefile)['Gene'])
        if (cgfile is not None):
            line_anno = misc.map_cancergene(line_anno, cgfile)
    else:
        line_anno = line_df

    #--------------------------------------------------------------------------#
    return circ_anno, line_anno, bp_pair



#------------------------------------------------------------------------------#
def plot_amplicon(circ_anno, line_anno = None, cn_amp = None,
                  genefile = None, bincnt = None, binsize = 10000,
                  ext = 0.3, fig_dir = 'plot', plot_n = 5,
                  fig_width = 15, fontsize = 12):
    """plot amplicon structure

    Parameters
    ----------
    circ_anno : circular amplicon dataframe
    line_anno : linear amplicon dataframe
    cn_amp : amplified segments dataframe
    genefile : gene annotation file
    bincnt : bin read counts per binsize genomic region
    ext : extent roi by 30 percent
    plot_dir : plot folder
    plot_n : only plot n structures for linear or cn
        Note : already plot all structures for circular
    """

    #--------------------------------------------------------------------------#
    # readin bed12 gene annotation file
    try:
        genebed = trackplot.read_bedf(genefile)
    except:
        print('cannot load')
    genebed = [g for g in genebed if (g.name != '')]

    #--------------------------------------------------------------------------#
    # plot amplified segments
    if (cn_amp is not None and not cn_amp.empty):
        # sort cn_amp by copy number gain
        cn_amp_list = trackplot.get_rois_with_score(cn_amp, binsize * 10)
        cn_amp_list = sorted(cn_amp_list, key=lambda x: x[3], reverse = True)

        # plot
        for tag, roi in enumerate(cn_amp_list[:plot_n]):
            # get sub dataframe
            df_sub = trackplot.get_sub_df(cn_amp, roi)

            # output pdf file name
            fig_out = os.path.join(fig_dir, 'ampseg_' + str(tag) + '.pdf')

            # prepare plot
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                df_sub, genebed, bincnt, binsize, ext, 'cn')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)

    #--------------------------------------------------------------------------#
    # plot circular amplicon
    if (not circ_anno.empty):
        for tag, circ_sub in circ_anno.groupby(['AmpliconID']):
            # output pdf file name
            fig_out = os.path.join(fig_dir, 'circular_' + tag + '.pdf')

            # prepare plot
            circ_sub_copy = circ_sub.copy()
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                circ_sub_copy, genebed, bincnt, binsize, ext, 'circular')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)

    #--------------------------------------------------------------------------#
    # plot linear amplicon
    if (line_anno is not None and not line_anno.empty):
        id2plot = list(
            line_anno.groupby('AmpliconID').CN.max().\
                sort_values(ascending = False)[0:plot_n].index)
        for tag in id2plot:
            # output pdf file name
            fig_out = os.path.join(fig_dir, 'Linear_' + tag + '.pdf')

            # get the sub dataframe
            line_sub = line_anno[line_anno['AmpliconID'] == tag]

            # prepare plot
            line_sub_copy = line_sub.copy()
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                line_sub_copy, genebed, bincnt, binsize, ext, 'linear')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)
