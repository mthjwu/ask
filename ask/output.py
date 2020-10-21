################################################################################
# output function including plotting
################################################################################
#------------------------------------------------------------------------------#
import os
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------#
def write_file_for_amplot(df_circ, bp_pair, dtype = 'circular',
                          plot_dir = 'plot/'):
    """
    write files for amplot
    """
    circ_bed, circ_arc, circ_roi = prepare_data_for_amplot(df_circ, bp_pair)
    
    file_arc = plot_dir + dtype + '_amplicons.arcs'
    file_bed = plot_dir + dtype + '_amplicons.bed'
    file_roi = plot_dir + dtype + '_amplicons_roi.bed'
        
    circ_bed.to_csv(file_bed, sep='\t', index=False, header=False)
    circ_arc.to_csv(file_arc, sep='\t', index=False, header=False)
    circ_roi.to_csv(file_roi, sep='\t', index=False, header=False)

    return circ_roi

#------------------------------------------------------------------------------#
def prepare_data_for_amplot(df_circ, bp_pair):
    """
    prepare segment bed file and arc file for plot

    df: dataframe of amplicon strutures
        circ_df or line_df
    """
    seg = df_circ.copy()
    link = bp_pair.copy()

    # make segment bed file
    seg['Name'] = list(seg.index)
    seg = seg[['Chrom', 'Start', 'End', 'Name', 'CN', 'Strand']].copy()

    # remove links not in amplicon structures
    keep = []
    for row in seg.itertuples():
        keep.append((row[1], row[2], 'L'))
        keep.append((row[1], row[3], 'R'))
    tf = [(row[1:4] in keep) and (row[4:7] in keep)
        for row in link.itertuples()]

    # make link file
    link = link[tf]
    if (not link.empty):
        link.insert(2, 'End1', link['Coord1'])
        link.insert(4, 'End2', link['Coord2'])
        link = link[['Chrom1', 'Coord1', 'End1', 'Chrom2',
                    'Coord2', 'End2', 'Count']].copy()
        # change to 0-based
        link['Coord1'] -= 1
        link['Coord2'] -= 1 
        link = link.drop_duplicates()

    # make the regions of interest
    op = []
    for gp1 in df_circ.groupby(['AmpliconID']):
        for gp in gp1[1].groupby(['Chrom']):
            op.append([gp[0], np.min(gp[1].iloc[:,1]),
                       np.max(gp[1].iloc[:,2])])
    roi = pd.DataFrame(op)

    # change to 0-based
    seg['Start'] -= 1
    roi.iloc[:,1] -= 1

    # sort bed file
    seg = seg.drop_duplicates().sort_values(['Chrom', 'Start', 'End'])
    roi = roi.drop_duplicates()

    return seg, link, roi


#------------------------------------------------------------------------------#
def make_amplot(rois, bed12 = None, trackfiles = [], dtype = 'circular',
                plot_dir = 'plot/'):
    """
    make amplicon plot
    """
    ## make plot
    for row in rois.itertuples():
        ext = 0.1 * (row[3] - row[2])
        roi = row[1] + ':' + str(int(row[2]-ext)) + '-' + str(int(row[3]+ext))
        file_arc = plot_dir + dtype + '_amplicons.arcs'
        file_bed = plot_dir + dtype + '_amplicons.bed'
        
        input_all = ' '.join(filter(None, trackfiles)) + ' '
        input_all += file_arc + ' '
        input_all += file_bed + ' '
        if (bed12 is not None): input_all += bed12 + ' '
        output_ini_default = plot_dir + '.' + dtype + '_tracks_' \
                             + str(int(row[0]+1)) + '.ini'
        output_ini = plot_dir + dtype + '_tracks_' + str(int(row[0]+1)) + '.ini'
        output_pdf = plot_dir + dtype + '_plot_' + str(int(row[0]+1)) + '.pdf'

        cmd1 = 'make_tracks_file --trackFiles ' + input_all + '-o ' \
                + output_ini_default
        # print(cmd1)

        cmd2 = "sed 's/labels = false/labels = True/g' " \
                + output_ini_default + " > " + output_ini
        # print(cmd2)

        cmd3 = 'pyGenomeTracks --tracks ' + output_ini + ' --region ' \
                + roi + ' --outFileName ' + output_pdf
        # print(cmd3)

        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)

#------------------------------------------------------------------------------#
def make_amplot_general(rois, trackfiles = [], dtype = 'cn',
                        plot_dir = 'plot/'):
    """
    make amplicon plot with as many input track files (need conda install tqdm)
    """
    ## make plot
    for row in rois.itertuples():
        ext = 0.05 * (row[3] - row[2])
        roi = row[1] + ':' + str(int(row[2]-ext)) + '-' + str(int(row[3]+ext))

        input_all = ' '.join(filter(None, trackfiles)) + ' '

        output_ini_default = plot_dir + '.' + dtype + '_tracks_' \
                             + str(int(row[0]+1)) + '.ini'
        output_ini = plot_dir + dtype + '_tracks_' \
                     + str(int(row[0]+1)) + '.ini'
        output_pdf = plot_dir + dtype + '_plot_' + str(int(row[0]+1)) + '.pdf'

        cmd1 = 'make_tracks_file --trackFiles ' + input_all \
                + '-o ' + output_ini_default
        # print(cmd1)

        cmd2 = "sed 's/labels = false/labels = True/g' " \
                + output_ini_default + " > " + output_ini
        # print(cmd2)

        cmd3 = 'pyGenomeTracks --tracks ' + output_ini + ' --region ' \
                + roi + ' --outFileName ' + output_pdf
        # print(cmd3)

        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
