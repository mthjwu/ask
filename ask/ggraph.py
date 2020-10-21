################################################################################
# GGraph class: search for circular amplicons
################################################################################
#------------------------------------------------------------------------------#
from collections import Counter, defaultdict
import numpy as np
import pandas as pd


#------------------------------------------------------------------------------#
import misc



#------------------------------------------------------------------------------#
# GGraph class: search for circular amplicons
#------------------------------------------------------------------------------#
class GGraph:

    # all nodes must be 0, 1, 2, 3, 4,...
    #--------------------------------------------------------------------------#
    def __init__(self):
        self.graph = defaultdict(list)
        self.vertex = defaultdict(list)
        self.edges = defaultdict(list)
        self.node_dict = {}
        self.single_segment_loop = []

    #--------------------------------------------------------------------------#
    def add_edge(self, a, b, weight):
        self.graph[a].append(b)
        self.graph[b].append(a)
        self.edges[(a, b)] = weight
        self.edges[(b, a)] = weight

    #--------------------------------------------------------------------------#
    def add_vertex(self, vertex_key, vertex_value):
        self.vertex[vertex_key] = vertex_value

    #--------------------------------------------------------------------------#
    def get_nodes(self):
        return list(self.graph.keys())

    #--------------------------------------------------------------------------#
    def dfs_search_circle(self, v, visited, parent, path, allpath):
        """
        depth-first search for circles in a graph

        Note: make visited and path have the same nodes
        """

        # mark current node as visited
        visited[v] = True

        # store current node to list
        path.append(v)

        # loop for adjacent nodes of the current node
        for v_adj in self.graph[v]:
            # edge type is different from the last one
            # one bp_pair edge and one segment edge
            if (parent == -1 \
                or self.edges[(parent, v)][0] != self.edges[(v, v_adj)][0]):
                if (visited[v_adj] == False): # recur if the node is not visited
                    self.dfs_search_circle(v_adj, visited, v, path, allpath)
                # cycle detected if the adjacent node is visited
                # and is not the parent of current node
                elif (parent != -1 and parent != v_adj): #and v_adj in path
                    path.append(v_adj) # add v_adj, it's the loop start
                    # store when loop detected
                    allpath.append(tuple(path)) # must use tuple! or can't run
                    path.pop() # remove v_adj from path
            else:
                if (len(self.graph[v]) == 1): # end of a path, linear amplicon
                    # store when to the branch end
                    allpath.append(tuple(path)) # must use tuple! or can't run

        # remove current edge from path and mark it as unvisited
        path.pop()
        visited[v] = False

        return allpath


    #--------------------------------------------------------------------------#
    def get_all_path_contain_circle(self):
        """
        get all paths containing circle
        also include paths not containing any circle
        """

        # init onjects
        path = []
        allpath = []
        n_node = len(self.get_nodes())

        # mark all nodes as not vsisited
        visited = [False]*(n_node)

        # detect circle in different subgraph
        # say until all nodes are visisted
        for i in range(n_node):
            if visited[i] == False:
                self.dfs_search_circle(i, visited, -1, path, allpath)

        return allpath


    #--------------------------------------------------------------------------#
    # amplicon processing functions
    #--------------------------------------------------------------------------#
    def path_unique(self, circ, type = 'circular'):
        """
        get unique paths of given list of paths

        type: circular or linear
        """
        # get all circle and unique
        if (type == 'circular'):
            circ = misc.unique([row[row.index(row[-1]):-1] for row in circ])
            circ = [x for x in circ if x] # remove empty path
        elif (type == 'linear'):
            circ = [row for row in circ if (row[-1] not in row[:-1])]

        # keep ones both start and end with inner
        circ = [row for row in circ \
            if (self.edges[row[:2]][0] == 'inner' and \
                self.edges[row[-2:]][0] == 'inner')]

        if (type == 'circular'):
            # rotate path to smallest in the first
            circ = [self.rotate_path(i) for i in circ]

            # invert path to make second one smaller than last one
            circ = misc.unique([
                self.invert_path(i) if (i[1] > i[-1]) else i for i in circ])
        elif (type == 'linear'):
            circ = misc.unique([
                i[::-1] if (i[0] > i[-1]) else i for i in circ])

        return circ

    #--------------------------------------------------------------------------#
    def path_longest(self, circ, type = 'circular'):
        """
        get the longest path

        if multiple longest path,
        choose the one with max edge weight

        if multiple max edge weight,
        choose the first one
        """
        # determine whether include the last to first loop
        if (type == 'circular'):
            last_idx = 1
        elif (type == 'linear'):
            last_idx = 0

        # get all longest
        cc = [len(row) for row in circ]
        idx = misc.all_max_index(cc)
        longest = [circ[i] for i in idx]

        # get the first one with max edge weight
        oop = []
        for row in longest:
            op = []
            for v in zip(row, row[1:] + row[:last_idx]):
                if (self.edges[v][0] == 'outer'):
                    op.append(self.edges[v][1])
            oop.append(np.sum(op))
        return longest[np.argmax(oop)]

    #--------------------------------------------------------------------------#
    def get_representitive_path(self, circ, type = 'circular'):
        """
        get representitive non-overlapping paths
        """
        paths = circ.copy()
        circ_repr = []

        while(paths):
            # search for longest path
            longest = self.path_longest(paths, type = type)
            circ_repr.append(longest)
            # search non-intersected path
            tf = [not bool(set(longest).intersection(row)) for row in paths]
            paths = [i for (i, v) in zip(paths, tf) if v]

        return(circ_repr)

    #--------------------------------------------------------------------------#
    def make_amplicon_df(self, circ, type = 'circular'):
        """
        convert ggraph to interpretable amplicon dataframe
        """
        # determine whether include the last to first loop
        if (type == 'circular'):
            last_idx = 1
            tag = 'circ_'
        elif (type == 'linear'):
            last_idx = 0
            tag = 'line_'

        op = []
        for idx in range(len(circ)):
            row = circ[idx]
            op_seg = []
            op_bpp = []
            for v in zip(row, row[1:] + row[:last_idx]):
                if (self.edges[v][0] == 'inner'):
                    op_seg.append(self.vertex[v[0]][0:3] \
                        + self.vertex[v[1]][0:3] + self.edges[v])
                else:
                    op_bpp.append(self.edges[v])

            # put zero count for the last one of the linear amplicon
            if (type == 'linear'):
                op_bpp.append((('outer', 0)))

            for i in range(len(op_seg)):
                row_seg = op_seg[i]
                row_bpp = op_bpp[i]
                
                if (row_seg[2] == 'L'):
                    op.append([row_seg[0], row_seg[1], row_seg[4], \
                        '+', row_bpp[1], row_seg[7], tag + str(idx)])
                else:
                    op.append([row_seg[0], row_seg[4], row_seg[1], \
                        '-', row_bpp[1], row_seg[7], tag + str(idx)])

        colnames = ['Chrom', 'Start', 'End', 'Strand', \
            'SplitCount', 'CN', 'AmpliconID']
        df = pd.DataFrame(op, columns = colnames)

        # also add single segment loops
        if (type == 'circular'):
            node_in_path = [v for path in circ for v in path]
            ssl = [row for row in self.single_segment_loop \
                if (self.node_dict[(row[0], row[1], 'L')] not in node_in_path)]

            # add index
            ssl_df = []
            idx = len(circ)
            for row in ssl:
                ssl_df.append(row + ['circ_' + str(idx)])
                idx += 1

            # make dataframe
            ssl_df = pd.DataFrame(ssl_df, columns = colnames)
            df = pd.concat([df, ssl_df])

        return df


    #--------------------------------------------------------------------------#
    def build_ggraph_from_bp(self, bp_pair, bp_fine, seg):
        """
        build ggraph from breakpoint data
        """

        for row in bp_fine.itertuples():
            self.add_vertex(row[0], row[1:])

        for row in bp_fine.itertuples():
            self.node_dict[row[1:4]] = row[0]

        for row in bp_pair.itertuples():
            a = self.node_dict[row[1:4]]
            b = self.node_dict[row[4:7]]
            w = ('outer', row[7])
            self.add_edge(a, b, w)

        for row in seg.itertuples():
            a = self.node_dict[(row[1], row[2], 'L')]
            b = self.node_dict[(row[1], row[3], 'R')]
            w = ('inner', row[4])
            # if already exist, it's a single segment loop
            if (a, b) in self.edges:
                self.single_segment_loop.append([row[1], row[2], row[3], \
                    '+', self.edges[(a, b)][1], row[4]])
                # don't add inner in this case
            else:
                self.add_edge(a, b, w)


    #--------------------------------------------------------------------------#
    @staticmethod
    def invert_path(path):
        path = path[::-1]
        return path[-1:] + path[:-1]

    @staticmethod
    def rotate_path(path):
        i = path.index(min(path))
        return path[i:]+path[:i]

    @staticmethod
    def is_new_path(path, paths):
        return not path in paths



    #--------------------------------------------------------------------------#
    # currently not in use
    #--------------------------------------------------------------------------#
    def dfs_connected(self, v, visited, visited_list):
        """
        depth-first search for connected nodes (not in use currently)
        """

        # mark current node as visited
        visited[v] = True

        # store current node to list
        visited_list.append(v)

        # recur for adjacent nodes of the current node
        for v_adj in self.graph[v]:
            if (visited[v_adj] == False): # add to list if not visited
                visited_list = self.dfs_connected(v_adj, visited, visited_list)

        return visited_list

    #--------------------------------------------------------------------------#
    def get_connected_subgraph(self):
        """
        get connected subgraph nodes (not in use currently)
        """

        # init objects
        n_node = len(self.get_nodes())
        visited = []
        op = []

        # mark all nodes as not vsisited
        visited = [False]*(n_node)

        # loop to search for all subgraphs
        for v in range(n_node):
            if (visited[v] == False):
                visited_list = []
                op.append(self.dfs_connected(v, visited, visited_list))
        return op
