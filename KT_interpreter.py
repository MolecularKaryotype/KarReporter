import copy

from Report_Genes import *
import re


class Aligned_Haplotype:
    id: int  # for debugging
    chrom: str
    alignment_score: int
    mt_aligned: [str]
    wt_aligned: [str]
    block_indices: {int: (int, int)}  # [block_id: (a, b)] where event block are [a, b)
    mt_blocks: {int: (str,)}  # insertion discordant blocks
    wt_blocks: {int: (str,)}  # deletion discordant blocks
    concordant_blocks: {int: (str,)}
    discordant_block_assignment: {int: str}  # block index: event name, exactly one event can be assigned to one block; only for discordant blocks
    size_dict: {str: int}  # seg_name: bp size
    mt_hap = []
    wt_hap = []

    def __init__(self, mt_aligned, wt_aligned, alignment_score, size_dict, id, chrom, mt_hap, wt_hap):
        self.id = id
        self.chrom = chrom
        self.mt_aligned = mt_aligned
        self.wt_aligned = wt_aligned
        self.alignment_score = alignment_score
        self.size_dict = size_dict
        self.mt_hap = mt_hap
        self.wt_hap = wt_hap

        ## get event_blocks
        event_block = []
        block_types = []  # [str], a block is either ins/del
        # we know len(wt_hap) == len(mt_hap) as they are global alignments
        hap_len = len(self.mt_aligned)
        seg_idx = 0
        while seg_idx < hap_len:
            mt_seg = self.mt_aligned[seg_idx]
            wt_seg = self.wt_aligned[seg_idx]
            if mt_seg == wt_seg:
                seg_idx += 1
                continue
            else:
                # this will be the first seg of (potentially) a section, that is indel
                if wt_seg == '-':
                    # insertion
                    # block ends whenever continuous section OR insertion section ends
                    continuous_end = continuous_extension(self.mt_aligned, seg_idx)
                    ins_end = seg_idx + 1
                    while ins_end < hap_len:
                        if self.wt_aligned[ins_end] == '-':
                            ins_end += 1
                        else:
                            break
                    block_end = min(continuous_end, ins_end)
                    block_types.append('ins')
                elif mt_seg == '-':
                    # deletion
                    continuous_end = continuous_extension(self.wt_aligned, seg_idx)
                    ins_end = seg_idx + 1
                    while ins_end < hap_len:
                        if self.mt_aligned[ins_end] == '-':
                            ins_end += 1
                        else:
                            break
                    block_end = min(continuous_end, ins_end)
                    block_types.append('del')
                else:
                    raise ValueError('mismatch not allowed')
                event_block.append((seg_idx, block_end))
                seg_idx = block_end

        ## extract concordant blocks, mt/wt blocks: first form blocks using INT-idx, then convert to STR
        discordant_blocks = {event_block_boundaries: block_types[event_block_idx] for event_block_idx, event_block_boundaries in enumerate(event_block)}
        # all gaps in discordant blocks are concordant blocks
        self.block_indices = {}
        self.wt_blocks = {}
        self.mt_blocks = {}
        self.concordant_blocks = {}
        self.discordant_block_assignment = {}
        block_id = 0
        p_endpoint = 0
        for discordant_block in discordant_blocks:
            c_startpoint = int(discordant_block[0])
            if p_endpoint < c_startpoint:
                # fill gap with concordant blocks
                # between two discordant blocks, the region belongs to at most 1 concordant block because the wt is always continuous
                # (i.e. will not split into 2 blocks)
                self.block_indices[block_id] = (p_endpoint, c_startpoint)
                self.concordant_blocks[block_id] = wt_aligned[p_endpoint: c_startpoint]
                block_id += 1
            self.block_indices[block_id] = discordant_block
            self.discordant_block_assignment[block_id] = ''
            if discordant_blocks[discordant_block] == 'del':
                self.wt_blocks[block_id] = wt_aligned[discordant_block[0]: discordant_block[1]]
            elif discordant_blocks[discordant_block] == 'ins':
                self.mt_blocks[block_id] = mt_aligned[discordant_block[0]: discordant_block[1]]
            block_id += 1
            p_endpoint = discordant_block[1]
        # add the last concordant block, if present
        if p_endpoint < len(self.wt_aligned):
            self.block_indices[block_id] = (p_endpoint, len(self.wt_aligned))
            self.concordant_blocks[block_id] = wt_aligned[p_endpoint: len(self.wt_aligned)]

        ## convert all block indices to STR
        def change_dict_key_type(input_dict):
            new_dict = {}
            for key, value in input_dict.items():
                new_dict[str(key)] = value
            return new_dict

        self.block_indices = change_dict_key_type(self.block_indices)
        self.discordant_block_assignment = change_dict_key_type(self.discordant_block_assignment)
        self.concordant_blocks = change_dict_key_type(self.concordant_blocks)
        self.mt_blocks = change_dict_key_type(self.mt_blocks)
        self.wt_blocks = change_dict_key_type(self.wt_blocks)

    def __str__(self):
        return "ID<{}>".format(self.id)

    def unique_segment_indices(self):
        return_set = set()
        for seg in self.wt_hap + self.mt_hap:
            seg_index = seg[:-1]
            return_set.add(seg_index)
        return return_set

    def get_block_segs_and_block_type(self, block_name):
        """
        given block_idx, get block segs and block type
        :param block_name:
        :return:
        """
        if block_name in self.wt_blocks:
            return self.wt_blocks[block_name], self.wt_blocks
        elif block_name in self.mt_blocks:
            return self.mt_blocks[block_name], self.mt_blocks
        elif block_name in self.concordant_blocks:
            return self.concordant_blocks[block_name], self.concordant_blocks
        else:
            return None, None

    def split_block(self, current_block_name, split_left_idx, split_right_idx):
        """
        split a block into [start, split_left_idx - 1], [split_left_idx, split_right_idx - 1], [split_right_idx, end]
        :param current_block_name:
        :param split_left_idx:
        :param split_right_idx:
        :return:
        """
        ## find the block type
        block_segs, blocks = self.get_block_segs_and_block_type(current_block_name)
        if split_left_idx == 0 and split_right_idx == len(block_segs):
            # no split required
            return current_block_name

        ## remove this block from all attributes
        if current_block_name in self.discordant_block_assignment:
            event_assignment = self.discordant_block_assignment[current_block_name]
            if len(event_assignment) > 0 and event_assignment != 'under investigation':
                raise RuntimeError()
        blocks.pop(current_block_name)
        block_indices = self.block_indices.pop(current_block_name)
        block_in_discordant_blocks = False
        if current_block_name in self.discordant_block_assignment:
            # will not have assigned value
            self.discordant_block_assignment.pop(current_block_name)
            block_in_discordant_blocks = True

        ## splitting
        middle_block_idx = 0
        sublist_segs = []
        sublist_seg_indicies = []
        seg_indices_shift = block_indices[0]
        seg_indices_end = block_indices[1]
        if split_left_idx != 0:
            # left-most
            sublist_segs.append(block_segs[:split_left_idx])
            sublist_seg_indicies.append((seg_indices_shift, seg_indices_shift + split_left_idx))
            middle_block_idx = 1
        # middle
        sublist_segs.append(block_segs[split_left_idx:split_right_idx])
        sublist_seg_indicies.append((seg_indices_shift + split_left_idx, seg_indices_shift + split_right_idx))
        if split_right_idx != len(block_segs):
            # right-most
            sublist_segs.append(block_segs[split_right_idx:])
            sublist_seg_indicies.append((seg_indices_shift + split_right_idx, seg_indices_end))

        ## re-assigning
        middle_block_name = ''
        for new_block_idx, (new_segs, new_seg_indices) in enumerate(zip(sublist_segs, sublist_seg_indicies)):
            new_block_name = current_block_name + chr(new_block_idx + 97)
            if new_block_idx == middle_block_idx:
                middle_block_name = new_block_name
            blocks[new_block_name] = new_segs
            self.block_indices[new_block_name] = new_seg_indices
            if block_in_discordant_blocks:
                self.discordant_block_assignment[new_block_name] = ''

        return middle_block_name

    def next_block(self, current_block_name):
        ### OPTIMIZE: add a variable in self to store the ordered list of block names, so no for loop needed
        current_block_value = block_value(current_block_name)
        min_value = float('inf')
        next_block = None
        for block in self.block_indices:
            value = block_value(block)
            if current_block_value < value < min_value:
                next_block = block
                min_value = value
        return next_block

    def previous_block(self, current_block_name):
        current_block_value = block_value(current_block_name)
        max_value = -1.0
        previous_block = None
        for block in self.block_indices:
            value = block_value(block)
            if max_value < value < current_block_value:
                previous_block = block
                max_value = value
        return previous_block

    def next_different_type_block(self, block_name):
        """
        :param block_name:
        :param different_type: True->search for next block of different type; False->search for next block of same type
        :return:
        """
        current_block_name = block_name
        _, original_block_dict = self.get_block_segs_and_block_type(current_block_name)
        skipped_distance = 0
        while True:
            current_block_name = self.next_block(current_block_name)
            if current_block_name is None:
                return None, None
            _, current_block_dict = self.get_block_segs_and_block_type(current_block_name)

            if current_block_dict != original_block_dict:
                return current_block_name, skipped_distance
            else:
                if self.discordant_block_assignment[current_block_name] in ['deletion', 'insertion']:
                    current_block_segs, _ = self.get_block_segs_and_block_type(current_block_name)
                    skipped_distance += section_size(current_block_segs, self.size_dict)

    def next_block_if_same_type(self, block_name):
        """
        :param block_name:
        :return: None if the next block is not of the same type, else, return next block's name
        """
        _, original_block_dict = self.get_block_segs_and_block_type(block_name)
        next_block_name = self.next_block(block_name)
        if next_block_name is None:
            return None
        else:
            _, next_block_dict = self.get_block_segs_and_block_type(next_block_name)
            if next_block_dict == original_block_dict:
                return next_block_name
            else:
                return None

    def previous_different_type_block(self, block_name):
        """
        :param block_name:
        :param different_type: True->search for next block of different type; False->search for next block of same type
        :return:
        """
        current_block_name = block_name
        _, original_block_dict = self.get_block_segs_and_block_type(current_block_name)
        skipped_distance = 0
        while True:
            current_block_name = self.previous_block(current_block_name)
            if current_block_name is None:
                return None, None
            _, current_block_dict = self.get_block_segs_and_block_type(current_block_name)

            if current_block_dict != original_block_dict:
                return current_block_name, skipped_distance
            else:
                # NOTE: only count distance of deletion and insertion
                if self.discordant_block_assignment[current_block_name] in ['deletion', 'insertion']:
                    current_block_segs, _ = self.get_block_segs_and_block_type(current_block_name)
                    skipped_distance += section_size(current_block_segs, self.size_dict)

    def previous_block_if_same_type(self, block_name):
        _, original_block_dict = self.get_block_segs_and_block_type(block_name)
        previous_block_name = self.previous_block(block_name)
        if previous_block_name is None:
            return None
        else:
            _, previous_block_dict = self.get_block_segs_and_block_type(previous_block_name)
            if previous_block_dict == original_block_dict:
                return previous_block_name
            else:
                return None

    def closest_block(self, current_block_name, direction, block_type):
        """
        get the closest left/right block of 'this' block;
        :param current_block_name: mt/wt/concordant
        :param direction: left/right
        :param block_type: 'this' block may not be in the block_type
        :return:
        """
        if block_type == 'mt':
            blocks = self.mt_blocks
        elif block_type == 'wt':
            blocks = self.wt_blocks
        elif block_type == 'concordant':
            blocks = self.concordant_blocks
        else:
            raise ValueError()
        if direction not in ['left', 'right']:
            raise ValueError()

        if direction == 'left':
            max_block = current_block_name
            c_value = block_value(current_block_name)
            max_value = -1.0
            for block in blocks:
                if max_value < block_value(block) < c_value:
                    max_block = block
                    max_value = block_value(block)
            if max_block == current_block_name:
                return None
            else:
                return max_block
        else:
            min_block = current_block_name
            c_value = block_value(current_block_name)
            min_value = float('inf')
            for block in blocks:
                if min_value > block_value(block) > c_value:
                    min_block = block
                    min_value = block_value(block)
            if min_block == current_block_name:
                return None
            else:
                return min_block

    def get_closest_real_seg(self, block_name, search_direction):
        """
        real block only exists in mt-block and concordant-block
        :param block_name:
        :param search_direction: left/right
        :return: list of indexed_segs corresponding to the closest block
        """
        if search_direction == 'left':
            left_mt = self.closest_block(block_name, 'left', 'mt')
            left_concordant = self.closest_block(block_name, 'left', 'concordant')
            if left_mt is None and left_concordant is None:
                return 'p-ter'
            else:
                if left_mt is None:
                    left_mt_val = -1.0
                else:
                    left_mt_val = block_value(left_mt)
                if left_concordant is None:
                    left_concordant_val = -1.0
                else:
                    left_concordant_val = block_value(left_concordant)
                if left_concordant_val > left_mt_val:
                    return self.concordant_blocks[left_concordant][-1]
                elif left_concordant_val < left_mt_val:
                    return self.mt_blocks[left_mt][-1]
                else:
                    raise RuntimeError()
        elif search_direction == 'right':
            right_mt = self.closest_block(block_name, 'right', 'mt')
            right_concordant = self.closest_block(block_name, 'right', 'concordant')
            if right_mt is None and right_concordant is None:
                return 'q-ter'
            else:
                if right_mt is None:
                    right_mt_val = float('inf')
                else:
                    right_mt_val = block_value(right_mt)
                if right_concordant is None:
                    right_concordant_val = float('inf')
                else:
                    right_concordant_val = block_value(right_concordant)
                if right_concordant_val < right_mt_val:
                    return self.concordant_blocks[right_concordant][0]
                elif right_concordant_val > right_mt_val:
                    return self.mt_blocks[right_mt][0]
                else:
                    raise RuntimeError()

    def report_SV(self, event_blocks, event_types):
        """
        :param event_blocks: cumulative dict, where value is a list of event_block, this helps to group paired-events that were on different chromosomes
        :param event_types: cumulative dict, where value is a list of event_type, same reason as above
        :return:
        """

        def get_block_boundaries(block_idx):
            if block_idx == 0:
                # first block
                left = 'p-ter'
            else:
                # previous_block, _ = self.get_block_segs_and_block_type(block_idx - 1)
                left = self.get_closest_real_seg(block_idx, 'left')
                # left = previous_block[-1]
            if block_idx == len(self.block_indices) - 1:
                # last block
                right = 'q-ter'
            else:
                # next_block, _ = self.get_block_segs_and_block_type(block_idx + 1)
                right = self.get_closest_real_seg(block_idx, 'right')
                # right = next_block[0]
            return left, right

        for mt_block_idx, mt_block in self.mt_blocks.items():
            mt_block_str = ','.join(list(mt_block))
            event_id = int(self.discordant_block_assignment[mt_block_idx].split(',')[-1])
            event_type = self.discordant_block_assignment[mt_block_idx].split(',')[0]
            left_boundary_seg, right_boundary_seg = get_block_boundaries(mt_block_idx)
            if event_id in event_blocks:
                event_blocks[event_id].append('{}.{}.mt({}).{}.{}'.format(self.id, mt_block_idx, mt_block_str, left_boundary_seg, right_boundary_seg))
                event_types[event_id].append(event_type)
            else:
                event_blocks[event_id] = ['{}.{}.mt({}).{}.{}'.format(self.id, mt_block_idx, mt_block_str, left_boundary_seg, right_boundary_seg)]
                event_types[event_id] = [event_type]
        for wt_block_idx, wt_block in self.wt_blocks.items():
            wt_block_str = ','.join(list(wt_block))
            event_id = int(self.discordant_block_assignment[wt_block_idx].split(',')[-1])
            event_type = self.discordant_block_assignment[wt_block_idx].split(',')[0]
            left_boundary_seg, right_boundary_seg = get_block_boundaries(wt_block_idx)
            if event_id in event_blocks:
                event_blocks[event_id].append('{}.{}.wt({}).{}.{}'.format(self.id, wt_block_idx, wt_block_str, left_boundary_seg, right_boundary_seg))
                event_types[event_id].append(event_type)
            else:
                event_blocks[event_id] = ['{}.{}.wt({}).{}.{}'.format(self.id, wt_block_idx, wt_block_str, left_boundary_seg, right_boundary_seg)]
                event_types[event_id] = [event_type]

        # sort each event_blocks by chr, then by block_idx
        for event_id, event_block_list in event_blocks.items():
            event_blocks[event_id] = sorted(event_block_list, key=lambda x: (int(x.split('.')[0]), x.split('.')[1]), reverse=False)

        return event_blocks, event_types

    def search_translocation_seed(self, query_section, query_block_name, seed_type, allow_inversion):
        """
        :param query_block_name:
        :param query_section:
        :param seed_type: 'ins' OR 'del', which type of seed to be searched of
        :param d:
        :param eps:
        :return: splitting information if found, None if not
        """
        ## assign the complementary reference blocks
        if seed_type == 'del':
            query_block = self.wt_blocks
        elif seed_type == 'ins':
            query_block = self.mt_blocks
        else:
            raise ValueError('invalid seed type')

        ## find biggest match
        max_size = -1
        max_size_ref_block = None
        max_size_query_sublist = None
        sublist_inverted = False
        for c_block_name, c_block in query_block.items():
            if len(self.discordant_block_assignment[c_block_name]) != 0:
                # already has assignment
                # TODO: can improve this by bipartite matching, instead of enumerative matching to greedily find the max size matching, one pair at a time
                continue
            overlap_query_sublist, _, _, overlap_size, c_sublist_inverted = max_size_overlap(c_block, query_section,
                                                                                             self.size_dict, -1,
                                                                                             allow_inversion=allow_inversion)
            if overlap_query_sublist != -1:
                if overlap_size > max_size:
                    max_size = overlap_size
                    max_size_ref_block = c_block_name
                    max_size_query_sublist = overlap_query_sublist
                    sublist_inverted = c_sublist_inverted

        ## split blocks for the matching
        if max_size_ref_block is None:
            return None
        ref_block_start_idx = sublist_idx(max_size_query_sublist, query_block[max_size_ref_block])
        ref_block_end_idx = ref_block_start_idx + len(max_size_query_sublist)
        query_block_start_idx = sublist_idx(max_size_query_sublist, query_section, inverted=sublist_inverted)
        query_block_end_idx = query_block_start_idx + len(max_size_query_sublist)
        return {'ref_block': {'max_size_ref_block': max_size_ref_block,
                              'ref_block_start_idx': ref_block_start_idx,
                              'ref_block_end_idx': ref_block_end_idx},
                'query_block': {'query_block_name': query_block_name,
                                'query_block_start_idx': query_block_start_idx,
                                'query_block_end_idx': query_block_end_idx}}


def block_value(block_name):
    # assumes there are at most 100 sub-blocks (in terms of ASCII)
    pattern = r'(\d+)(\D*)'
    match = re.match(pattern, block_name)
    if len(match.group(2)) == 0:
        char_value = 0.0
    else:
        char_value = 0.0
        multiplier = 0.01
        for character in match.group(2):
            char_value += multiplier * (ord(character) - 96)
            multiplier *= 0.01
    return int(match.group(1)) + char_value


def lcs(list1, list2, size_dict):
    """
    longest common subsequence finding
    :param list1:
    :param list2:
    :param size_dict:
    :return:
    """
    indel_penalty_per_nt = -1
    alignment_1 = []
    alignment_2 = []
    scoring_matrix = [[0 for i in range(0, len(list2) + 1)] for j in range(0, len(list1) + 1)]
    backtrack_matrix = [["" for i in range(0, len(list2) + 1)] for j in range(0, len(list1) + 1)]

    # initialize starting grid
    scoring_matrix[0][0] = 0
    for row_index in range(1, len(list1) + 1):
        current_segment = list1[row_index - 1]
        len_current_segment = size_dict[current_segment[:-1]]  # remove sign
        scoring_matrix[row_index][0] = scoring_matrix[row_index - 1][0] + len_current_segment * indel_penalty_per_nt
        backtrack_matrix[row_index][0] = "down"

    for col_index in range(1, len(list2) + 1):
        current_segment = list2[col_index - 1]
        len_current_segment = size_dict[current_segment[:-1]]  # remove sign
        scoring_matrix[0][col_index] = scoring_matrix[0][col_index - 1] + len_current_segment * indel_penalty_per_nt
        backtrack_matrix[0][col_index] = "rigt"

    # DP recursion
    for row_index in range(1, len(list1) + 1):
        for col_index in range(1, len(list2) + 1):
            len_down_segment = size_dict[list1[row_index - 1][:-1]]
            down_value = scoring_matrix[row_index - 1][col_index] + len_down_segment * indel_penalty_per_nt

            len_right_segment = size_dict[list2[col_index - 1][:-1]]
            right_value = scoring_matrix[row_index][col_index - 1] + len_right_segment * indel_penalty_per_nt

            if list1[row_index - 1] != list2[col_index - 1]:
                # mismatch: not allowed
                diagonal_value = float('-inf')
            else:
                # match
                diagonal_value = scoring_matrix[row_index - 1][col_index - 1]

            if diagonal_value >= down_value and diagonal_value >= right_value:
                scoring_matrix[row_index][col_index] = diagonal_value
                backtrack_matrix[row_index][col_index] = "diag"
            elif down_value >= right_value:
                scoring_matrix[row_index][col_index] = down_value
                backtrack_matrix[row_index][col_index] = "down"
            else:
                scoring_matrix[row_index][col_index] = right_value
                backtrack_matrix[row_index][col_index] = "rigt"

    # backtracking
    final_score = scoring_matrix[len(list1)][len(list2)]
    current_row = len(list1)
    current_col = len(list2)

    while True:
        if current_row == 0 and current_col == 0:
            break
        if backtrack_matrix[current_row][current_col] == "diag":
            alignment_1.insert(0, list1[current_row - 1])
            alignment_2.insert(0, list2[current_col - 1])
            current_col -= 1
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "down":
            alignment_1.insert(0, list1[current_row - 1])
            alignment_2.insert(0, "-")
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "rigt":
            alignment_1.insert(0, "-")
            alignment_2.insert(0, list2[current_col - 1])
            current_col -= 1
        else:
            raise RuntimeError("error in backtrack matrix")

    return final_score, alignment_1, alignment_2


def interpret_haplotypes(mt_hap_list: [[str]], wt_hap_list: [[str]], chrom_identities: [str], segment_size_dict: {str: int}, d=500000, eps=500000):
    """
    Input a haplotype and a WT, report SVs; hap from mt list must correspond to hap from wt list
    :param segment_size_dict: {str: int} segment mapped to size of the segment (eg. 11: 12,345)
    :param mt_hap_list:
    :param wt_hap_list:
    :param chrom_identities:
    :return:
    """
    event_id = 0  # static variable, keep different events separate, keep paired events under same id (eg. balanced trans)
    hap_id = 0  # used for naming haplotypes, preserves the order they come-in in mt_hap/wt_hap lists
    aligned_haplotypes = []
    for idx, mt_hap in enumerate(mt_hap_list):
        wt_hap = wt_hap_list[idx]
        c_chrom = chrom_identities[idx]
        a1, a2, a3 = lcs(mt_hap, wt_hap, segment_size_dict)
        aligned_haplotypes.append(Aligned_Haplotype(a2, a3, a1, segment_size_dict, hap_id, c_chrom, mt_hap, wt_hap))
        hap_id += 1

    def genomewide_seed_search_for_translocation(query_section, query_block_name_presplit, query_hap_idx, query_type):
        """
        :param query_section: (segs) to be searched for
        :param query_block_name_presplit: needed for writing discordant_block_assignment if seed found
        :param query_hap_idx: needed to prioritize searching inter-chr first, then intra-chr if inter-chr not found
        :param query_type:
        :return:
        """
        for aligned_hap_idx1, aligned_hap1 in enumerate(aligned_haplotypes):
            if aligned_hap_idx1 == query_hap_idx:
                # first get inter-chr event
                continue
            splitting_info = aligned_hap1.search_translocation_seed(query_section, query_block_name_presplit, query_type, True)
            if splitting_info is not None:
                ref_hap = aligned_haplotypes[aligned_hap_idx1]
                query_hap = aligned_haplotypes[query_hap_idx]
                ref_block_name_postsplit = ref_hap.split_block(splitting_info['ref_block']['max_size_ref_block'],
                                                               splitting_info['ref_block']['ref_block_start_idx'],
                                                               splitting_info['ref_block']['ref_block_end_idx'])
                query_block_name_postsplit = query_hap.split_block(splitting_info['query_block']['query_block_name'],
                                                                   splitting_info['query_block']['query_block_start_idx'],
                                                                   splitting_info['query_block']['query_block_end_idx'])
                aligned_hap1.discordant_block_assignment[ref_block_name_postsplit] = 'balanced_translocation,{}'.format(event_id)
                aligned_haplotypes[query_hap_idx].discordant_block_assignment[query_block_name_postsplit] = 'balanced_translocation,{}'.format(event_id)
                seed_found = True
                if ref_block_name_postsplit != splitting_info['ref_block']['max_size_ref_block'] or \
                        query_block_name_postsplit != splitting_info['query_block']['query_block_name']:
                    block_splitted = True
                else:
                    block_splitted = False
                return seed_found, block_splitted
        # at last, check intra-chr event
        # do not allow intra-chr reciprocal events to have inversion, otherwise, normal inversions will be marked
        splitting_info = aligned_haplotypes[query_hap_idx].search_translocation_seed(query_section, query_block_name_presplit, query_type, False)
        if splitting_info is not None:
            ref_hap = aligned_haplotypes[query_hap_idx]
            query_hap = aligned_haplotypes[query_hap_idx]
            ref_block_name_postsplit = ref_hap.split_block(splitting_info['ref_block']['max_size_ref_block'],
                                                           splitting_info['ref_block']['ref_block_start_idx'],
                                                           splitting_info['ref_block']['ref_block_end_idx'])
            query_block_name_postsplit = query_hap.split_block(splitting_info['query_block']['query_block_name'],
                                                               splitting_info['query_block']['query_block_start_idx'],
                                                               splitting_info['query_block']['query_block_end_idx'])
            aligned_haplotypes[query_hap_idx].discordant_block_assignment[ref_block_name_postsplit] = 'balanced_translocation,{}'.format(event_id)
            aligned_haplotypes[query_hap_idx].discordant_block_assignment[query_block_name_postsplit] = 'balanced_translocation,{}'.format(event_id)

            seed_found = True
            if ref_block_name_postsplit != splitting_info['ref_block']['max_size_ref_block'] or \
                    query_block_name_postsplit != splitting_info['query_block']['query_block_name']:
                block_splitted = True
            else:
                block_splitted = False
            return seed_found, block_splitted

        return False, False

    ### Order of resolving: balanced-translocation, inv, dup-inv, tandem-dup, del, unbalanced-translocation
    ## resolve all balanced translocations
    # re-run all haps when splitting occurs, as more balanced translocation can be found
    run = True
    while run:
        splitted = False
        for aligned_hap_idx, aligned_hap in enumerate(aligned_haplotypes):
            for c_wt_block_name, c_wt_block in aligned_hap.wt_blocks.items():
                if len(aligned_hap.discordant_block_assignment[c_wt_block_name]) > 0:
                    # already has assignment
                    continue
                # block it from being matched in intra-chr event; maybe not necessary
                aligned_hap.discordant_block_assignment[c_wt_block_name] = 'under investigation'
                balanced_translocation_found, c_splitted = genomewide_seed_search_for_translocation(c_wt_block, c_wt_block_name, aligned_hap_idx, 'ins')
                if balanced_translocation_found:
                    event_id += 1
                    if c_splitted:
                        splitted = True
                    break
                else:
                    # reset block assignment for query
                    aligned_hap.discordant_block_assignment[c_wt_block_name] = ''
            for c_mt_block_name, c_mt_block in aligned_hap.mt_blocks.items():
                if len(aligned_hap.discordant_block_assignment[c_mt_block_name]) > 0:
                    # already has assignment
                    continue
                aligned_hap.discordant_block_assignment[c_mt_block_name] = 'under investigation'  # block it from being matched in intra-chr event
                balanced_translocation_found, c_splitted = genomewide_seed_search_for_translocation(c_mt_block, c_mt_block_name, aligned_hap_idx, 'del')
                if balanced_translocation_found:
                    event_id += 1
                    if c_splitted:
                        splitted = True
                    break
                else:
                    # reset block assignment for query
                    aligned_hap.discordant_block_assignment[c_mt_block_name] = ''
        if splitted is False:
            run = False

    ## inversion: mt{- k-} wt{k+ -} OR mt{k- -} wt{- k+}
    for aligned_hap in aligned_haplotypes:
        # case: mt{k- -} wt{- k+}, ins of inverted-seg + del of uninverted-seg
        for c_mt_block_name, c_mt_block in aligned_hap.mt_blocks.items():
            if c_mt_block[0][-1] == "+":
                # not inverted
                continue
            if len(aligned_hap.discordant_block_assignment[c_mt_block_name]) > 0:
                continue
            if aligned_hap.next_block(c_mt_block_name) is None or \
                    aligned_hap.next_block(c_mt_block_name) not in aligned_hap.wt_blocks or \
                    len(aligned_hap.discordant_block_assignment[aligned_hap.next_block(c_mt_block_name)]) > 0:
                # this is the last block, next block is not a del-block, OR next del-block has assignment
                continue
            uninverted_segs = [seg[:-1] + '+' for seg in c_mt_block[::-1]]
            next_del_block = aligned_hap.wt_blocks[aligned_hap.next_block(c_mt_block_name)]
            seed_start, seed_end = is_seeded(next_del_block, uninverted_segs, segment_size_dict, indel_direction='both', d=d, eps=eps)
            if seed_start != -1:
                aligned_hap.discordant_block_assignment[c_mt_block_name] = 'inversion,{}'.format(event_id)
                aligned_hap.discordant_block_assignment[aligned_hap.next_block(c_mt_block_name)] = 'inversion,{}'.format(event_id)
                event_id += 1
        # case: mt{- k-} wt{k+ -}, del of uninverted-seg + ins of inverted-seg
        for c_wt_block_name, c_wt_block in aligned_hap.wt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_wt_block_name]) > 0:
                continue
            if aligned_hap.next_block(c_wt_block_name) is None or \
                    aligned_hap.next_block(c_wt_block_name) not in aligned_hap.mt_blocks or \
                    len(aligned_hap.discordant_block_assignment[aligned_hap.next_block(c_wt_block_name)]) > 0:
                # this is the last block, next block is not a ins-block, OR next ins-block has assignment
                continue
            inverted_segs = [seg[:-1] + '-' for seg in c_wt_block[::-1]]
            next_ins_block = aligned_hap.mt_blocks[aligned_hap.next_block(c_wt_block_name)]
            seed_start, seed_end = is_seeded(next_ins_block, inverted_segs, segment_size_dict, indel_direction='both', d=d, eps=eps)
            if seed_start != -1:
                aligned_hap.discordant_block_assignment[c_wt_block_name] = 'inversion,{}'.format(event_id)
                aligned_hap.discordant_block_assignment[aligned_hap.next_block(c_wt_block_name)] = 'inversion,{}'.format(event_id)
                event_id += 1

    ## duplication inversion
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_name, c_mt_block in aligned_hap.mt_blocks.items():
            # left-dup-inv: ins of inverted-seg + concordant block w/ uninverted-seg as seed
            if c_mt_block[0][-1] == "+":
                # not inverted
                continue
            if len(aligned_hap.discordant_block_assignment[c_mt_block_name]) > 0:
                continue
            uninverted_segs = [seg[:-1] + '+' for seg in c_mt_block[::-1]]
            # this is not the last block, AND next block is concordant
            if aligned_hap.next_block(c_mt_block_name) is not None and \
                    aligned_hap.next_block(c_mt_block_name) in aligned_hap.concordant_blocks:
                next_concordant_block = aligned_hap.concordant_blocks[aligned_hap.next_block(c_mt_block_name)]
                seed_start, seed_end = is_seeded(next_concordant_block, uninverted_segs, segment_size_dict, indel_direction='left', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_name] = 'left_duplication_inversion,{}'.format(event_id)
                    event_id += 1
                    continue
            # right-dup-inv: concordant block w/ uninverted-seg as seed + ins of inverted-seg
            # this is not the first block, AND previous block is concordant
            if aligned_hap.previous_block(c_mt_block_name) is not None and \
                    aligned_hap.previous_block(c_mt_block_name) in aligned_hap.concordant_blocks:
                previous_concordant_block = aligned_hap.concordant_blocks[aligned_hap.previous_block(c_mt_block_name)]
                seed_start, seed_end = is_seeded(previous_concordant_block, uninverted_segs, segment_size_dict, indel_direction='right', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_name] = 'right_duplication_inversion,{}'.format(event_id)
                    event_id += 1

    ## tandem-dup: mt{k+ k+}, wt{- k+} OR wt{k+ -}
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_name, c_mt_block in aligned_hap.mt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_mt_block_name]) > 0:
                continue
            block_segs = list(c_mt_block)
            # case: mt{k+ k+}, wt{- k+}, ins seg + concordant block w/ ins-seg as seed
            # this is not the last block, AND next block is concordant
            if aligned_hap.next_block(c_mt_block_name) is not None and \
                    aligned_hap.next_block(c_mt_block_name) in aligned_hap.concordant_blocks:
                next_concordant_block = aligned_hap.concordant_blocks[aligned_hap.next_block(c_mt_block_name)]
                seed_start, seed_end = is_seeded(next_concordant_block, block_segs, segment_size_dict, indel_direction='left', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_name] = 'tandem_duplication,{}'.format(event_id)
                    event_id += 1
                    continue
            # case: mt{k+ k+}, wt{k+ -}, concordant block w/ ins-seg as seed + ins seg
            # this is not the first block, AND previous block is concordant
            if aligned_hap.previous_block(c_mt_block_name) is not None and \
                    aligned_hap.previous_block(c_mt_block_name) in aligned_hap.concordant_blocks:
                previous_concordant_block = aligned_hap.concordant_blocks[aligned_hap.previous_block(c_mt_block_name)]
                seed_start, seed_end = is_seeded(previous_concordant_block, block_segs, segment_size_dict, indel_direction='right', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_name] = 'tandem_duplication,{}'.format(event_id)
                    event_id += 1

    ## deletion: all remaining wt_blocks are deletions
    for aligned_hap in aligned_haplotypes:
        for c_wt_block_name, c_wt_block in aligned_hap.wt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_wt_block_name]) > 0:
                continue
            else:
                aligned_hap.discordant_block_assignment[c_wt_block_name] = 'deletion,{}'.format(event_id)
                event_id += 1

    ## duplicated insertion: all remaining mt_blocks are insertions (i.e. unbalanced translocations)
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_name, c_mt_block in aligned_hap.mt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_mt_block_name]) > 0:
                continue
            else:
                aligned_hap.discordant_block_assignment[c_mt_block_name] = 'insertion,{}'.format(event_id)
                event_id += 1

    ## congregate events for report
    event_blocks = {}
    event_types = {}
    for aligned_hap in aligned_haplotypes:
        event_blocks, event_types = aligned_hap.report_SV(event_blocks, event_types)
    sorted_event_id = sorted(list(event_blocks.keys()))

    # if for the same event ID, we have different event types, we have an issue, otherwise, name the event as the singular name
    conglomerated_event_types = {}
    for c_id, event_type_list in event_types.items():
        c_event_type = event_type_list[0]
        for event_type in event_type_list:
            if event_type != c_event_type:
                raise ValueError('same ID, multiple event types')
        conglomerated_event_types[c_id] = c_event_type

    ## consolidate neighboring blocks for events with paired blocks
    ## balanced translocations will enumerate all possible cycles and mark a set of min-distance cycles as the set of reciprocal translocations
    reciprocal_cycles = []

    def find_reciprocal_cycles(current_path_idx, next_block_name,
                               events_visited, blocks_visited, current_distance):
        # Note: current_block_name and next_block_name are both on current_path_idx
        if next_block_name is None:
            ### not a cycle and ended
            return
        next_event_type = aligned_haplotypes[current_path_idx].discordant_block_assignment[next_block_name].split(',')[0]
        if next_event_type != 'balanced_translocation':
            # cycle exhausted
            return
        next_event_id = int(aligned_haplotypes[current_path_idx].discordant_block_assignment[next_block_name].split(',')[1])
        if next_event_id == events_visited[0]:
            ### cycle closed: next block's event ID returns to the original event ID
            reciprocal_cycles.append({'blocks': blocks_visited,
                                      'events': events_visited,
                                      'distance': current_distance,
                                      'is_cycle': True})
        else:
            if next_event_id in events_visited:
                # cycle is convoluted, traced back to itself in the middle (like small epsilon)
                return
            ### continue search
            ## move across the event ID
            next_event_pair = event_blocks[next_event_id]
            # either the first/second is the paired block
            path_idx1 = int(next_event_pair[0].split('.')[0])
            block_idx1 = next_event_pair[0].split('.')[1]
            path_idx2 = int(next_event_pair[1].split('.')[0])
            block_idx2 = next_event_pair[1].split('.')[1]
            if path_idx1 == current_path_idx and block_idx1 == next_block_name:
                paired_path_idx = path_idx2
                paired_block_name = block_idx2
            elif path_idx2 == current_path_idx and block_idx2 == next_block_name:
                paired_path_idx = path_idx1
                paired_block_name = block_idx1
            else:
                raise RuntimeError('bug in event finding, block pair located was not found')
            ## search left and right: skipping over all blocks of current-type, iterate over all blocks of next type
            left_block, left_skipped_distance = aligned_haplotypes[paired_path_idx].previous_different_type_block(paired_block_name)
            right_block, right_skipped_distance = aligned_haplotypes[paired_path_idx].next_different_type_block(paired_block_name)
            # we opt to do not inherit the skip distance of the same type
            # left_skipped_distance = 0
            # right_skipped_distance = 0
            if left_block not in aligned_haplotypes[paired_path_idx].discordant_block_assignment:
                left_block = None
            if right_block not in aligned_haplotypes[paired_path_idx].discordant_block_assignment:
                right_block = None
            while left_block is not None:
                find_reciprocal_cycles(paired_path_idx, left_block,
                                       events_visited + [next_event_id],
                                       blocks_visited + [f"{current_path_idx}.{next_block_name}", f"{paired_path_idx}.{paired_block_name}"],
                                       current_distance + left_skipped_distance)
                # venture leftward for the same block type
                left_block_segs, _ = aligned_haplotypes[paired_path_idx].get_block_segs_and_block_type(left_block)
                left_skipped_distance += section_size(left_block_segs, segment_size_dict)
                left_block = aligned_haplotypes[paired_path_idx].previous_block_if_same_type(left_block)
            while right_block is not None:
                find_reciprocal_cycles(paired_path_idx, right_block,
                                       events_visited + [next_event_id],
                                       blocks_visited + [f"{current_path_idx}.{next_block_name}", f"{paired_path_idx}.{paired_block_name}"],
                                       current_distance + right_skipped_distance)
                # venture rightward for the same block type
                right_block_segs, _ = aligned_haplotypes[paired_path_idx].get_block_segs_and_block_type(right_block)
                right_skipped_distance += section_size(right_block_segs, segment_size_dict)
                right_block = aligned_haplotypes[paired_path_idx].next_block_if_same_type(right_block)

    for event_id, event_type in event_types.items():
        event_type = event_type[0]
        event_block_list = event_blocks[event_id]
        if event_type == 'inversion':
            if len(event_block_list) != 2:
                raise RuntimeError('inversion needs to have 2 consecutive discordant blocks')
            block1_segs = event_block_list[0].split('.')[2].split('(')[1].split(')')[0].split(',')
            block2_segs = event_block_list[1].split('.')[2].split('(')[1].split(')')[0].split(',')
            block1_size = section_size(block1_segs, segment_size_dict)
            block2_size = section_size(block2_segs, segment_size_dict)
            bp1 = event_block_list[0].split('.')[3]
            bp2 = event_block_list[1].split('.')[4]
            if block1_segs[0][-1] == '-':
                # earlier is ins of inverted-segs
                blocks1_uninverted_segs = [seg[:-1] + '+' for seg in block1_segs[::-1]]
                block2_uninverted_segs = block2_segs
            elif block2_segs[0][-1] == '-':
                # later is ins of inverted-segs
                blocks1_uninverted_segs = block1_segs
                block2_uninverted_segs = [seg[:-1] + '+' for seg in block2_segs[::-1]]
            else:
                raise RuntimeError('inversion not found')
            if block1_size >= block2_size:
                c_block_segs = blocks1_uninverted_segs
            else:
                c_block_segs = block2_uninverted_segs
            c_block_str = '({})'.format(','.join(c_block_segs))
            new_block_list = ['{}.{}.{}.{}.{}'.format(event_block_list[0].split('.')[0],
                                                      event_block_list[0].split('.')[1] + ',' + event_block_list[1].split('.')[1],
                                                      c_block_str + ',' + event_block_list[0].split('.')[2] + ',' + event_block_list[1].split('.')[2],
                                                      bp1,
                                                      bp2)]
            event_blocks[event_id] = new_block_list
        elif event_type == 'balanced_translocation':
            current_path_idx = int(event_block_list[1].split('.')[0])
            current_block_name = event_block_list[1].split('.')[1]
            event_pair = event_blocks[event_id]
            # either the first/second is the paired block
            path_idx1 = int(event_pair[0].split('.')[0])
            block_idx1 = event_pair[0].split('.')[1]
            path_idx2 = int(event_pair[1].split('.')[0])
            block_idx2 = event_pair[1].split('.')[1]
            if path_idx1 == current_path_idx and block_idx1 == current_block_name:
                paired_path_idx = path_idx2
                paired_block_name = block_idx2
            elif path_idx2 == current_path_idx and block_idx2 == current_block_name:
                paired_path_idx = path_idx1
                paired_block_name = block_idx1
            else:
                raise RuntimeError('bug in event finding, block pair located was not found')
            events_visited = [event_id]
            blocks_visited = [f"{current_path_idx}.{current_block_name}", f"{paired_path_idx}.{paired_block_name}"]

            ## search left and right: skipping over all blocks of current-type, iterate over all blocks of next type
            left_block, left_skipped_distance = aligned_haplotypes[paired_path_idx].previous_different_type_block(paired_block_name)
            right_block, right_skipped_distance = aligned_haplotypes[paired_path_idx].next_different_type_block(paired_block_name)
            # we opt to do not inherit the skip distance of the same type
            # left_skipped_distance = 0
            # right_skipped_distance = 0
            if left_block not in aligned_haplotypes[paired_path_idx].discordant_block_assignment:
                left_block = None
            if right_block not in aligned_haplotypes[paired_path_idx].discordant_block_assignment:
                right_block = None
            while left_block is not None:
                find_reciprocal_cycles(paired_path_idx, left_block,
                                       events_visited,
                                       blocks_visited,
                                       left_skipped_distance)
                # venture leftward for the same block type
                left_block_segs, _ = aligned_haplotypes[paired_path_idx].get_block_segs_and_block_type(left_block)
                left_skipped_distance += section_size(left_block_segs, segment_size_dict)
                left_block = aligned_haplotypes[paired_path_idx].previous_block_if_same_type(left_block)
            while right_block is not None:
                find_reciprocal_cycles(paired_path_idx, right_block,
                                       events_visited,
                                       blocks_visited,
                                       right_skipped_distance)
                # venture rightward for the same block type
                right_block_segs, _ = aligned_haplotypes[paired_path_idx].get_block_segs_and_block_type(right_block)
                right_skipped_distance += section_size(right_block_segs, segment_size_dict)
                right_block = aligned_haplotypes[paired_path_idx].next_block_if_same_type(right_block)

    # remove cycles with large distance
    cycle_with_large_distance = []
    for cycle_idx, cycle in enumerate(reciprocal_cycles):
        if cycle['distance'] >= 7000000:
            cycle_with_large_distance.append(cycle_idx)
    reciprocal_cycles = [cycle for cycle_idx, cycle in enumerate(reciprocal_cycles) if cycle_idx not in cycle_with_large_distance]
    # remove identical cycles after removing cycle with large distance, as different directions can result in different distances
    reciprocal_cycles = remove_identical_cycles(reciprocal_cycles)
    selected_cycles = []  # for later sorting purpose

    ## iteratively resolve the most represented event
    while len(reciprocal_cycles) != 0:
        # get most represented event
        event_occurrence = {}
        for cycle in reciprocal_cycles:
            for event_name in cycle['events']:
                if event_name in event_occurrence:
                    event_occurrence[event_name] += 1
                else:
                    event_occurrence[event_name] = 1
        most_represented_event_id = -1
        max_occurrence = -1
        for event_id, occurrence in event_occurrence.items():
            if occurrence > max_occurrence:
                most_represented_event_id = event_id
                max_occurrence = occurrence

        # get cycle with min distance that has this id
        min_distance = -1
        min_distance_cycle = -1
        for cycle_idx, cycle in enumerate(reciprocal_cycles):
            if most_represented_event_id in cycle['events']:
                if min_distance == -1 or cycle['distance'] < min_distance:
                    min_distance = cycle['distance']
                    min_distance_cycle = cycle

        # resolve this cycle
        selected_cycles.append(min_distance_cycle['events'])
        for event_id_idx, event_id_itr in enumerate(min_distance_cycle['events']):
            if event_id_idx >= len(min_distance_cycle['events']) - 1:
                # loop back to the first event
                next_event_id_idx = 0
            else:
                next_event_id_idx = event_id_idx + 1
            if conglomerated_event_types[event_id_itr] != 'balanced_translocation':
                raise RuntimeError('illegal labeling of balanced translocation')
            conglomerated_event_types[event_id_itr] = 'balanced_translocation_associated' + '<{}>'.format(min_distance_cycle['events'][next_event_id_idx])

        # remove all cycles with intersection with the current cycle's events
        resolved_events = set(min_distance_cycle['events'])
        conflicting_cycle_idx = []
        for cycle_idx, cycle in enumerate(reciprocal_cycles):
            current_events = set(cycle['events'])
            intersecting_events = current_events.intersection(resolved_events)
            if len(intersecting_events) > 0:
                conflicting_cycle_idx.append(cycle_idx)
        reciprocal_cycles = [cycle for cycle_idx, cycle in enumerate(reciprocal_cycles) if cycle_idx not in conflicting_cycle_idx]

    ## all remaining balanced translocation are unreciprocal (they do not belong to a cycle)
    nonreciporcal_translocation_events = []
    for event_id, event_type in conglomerated_event_types.items():
        if event_type == 'balanced_translocation':
            nonreciporcal_translocation_events.append(event_id)
            conglomerated_event_types[event_id] = 'balanced_translocation_unassociated'

    output_list = []
    marked_events = []
    ## event order: reciprocal translocation, nonreciprocal translocation, other events
    for cycle in selected_cycles:
        marked_events += cycle
        for event_id in cycle:
            print('event<{}>,type<{}>,blocks<{}>'.format(event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
            output_list.append((event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
    for event_id in nonreciporcal_translocation_events:
        marked_events.append(event_id)
        print('event<{}>,type<{}>,blocks<{}>'.format(event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
        output_list.append((event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
    for event_id in sorted_event_id:
        if event_id not in marked_events:
            print('event<{}>,type<{}>,blocks<{}>'.format(event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
            output_list.append((event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
    return output_list, aligned_haplotypes


def remove_identical_cycles(input_cycles):
    current_idx = 0
    modified_cycles = copy.deepcopy(input_cycles)
    while current_idx != len(modified_cycles):
        current_cycle = modified_cycles[current_idx]['blocks']
        # search for all later cycles if they are the same as itself
        duplicate_cycle_idx = []
        for other_cycle_idx, other_cycle_obj in enumerate(modified_cycles[current_idx + 1:]):
            other_cycle_idx += current_idx + 1
            other_cycle = other_cycle_obj['blocks']
            if current_cycle[0] not in other_cycle or len(current_cycle) != len(other_cycle):
                # cannot be the same
                continue
            other_start_idx = other_cycle.index(current_cycle[0])
            forward_identical = True
            for i in range(len(current_cycle)):
                if current_cycle[i] != other_cycle[(i + other_start_idx) % len(current_cycle)]:
                    forward_identical = False
                    break
            backward_identical = True
            for i in range(len(current_cycle)):
                if current_cycle[i] != other_cycle[(other_start_idx - i) % len(current_cycle)]:
                    backward_identical = False
                    break
            if forward_identical or backward_identical:
                duplicate_cycle_idx.append(other_cycle_idx)
        modified_cycles = [keep for keep_idx, keep in enumerate(modified_cycles) if keep_idx not in duplicate_cycle_idx]
        current_idx += 1
    return modified_cycles


def continuous_extension(input_hap, idx_ptr):
    """
    start from the idx_ptr location, moving leftward until no longer can form a continuous section
    :param idx_ptr:
    :param input_hap:
    :return: final_ptr, where [idx_ptr, final_ptr) is a continuous section
    """
    final_ptr = idx_ptr + 1
    hap_len = len(input_hap)
    if idx_ptr >= hap_len:
        raise ValueError('idx_ptr out of bound')

    current_seg = input_hap[idx_ptr]
    section_sign = current_seg[-1]
    while True:
        if final_ptr == hap_len:
            break

        next_seg = input_hap[final_ptr]
        if next_seg == '-':
            break
        if section_sign == '+':
            if next_seg[-1] != '+':
                break
            elif int(next_seg[:-1]) != int(current_seg[:-1]) + 1:
                break
        elif section_sign == '-':
            if next_seg[-1] != '-':
                break
            elif int(next_seg[:-1]) != int(current_seg[:-1]) - 1:
                break
        else:
            raise ValueError('sign error')
        current_seg = input_hap[final_ptr]
        final_ptr += 1

    return final_ptr


def sublist_idx(list1, list2, inverted=False):
    # list 1 is small, list 2 is large (presumed superlist)
    if inverted:
        copy_list1 = invert_seg_list(list1)
    else:
        copy_list1 = copy.deepcopy(list1)
    for i in range(len(list2) - len(copy_list1) + 1):
        if list2[i:i + len(copy_list1)] == copy_list1:
            return i
    return None


def invert_seg_list(input_seg_list):
    inverted_seg_list = []
    for seg in input_seg_list[::-1]:
        if seg[-1] == '-':
            inverted_seg_list.append(seg[:-1] + '+')
        elif seg[-1] == '+':
            inverted_seg_list.append(seg[:-1] + '-')
        else:
            raise RuntimeError()
    return inverted_seg_list


def max_size_overlap(supergroup_section, query_section, size_dict, d, allow_inversion=False):
    """
    :param allow_inversion: whether the match in query_section can be inverted
    :param supergroup_section:
    :param query_section:
    :param size_dict:
    :param d: min-size accepted; -1 if size does not matter
    :return:
    """
    all_sublists = [query_section[i:j + 1] for i in range(len(query_section)) for j in range(i, len(query_section))]
    sublist_sizes = [section_size(sublist, size_dict) for sublist in all_sublists]
    inverted_sublist_start_idx = -1
    if allow_inversion and query_section[0][-1] == "+":
        # we only try to flip the query section when it is the positive orientation (the negative orientation will have a counterpart to come, which is positive
        inverted_sublists = []
        for sublist in all_sublists:
            inverted_sublists.append(invert_seg_list(sublist))
        inverted_sublist_start_idx = len(all_sublists)
        all_sublists += inverted_sublists
        sublist_sizes += sublist_sizes

    sublists_found = []
    seed_start_location_in_supergroup = []
    max_size = -1
    max_size_sublist_idx = -1
    for query_sublist_idx, query_sublist in enumerate(all_sublists):
        current_size = sublist_sizes[query_sublist_idx]
        if d != -1 and current_size < d:
            seed_start_location_in_supergroup.append(-1)
            continue
        current_start_location = sublist_idx(query_sublist, supergroup_section)
        if current_start_location is not None:
            sublists_found.append(query_sublist)
            seed_start_location_in_supergroup.append(current_start_location)
            if current_size > max_size:
                max_size = current_size
                max_size_sublist_idx = query_sublist_idx
        else:
            seed_start_location_in_supergroup.append(-1)

    if max_size_sublist_idx == -1:
        # not found
        return -1, -1, -1, -1, -1
    max_size_sublist = all_sublists[max_size_sublist_idx]
    max_size_sublist_start_location = seed_start_location_in_supergroup[max_size_sublist_idx]
    max_size_sublist_end_location = seed_start_location_in_supergroup[max_size_sublist_idx] + len(max_size_sublist)
    if inverted_sublist_start_idx == -1:
        sublist_inverted = False
    else:
        if max_size_sublist_idx >= inverted_sublist_start_idx:
            sublist_inverted = True
        else:
            sublist_inverted = False
    return max_size_sublist, max_size_sublist_start_location, max_size_sublist_end_location, max_size, sublist_inverted


def is_seeded(supergroup_section, cont_section, size_dict, indel_direction, d, eps):
    """
    whether cont_section is seeded in input_hap, with significant size (>d)
    :param supergroup_section:
    :param cont_section:
    :param size_dict:
    :param indel_direction: 'both', 'left', OR 'right': which direction to compute indels size on
    :param d: required seed overlap size; -1 to ignore this filtering
    :param eps: epsilon, allowed 1/2 indel size (i.e. 2 * eps >= indel size)
    :return: [start, end) indices in the supergroup
    """
    max_size_sublist, max_size_sublist_start_location, max_size_sublist_end_location, _, _ = max_size_overlap(supergroup_section, cont_section, size_dict, d)
    if max_size_sublist == -1:
        return -1, -1

    del_size = section_size(cont_section, size_dict) - section_size(max_size_sublist, size_dict)
    left_ins_size = section_size(supergroup_section[:max_size_sublist_start_location], size_dict)
    right_ins_size = section_size(supergroup_section[max_size_sublist_end_location:], size_dict)

    if indel_direction == 'both':
        indel_size = del_size + left_ins_size + right_ins_size
        if indel_size <= 2 * eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    elif indel_direction == 'left':
        indel_size = del_size + left_ins_size
        if indel_size <= eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    elif indel_direction == 'right':
        indel_size = del_size + right_ins_size
        if indel_size <= eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    else:
        raise ValueError('indel_direction parameter input illegal')

    # indel size too large
    return -1, -1


def section_size(input_section, size_dict):
    tot_size = 0
    for seg in input_section:
        tot_size += size_dict[seg[:-1]]
    return tot_size


def populate_wt_indexed_lists(mt_path_chrs, wt_path_dict):
    wt_indexed_lists = []
    for path_chr in mt_path_chrs:
        wt_indexed_lists.append(wt_path_dict[path_chr])
    return wt_indexed_lists


def form_dependent_clusters(input_events, i_aligned_haplotypes, index_to_segment_dict):
    """
    This code can be optimized
    :param input_events:
    :return: list of list, inner-list contains the path_id of the cluster, outer-list sorted by min-value path_id in the list
    """

    def hap_origins(aligned_haplotype):
        ## returns the set of chromosomal origins this MT hap contains
        return_set = set()
        for seg in aligned_haplotype.mt_hap:
            indexed_hap = int(seg[:-1])
            return_set.add(index_to_segment_dict[indexed_hap].chr_name)
        return return_set

    def get_aligned_hap(i_id):
        ## find the aligned hap with the id of interest
        for aligned_hap in i_aligned_haplotypes:
            if aligned_hap.id == i_id:
                return aligned_hap
        raise RuntimeError('no aligned hap of this ID is found')

    def get_event_chr(event_info, seg_pattern=r'\(.*?\)'):
        event_scripts = event_info[2]
        event_paths = set()
        event_segs = set()
        for block in event_scripts:
            event_paths.add(int(block.split('.')[0]))
            seg_string = block.split('.')[2]
            matches = re.findall(seg_pattern, seg_string)
            for match in matches:
                match = match.replace('(', '')
                match = match.replace(')', '')
                match = match.replace('+', '')
                match = match.replace('-', '')
                match = match.split(',')
                for indexed_seg in match:
                    event_segs.add(index_to_segment_dict[int(indexed_seg)])

        c_involved_chr = set()
        for event_path in event_paths:
            c_involved_chr = c_involved_chr.union(hap_origins(get_aligned_hap(event_path)))
        for event_seg in event_segs:
            c_involved_chr.add(event_seg.chr_name)
        return c_involved_chr

    ## cluster the chromosomes
    cluster_chrs = []
    for event_info_itr in input_events:
        involved_chr = get_event_chr(event_info_itr)

        ## append/merge cluster if exists, otherwise create new cluster
        cluster_found = []
        for cluster_idx, cluster in enumerate(cluster_chrs):
            # search in all clusters
            chr_intersection = involved_chr.intersection(cluster)
            if len(chr_intersection) > 0:
                cluster_found.append(cluster_idx)

        if len(cluster_found) == 1:
            cluster_chrs[cluster_found[0]] = involved_chr.union(cluster_chrs[cluster_found[0]])
        elif len(cluster_found) > 0:
            # union all chrs found
            new_cluster = set()
            for cluster_idx in cluster_found:
                new_cluster = new_cluster.union(cluster_chrs[cluster_idx])
            new_cluster = new_cluster.union(involved_chr)
            # remove these sets from the previous cluster list
            new_cluster_chrs = [cluster for cluster_idx, cluster in enumerate(cluster_chrs) if cluster_idx not in cluster_found]
            cluster_chrs = new_cluster_chrs
            cluster_chrs.append(new_cluster)
        else:
            # i.e. len(cluster_found) == 0
            cluster_chrs.append(involved_chr)

    dependent_clusters = [[] for _ in cluster_chrs]
    cluster_events = [[] for _ in cluster_chrs]
    ## for each aligned_hap, it should belong to exactly one bin
    for aligned_hap in i_aligned_haplotypes:
        if len(aligned_hap.discordant_block_assignment) == 0:
            # non-event hap
            continue
        hap_chr = hap_origins(aligned_hap)
        cluster_found = []
        for cluster_idx, cluster in enumerate(cluster_chrs):
            if len(cluster.intersection(hap_chr)) > 0:
                cluster_found.append(cluster_idx)
        if len(cluster_found) != 1:
            raise RuntimeError(cluster_found)
        dependent_clusters[cluster_found[0]].append(aligned_hap.id)

    ## for each event, it should belong to exactly one bin
    for event_info_itr in input_events:
        event_chr = get_event_chr(event_info_itr)
        cluster_found = []
        for cluster_idx, cluster in enumerate(cluster_chrs):
            if len(cluster.intersection(event_chr)) > 0:
                cluster_found.append(cluster_idx)
        if len(cluster_found) != 1:
            raise RuntimeError()
        cluster_events[cluster_found[0]].append(event_info_itr)

    ## sort clusters by min_element
    min_values = []
    for bin_idx, bin_itr in enumerate(cluster_chrs):
        chr_numbers = []
        for chr_itr in bin_itr:
            c_chr_number = chr_itr.replace('Chr', '')
            if c_chr_number.upper() == 'X':
                chr_numbers.append(23)
            elif c_chr_number.upper() == 'Y':
                chr_numbers.append(24)
            else:
                chr_numbers.append(int(c_chr_number))
        min_values.append(min(chr_numbers))

    zipped_clusters = zip(min_values, cluster_chrs, dependent_clusters, cluster_events)
    sorted_zipped_clusters = sorted(zipped_clusters, key=lambda x: x[0])
    sorted_min_values, sorted_cluster_chrs, sorted_dependent_clusters, sorted_cluster_events = zip(*sorted_zipped_clusters)

    return sorted_dependent_clusters, sorted_cluster_events


########################PERIPHERAL UTILS#######################
def gather_breakpoints(breakpoints: []):
    all_breakpoints = []
    for c_breakpoint in breakpoints:
        if c_breakpoint[1] is None:
            continue
        all_breakpoints.append((c_breakpoint[0], c_breakpoint[1]))
    return all_breakpoints


def get_genes_near_breakpoints(breakpoints: [(str, int)], proximity=50000):
    breakpoint_ranges = []
    for c_breakpoint in breakpoints:
        breakpoint_ranges.append((c_breakpoint[0],
                                  c_breakpoint[1] - proximity,
                                  c_breakpoint[1] + proximity))
    genes_in_regions = set()
    for breakpoint_range in breakpoint_ranges:
        genes_in_regions = genes_in_regions.union(get_genes_in_region(*breakpoint_range))
    return list(genes_in_regions)


def report_on_genes_based_on_breakpoints(breakpoints):
    breakpoints = gather_breakpoints(breakpoints)
    genes_near_bp = get_genes_near_breakpoints(breakpoints)
    DDG_df = get_DDG_overlapped_genes(genes_near_bp)
    DDG_gene_list, DDG_disease_list = tostring_gene_disease_omim(DDG_df)
    return genes_near_bp, list(zip(DDG_gene_list, DDG_disease_list))


def report_cnv_genes_on_region(chrom, start, end, proximity=50000):
    if start > end:
        temp = start
        start = end
        end = temp
    start = max(0, start - proximity)
    end = end + proximity
    genes = get_genes_in_region('chr' + chrom, start, end)
    DDG_df = get_DDG_overlapped_genes(genes)
    DDG_gene_list, DDG_disease_list = tostring_gene_disease_omim(DDG_df)
    return genes, list(zip(DDG_gene_list, DDG_disease_list))


def conglomerate_cn_genes(input_genes_reports):
    """
    if a gene is +1 and -1 CN, then it is reported in neither;
    :param input_genes_reports:
    :return:
    """
    cnv_genes = {}
    responsible_events = {}
    cnv_highlighted_genes = {}
    for gene_report_idx, gene_report in enumerate(input_genes_reports):
        one_based_index = gene_report_idx + 1
        for gene in gene_report['cnv_genes']:
            if gene in cnv_genes:
                cnv_genes[gene] += gene_report['cnv']
                responsible_events[gene].append(one_based_index)
            else:
                cnv_genes[gene] = gene_report['cnv']
                responsible_events[gene] = [one_based_index]
        for gene_info in gene_report['cnv_genes_highlight']:
            gene_tuple = gene_info[0]
            cnv_highlighted_genes[gene_tuple[0]] = gene_info
    new_cnv_genes = {gene: cnv for gene, cnv in cnv_genes.items() if cnv != 0}
    new_responsible_events = {gene: events for gene, events in responsible_events.items() if cnv_genes[gene] != 0}
    return new_cnv_genes, new_responsible_events, cnv_highlighted_genes


def get_ucsc_url(chrom, start_pos, end_pos, db='hg38'):
    prefix = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db={}' \
             '&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position='.format(db)
    suffix = '{}%3A{}%2D{}'.format(chrom.lower(), start_pos, end_pos)
    return prefix + suffix


def sort_events(input_events):
    def sort_key(event):
        event_type = event[1]
        if 'translocation' in event_type:
            return -1.0
        else:
            path_cluster_idx = float(event[2][0].split('.')[0]) * 100 + block_value(event[2][0].split('.')[1])
            return float(path_cluster_idx)

    return sorted(input_events, key=sort_key)


def format_report(interpreted_events, aligned_haplotypes, index_to_segment_dict, debug=False):
    iscn_events = []
    gene_reports = []
    associated_event_already_reported = []
    for event in interpreted_events:
        event_id = event[0]
        event_type = event[1]
        if event_id in associated_event_already_reported:
            continue
        cn_signature = 0
        cn_changed_genes = []
        cn_changed_genes_highlight = []
        if not event_type.startswith('balanced_translocation'):
            # then only 1 block for each event
            if len(event[2]) != 1:
                raise RuntimeError('not exactly 1 block notation')
            path_idx = int(event[2][0].split('.')[0])
            path_chr = aligned_haplotypes[path_idx].chrom[3:]
            event_segs = event[2][0].split('.')[2].split(')')[0].split('(')[1].split(',')
            left_event_seg = index_to_segment_dict[int(event_segs[0][:-1])]
            right_event_seg = index_to_segment_dict[int(event_segs[-1][:-1])]
            left_event_seg_dir = event_segs[0][-1]
            right_event_seg_dir = event_segs[-1][-1]
            if left_event_seg_dir == '-':
                # we assume the event segments are continuous
                if right_event_seg_dir != '-':
                    raise RuntimeError('event segs not in the same direction')
                # only maintain the directionality if it is an INS()
                if event_type != 'insertion':
                    bp2 = right_event_seg.start
                    bp3 = left_event_seg.end
                    bp2_chr = right_event_seg.chr_name
                    bp3_chr = left_event_seg.chr_name
                else:
                    bp2 = left_event_seg.end
                    bp3 = right_event_seg.start
                    bp2_chr = left_event_seg.chr_name
                    bp3_chr = right_event_seg.chr_name
            else:
                bp2 = left_event_seg.start
                bp3 = right_event_seg.end
                bp2_chr = left_event_seg.chr_name
                bp3_chr = right_event_seg.chr_name
            if event[2][0].split('.')[3] == 'p-ter':
                bp1 = None
                bp1_chr = None
            else:
                bp1_seg = index_to_segment_dict[int(event[2][0].split('.')[3][:-1])]
                bp1_chr = bp1_seg.chr_name
                if event[2][0].split('.')[3][-1] == '+':
                    bp1 = bp1_seg.start
                else:
                    bp1 = bp1_seg.end
            if event[2][0].split('.')[4] == 'q-ter':
                bp4 = None
                bp4_chr = None
            else:
                bp4_seg = index_to_segment_dict[int(event[2][0].split('.')[4][:-1])]
                bp4_chr = bp4_seg.chr_name
                if event[2][0].split('.')[4][-1] == '+':
                    bp4 = bp4_seg.start
                else:
                    bp4 = bp4_seg.end
            if bp1 is not None:
                bp1_band = get_band_location(bp1_chr, bp1)
            else:
                bp1_band = 'pter'
            bp2_band = get_band_location(bp2_chr, bp2)
            bp3_band = get_band_location(bp3_chr, bp3)
            if bp4 is not None:
                bp4_band = get_band_location(bp4_chr, bp4)
            else:
                bp4_band = 'qter'

            if event_type == 'deletion':
                if bp2_band != bp3_band:
                    main_str = 'del({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'del({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'deletion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = -1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(path_chr, bp2, bp3)
            elif event_type == 'inversion':
                if bp2_band != bp3_band:
                    main_str = 'inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'inversion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
            elif event_type == 'tandem_duplication':
                if bp2_band != bp3_band:
                    main_str = 'dup({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'dup({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'tandem duplication on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = 1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(path_chr, bp2, bp3)
            elif event_type == 'left_duplication_inversion':
                if bp2_band != bp3_band:
                    main_str = 'left-dup-inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'left-dup-inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'left duplication inversion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = 1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(path_chr, bp2, bp3)
            elif event_type == 'right_duplication_inversion':
                if bp2_band != bp3_band:
                    main_str = 'right-dup-inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'right-dup-inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'right duplication inversion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = 1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(path_chr, bp2, bp3)
            elif event_type == 'insertion':
                # different report format if insertion is from different chr
                if 'Chr' + path_chr == bp2_chr:
                    # TODO: check ISCN syntax if bp2_band == bp3_band
                    main_str = 'ins({})({}{}{})'.format(path_chr, bp1_band, bp2_band, bp3_band)
                else:
                    main_str = 'ins({};{})({};{}{})'.format(path_chr, bp2_chr, bp1_band, bp2_band, bp3_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                if bp1 is None:
                    bp1_text = bp1_band  # TODO: use 0/len(Chr) for pter/qter
                else:
                    bp1_text = format(bp1, ',d')
                iscn_interpretation = 'duplicated-insertion of Chr{}: {} into Chr{}: {} ({})'. \
                    format(bp2_chr[3:], chr_range, path_chr, bp1_text, bp1_band)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = 1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(bp2_chr[3:], bp2, bp3)
            else:
                # continue
                raise RuntimeError('event not in allowed list')
        elif event_type.startswith("balanced_translocation_associated"):
            # TODO: only works with 2-break reciprocal balanced translocation
            o_event_id = int(event_type.split('<')[1].split('>')[0])
            # get extract the other event
            co_event_idx = -1
            for o_event_idx, event_itr in enumerate(interpreted_events):
                co_event_id = event_itr[0]
                if co_event_id == o_event_id:
                    co_event_idx = o_event_idx
                    break
            # check if o_event associate back
            o_event = interpreted_events[co_event_idx]
            if int(o_event[1].split('<')[1].split('>')[0]) != event_id:
                print('more than 2-break reciprocal translocation detected')
                raise RuntimeError('more than 2-breaks detected')
            # get breakpoints and determine the swaps by number of qter/pter
            c_event_info = event[2]
            o_event_info = o_event[2]
            event_bps = c_event_info[0].split('.')[3:5] + c_event_info[1].split('.')[3:5] + o_event_info[0].split('.')[3:5] + o_event_info[1].split('.')[3:5]
            pter_idx = []
            qter_idx = []
            for event_bp_idx, event_bp_itr in enumerate(event_bps):
                if event_bp_itr == 'p-ter':
                    pter_idx.append(event_bp_idx)
                elif event_bp_itr == 'q-ter':
                    qter_idx.append(event_bp_idx)
            if len(pter_idx) + len(qter_idx) != 2 or len(pter_idx) == len(qter_idx):
                pass
                # raise RuntimeError('non-terminal 2-break reciprocal translocation detected')
            # locate endpoint of event segments, if p-side, choose right, if q-side, choose left
            indexed_event_segs1 = c_event_info[0].split('.')[2].split(')')[0].split('(')[1].split(',')
            indexed_event_segs2 = o_event_info[0].split('.')[2].split(')')[0].split('(')[1].split(',')
            typed_event_segs1 = []
            typed_event_segs2 = []
            for indexed_event_seg in indexed_event_segs1:
                typed_event_segs1.append(index_to_segment_dict[int(indexed_event_seg[:-1])])
            for indexed_event_seg in indexed_event_segs2:
                typed_event_segs2.append(index_to_segment_dict[int(indexed_event_seg[:-1])])
            seg1_left_bp = typed_event_segs1[0].start
            seg1_right_bp = typed_event_segs1[-1].end
            seg2_left_bp = typed_event_segs2[0].start
            seg2_right_bp = typed_event_segs2[-1].end
            seg1_left_band = get_band_location(typed_event_segs1[0].chr_name, seg1_left_bp)
            seg1_right_band = get_band_location(typed_event_segs1[-1].chr_name, seg1_right_bp)
            seg2_left_band = get_band_location(typed_event_segs2[0].chr_name, seg2_left_bp)
            seg2_right_band = get_band_location(typed_event_segs2[-1].chr_name, seg2_right_bp)
            if seg1_left_band[0] == 'p':
                is_pside = True
            elif seg1_left_band[0] == 'q':
                is_pside = False
            else:
                raise RuntimeError('illegal band name')
            chr1 = typed_event_segs1[0].chr_name[3:]  # assumes the segments have the same chr
            chr2 = typed_event_segs2[0].chr_name[3:]
            if is_pside:
                bp1 = typed_event_segs1[-1].end
                bp1_band = seg1_right_band
                bp2 = typed_event_segs2[-1].end
                bp2_band = seg2_right_band
            else:
                bp1 = typed_event_segs1[-1].start
                bp1_band = seg1_left_band
                bp2 = typed_event_segs2[-1].start
                bp2_band = seg2_left_band

            # if there is a sex chr, place it first
            flip_order = False
            if chr2.lower() == 'x':
                if chr1.lower() != 'x':
                    # temp_chr1, temp_bp1, temp_bp1_band = chr1, bp1, bp1_band
                    # chr1, bp1, bp1_band = chr2, bp2, bp2_band
                    # chr2, bp2, bp2_band = temp_chr1, temp_bp1, temp_bp1_band
                    flip_order = True
            elif chr2.lower() == 'y':
                if chr1.lower() != 'x' and chr1.lower() != 'y':
                    # temp_chr1, temp_bp1, temp_bp1_band = chr1, bp1, bp1_band
                    # chr1, bp1, bp1_band = chr2, bp2, bp2_band
                    # chr2, bp2, bp2_band = temp_chr1, temp_bp1, temp_bp1_band
                    flip_order = True
            elif chr1.lower() not in ['x', 'y'] and int(chr1) > int(chr2):
                flip_order = True
            if not flip_order:
                main_str = 't({};{})({};{})'.format(chr1, chr2, bp1_band, bp2_band)
                chr_range1 = chr_range_tostr(seg1_left_bp, seg1_right_bp, seg1_left_band, seg1_right_band)
                chr_range2 = chr_range_tostr(seg2_left_bp, seg2_right_bp, seg2_left_band, seg2_right_band)
                iscn_interpretation = 'balanced translocation between Chr{} and Chr{}, ' \
                                      'between segments Chr{}: {} and Chr{}: {}'. \
                    format(chr1, chr2,
                           chr1, chr_range1,
                           chr2, chr_range2)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(chr1, seg1_left_bp), (chr1, seg1_right_bp),
                                                                                     (chr2, seg2_left_bp), (chr2, seg2_right_bp)])
            else:
                main_str = 't({};{})({};{})'.format(chr2, chr1, bp2_band, bp1_band)
                chr_range1 = chr_range_tostr(seg2_left_bp, seg2_right_bp, seg2_left_band, seg2_right_band)
                chr_range2 = chr_range_tostr(seg1_left_bp, seg1_right_bp, seg1_left_band, seg1_right_band)
                iscn_interpretation = 'balanced translocation between Chr{} and Chr{}, ' \
                                      'between segments Chr{}: {} and Chr{}: {}'. \
                    format(chr2, chr1,
                           chr2, chr_range1,
                           chr1, chr_range2)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(chr2, seg2_left_bp), (chr2, seg2_right_bp),
                                                                                     (chr1, seg1_left_bp), (chr1, seg1_right_bp)])
            associated_event_already_reported.append(o_event_id)
        elif event_type.startswith("balanced_translocation_unassociated"):
            event_info = event[2]
            if event_info[0].split('.')[2].startswith('mt'):
                del_idx = 0
                ins_idx = 1
            else:
                del_idx = 1
                ins_idx = 0
            ins_path_idx = int(event_info[ins_idx].split('.')[0])
            ins_chr = aligned_haplotypes[ins_path_idx].chrom[3:]
            indexed_event_segs = event_info[0].split('.')[2].split(')')[0].split('(')[1].split(',')
            typed_event_segs = []
            for indexed_event_seg in indexed_event_segs:
                typed_event_segs.append(index_to_segment_dict[int(indexed_event_seg[:-1])])
            event_seg_left_bp = typed_event_segs[0].start
            event_seg_right_bp = typed_event_segs[-1].end
            event_seg_left_band = get_band_location(typed_event_segs[0].chr_name, event_seg_left_bp)
            event_seg_right_band = get_band_location(typed_event_segs[-1].chr_name, event_seg_right_bp)
            if typed_event_segs[0].chr_name != typed_event_segs[-1].chr_name:
                raise RuntimeError('diff chr in event segs')
            else:
                event_seg_chr = typed_event_segs[0].chr_name[3:]
            ins_site_left_seg = event_info[ins_idx].split('.')[3]
            if ins_site_left_seg == 'p-ter':
                ins_site_left_bp = 0
                ins_site_left_band = 'p-ter'
            else:
                ins_site_left_bp = index_to_segment_dict[int(ins_site_left_seg[:-1])].start
                ins_site_left_band = get_band_location('chr' + ins_chr, ins_site_left_bp)
            main_str = 'ins-t({};{})({};{}{})'.format(ins_chr, event_seg_chr, ins_site_left_band, event_seg_left_band, event_seg_right_band)
            chr_range = chr_range_tostr(event_seg_left_bp, event_seg_right_bp, event_seg_left_band, event_seg_right_band)
            iscn_interpretation = 'balanced non-reciprocal translocation of Chr{}: {} into Chr{}: {}({})' \
                .format(event_seg_chr, chr_range, ins_chr, ins_site_left_bp, ins_site_left_band)
            bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(event_seg_chr, event_seg_left_bp),
                                                                                 (event_seg_chr, event_seg_left_bp),
                                                                                 (ins_chr, ins_site_left_bp)])
        else:
            raise RuntimeError('illegal type assigned')
        if len(main_str) == 0 or len(iscn_interpretation) == 0:
            raise RuntimeError('missed interpretation')
        if debug:
            iscn_interpretation += '\t {}'.format(event)
        iscn_events.append([main_str, iscn_interpretation])
        gene_reports.append({'bp_genes': bp_genes,
                             'bp_genes_highlight': bp_genes_highlight,
                             'cnv': cn_signature,
                             'cnv_genes': cn_changed_genes,
                             'cnv_genes_highlight': cn_changed_genes_highlight})
    if len(iscn_events) != len(gene_reports):
        raise RuntimeError('unmatched reports')
    return iscn_events, gene_reports


def format_genes_report(genes_report):
    conglomerated_cnv_genes, cnv_genes_events, cnv_highlighted_genes = conglomerate_cn_genes(genes_report)
    formatted_genes_report = []
    for event_idx, report_dict in enumerate(genes_report):
        one_based_idx = event_idx + 1
        for entry in report_dict['bp_genes_highlight']:
            gene_entry = entry[0]
            disease_entry = entry[1]
            gene_name = gene_entry[0]
            gene_omim = gene_entry[1]
            disease_names = []
            disease_omim = []
            for disease in disease_entry:
                disease_names.append(disease[0])
                disease_omim.append(disease[1])
            formatted_genes_report.append({'SV': one_based_idx,
                                           'rationale': 'breakpoint proximal',
                                           'gene name': gene_name,
                                           'gene omim': gene_omim,
                                           'diseases': disease_names,
                                           'disease omims': disease_omim})
    for gene_name, cnv in conglomerated_cnv_genes.items():
        if gene_name in cnv_highlighted_genes:
            related_events = cnv_genes_events[gene_name]
            gene_omim = cnv_highlighted_genes[gene_name][0][1]
            disease_names = []
            disease_omim = []
            for disease in cnv_highlighted_genes[gene_name][1]:
                disease_names.append(disease[0])
                disease_omim.append(disease[1])
            if cnv > 0:
                cnv_str = "+" + str(cnv)
            else:
                cnv_str = str(cnv)
            formatted_genes_report.append({'SV': ','.join([str(i) for i in sorted(related_events)]),
                                           'rationale': 'CN' + cnv_str,
                                           'gene name': gene_name,
                                           'gene omim': gene_omim,
                                           'diseases': disease_names,
                                           'disease omims': disease_omim})
    return formatted_genes_report


def chr_range_tostr(bpa, bpb, bpa_band, bpb_band):
    return "{}-{} ({} - {})".format(format(bpa, ',d'), format(bpb, ',d'), bpa_band, bpb_band)


def test_interpreter():
    i_mt_hap1 = ['1+', '1+', '2+', '3-', '4+', '9+', '10+']
    i_mt_hap2 = ['7+', '8+', '5+', '6+']
    i_mt_hap3 = ['11+', '12+', '11-']
    i_mt_hap4 = ['14-', '13+', '14+']
    i_mt_hap5 = ['15+', '16+', '17-', '16-', '18+']
    i_wt_hap1 = ['1+', '2+', '3+', '4+', '5+', '6+']
    i_wt_hap2 = ['7+', '8+', '9+', '10+']
    i_wt_hap3 = ['11+', '12+']
    i_wt_hap4 = ['13+', '14+']
    i_wt_hap5 = ['15+', '16+', '17+', '18+']
    i_mt_list = [i_mt_hap1, i_mt_hap2, i_mt_hap3, i_mt_hap4, i_mt_hap5]
    i_wt_list = [i_wt_hap1, i_wt_hap2, i_wt_hap3, i_wt_hap4, i_wt_hap5]
    i_size_dict = {str(i): 1 for i in range(19)}
    # print(lcs(i_mt_hap, i_wt_hap, i_size_dict))
    out = interpret_haplotypes(i_mt_list, i_wt_list, i_size_dict, 1, 1)
    print(out)


def test_3break_qter():
    i_wt_hap0 = ['1+', '2+', '3+']
    i_wt_hap1 = ['4+', '5+', '6+']
    i_wt_hap2 = ['7+', '8+', '9+']
    i_mt_hap0 = ['1+', '2+', '9+']
    i_mt_hap1 = ['4+', '5+', '3+']
    i_mt_hap2 = ['7+', '8+', '6+']
    chrom_id = ['Chr1', 'Chr2', 'Chr3']
    i_mt_list = [i_mt_hap0, i_mt_hap1, i_mt_hap2]
    i_wt_list = [i_wt_hap0, i_wt_hap1, i_wt_hap2]
    i_size_dict = {str(i): 1 for i in range(19)}
    out = interpret_haplotypes(i_mt_list, i_wt_list, chrom_id, i_size_dict, 1, 1)
    print(out)


def test_reciprocal_trans():
    i_wt_hap0 = ['1+', '2+', '3+']
    i_wt_hap1 = ['4+', '5+', '6+']
    i_mt_hap0 = ['1+', '2+', '6+']
    i_mt_hap1 = ['3+', '4+', '5+']
    chrom_id = ['Chr1', 'Chr2']
    i_mt_list = [i_mt_hap0, i_mt_hap1]
    i_wt_list = [i_wt_hap0, i_wt_hap1]
    i_size_dict = {str(i): 1 for i in range(19)}
    out = interpret_haplotypes(i_mt_list, i_wt_list, chrom_id, i_size_dict, 1, 1)
    print(out)


def test_ucsc_url():
    x = get_ucsc_url('Chr2', 130271298, 131386307)
    print(x)


# def test_segs_union():
#     l1 = ['7+', '8+']
#     l2 = ['7+']
#     l3 = segs_union(l1, l2)
#     print(l3)

if __name__ == "__main__":
    test_ucsc_url()
