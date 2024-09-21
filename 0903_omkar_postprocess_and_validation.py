import os, sys

from KarUtils import *

### test_validation
# test_negative_dir = "example_input/b17_karsim/"
# test_negative_dir = "D:\\dremsek_OMKar_output_paths\\"
# test_positive_dir = "example_input/validation_positive_tests/"
# validate_OMKar_output(test_negative_dir)


# seg1 = Segment('chr1', 1, 2)
# seg2 = Segment('chr1', 2, 3)
# seg3 = Segment('chr1', 1, 2)
# seg4 = Segment('chr1', 2, 3)
# seg2.invert()
# seg2.invert()
#
# x = [(-2, [seg1, seg2])]
# print((-2, [seg3, seg4] )in x)


### heapq
# import heapq
#
# heap = []
# x = 'abcd'
# for start_idx in range(len(x)):
#     for end_idx in range(start_idx + 1, len(x)):
#         lst = x[start_idx:end_idx + 1]
#         heapq.heappush(heap, (-len(lst), lst))
# print(heap)


# def find_longest_sublist(lists):
#     def find_increasing_sublists(lst):
#         """Helper function to find all strictly increasing sublists by +1 in a list."""
#         sublists = []
#         start = 0
#
#         for i in range(1, len(lst)):
#             if lst[i] != lst[i - 1] + 1:
#                 if i - start > 1:
#                     sublists.append(lst[start:i])
#                 start = i
#
#         # Add the last found sublist if it exists
#         if len(lst) - start > 1:
#             sublists.append(lst[start:])
#
#         return sublists
#
#     def find_common_sublists(list1_sublists, list2_sublists):
#         """Helper function to find common sublists across two sets of sublists."""
#         common_sublists = []
#         for sublist1 in list1_sublists:
#             for sublist2 in list2_sublists:
#                 # Check if the sublists are identical
#                 if sublist1 == sublist2:
#                     common_sublists.append(sublist1)
#         return common_sublists
#
#     # Find strictly increasing sublists for each ordered list
#     all_sublists = [find_increasing_sublists(lst) for lst in lists]
#
#     # Compare all sublists between the first and all others
#     common_sublists = all_sublists[0]
#     for i in range(1, len(all_sublists)):
#         common_sublists = find_common_sublists(common_sublists, all_sublists[i])
#
#     # If multiple common sublists are found, return the longest one
#     if common_sublists:
#         return max(common_sublists, key=len)
#     else:
#         return []
#
#
# # Test the function with the examples
# input1 = [[1, 2, 3, 4], [1, 2, 2, 3, 4]]
# input2 = [[1, 2, 3, 4], [1, 2, 3, 4, 3]]
#
# print(find_longest_sublist(input1))  # Output: [3, 4]
# print(find_longest_sublist(input2))  # Output: [1, 2]


### MK's I/O
# input_dir = '0903_test_files/input/'
# output_dir = '0903_test_files/output/'
# for file_name in os.listdir(input_dir):
#     print(file_name)
#     input_filepath = os.path.join(input_dir, file_name)
#     output_filepath = os.path.join(output_dir, file_name)
#     path_list, segment_dict = read_OMKar_output(input_filepath, return_segment_dict=True)
#     # segment_obj_to_idx_dict = reverse_dict(segment_dict)
#     processed_path_list, segment_obj_to_idx_dict = post_process_OMKar_output(path_list)
#     write_MK_file(output_filepath, processed_path_list, segment_obj_to_idx_dict)


### debug post-processing
input_filepath = "D:\\keyhole_0717_paths\\2254.txt"
output_filepath = "0903_test_files/output/2254.txt"
path_list, segment_dict = read_OMKar_output(input_filepath, return_segment_dict=True)
processed_path_list, segment_obj_to_idx_dict = post_process_OMKar_output(path_list)
write_MK_file(output_filepath, processed_path_list, segment_obj_to_idx_dict)


# batch_post_process_OMKar_output('0903_test_files/input/', '0903_test_files/test2/')

