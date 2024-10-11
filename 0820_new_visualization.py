from KT_visualizer import *
from generate_content import *


omkar_file_path = 'example_input/b17_karsim/23X_1q21_recurrent_microduplication_r1.1.txt'
mt_indexed_lists, mt_path_chrs, segment_to_index_dict, segment_size_dict = read_OMKar_to_indexed_list(omkar_file_path)
mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
wt_path_dict = generate_wt_from_OMKar_output(segment_to_index_dict)
wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)

c_vis_input = generate_cytoband_visualizer_input(events, aligned_haplotypes, segment_to_index_dict)
vis_input_used = [c_vis_input[0], c_vis_input[1], c_vis_input[2], c_vis_input[3]]
make_image(vis_input_used, max_chr_length(vis_input_used), 'simulation_report/test_image/cytoband', IMG_LENGTH_SCALE_VERTICAL_SPLIT)

c_vis_input = generate_segment_visualizer_input(events, aligned_haplotypes, segment_to_index_dict)
vis_input_used = [c_vis_input[0], c_vis_input[1], c_vis_input[2], c_vis_input[3]]
make_image(vis_input_used, max_chr_length(vis_input_used), 'simulation_report/test_image/segmentview', IMG_LENGTH_SCALE_VERTICAL_SPLIT)