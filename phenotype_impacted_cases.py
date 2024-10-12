import os

from Report_Genes import *
from KT_visualizer import *
from KarInterpreter import *
from KarUtils import *


def get_boundary_coordinate(boundary_type, boundary_segment, segment_index_to_type_dict):
    if boundary_type not in ['left', 'right']:
        raise ValueError()

    boundary_segment_index = boundary_segment[:-1]
    boundary_segment_orientation = boundary_segment[-1]
    boundary_segment_object = segment_index_to_type_dict[int(boundary_segment_index)]
    output_chrom = boundary_segment_object.chr_name

    if boundary_segment_orientation not in ['-', '+']:
        raise RuntimeError()

    if boundary_type == 'left':
        if boundary_segment_orientation == '-':
            output_pos = boundary_segment_object.end
        else:
            output_pos = boundary_segment_object.start
    else:
        if boundary_segment_orientation == '-':
            output_pos = boundary_segment_object.start
        else:
            output_pos = boundary_segment_object.end
    return {'chrom': output_chrom, 'pos': output_pos}


def report_impacted_genes_by_boundaries(i_left_boundary_seg, i_right_boundary_seg, segment_index_to_type_dict, position_estimate):
    left_boundary_coordinate = get_boundary_coordinate('left', i_left_boundary_seg, segment_index_to_type_dict)
    right_boundary_coordinate = get_boundary_coordinate('right', i_left_boundary_seg, segment_index_to_type_dict)
    bps = [(left_boundary_coordinate['chrom'], left_boundary_coordinate['pos']), (right_boundary_coordinate['chrom'], right_boundary_coordinate['pos'])]
    genes_near_bp = get_genes_near_breakpoints(bps, proximity=position_estimate)
    DDG_df = get_DDG_overlapped_genes(genes_near_bp)
    DDG_df = DDG_df[DDG_df['allelic requirement'] == 'monoallelic_autosomal' & DDG_df['mutation consequence'].isin('absent gene product')]


def report_cnv_minus_one_impacted_genes(cnv_segments_dict, segment_index_to_type_dict):
    pass


forbidden_region_file = get_metadata_file_path('acrocentric_telo_cen.bed')
ddg2p_file = get_metadata_file_path('DDG2P_14_11_2023.csv')
omkar_output_dir = '/Volumes/Data/keyhole_0926_paths/'
file_of_interest = []
skip = []
interpretations = {}

files = [file for file in os.listdir(omkar_output_dir)]
# files = sorted(files, key=int_file_keys)
for file in files:
    if file.startswith('.'):
        # macOS trashfiles
        continue
    file_event_type_reports = []
    if file_of_interest:
        if file.split('.')[0] not in file_of_interest:
            continue
    if skip:
        if file.split('.')[0] in skip:
            continue
    filename = file.split('.')[0]
    file_path = omkar_output_dir + file
    print(file)
    mt_indexed_lists, mt_path_chrs, segment_to_index_dict, segment_size_dict = read_OMKar_to_indexed_list(file_path, forbidden_region_file)
    index_to_segment_dict = reverse_dict(segment_to_index_dict)
    mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
    wt_path_dict = generate_wt_from_OMKar_output(segment_to_index_dict)
    wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
    events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)
    phenotype_report = []
    cnv_segments = {}
    for event_info in events:
        c_phenotype_report = {}
        event_type = event_info[1]
        if 'translocation' in event_type:
            event_cnv = 0
            annotated_event_type = 'balanced translocation/transposition'
        elif event_type == 'inversion':
            event_cnv = 0
            annotated_event_type = 'inversion'
        elif event_type in ['deletion']:
            event_cnv = -1
            annotated_event_type = event_type
        elif event_type in ['insertion', 'tandem_duplication', 'left_duplication_inversion', 'right_duplication_inversion']:
            event_cnv = 1
            annotated_event_type = event_type
        else:
            print(event_info)
            raise RuntimeError
        c_phenotype_report['CNV'] = event_cnv
        c_phenotype_report['event_type'] = annotated_event_type
        c_phenotype_report['event_info'] = event_info[2]
        phenotype_report.append(c_phenotype_report)
        if annotated_event_type == 'balanced translocation/transposition':
            # take boundaries from the wt
            if 'wt' in event_info[2][0]:
                wt_segments = event_info[2][0].split('(')[1].split(')')[0].split(',')
            else:
                wt_segments = event_info[2][1].split('(')[1].split(')')[0].split(',')
            left_boundary_seg = wt_segments[0]
            right_boundary_seg = wt_segments[-1]
            report_impacted_genes_by_boundaries(left_boundary_seg, right_boundary_seg, index_to_segment_dict, 20000)
        elif annotated_event_type == 'inversion':
            inverted_segments = event_info[2][0].split('(')[1].split(')')[0].split(',')
            left_boundary_seg = inverted_segments[0]
            right_boundary_seg = inverted_segments[-1]
            report_impacted_genes_by_boundaries(left_boundary_seg, right_boundary_seg, index_to_segment_dict, 20000)
        elif annotated_event_type == 'deletion':
            deleted_segments = event_info[2][0].split('(')[1].split(')')[0].split(',')
            for segment in deleted_segments:
                seg_index = int(segment[:-1])
                if seg_index in cnv_segments:
                    cnv_segments[seg_index] -= 1
                else:
                    cnv_segments[seg_index] = -1
        elif annotated_event_type in ['insertion', 'tandem_duplication', 'left_duplication_inversion', 'right_duplication_inversion']:
            amplified_segments = event_info[2][0].split('(')[1].split(')')[0].split(',')
            for segment in amplified_segments:
                seg_index = int(segment[:-1])
                if seg_index in cnv_segments:
                    cnv_segments[seg_index] += 1
                else:
                    cnv_segments[seg_index] = 1
        else:
            continue




    interpretations[file] = {'events': events, 'aligned_haplotypes': aligned_haplotypes, 'phenotypes': phenotype_report}

