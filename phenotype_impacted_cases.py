import os
from collections import Counter

from numpy import number

from Report_Genes import *
from KT_visualizer import *
from KarInterpreter import *
from KarUtils import *


def count_number_of_chrx(i_mt_path_chrs):
    chrx_count = 0
    for i in i_mt_path_chrs:
        if i == 'ChrX':
            chrx_count += 1
    return chrx_count


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


def report_impacted_genes_by_boundaries(i_left_boundary_seg, i_right_boundary_seg, segment_index_to_type_dict, position_estimate, number_of_chrx):
    left_boundary_coordinate = get_boundary_coordinate('left', i_left_boundary_seg, segment_index_to_type_dict)
    right_boundary_coordinate = get_boundary_coordinate('right', i_right_boundary_seg, segment_index_to_type_dict)
    bps = [(left_boundary_coordinate['chrom'], left_boundary_coordinate['pos']), (right_boundary_coordinate['chrom'], right_boundary_coordinate['pos'])]
    genes_near_bp = get_genes_near_breakpoints(bps, proximity=position_estimate)
    once_interruption_gene = [gene for gene, cnt in Counter(genes_near_bp).items() if cnt == 1]
    multiple_interruption_gene = [gene for gene, cnt in Counter(genes_near_bp).items() if cnt > 1]
    DDG_df_once = get_DDG_overlapped_genes(once_interruption_gene)
    DDG_df_multiple = get_DDG_overlapped_genes(multiple_interruption_gene)
    if number_of_chrx == 1:
        DDG_df_once = DDG_df_once[
            (DDG_df_once['allelic requirement'].isin(['monoallelic_autosomal', 'monoallelic_X_hem'])) & (DDG_df_once['mutation consequence'].isin(['absent gene product']))]
    elif number_of_chrx > 1:
        DDG_df_once = DDG_df_once[
            (DDG_df_once['allelic requirement'].isin(['monoallelic_autosomal', 'monoallelic_X_het'])) & (
                DDG_df_once['mutation consequence'].isin(['absent gene product']))]
    else:
        raise RuntimeError()
    DDG_df_multiple = DDG_df_multiple[
        (DDG_df_multiple['allelic requirement'].isin(['biallelic_autosomal'])) & (DDG_df_multiple['mutation consequence'].isin(['absent gene product']))]
    result = DDG_df_once[['gene symbol', 'disease name', 'mutation consequence', 'allelic requirement', 'confidence category']].to_dict(orient='records')
    result += DDG_df_multiple[['gene symbol', 'disease name', 'mutation consequence', 'allelic requirement', 'confidence category']].to_dict(orient='records')
    return result


def report_impacted_genes_by_CN(cnv_segments_dict, segment_index_to_type_dict, position_estimate, number_of_chrx):
    def check_presence_in_result(gene, disease):
        for itr in result:
            if itr['gene symbol'] == gene and itr['disease name'] == disease:
                return True

    result = []
    for c_segment, CNV in cnv_segments_dict.items():
        left_boundary_pos = segment_index_to_type_dict[c_segment].start
        right_boundary_pos = segment_index_to_type_dict[c_segment].end
        chrom = segment_index_to_type_dict[c_segment].chr_name[3:]
        start_pos = max(0, left_boundary_pos - position_estimate)
        end_pos = right_boundary_pos + position_estimate
        genes = get_genes_in_region('chr' + chrom, start_pos, end_pos)
        DDG_df = get_DDG_overlapped_genes(genes)
        if CNV == -1:
            if number_of_chrx == 1:
                DDG_df = DDG_df[
                    (DDG_df['allelic requirement'].isin(['monoallelic_autosomal', 'monoallelic_X_hem'])) & (
                        DDG_df['mutation consequence'].isin(['absent gene product']))]
            elif number_of_chrx > 1:
                DDG_df = DDG_df[
                    (DDG_df['allelic requirement'].isin(['monoallelic_autosomal', 'monoallelic_X_het'])) & (
                        DDG_df['mutation consequence'].isin(['absent gene product']))]
            else:
                raise RuntimeError()
        elif CNV == -2:
            DDG_df = DDG_df[
                (DDG_df['allelic requirement'].isin(['monoallelic_autosomal', 'biallelic_autosomal', 'monoallelic_X_hem', 'monoallelic_X_het'])) & (
                    DDG_df['mutation consequence'].isin(['absent gene product']))]
        else:
            # almost no gain in gene product causal in DDG2P
            continue
        c_result = DDG_df[['gene symbol', 'disease name', 'mutation consequence', 'allelic requirement', 'confidence category']].to_dict(orient='records')
        for entry in c_result:
            entry_gene = entry['gene symbol']
            entry_disease = entry['disease name']
            if not check_presence_in_result(entry_gene, entry_disease):
                # don't add duplicates, they are only possible through the boundary-allowance expansion
                result.append(entry)
    return result


forbidden_region_file = get_metadata_file_path('acrocentric_telo_cen.bed')
ddg2p_file = get_metadata_file_path('DDG2P_14_11_2023.csv')
# omkar_output_dir = '/Volumes/Data/keyhole_0926_paths/'
# omkar_output_dir = '/media/zhaoyang-new/workspace/keyhole/0926_paths/'
# omkar_output_dir = '/media/zhaoyang-new/workspace/keyhole/0729_paths/'
omkar_output_dir = '/media/zhaoyang-new/workspace/sunnyside/0926_paths/'
keyhole_previously_unexplained_G2P = [889, 1251, 1294, 1301, 1404, 1544, 1562, 2081, 2085, 2243, 2245, 2249, 2254, 2276, 2280, 2281, 2282, 2283, 2284, 2327, 2335]
keyhole_newly_found_G2P = [2081, 2276, 2280, 2281, 2282]
mosaic_or_polyploidy = [510, 2173, 1530, 1996]
debug = [2243]

file_of_interest = []
skip = mosaic_or_polyploidy
interpretations = {}

files = [file for file in os.listdir(omkar_output_dir)]
# files = sorted(files, key=int_file_keys)
for file in files:
    if file.startswith('.'):
        # macOS trashfiles
        continue
    file_event_type_reports = []
    if file_of_interest:
        if int(file.split('.')[0]) not in file_of_interest:
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
    count_of_xchr = count_number_of_chrx(mt_path_chrs)
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
            c_phenotype_report['boundaries_genes'] = report_impacted_genes_by_boundaries(left_boundary_seg, right_boundary_seg, index_to_segment_dict, 20000, count_of_xchr)
        elif annotated_event_type == 'inversion':
            inverted_segments = event_info[2][0].split('(')[1].split(')')[0].split(',')
            left_boundary_seg = inverted_segments[0]
            right_boundary_seg = inverted_segments[-1]
            c_phenotype_report['boundaries_genes'] = report_impacted_genes_by_boundaries(left_boundary_seg, right_boundary_seg, index_to_segment_dict, 20000, count_of_xchr)
        elif annotated_event_type == 'deletion':
            deleted_segments = event_info[2][0].split('(')[1].split(')')[0].split(',')
            for segment in deleted_segments:
                seg_index = int(segment[:-1])
                if seg_index in cnv_segments:
                    cnv_segments[seg_index] -= 1
                else:
                    cnv_segments[seg_index] = -1
            c_phenotype_report['boundaries_genes'] = []
        elif annotated_event_type in ['insertion', 'tandem_duplication', 'left_duplication_inversion', 'right_duplication_inversion']:
            amplified_segments = event_info[2][0].split('(')[1].split(')')[0].split(',')
            for segment in amplified_segments:
                seg_index = int(segment[:-1])
                if seg_index in cnv_segments:
                    cnv_segments[seg_index] += 1
                else:
                    cnv_segments[seg_index] = 1
            c_phenotype_report['boundaries_genes'] = []
        else:
            continue

    interpretations[file] = {'events': events,
                             'aligned_haplotypes': aligned_haplotypes,
                             'phenotypes': phenotype_report,
                             'CN_genes': report_impacted_genes_by_CN(cnv_segments, index_to_segment_dict, 20000, count_of_xchr)}


def format_dicts_list(dicts_list):
    formatted_lines = []
    for d in dicts_list:
        formatted_lines.append(str(d))
    return "\n\t".join(formatted_lines)


summary_event_counts = {}
file_with_balanced_boundary_G2P = []
file_with_unbalanced_G2P = []
for file, interpret in interpretations.items():
    if interpret['CN_genes']:
        formated_string = format_dicts_list(interpret['CN_genes'])
        print(f"{file}, CN: {formated_string}")
        file_with_unbalanced_G2P.append(file)
    phenotype = interpret['phenotypes']
    for phenotype_report in phenotype:
        boundary_gene_report = phenotype_report['boundaries_genes']
        event_type = phenotype_report['event_type']
        if boundary_gene_report:
            file_with_balanced_boundary_G2P.append(file)
            print(f"{file}, boundary: {boundary_gene_report}")
            break
        if event_type not in summary_event_counts:
            summary_event_counts[event_type] = 1
        else:
            summary_event_counts[event_type] += 1

balanced_event_count = 0
unbalanced_event_count = 0
for event_type_name, count in summary_event_counts.items():
    if event_type_name in ['balanced translocation/transposition', 'inversion']:
        balanced_event_count += count
    elif event_type_name in ['deletion', 'insertion', 'tandem_duplication', 'left_duplication_inversion', 'right_duplication_inversion']:
        unbalanced_event_count += count
    else:
        print(event_type_name)
        raise RuntimeError()

print(f"summary_event_counts: {summary_event_counts}")
print(f"file_with_balanced_boundary_G2P ({len(file_with_balanced_boundary_G2P)}): {sorted(file_with_balanced_boundary_G2P)}")
print(f"file_with_unbalanced_G2P ({len(file_with_unbalanced_G2P)}): {sorted(file_with_unbalanced_G2P)}")
print(f"balanced_event_count: {balanced_event_count}")
print(f"unbalanced_event_count: {unbalanced_event_count}")

