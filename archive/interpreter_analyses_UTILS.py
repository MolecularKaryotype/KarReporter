from KarUtils.read_OMKar_output import *
from generate_content import *
from KarUtils.utils import *

import os

# data_dir = '/media/zhaoyang-new/workspace/sunnyside/OMKar_output_paths/'
# data_dir = '/media/zhaoyang-new/workspace/keyhole/OMKar_output_paths/'
data_dir = '/media/zhaoyang-new/workspace/paul_dremsek/omkar_output/'
# data_dir = '/Users/zhaoyangjia/Library/CloudStorage/OneDrive-UCSanDiego/Bafna_Lab/Paul_Dremsek_OMKar_output/'
forbidden_region_file = "Metadata/acrocentric_telo_cen.bed"


def batch_interpret(omkar_output_dir):
    files = [file for file in os.listdir(omkar_output_dir)]
    files = sorted(files, key=dremsek_file_keys)

    for file in files:
        if file == '17.txt':
        # if True:
            file_path = omkar_output_dir + file
            print(file)
            mt_indexed_lists, mt_path_chrs, segment_dict, segment_size_dict = read_OMKar_to_indexed_list(file_path, forbidden_region_file)
            mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
            wt_path_dict = generate_wt_from_OMKar_output(segment_dict)
            wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
            events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)
            main_bullets, sub_bullets = format_report(events, aligned_haplotypes, reverse_dict(segment_dict))
            for bullet_idx, main_bullet in enumerate(main_bullets):
                sub_bullet = sub_bullets[bullet_idx]
                print(main_bullet)
                print(sub_bullet)
            # generate_pdf_report('report_PDF/', file.split('.')[0] + '.report.pdf', main_bullets, sub_bullets)
            print()


def populate_wt_indexed_lists(mt_path_chrs, wt_path_dict):
    wt_indexed_lists = []
    for path_chr in mt_path_chrs:
        wt_indexed_lists.append(wt_path_dict[path_chr])
    return wt_indexed_lists


def dremsek_file_keys(f):
    return int(f.split('.')[0])


if __name__ == "__main__":
    batch_interpret(data_dir)
