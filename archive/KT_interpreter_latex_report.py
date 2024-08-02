import re
import datetime
import os

from KarUtils import *
from KT_visualizer import *


def generate_latex_frontpage(title,
                             genefile_name,
                             breakpoint_reporting_proximity,
                             interpretation_insertion_threshold,
                             interpretation_deletion_threshold,
                             omkar_version='xxx'):
    output_str = ''

    def append_action(input_str, o):
        o = o + input_str + "\n"
        return o

    output_str += "\\catcode`\_=12\n"
    output_str += "\\documentclass[12pt]{article}\n"
    output_str += "\\input{macros}\n"
    output_str += "\\usepackage[letterpaper, margin=0.75in]{geometry}\n"
    output_str += "\\setcounter{secnumdepth}{0}\n"
    output_str += "\\usepackage{graphicx}\n"
    output_str += "\\usepackage{setspace}\n"
    output_str += "\\usepackage{titling}\n"
    output_str += "\\usepackage{enumitem}\n"
    # output_str += "\\usepackage[T1]{fontenc}\n"
    output_str += "\\renewcommand\\maketitlehookc{\\vspace{-10ex}}\n"
    # output_str += "\\graphicspath{ {./images/} }\n"

    output_str += "\n"
    output_str += "\\begin{document}\n"
    output_str += "\n"
    output_str += "\\title{" + title + "}\n"
    today = datetime.today()
    formatted_date = today.strftime("%b.%dth, %Y")
    output_str += "\\date{" + str(formatted_date) + "}\n"
    # output_str += "\\maketitle\n"
    output_str += "\n"
    # output_str += "\\textbf{{Samples Included: }} "
    # output_str += ', '.join(sample_names) + "\n"
    output_str += "\n"

    output_str = append_action('\\paragraph{Parameters}', output_str)
    output_str = append_action('\\begin{packed_itemize}', output_str)
    output_str = append_action('\\item OMKar version: {}'.format(omkar_version), output_str)
    output_str = append_action('\\item Genome assembly: hg38', output_str)
    output_str = append_action('\\item Genes: protein coding', output_str)
    output_str = append_action('\\item Breakpoint Gene Reporting Proximity: {}'.format(breakpoint_reporting_proximity), output_str)
    output_str = append_action('\\item Threshold for event insertion size: {}'.format(interpretation_insertion_threshold), output_str)
    output_str = append_action('\\item Threshold for event deletion size: {}'.format(interpretation_deletion_threshold), output_str)
    output_str = append_action('\\item Supported SV types:', output_str)
    output_str = append_action('\\begin{packed_itemize}', output_str)
    output_str = append_action('\\item Deletion', output_str)
    output_str = append_action('\\item Inversion', output_str)
    output_str = append_action('\\item Single/repeated Tandem-duplication', output_str)
    output_str = append_action('\\item Left/right Duplication-inversion', output_str)
    output_str = append_action('\\item 2/multi-break Reciprocal-balanced-translocation', output_str)
    output_str = append_action('\\item Nonreciprocal-balanced-translocation', output_str)
    output_str = append_action('\\item Duplicated-insertion', output_str)
    output_str = append_action('\\end{packed_itemize}', output_str)
    output_str = append_action('\\end{packed_itemize}', output_str)
    output_str += "\n"
    output_str += "\\newpage\n"

    return output_str


def batch_generate_latex_case_str(omkar_output_dir, image_dir, compile_image=False):
    """
    :param omkar_output_dir: assume files are named as (int).txt
    :param image_dir: relative path to the latex DIR, assume exists an image file with the same name, but (int).pdf
    :return:
    """
    final_str = ""
    cases_with_events = []
    files = [file for file in os.listdir(omkar_output_dir)]
    # files = sorted(files, key=int_file_keys)  # TODO: turn back on for files with int ID

    # highlight_files = ['23X_1q21_recurrent_microduplication_r1',
    #                    '23X_22q11_duplication_r2',
    #                    '23X_Angelman_r1',
    #                    '23Y_Cri_du_Chat_r1',
    #                    '23Y_WAGR_11p13_deletion_r2']
    highlight_files = []
    highlight_files = [i + '.1.txt' for i in highlight_files]
    files = [file for file in files if file not in highlight_files]
    files = highlight_files + files

    for file in files:
        # if True:
        # if file in ['3.txt', '39.txt', '49.txt', '12.txt', '45.txt']:
        if file == 'CMT1A_example.txt':
            filename = file.split('.')[0]
            file_path = omkar_output_dir + file
            print(file)
            mt_indexed_lists, mt_path_chrs, segment_dict, segment_size_dict = read_OMKar_to_indexed_list(file_path, forbidden_region_file)
            mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
            wt_path_dict = generate_wt_from_OMKar_output(segment_dict)
            wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
            events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)
            index_to_segment_dict = reverse_dict(segment_dict)
            if len(events) == 0:
                continue
            dependent_clusters, cluster_events = form_dependent_clusters(events, aligned_haplotypes, index_to_segment_dict)
            print(dependent_clusters)
            final_str += "\\section{{Sample Id: {}}}\n".format(filename)
            final_str += "\n"
            ## iterate over all clusters
            n_clusters = len(dependent_clusters)
            for image_cluster_idx, (c_cluster, c_events) in enumerate(zip(dependent_clusters, cluster_events)):
                ## include all homologues
                event_chr = set()
                for cluster_idx in c_cluster:
                    event_chr.add(aligned_haplotypes[cluster_idx].chrom)
                hap_idx_to_plot = []
                for hap_idx, hap in enumerate(aligned_haplotypes):
                    if hap.chrom in event_chr:
                        hap_idx_to_plot.append(hap_idx)

                # if len(hap_idx_to_plot) > 4:
                #     raise RuntimeError('more than 4 chrom selected')
                c_aligned_haplotypes = [aligned_haplotypes[i] for i in hap_idx_to_plot]

                ## generate report text
                c_events = sort_events(c_events)
                iscn_events, genes_report = format_report(c_events, aligned_haplotypes, reverse_dict(segment_dict))
                ## generate image
                c_vis_input = generate_visualizer_input(c_events, c_aligned_haplotypes, segment_dict)

                def vis_key(input_vis):
                    chr_val = input_vis['chr'][3:]
                    if chr_val == "X":
                        return_val = 23.0
                    elif chr_val == "Y":
                        return_val = 24.0
                    else:
                        return_val = float(chr_val)
                    if input_vis['highlight']:
                        return_val += 0.5  # highlight always later
                    return return_val

                c_vis_input = sorted(c_vis_input, key=vis_key)

                image_prefix = "{}/{}_imagecluster{}".format(image_dir, filename, image_cluster_idx)
                image_path = image_prefix + '_rotated.png'
                overleaf_relative_image_path = image_dir.replace('latex_reports/', '') + image_path.split('/')[-1]
                pycharm_relative_image_path = image_path
                if compile_image:
                    if len(c_vis_input) <= 4:
                        make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_VERTICAL_SPLIT)
                    else:
                        make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_HORIZONTAL_SPLIT)

                final_str += "\\subsection{{Event Cluster {} (of {})}}\n".format(image_cluster_idx + 1, n_clusters)
                final_str += "\n"
                # final_str += "\\noindent\n"
                # final_str += "\\begin{wrapfigure}{r}{0.5\\textwidth}\n"

                # final_str += "\\end{wrapfigure}\n"

                # Iterate through all events
                final_str += "\\begin{minipage}[t][][t]{0.5\\textwidth}\n"
                final_str += "\\vspace*{0pt}\n"
                final_str += "\\paragraph{SVs}\n"
                final_str += "\\medskip\n"
                final_str += "\\begin{flushleft}\n"
                SV_counter = 1
                for bullet_idx, (main_str, iscn_interpretation) in enumerate(iscn_events):
                    iscn_interpretation = hyperlink_iscn_interpretation(iscn_interpretation)
                    # format space-symbol better
                    iscn_interpretation = iscn_interpretation.replace(' on ', '&^&^&^&')
                    iscn_interpretation = iscn_interpretation.replace(' between ', '!@!@!@!@!')
                    iscn_interpretation = iscn_interpretation.replace(' ', '\\,')
                    iscn_interpretation = iscn_interpretation.replace('&^&^&^&', '\\,on ')
                    iscn_interpretation = iscn_interpretation.replace('!@!@!@!@!', ' between ')

                    # final_str += "\\item \\textbf{{{}}}.\\,{}\n".format(main_str, iscn_interpretation)
                    final_str += "\\textbf{{{}.\\;{}}}. {}\n\n".format(SV_counter, main_str, iscn_interpretation)
                    SV_counter += 1
                final_str += "\\end{flushleft}\n"
                final_str += "\n"
                final_str += "\\paragraph{Impacted genes in DDG2P}\\;\n"
                final_str += latex_gene_table(genes_report)
                final_str += "\\end{minipage}%\n"
                final_str += "\\hfill\n"
                final_str += "\\begin{minipage}[t][][t]{0.5\\textwidth}\n"
                final_str += "\\vspace*{0pt}\n"
                final_str += "\\centering\n"
                # final_str += "\\fbox{{\\includegraphics[width=\\linewidth]{{{}}}}}\n".format(overleaf_relative_image_path)
                final_str += "\\includegraphics[width=\\linewidth]{{{}}}\n".format(overleaf_relative_image_path)
                final_str += "\\captionof{{figure}}{{sample {}, event cluster {}}}\n".format(filename, image_cluster_idx + 1)
                final_str += "\\end{minipage}\n"
                final_str += "\\newpage\n\n"

    final_str += "\n"
    final_str += "\\end{document}\n"

    #         image_path = "{}/{}".format(dremsek_images, str(filename).zfill(3) + image_suffix)
    #         # image_path = "{}/{}".format(dremsek_images, str(filename) + image_suffix)
    #
    #         ## insert image
    #         if os.path.exists('latex_reports/' + image_path):
    #             # make sure the image file exists
    #             final_str += "\\section{{Sample Id: {}}}\n".format(filename)
    #             final_str += "\\begin{figure}[h!]\n"
    #             final_str += "\\centering\n"
    #             final_str += "\\includegraphics[width=3in]{{{}}}\n".format(image_path)
    #             final_str += "\\caption{{\\footnotesize Chromosomes with aberrant karyotypes}}\n"
    #             final_str += "\\label{{fig:karyotype_id{}}}\n".format(filename)
    #             final_str += "\\end{figure}\n"
    #             print('latex_reports/' + image_path, 'found')
    #         else:
    #             final_str += "\\section{{Sample Id: {}}}\n".format(filename)
    #             final_str += "\n"
    #             print('latex_reports/' + image_path, 'not found')
    #
    #         ## Iterate through all events
    #         final_str += "\\paragraph{Events}\n"
    #         final_str += "\\begin{packed_enum}\n"
    #         for bullet_idx, (main_str, iscn_interpretation) in enumerate(iscn_events):
    #             iscn_interpretation = hyperlink_iscn_interpretation(iscn_interpretation)
    #             final_str += "\\item {{\\bf {}}}. {}\n".format(main_str, iscn_interpretation)
    #         final_str += "\\end{packed_enum}\n"
    #
    #         final_str += "\n"
    #         final_str += "\\paragraph{Impacted genes in DDG2P}$\\;$\\\\\\\\\n"
    #         final_str += latex_gene_table(genes_report)
    #
    #         final_str += "\n"
    #         final_str += "\\newpage\n"
    # final_str += "\n"
    # final_str += "\\end{document}\n"
    return final_str, cases_with_events


def generate_latex_report(output_filename_prefix, front_page_str, batch_case_str):
    directory_path = os.path.dirname(output_filename_prefix) + '/'
    latex_path = output_filename_prefix + '.tex'
    with open(latex_path, 'w') as fp_write:
        fp_write.write(front_page_str)
        fp_write.write(batch_case_str)
    os.chdir('/'.join(output_filename_prefix.split('/')[:-1]))
    relative_latex_path = latex_path.split('/')[-1]
    subprocess.run(['pdflatex', relative_latex_path])


def int_file_keys(f):
    return int(f.split('.')[0])


def latex_gene_table(genes_report):
    ## format genes to report
    formatted_genes_report = format_genes_report(genes_report)

    ## form latex table
    if len(formatted_genes_report) == 0:
        return '\\\\\\quad None\n\n'
    return_str = "\\medskip\n"
    return_str += "{\\\\\\scriptsize\n"
    # return_str += "\\begin{flushleft}\n"
    return_str += "\\begin{tabular}{|llll|}\\hline\n"
    return_str += "SV & Rationale & Gene Name & Gene Omim  \\\\\\hline\n"
    for entry_idx, entry in enumerate(formatted_genes_report):
        new_line = '{SV} & {Rationale} & {Gene_Name} & {Gene_Omim} \\\\'. \
            format(SV=entry['SV'],
                   Rationale=entry['rationale'],
                   Gene_Name=entry['gene name'],
                   Gene_Omim=entry['gene omim'])
        if entry_idx == len(formatted_genes_report) - 1:
            return_str += new_line + "\\hline\n"
        else:
            return_str += new_line + "\n"
    return_str += "\\end{tabular}\n"
    return_str += "}\n"
    # return_str += "\\end{flushleft}\n"
    return return_str


def latex_hyperlink_coordinates(input_str, proximity=50000):
    return_dict = {}  # {replacement_string: hyperlinked_string}

    pattern = r'Chr(\d+|X|Y): (\d{1,3}(?:,\d{3})*)-(\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        start_pos = int(match.group(2).replace(',', ''))
        end_pos = int(match.group(3).replace(',', ''))
        ucsc_url = get_ucsc_url('chr' + match.group(1), start_pos, end_pos)
        hyperlinked_str = '\\href{{{}}}{{{}}}'.format(ucsc_url, replacement_str)
        return_dict[replacement_str] = hyperlinked_str

    pattern = r'Chr(\d+): (\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        chrom = 'chr' + match.group(1)
        pos = int(match.group(2).replace(',', ''))
        c_chr_length = get_chr_length_from_forbidden_file(chrom)
        ucsc_url = get_ucsc_url(chrom, max(0, pos - proximity), min(c_chr_length, pos + proximity))
        hyperlinked_str = '\\href{{{}}}{{{}}}'.format(ucsc_url, replacement_str)
        return_dict[replacement_str] = hyperlinked_str

    return return_dict


def hyperlink_iscn_interpretation(input_str):
    hyperlinked_mapping = latex_hyperlink_coordinates(input_str)
    for replacement_str, hyperlinked_str in hyperlinked_mapping.items():
        input_str = input_str.replace(replacement_str, hyperlinked_str)
    return input_str


def test_latex(output_name, compile_image):
    output_path = 'latex_reports/{}'.format(output_name)
    batch_case_str, cases_in_report = batch_generate_latex_case_str(data_dir, image_dir, compile_image=compile_image)
    front_str = generate_latex_frontpage('{} Data'.format(' '.join(output_name.split('_'))),
                                         'hg38 all coding genes',
                                         50,
                                         200,
                                         200)
    generate_latex_report(output_path, front_str, batch_case_str)


if __name__ == "__main__":
    forbidden_region_file = "Metadata/acrocentric_telo_cen.bed"
    # test_interpreter()
    # test_segs_union()
    # test_reciprocal_trans()
    # c_output_name, data_dir, dremsek_images = 'Dremsek', 'real_case_data/dremsek_OMKar_output_paths/', 'latex_reports/paul_dremsek_plots_new/'
    # c_output_name, data_dir, dremsek_images = 'Keyhole', 'real_case_data/keyhole_OMKar_output_paths/', 'latex_reports/keyhole_plots_new/'
    # c_output_name, data_dir, dremsek_images = 'Sunnyside', 'real_case_data/sunnyside_OMKar_output_paths/', 'latex_reports/sunnyside_plots_new/'
    c_output_name, data_dir, image_dir = 'ACC_Simulation', 'omkar_analyses_pipeline/builds/b14/omkar_paths/', 'latex_reports/ACC_simulation_plots/'
    # batch_case_str = batch_generate_latex_case_str(data_dir, dremsek_images)
    os.makedirs(image_dir, exist_ok=True)
    test_latex(c_output_name, compile_image=True)
