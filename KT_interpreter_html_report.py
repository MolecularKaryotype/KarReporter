from jinja2 import Environment, FileSystemLoader
import os
import base64
import re

from KarInterpreter import *
from format_report import *
from KarUtils import *

def image_to_base64(image_path):
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except FileNotFoundError:
        print(f"Error: File {image_path} not found.")
        return ""


def int_file_keys(filename):
    return int(filename.split('_')[0])


def html_hyperlink_coordinates(input_str, proximity=50000):
    return_dict = {}  # {replacement_string: hyperlinked_string}

    pattern = r'Chr(\d+|X|Y): (\d{1,3}(?:,\d{3})*)-(\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        start_pos = int(match.group(2).replace(',', ''))
        end_pos = int(match.group(3).replace(',', ''))
        ucsc_url = get_ucsc_url('chr' + match.group(1), start_pos, end_pos)
        hyperlinked_str = f'<a href="{ucsc_url}">{replacement_str}</a>'
        return_dict[replacement_str] = hyperlinked_str

    pattern = r'Chr(\d+): (\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        chrom = 'chr' + match.group(1)
        pos = int(match.group(2).replace(',', ''))
        c_chr_length = get_chr_length_from_forbidden_file(chrom)
        ucsc_url = get_ucsc_url(chrom, max(0, pos - proximity), min(c_chr_length, pos + proximity))
        hyperlinked_str = f'<a href="{ucsc_url}">{replacement_str}</a>'
        return_dict[replacement_str] = hyperlinked_str

    return return_dict


def hyperlink_iscn_interpretation(input_str):
    hyperlinked_mapping = html_hyperlink_coordinates(input_str)
    for replacement_str, hyperlinked_str in hyperlinked_mapping.items():
        input_str = input_str.replace(replacement_str, hyperlinked_str)
    return input_str


def test(compile_image, cases_of_interest, title, omkar_output_dir, image_output_dir, output_file, debug=False, skip=None):
    os.makedirs(image_output_dir, exist_ok=True)

    # one tuple per cluster (event)
    headers, cases_with_events, image_paths, iscn_reports, genes_reports, debug_outputs = batch_populate_contents(omkar_output_dir, image_output_dir,
                                                                                          file_of_interest=cases_of_interest, compile_image=compile_image, debug=debug, skip=skip)
    images_base64 = [image_to_base64(img) for img in image_paths]

    formatted_genes_reports = [format_genes_report(genes_report) for genes_report in genes_reports]
    columns_order = ['SV', 'rationale', 'gene name', 'gene omim']

    ## hyperlinking
    for iscn_report in iscn_reports:
        for iscn_report_idx, (iscn, sv_interpretation) in enumerate(iscn_report):
            hyperlinked_sv_interpretation = hyperlink_iscn_interpretation(sv_interpretation)
            iscn_report[iscn_report_idx][1] = hyperlinked_sv_interpretation

    content = [(header, text, image, table_content, debug) for header, text, image, table_content, debug in zip(headers, iscn_reports, images_base64, formatted_genes_reports, debug_outputs)]

    env = Environment(loader=FileSystemLoader('html_reports/'))
    template = env.get_template('template.html')
    rendered_html = template.render(title=title, content=content, columns_order=columns_order, debug=debug)

    with open(output_file, 'w') as f:
        f.write(rendered_html)
    print(f"HTML file generated: {os.path.abspath(output_file)}")


def manual_test():
    # Define the data
    title = "My Text and Images"
    texts = ['text1',
            'text2',
            'text3']

    ## ZJ: image paths need to be relative path
    images = ['/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster0_rotated.png',
              '/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster1_rotated.png',
              '/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster2_rotated.png']
    images_base64 = []
    for img in images:
        images_base64.append(image_to_base64(img))


    table_contents = [
        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]],

        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]],

        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]]
    ]

    content = [(text, image, table_content) for text, image, table_content in zip(texts, images_base64, table_contents)]

    # Create an environment for Jinja2
    env = Environment(loader=FileSystemLoader('html_reports/'))
    template = env.get_template('template.html')

    # Render the template with the data
    rendered_html = template.render(title=title, content=content)

    # Write the rendered HTML to a file
    output_file = 'html_reports/test.html'
    with open(output_file, 'w') as f:
        f.write(rendered_html)

    print(f"HTML file generated: {os.path.abspath(output_file)}")


if __name__ == "__main__":
    forbidden_region_file = "Metadata/acrocentric_telo_cen.bed"
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'keyhole', 'real_case_data/keyhole_OMKar_output_paths/', 'html_reports/keyhole_plots/', 'html_reports/keyhole.html'
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'sunnyside', 'real_case_data/sunnyside_OMKar_output_paths/', 'html_reports/sunnyside_plots/', 'html_reports/sunnyside.html'
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'dremsek', 'real_case_data/dremsek_OMKar_output_paths/', 'html_reports/dremsek_plots/', 'html_reports/dremsek.html'
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'Dremsek_b14', 'real_case_data/dremsek_b14_paths/', 'html_reports/dremsek_plots_b14/', 'html_reports/dremsek_b14.html'
    i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'karsim', 'omkar_analyses_pipeline/builds/b14/omkar_paths/', 'html_reports/karsim_plots/', 'html_reports/karsim.html'
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'test', 'omkar_analyses_pipeline/builds/b14/omkar_paths/', 'html_reports/test_plots/', 'html_reports/test3.html'
    # test(True, ['23X_22q11_duplication_r2.1.txt'], i_title, i_omkar_output_dir, i_image_output_dir, i_output_file, debug=True)
    # test(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_file, debug=True, skip=['23X_Xp11_22_Microduplication_r2'])
    test(True, ['23X_22q11_duplication_r2'], i_title, i_omkar_output_dir, i_image_output_dir, i_output_file, debug=True)
    # test(True, ['2280'], i_title, i_omkar_output_dir, i_image_output_dir, 'html_reports/2280.html', debug=True)
    #
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'sunnyside', 'real_case_data/sunnyside_OMKar_output_paths/', 'html_reports/sunnyside_plots/', 'html_reports/sunnyside.html'
    # test(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_file, debug=True)
