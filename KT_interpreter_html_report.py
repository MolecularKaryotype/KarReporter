from random import sample

from jinja2 import Environment, FileSystemLoader
import base64
import os
import shutil
from generate_content import *
from KarUtils import *
import copy


def image_to_base64(image_path):
    """
    Convert an image file into a base64 embedding
    """
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except FileNotFoundError:
        print(f"Error: File {image_path} not found.")
        return ""


def int_file_keys(filename):
    return int(filename.split('_')[0])


def html_hyperlink_coordinates(input_str, proximity=50000):
    """
    given an input string, which may contain a genomic coordinate substring, generate a mapping for each genomic
    coordinate to the UCSC Genome Browser Hyperlink
    @param input_str:
    @param proximity: If the genomic coordinate is a single position (instead of two, forming a range),
                      we generate a range of position +/- proximity
    @return: mapping dict
    """
    return_dict = {}  # {replacement_string: hyperlinked_string}

    pattern = r'Chr(\d+|X|Y): (\d{1,3}(?:,\d{3})*)-(\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        start_pos = int(match.group(2).replace(',', ''))
        end_pos = int(match.group(3).replace(',', ''))
        ucsc_url = get_ucsc_url('chr' + match.group(1), start_pos, end_pos)
        hyperlinked_str = f'<a href="{ucsc_url}" target="_blank">{replacement_str}</a>'
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


def generate_html_report(compile_image, cases_of_interest, title, data_dir, image_output_dir, output_dir, debug=False, skip=None):
    """

    @param compile_image: bool, whether the images were already compiled/not needing update
    @param cases_of_interest: list, if None, report all cases; otherwise, report only cases in the list
    @param title: title of the report
    @param data_dir:
    @param image_output_dir: where to store the image
    @param output_dir: path to the output files
    @param debug: bool, whether to include debug info in the report
    @param skip: list, cases to skip
    @return: N/a
    """
    os.makedirs(image_output_dir, exist_ok=True)

    # one tuple per cluster (event), where the output of batch_populate_html_contents is all the clusters (from all case files)
    (filenames, clusters, headers, cases_with_events,
     image1_paths, image2_paths,
     iscn_reports, genes_reports,
     case_event_type_reports, case_complexities, DDG2P_interruptions, DDG2P_CNV, summary_image_paths, summary_preview_image_paths,
     bed_header, bed_rows,
     debug_outputs) = batch_populate_html_contents(data_dir, image_output_dir, file_of_interest=cases_of_interest, compile_image=compile_image, debug=debug, skip=skip)

    ## summarizing number of events
    summary_event_counts = {}
    for sample_str in case_event_type_reports:
        parsed_str = copy.deepcopy(sample_str)
        parsed_str = parsed_str.replace('<b>', '').replace('</b>', '').split(', ')
        if not sample_str:
            # case without interpreted event
            continue
        for info in parsed_str:
            info = info.split(': ')
            event_type = info[0]
            event_count = int(info[1])
            if event_type not in summary_event_counts:
                summary_event_counts[event_type] = event_count
            else:
                summary_event_counts[event_type] += event_count
    print({key: summary_event_counts[key] for key in sorted(summary_event_counts)})

    images1_base64 = [image_to_base64(img) for img in image1_paths]
    images2_base64 = [image_to_base64(img) for img in image2_paths]
    os.makedirs(output_dir + "/karyotypes/", exist_ok=True)
    for img in summary_image_paths:
        shutil.copy(img, output_dir + "/karyotypes/")
    for img in summary_preview_image_paths:
        shutil.copy(img, output_dir + "/karyotypes/")
    summary_image_names = [img.split('/')[-1] for img in summary_image_paths]
    summary_preview_image_names = [img.split('/')[-1] for img in summary_preview_image_paths]
    with open("bootstrap/static/assets/pics/magnifying_glass_icon_reflected.txt") as fp_read:
        magnifying_glass_icon = fp_read.readline().strip()

    formatted_genes_reports = [format_genes_report(genes_report) for genes_report in genes_reports]
    columns_order = ['SV', 'gene name', 'gene omim', 'rationale']

    ## hyperlinking
    for iscn_report in iscn_reports:
        for iscn_report_idx, (iscn, sv_interpretation) in enumerate(iscn_report):
            hyperlinked_sv_interpretation = hyperlink_iscn_interpretation(sv_interpretation)
            iscn_report[iscn_report_idx][1] = hyperlinked_sv_interpretation

    ## for old report
    # content = [(header, text, image, table_content, debug_info) for header, text, image, table_content, debug_info in
    #            zip(headers, iscn_reports, images1_base64, formatted_genes_reports, debug_outputs)]
    # env = Environment(loader=FileSystemLoader('./'))
    # template = env.get_template('template.html')
    # rendered_html = template.render(title=title, content=content, columns_order=columns_order, debug=debug)
    dashboard = [(filename, cluster, case_event_type_report, case_complexity, DDG2P_interruptions, DDG2P_CNV, summary_image, summary_preview_image)
                 for filename, cluster, case_event_type_report, case_complexity, DDG2P_interruptions, DDG2P_CNV, summary_image, summary_preview_image in
                 zip(filenames, clusters, case_event_type_reports, case_complexities, DDG2P_interruptions, DDG2P_CNV, summary_image_names, summary_preview_image_names)]
    env1 = Environment(loader=FileSystemLoader('./bootstrap'))
    newtemplate = env1.get_template('dashboard.html')
    newrendered_html = newtemplate.render(title=title, content=dashboard, debug=debug)

    #create output folder if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir,exist_ok=True)
        #Copy the static folder into output folder and overwrite
    if os.path.exists(output_dir+"/static"):
        shutil.rmtree(output_dir+"/static")
    shutil.copytree("bootstrap/static", output_dir+"/static", dirs_exist_ok=True)
    ###

    with open(output_dir+"/"+"dashboard.html", 'w') as f:
        f.write(newrendered_html)

    # with open(output_dir+"/oldreport.html", 'w') as f:
    #     f.write(rendered_html)

    start = 0
    index = 0
    reporttemplate = env1.get_template('report.html')
    # os.makedirs(output_dir+"/cases_html", exist_ok=True)
    for cluster in clusters:
        filtered_headers = headers[start:start+cluster]
        filtered_images1 = images1_base64[start:start+cluster]
        filtered_images2 = images2_base64[start:start+cluster]
        filtered_iscn = iscn_reports[start:start+cluster]
        filtered_gene_reports = formatted_genes_reports[start:start+cluster]
        filtered_debug = debug_outputs[start:start+cluster]
        filtered_bed_rows = bed_rows[start:start+cluster]
        start = start+cluster
        report_title = filenames[index]
        index+=1

        # this is the content for each cluster
        filtered_content = [(header, image1, image2, iscn, gene_report, debug_info, bed_row) for header, image1, image2, iscn, gene_report, debug_info, bed_row in
               zip(filtered_headers, filtered_images1, filtered_images2, filtered_iscn, filtered_gene_reports, filtered_debug, filtered_bed_rows)]
        filtered_report = reporttemplate.render(title=report_title, content=filtered_content, columns_order=columns_order, debug=debug, mag_icon=magnifying_glass_icon, bed_header=bed_header)
        with open(f"{output_dir}/{report_title}.html", 'w') as f:
            f.write(filtered_report)

    print(f"HTML file generated")


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
    output_file = 'html_reports/generate_html_report.html'
    with open(output_file, 'w') as f:
        f.write(rendered_html)

    print(f"HTML file generated: {os.path.abspath(output_file)}")


if __name__ == "__main__":
    forbidden_region_file = "KarUtils/Metadata/acrocentric_telo_cen.bed"
    i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'example_run', 'example_input/b17_karsim/', 'output/plots/', 'output'
    generate_html_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir, debug=True)
