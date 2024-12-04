from jinja2 import Environment, FileSystemLoader
import base64
import os
import shutil
import copy


from .generate_content import *
from .KarUtils import *


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


def generate_html_report(compile_image, cases_of_interest, title, data_dir, image_output_dir, output_dir, omkar_input_data_dir, debug=False, skip=None):
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
    @param omkar_input_data_dir:
    """
    os.makedirs(image_output_dir, exist_ok=True)
    os.makedirs(output_dir + "/full_karyotype_images/", exist_ok=True)
    kr_dir = os.path.dirname(os.path.abspath(__file__))

    # contains sample-level and cluster-level info
    # each is a list, split either by sample or by cluster
    # sample-level is for report summary, and cluster-level is for each subpage in the report
    (filenames, n_clusters, cluster_headers, cases_with_events,
        image1_paths, image2_paths,
        iscn_reports_full, iscn_reports_partial, genes_reports_full, genes_reports_partial,
        case_event_type_reports_full, case_event_type_reports_partial, case_complexities_full, case_complexities_partial,
        DDG2P_interruptions_full, DDG2P_interruptions_partial, DDG2P_CNV_full, DDG2P_CNV_partial,
        summary_image_paths, summary_preview_image_paths,
        bed_header, bed_rows,
        debug_outputs) = batch_populate_html_contents(data_dir, image_output_dir, omkar_input_data_dir, file_of_interest=cases_of_interest, compile_image=compile_image, debug=debug, skip=skip)

    ## content for report summary
    # format number of events, also log number of events for summary count
    event_order = ['reciprocal_translocation',
                   'nonreciprocal_translocation',
                   'duplication_inversion',
                   'inversion',
                   'duplicated_insertion',
                   'tandem_duplication',
                   'deletion']
    case_event_type_report_strs = []
    summary_event_counts = {}
    for (event_tally_full, event_tally_partial) in zip(case_event_type_reports_full, case_event_type_reports_partial):
        c_str = []
        for e in event_order:
            full_count = 0 if not e in event_tally_full else event_tally_full[e]
            partial_count = 0 if not e in event_tally_partial else event_tally_partial[e]
            total_count = full_count + partial_count
            if total_count != 0:
                c_str.append(f"<b>{e.replace('_', ' ')}: {total_count}</b>&#8202;({partial_count})")
                if e not in summary_event_counts:
                    summary_event_counts[e] = total_count
                else:
                    summary_event_counts[e] += total_count
        case_event_type_report_strs.append(', '.join(c_str))
    print({key: summary_event_counts[key] for key in sorted(summary_event_counts)})

    # case complexities
    case_complexity_strs = []
    for (complexity_full, complexity_partial) in zip(case_complexities_full, case_complexities_partial):
        c_str = f"{complexity_full + complexity_partial}&#8202;({complexity_partial})"
        case_complexity_strs.append(c_str)

    # DDG2P interruptions
    DDG2P_interruption_strs = []
    for (interruption_full, interruption_partial) in zip(DDG2P_interruptions_full, DDG2P_interruptions_partial):
        c_str = f"{interruption_full + interruption_partial}&#8202;({interruption_partial})"
        DDG2P_interruption_strs.append(c_str)

    # DDG2P CNV
    DDG2P_cnv_strs = []
    for (cnv_full, cnv_partial) in zip(DDG2P_CNV_full, DDG2P_CNV_partial):
        c_str = f"{cnv_full + cnv_partial}&#8202;({cnv_partial})"
        DDG2P_cnv_strs.append(c_str)

    # summary images
    summary_image_names = [img.split('/')[-1] for img in summary_image_paths]
    summary_preview_image_names = [img.split('/')[-1] for img in summary_preview_image_paths]
    with open(f"{kr_dir}/bootstrap/static/assets/pics/magnifying_glass_icon_reflected.txt") as fp_read:
        magnifying_glass_icon = fp_read.readline().strip()

    summary_report = [(filename, cluster, case_event_type_report, case_complexity, DDG2P_interruptions, DDG2P_CNV, summary_image, summary_preview_image)
                      for filename, cluster, case_event_type_report, case_complexity, DDG2P_interruptions, DDG2P_CNV, summary_image, summary_preview_image in
                      zip(filenames, n_clusters, case_event_type_report_strs, case_complexity_strs, DDG2P_interruption_strs, DDG2P_cnv_strs, summary_image_names, summary_preview_image_names)]
    env1 = Environment(loader=FileSystemLoader(f'{kr_dir}/bootstrap'))
    newtemplate = env1.get_template('report_summary.html')
    newrendered_html = newtemplate.render(title=title, content=summary_report, debug=debug)

    ## copying html formatting assets
    if not os.path.exists(output_dir):
        os.makedirs(output_dir,exist_ok=True)
    if os.path.exists(output_dir+"/static"):
        shutil.rmtree(output_dir+"/static")

    ## replaced the dirs_exist_ok=True; which does not exist before py3.8
    source = f"{kr_dir}/bootstrap/static"
    destination = output_dir + "/static"
    # Check if destination exists
    if os.path.exists(destination):
        # Copy the contents of the source directory into the destination
        for item in os.listdir(source):
            src_path = os.path.join(source, item)
            dst_path = os.path.join(destination, item)

            if os.path.isdir(src_path):
                # Copy directory recursively
                shutil.copytree(src_path, dst_path)
            else:
                # Copy individual files
                shutil.copy2(src_path, dst_path)
    else:
        # If the destination doesn't exist, use copytree as usual
        shutil.copytree(source, destination)

    with open(output_dir+"/"+"report_summary.html", 'w') as f:
        f.write(newrendered_html)

    ## content for each cluster subpage
    # hyperlinking
    for ir in iscn_reports_full:
        for iscn_report_idx, (iscn, sv_interpretation) in enumerate(ir):
            hyperlinked_sv_interpretation = hyperlink_iscn_interpretation(sv_interpretation)
            ir[iscn_report_idx][1] = hyperlinked_sv_interpretation
    for ir in iscn_reports_partial:
        for iscn_report_idx, (iscn, sv_interpretation) in enumerate(ir):
            hyperlinked_sv_interpretation = hyperlink_iscn_interpretation(sv_interpretation)
            ir[iscn_report_idx][1] = hyperlinked_sv_interpretation

    # cluster images
    image1_names = [img.split('/')[-1] for img in image1_paths]
    image2_names = [img.split('/')[-1] for img in image2_paths]

    # genes report
    formatted_genes_reports_full = [format_genes_report(genes_report) for genes_report in genes_reports_full]
    formatted_genes_reports_partial = [format_genes_report(genes_report) for genes_report in genes_reports_partial]
    columns_order = ['SV', 'gene', 'rationale', 'allelic req.', 'mutatation req.', 'confidence', 'organ']
    start = 0
    index = 0
    reporttemplate = env1.get_template('report.html')
    # os.makedirs(output_dir+"/cases_html", exist_ok=True)
    for n_cluster in n_clusters:
        ## each sample page is generated with n_cluster subpages, where each input data is spliced for the n_cluster's info
        filtered_headers = cluster_headers[start:start+n_cluster]
        filtered_images1 = image1_names[start:start+n_cluster]
        filtered_images2 = image2_names[start:start+n_cluster]
        filtered_iscn_full = iscn_reports_full[start:start+n_cluster]
        filtered_iscn_partial = iscn_reports_partial[start:start + n_cluster]
        filtered_gene_reports_full = formatted_genes_reports_full[start:start+n_cluster]
        filtered_gene_reports_partial = formatted_genes_reports_partial[start:start + n_cluster]
        filtered_debug = debug_outputs[start:start+n_cluster]
        filtered_bed_rows = bed_rows[start:start+n_cluster]
        start = start+n_cluster
        report_title = filenames[index]
        index += 1

        # this is the content for each cluster
        filtered_content = [(header, image1, image2, iscn_full, iscn_partial, gene_report_full, gene_report_partial, debug_info, bed_row) for header, image1, image2, iscn_full, iscn_partial, gene_report_full, gene_report_partial, debug_info, bed_row in
               zip(filtered_headers, filtered_images1, filtered_images2, filtered_iscn_full, filtered_iscn_partial, filtered_gene_reports_full, filtered_gene_reports_partial, filtered_debug, filtered_bed_rows)]
        filtered_report = reporttemplate.render(title=report_title, content=filtered_content, columns_order=columns_order, debug=debug, mag_icon=magnifying_glass_icon, bed_header=bed_header)
        with open(f"{output_dir}/{report_title}.html", 'w') as f:
            f.write(filtered_report)

    print(f"HTML file generated")

if __name__ == "__main__":
    forbidden_region_file = "KarUtils/Metadata/acrocentric_telo_cen.bed"
    i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'example_run', 'example_input/b17_karsim/', 'output/plots/', 'output'
    generate_html_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir, debug=True)
