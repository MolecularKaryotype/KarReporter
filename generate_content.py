import os

from .Report_Genes import *
from .KT_visualizer import *
from .KarInterpreter import *
from .KarUtils import *


def format_report(interpreted_events_full, interpreted_events_partial, aligned_haplotypes, index_to_segment_dict, debug=False):
    """
    Given the output of KarInterpreter, generate the event string and genes of interest for downstream report
    @param interpreted_events_full: fully supported by SV and CNV calls
    @param interpreted_events_partial: partially supported by SV and CNV calls
    @param aligned_haplotypes:
    @param index_to_segment_dict:
    @param debug:
    @return:
    """

    def generate_parsed_report(interpreted_events, iscn_events, gene_reports, event_type_reports):
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
                    event_type_reports['deletion'] += 1
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
                    event_type_reports['inversion'] += 1
                    if bp2_band != bp3_band:
                        main_str = 'inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                    else:
                        main_str = 'inv({})({})'.format(path_chr, bp2_band)
                    chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                    iscn_interpretation = 'inversion on Chr{}: {}'.format(path_chr, chr_range)
                    bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                elif event_type == 'tandem_duplication':
                    event_type_reports['tandem_duplication'] += 1
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
                    event_type_reports['left_duplication_inversion'] += 1
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
                    event_type_reports['right_duplication_inversion'] += 1
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
                    event_type_reports['insertion'] += 1
                    # different report format if insertion is from different chr
                    if 'Chr' + path_chr == bp2_chr:
                        # TODO: check ISCN syntax if bp2_band == bp3_band
                        main_str = 'ins({})({}{}{})'.format(path_chr, bp1_band, bp2_band, bp3_band)
                    else:
                        main_str = 'ins({};{})({};{}{})'.format(path_chr, bp2_chr[3:], bp1_band, bp2_band, bp3_band)
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
                event_type_reports['reciprocal_translocation'] += 1
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
                    # raise RuntimeError('more than 2-breaks detected')
                # get breakpoints and determine the swaps by number of qter/pter
                c_event_info = event[2]
                o_event_info = o_event[2]
                event_bps = c_event_info[0].split('.')[3:5] + c_event_info[1].split('.')[3:5] + o_event_info[0].split('.')[3:5] + o_event_info[1].split('.')[
                                                                                                                                  3:5]
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
                event_type_reports['nonreciprocal_translocation'] += 1
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
                iscn_interpretation = 'balanced non-reciprocal translocation of Chr{}: {} into Chr{}: {} ({})' \
                    .format(event_seg_chr, chr_range, ins_chr, f"{ins_site_left_bp:,}", ins_site_left_band)
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

    iscn_events_full = []
    iscn_events_partial = []
    gene_reports_full = []
    gene_reports_partial = []
    event_type_reports_full = {
        "deletion":0,
        "inversion":0,
        "tandem_duplication":0,
        "left_duplication_inversion":0,
        "right_duplication_inversion":0,
        "insertion":0,
        'nonreciprocal_translocation':0,
        'reciprocal_translocation':0
    }
    event_type_reports_partial = {
        "deletion": 0,
        "inversion": 0,
        "tandem_duplication": 0,
        "left_duplication_inversion": 0,
        "right_duplication_inversion": 0,
        "insertion": 0,
        'nonreciprocal_translocation': 0,
        'reciprocal_translocation': 0
    }
    generate_parsed_report(interpreted_events_full, iscn_events_full, gene_reports_full, event_type_reports_full)
    generate_parsed_report(interpreted_events_partial, iscn_events_partial, gene_reports_partial, event_type_reports_partial)

    return iscn_events_full, iscn_events_partial, gene_reports_full, gene_reports_partial, event_type_reports_full, event_type_reports_partial


def chr_range_tostr(bpa, bpb, bpa_band, bpb_band):
    return "{}-{} ({} - {})".format(format(bpa, ',d'), format(bpb, ',d'), bpa_band, bpb_band)


def get_segment_multiplicity(mt_list):
    mult_dict = defaultdict(int)
    for mt in mt_list:
        segs = [int(x[:-1]) for x in mt]
        for seg in segs:
            mult_dict[seg] += 1
    return mult_dict


def batch_populate_html_contents(omkar_output_dir, image_dir, omkar_input_data_dir, file_of_interest=None, compile_image=False, debug=False, skip=None, forbidden_region_file=get_metadata_file_path('acrocentric_telo_cen.bed')):
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

    cluster_headers = []
    filenames = []
    n_clusters = []
    cases_with_events = []
    image1_paths = []
    image2_paths = []
    summary_image_paths = []
    summary_preview_image_paths = []
    debug_outputs = []  # list of dicts [{'segs': [], 'mt_haps': [], 'wt_haps': []}]
    folders = [directory for directory in os.listdir(omkar_output_dir)]
    bed_header = []
    bed_rows = []

    # vars with distinction between full/partially captured SVs
    iscn_reports_full, genes_reports_full, case_event_type_reports_full, DDG2P_interruptions_full, DDG2P_CNV_full, case_complexities_full = [], [], [], [], [], []
    iscn_reports_partial, genes_reports_partial, case_event_type_reports_partial, DDG2P_interruptions_partial, DDG2P_CNV_partial, case_complexities_partial = [], [], [], [], [], []

    for folder in folders:
        if folder == '.DS_Store':
            continue
        file = f"{folder}/{folder}.txt"
        bed_filepath = f"{omkar_output_dir}/{folder}/{folder}_SV.bed"

        ## parse bed file
        bed_df = pd.read_csv(bed_filepath, sep='\t')
        if not bed_header:
            bed_header = bed_df.columns.tolist()
        else:
            # check all bed has same col-headers
            c_bed_header = bed_df.columns.tolist()
            if c_bed_header != bed_header:
                raise RuntimeError('inconsistent bed headers across outputs')

        file_event_type_reports_full = []
        file_event_type_reports_partial = []
        file_DDG2P_interruptions_full = 0
        file_DDG2P_interruptions_partial = 0
        file_DDG2P_CNV_full = 0
        file_DDG2P_CNV_partial = 0
        if file_of_interest is not None:
            if file.split('.')[0] not in file_of_interest:
                continue
        if skip is not None:
            if file.split('.')[0] in skip:
                continue
        filename = file.split('/')[-1].split('.')[0]
        filenames.append(filename)
        file_path = omkar_output_dir + file
        print(file)
        mt_indexed_lists, mt_path_chrs, segment_to_index_dict, segment_size_dict = read_OMKar_to_indexed_list(file_path, forbidden_region_file)
        index_to_segment_dict = reverse_dict(segment_to_index_dict)
        mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
        wt_path_dict = generate_wt_from_OMKar_output(segment_to_index_dict)
        wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
        events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)
        if len(events) == 0:
            continue
        else:
            cases_with_events.append(filename)

        ## locate smap and cnv calls
        smap_filepath1 = f"{omkar_input_data_dir}/{filename}/exp_refineFinal1_merged_filter_inversions.smap"
        smap_filepath2 = f"{omkar_input_data_dir}/{filename}/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_filter_inversions.smap"
        cnv_filepath1 = f"{omkar_input_data_dir}/{filename}/cnv_calls_exp.txt"
        cnv_filepath2 = f"{omkar_input_data_dir}/{filename}/contigs/alignmolvref/copynumber/cnv_calls_exp.txt"
        smap_filepath = smap_filepath1 if os.path.isfile(smap_filepath1) else smap_filepath2
        cnv_filepath = cnv_filepath1 if os.path.isfile(cnv_filepath1) else cnv_filepath2

        # summary view image
        events = sort_events(events)
        # this is a stable separation, keeping the sorted order
        segment_multiplicity = get_segment_multiplicity(mt_indexed_lists)
        c_events_full, c_events_partial = separate_full_vs_partial_called_events(events, smap_filepath, cnv_filepath, index_to_segment_dict, segment_multiplicity, debug=True)
        summary_image_prefix = "{}/../full_karyotype_images/{}_cytoband_summary".format(image_dir, filename)
        image_path = summary_image_prefix + "_merged_rotated.png"
        preview_image_path = summary_image_prefix + "_merged_rotated_preview.png"
        relative_image_path = image_dir.replace('latex_reports/', '') + image_path.split('/')[-1]
        relative_preview_image_path = image_dir.replace('latex_reports/', '') + preview_image_path.split('/')[-1]
        summary_image_paths.append(relative_image_path)
        summary_preview_image_paths.append(relative_preview_image_path)
        if compile_image:
            summary_vis_input = generate_cytoband_visualizer_input(c_events_full, c_events_partial, aligned_haplotypes, segment_to_index_dict)
            summary_vis_input = sorted(summary_vis_input, key=vis_key)
            make_summary_image(filename, summary_vis_input, summary_image_prefix)

        dependent_clusters, cluster_events = form_dependent_clusters(events, aligned_haplotypes, index_to_segment_dict)
        print(dependent_clusters)
        ## iterate over all clusters
        n_cluster = len(dependent_clusters)
        n_clusters.append(n_cluster)
        for image_cluster_idx, (c_cluster, c_events) in enumerate(zip(dependent_clusters, cluster_events)):
            # to remove all later file names, check cluster_idx != 0
            
            cluster_headers.append('Chromosomal Cluster {} (of {})'.format(image_cluster_idx + 1, n_cluster))
            ## include all homologues
            event_chr = set()
            for cluster_idx in c_cluster:
                event_chr.add(aligned_haplotypes[cluster_idx].chrom)
            hap_idx_to_plot = []
            for hap_idx, hap in enumerate(aligned_haplotypes):
                if hap.chrom in event_chr:
                    hap_idx_to_plot.append(hap_idx)

            c_aligned_haplotypes = [aligned_haplotypes[i] for i in hap_idx_to_plot]

            ## filter bed table to only keep the rows with the corrsponding Chr
            chr_numbers = []
            for chr_itr in event_chr:
                if chr_itr[3:].lower() == 'x':
                    chr_numbers.append(23)
                elif chr_itr[3:].lower() == 'y':
                    chr_numbers.append(24)
                else:
                    chr_numbers.append(int(chr_itr[3:]))
            filtered_bed_df = bed_df[bed_df['chromosome'].isin(chr_numbers)]
            bed_rows.append(filtered_bed_df.values.tolist())

            ## generate report text
            c_events = sort_events(c_events)
            # this is a stable separation, keeping the sorted order
            c_events_full, c_events_partial = separate_full_vs_partial_called_events(c_events, smap_filepath, cnv_filepath, index_to_segment_dict, segment_multiplicity)
            iscn_events_full, iscn_events_partial, genes_report_full, genes_report_partial, cluster_event_type_reports_full, cluster_event_type_reports_partial = format_report(c_events_full, c_events_partial, aligned_haplotypes, index_to_segment_dict, debug=debug)

            ## generate images
            c_vis_input = generate_cytoband_visualizer_input(c_events_full, c_events_partial, c_aligned_haplotypes, segment_to_index_dict)
            c_vis_input = sorted(c_vis_input, key=vis_key)
            image_prefix = "{}/{}_cytoband_imagecluster{}".format(image_dir, filename, image_cluster_idx)
            image_path = image_prefix + '_rotated.png'
            relative_image_path = image_dir.replace('latex_reports/', '') + image_path.split('/')[-1]
            if compile_image:
                if len(c_vis_input) <= 4:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_VERTICAL_SPLIT)
                else:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_HORIZONTAL_SPLIT)
            image1_paths.append(relative_image_path)

            c_vis_input = generate_segment_visualizer_input(c_events_full, c_events_partial, c_aligned_haplotypes, segment_to_index_dict, label_centromere=True)
            c_vis_input = sorted(c_vis_input, key=vis_key)
            image_prefix = "{}/{}_segmentview_imagecluster{}".format(image_dir, filename, image_cluster_idx)
            image_path = image_prefix + '_rotated.png'
            relative_image_path = image_dir.replace('latex_reports/', '') + image_path.split('/')[-1]
            if compile_image:
                if len(c_vis_input) <= 4:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_VERTICAL_SPLIT)
                else:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix,
                               IMG_LENGTH_SCALE_HORIZONTAL_SPLIT)
            image2_paths.append(relative_image_path)

            iscn_reports_full.append(iscn_events_full)
            genes_reports_full.append(genes_report_full)
            file_event_type_reports_full.append(cluster_event_type_reports_full)
            iscn_reports_partial.append(iscn_events_partial)
            genes_reports_partial.append(genes_report_partial)
            file_event_type_reports_partial.append(cluster_event_type_reports_partial)
            ## count the number of DDG2P genes interruption in each category
            for sv_gene_report in genes_report_full:
                file_DDG2P_interruptions_full += len(sv_gene_report['bp_genes_highlight'])
                file_DDG2P_CNV_full += len(sv_gene_report['cnv_genes_highlight'])
            for sv_gene_report in genes_report_partial:
                file_DDG2P_interruptions_partial += len(sv_gene_report['bp_genes_highlight'])
                file_DDG2P_CNV_partial += len(sv_gene_report['cnv_genes_highlight'])

            ## generate debug output
            debug_segs = set()
            debug_mt_haps = []
            debug_wt_haps = []
            debug_hap_ids = []
            debug_mt_aligned = []
            debug_wt_aligned = []
            for aligned_haplotype in c_aligned_haplotypes:
                unique_segs = aligned_haplotype.unique_segment_indices()
                for seg in unique_segs:
                    seg_object = index_to_segment_dict[int(seg)]
                    debug_segs.add((seg, seg_object.chr_name, f"{seg_object.start:,}", f"{seg_object.end:,}", f"{len(seg_object):,}"))
            debug_segs = list(debug_segs)
            debug_segs = sorted(debug_segs, key=lambda x: int(x[0]))

            for c_vis in c_vis_input:
                debug_hap_ids.append(c_vis['hap_id'])
                # find the correct aligned_haplotype for mt/wt
                hap_found = False
                for aligned_haplotype in aligned_haplotypes:
                    if aligned_haplotype.id == c_vis['hap_id']:
                        debug_mt_haps.append(str(aligned_haplotype.mt_hap).replace("'", ''))
                        debug_wt_haps.append(str(aligned_haplotype.wt_hap).replace("'", ''))
                        debug_mt_aligned.append(aligned_haplotype.mt_aligned)
                        debug_wt_aligned.append(aligned_haplotype.wt_aligned)
                        hap_found = True
                        break
                if not hap_found:
                    raise RuntimeError('hap not found')
            debug_outputs.append({'segs': debug_segs, 'mt_haps': debug_mt_haps, 'wt_haps': debug_wt_haps, 'IDs': debug_hap_ids,
                                  'mt_aligned': debug_mt_aligned, 'wt_aligned': debug_wt_aligned})

        # Add the new dictionary into the case_event_type_reports
        event_multiplicity_str_full, case_complexity_full = parse_event_multiplicities(file_event_type_reports_full)
        event_multiplicity_str_partial, case_complexity_partial = parse_event_multiplicities(file_event_type_reports_partial)
        case_event_type_reports_full.append(event_multiplicity_str_full)
        case_complexities_full.append(case_complexity_full)
        case_event_type_reports_partial.append(event_multiplicity_str_partial)
        case_complexities_partial.append(case_complexity_partial)

        # DDG2P case summary
        DDG2P_interruptions_full.append(file_DDG2P_interruptions_full)
        DDG2P_CNV_full.append(file_DDG2P_CNV_full)
        DDG2P_interruptions_partial.append(file_DDG2P_interruptions_partial)
        DDG2P_CNV_partial.append(file_DDG2P_CNV_partial)
        
    return (filenames, n_clusters, cluster_headers, cases_with_events,
            image1_paths, image2_paths,
            iscn_reports_full, iscn_reports_partial, genes_reports_full, genes_reports_partial,
            case_event_type_reports_full, case_event_type_reports_partial, case_complexities_full, case_complexities_partial,
            DDG2P_interruptions_full, DDG2P_interruptions_partial, DDG2P_CNV_full, DDG2P_CNV_partial,
            summary_image_paths, summary_preview_image_paths,
            bed_header, bed_rows,
            debug_outputs)


def contains_seg_with_cnv(segs, index_to_segment_dict, segment_multiplicity, df, cnv_type, cnv_threshold=0.2, size_threshold=5000):
    """
    given a list of segments (with orientation), and a dict of full genome multiplicity, determine if exists a seg with cnv deviation from expected_cn
    :param size_threshold:
    :param cnv_threshold:
    :param cnv_type: in ['gain', 'loss']
    :param df:
    :param index_to_segment_dict:
    :param segs:
    :param segment_multiplicity:
    :return:
    """
    stripped_segs = [int(s[:-1]) for s in segs]
    cnv_size = 0
    for seg in stripped_segs:
        seg_obj = index_to_segment_dict[seg]
        chrom = convert_chrom(seg_obj.chr_name.replace('Chr', ''))
        avg_cn, expected_cn = region_reported_cn(df, chrom, seg_obj.start, seg_obj.end)
        if expected_cn == segment_multiplicity[seg]:
            # complementary event elsewhere
            continue
        if cnv_type == 'gain' and avg_cn - expected_cn < cnv_threshold:
            cnv_size += len(seg_obj)
        elif cnv_type == 'loss' and expected_cn - avg_cn < cnv_threshold:
            cnv_size += len(seg_obj)
    return True if cnv_size >= size_threshold else False

def get_sv_edge(sv_type, index_to_segment_dict, indel_segs):
    """
    :param sv_type: currently supported types in ['tandem duplication', 'left/right duplication inversion', 'deletion']
    :param index_to_segment_dict:
    :param indel_segs:
    :return:
    """
    segs = [int(s[:-1]) for s in indel_segs]
    if sv_type == 'deletion':
        chrom = index_to_segment_dict[segs[0]].chr_name
        p1 = index_to_segment_dict[segs[0]].start
        p2 = index_to_segment_dict[segs[-1]].end
        return [(sv_type, chrom, p1, chrom, p2)]
    elif sv_type == 'tandem_duplication':
        chrom = index_to_segment_dict[segs[0]].chr_name
        p1 = index_to_segment_dict[segs[0]].start
        p2 = index_to_segment_dict[segs[-1]].end
        return [(sv_type, chrom, p2, chrom, p1)]
    elif sv_type == 'left_duplication_inversion':
        chrom = index_to_segment_dict[segs[0]].chr_name
        p1 = index_to_segment_dict[segs[-1]].start
        p2 = index_to_segment_dict[segs[0]].end
        return [(sv_type, chrom, p1, chrom, p2), (sv_type, chrom, p1, chrom, p1)]
    elif sv_type == 'right_duplication_inversion':
        chrom = index_to_segment_dict[segs[0]].chr_name
        p1 = index_to_segment_dict[segs[-1]].start
        p2 = index_to_segment_dict[segs[0]].end
        return [(sv_type, chrom, p1, chrom, p2), (sv_type, chrom, p1, chrom, p1)]


def separate_full_vs_partial_called_events(sorted_events, smap_filepath, cnv_filepath, index_to_segment_dict, segment_multiplicity, debug=False):
    smap_df = smap_to_df(smap_filepath)
    cnv_df = cnv_to_df(cnv_filepath)
    full_event = []
    partial_event = []
    sv_missed_event = []  # for stats collection
    cnv_missed_event = []  # for stats collection
    sv_missed_svedge = []
    cnv_missed_svedge = []
    for e in sorted_events:
        if e[1] == 'deletion':
            deleted_segments = e[2][0].split('.')[2].replace('wt(', '').replace(')', '').split(',')
            left_segment = int(deleted_segments[0][:-1])  # strip orientation, wt has to be +
            right_segment = int(deleted_segments[-1][:-1])
            left_boundary = index_to_segment_dict[left_segment].start
            right_boundary = index_to_segment_dict[right_segment].end
            chrom = convert_chrom(index_to_segment_dict[left_segment].chr_name.replace('Chr', ''))
            # SV-missed
            smap_status = edge_in_smap(smap_df, chrom, chrom, left_boundary, right_boundary, [], 0.0, 50000)
            if not smap_status:
                partial_event.append(e)
                sv_missed_event.append(e)
                sv_missed_svedge.append(get_sv_edge(e[1], index_to_segment_dict, deleted_segments))
                continue
            # CNV-missed
            # debug use
            avg_cn, expected_cn = region_reported_cn(cnv_df, chrom, left_boundary, right_boundary)
            print(f"1223debug: {e}, {expected_cn}, {avg_cn}, {chrom}, {abs(right_boundary - left_boundary)}")
            if contains_seg_with_cnv(deleted_segments, index_to_segment_dict, segment_multiplicity, cnv_df, 'loss'):
                partial_event.append(e)
                cnv_missed_event.append(e)
                cnv_missed_svedge.append(get_sv_edge(e[1], index_to_segment_dict, deleted_segments))
                continue
        elif e[1] == 'tandem_duplication':
            duplicated_segments = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')
            left_segment = int(duplicated_segments[0][:-1])
            right_segment = int(duplicated_segments[-1][:-1])
            left_boundary = index_to_segment_dict[left_segment].start
            right_boundary = index_to_segment_dict[right_segment].end
            chrom = convert_chrom(index_to_segment_dict[left_segment].chr_name.replace('Chr', ''))
            # SV-missed
            smap_status = edge_in_smap(smap_df, chrom, chrom, left_boundary, right_boundary, [], 0.0, 50000)
            if not smap_status:
                partial_event.append(e)
                sv_missed_event.append(e)
                sv_missed_svedge.append(get_sv_edge(e[1], index_to_segment_dict, duplicated_segments))
                continue
            # CNV-missed
            avg_cn, expected_cn = region_reported_cn(cnv_df, chrom, left_boundary, right_boundary)
            print(f"1223debug: {e}, {expected_cn}, {avg_cn}")
            if contains_seg_with_cnv(duplicated_segments, index_to_segment_dict, segment_multiplicity, cnv_df, 'gain'):
                partial_event.append(e)
                cnv_missed_event.append(e)
                cnv_missed_svedge.append(get_sv_edge(e[1], index_to_segment_dict, duplicated_segments))
                continue
        elif e[1] == 'insertion':
            seg1 = e[2][0].split('.')[3]
            seg2 = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')[0]
            seg3 = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')[-1]
            seg4 = e[2][0].split('.')[4]
            duplicated_segments = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')
            if seg1 == 'p-ter':
                bp1 = -99
                chrom1 = -1
            else:
                bp1 = index_to_segment_dict[int(seg1[:-1])].end if seg1[-1] == '+' else index_to_segment_dict[int(seg1[:-1])].start
                chrom1 = convert_chrom(index_to_segment_dict[int(seg1[:-1])].chr_name.replace('Chr', ''))
            bp2 = index_to_segment_dict[int(seg2[:-1])].start if seg2[-1] == '+' else index_to_segment_dict[int(seg2[:-1])].end
            bp3 = index_to_segment_dict[int(seg3[:-1])].end if seg3[-1] == '+' else index_to_segment_dict[int(seg3[:-1])].start
            chrom2 = convert_chrom(index_to_segment_dict[int(seg2[:-1])].chr_name.replace('Chr', ''))
            chrom3 = convert_chrom(index_to_segment_dict[int(seg3[:-1])].chr_name.replace('Chr', ''))
            if seg4 == 'q-ter':
                bp4 = -99
                chrom4 = -1
            else:
                bp4 = index_to_segment_dict[int(seg4[:-1])].start if seg4[-1] == '+' else index_to_segment_dict[int(seg4[:-1])].end
                chrom4 = convert_chrom(index_to_segment_dict[int(seg4[:-1])].chr_name.replace('Chr', ''))
            # SV-missed
            # if terminal, we will require one fewer edge
            smap_status12 = edge_in_smap(smap_df, chrom1, chrom2, bp1, bp2, [], 0.0, 50000) if bp1 != -99 else True
            smap_status34 = edge_in_smap(smap_df, chrom3, chrom4, bp3, bp4, [], 0.0, 50000) if bp4 != -99 else True
            if (not smap_status12) or (not smap_status34):
                partial_event.append(e)
                sv_missed_event.append(e)
                svedge = []
                if chrom1 != -1:
                    svedge.append(('insertion', inverse_convert_chrom(chrom1), bp1, inverse_convert_chrom(chrom2), bp2))
                if chrom4 != -1:
                    svedge.append(('insertion', inverse_convert_chrom(chrom3), bp3, inverse_convert_chrom(chrom4), bp4))
                sv_missed_svedge.append(svedge)
                continue
            # CNV-missed
            avg_cn, expected_cn = region_reported_cn(cnv_df, chrom2, bp2, bp3)
            print(f"1223debug: {e}, {expected_cn}, {avg_cn}")
            if contains_seg_with_cnv(duplicated_segments, index_to_segment_dict, segment_multiplicity, cnv_df, 'gain'):
                partial_event.append(e)
                cnv_missed_event.append(e)
                svedge = []
                if chrom1 != -1:
                    svedge.append(('insertion', inverse_convert_chrom(chrom1), bp1, inverse_convert_chrom(chrom2), bp2))
                if chrom4 != -1:
                    svedge.append(('insertion', inverse_convert_chrom(chrom3), bp3, inverse_convert_chrom(chrom4), bp4))
                cnv_missed_svedge.append(svedge)
                continue
        elif e[1] == 'left_duplication_inversion':
            seg1 = e[2][0].split('.')[3]
            seg2 = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')[0]
            seg3 = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')[-1]
            seg4 = e[2][0].split('.')[4]
            duplicated_segments = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')
            if seg1 == 'p-ter':
                bp1 = -99
            else:
                bp1 = index_to_segment_dict[int(seg1[:-1])].end if seg1[-1] == '+' else index_to_segment_dict[int(seg1[:-1])].start
            bp2 = index_to_segment_dict[int(seg2[:-1])].start if seg2[-1] == '+' else index_to_segment_dict[int(seg2[:-1])].end
            bp3 = index_to_segment_dict[int(seg3[:-1])].end if seg3[-1] == '+' else index_to_segment_dict[int(seg3[:-1])].start
            bp4 = bp3  # self-edge
            chrom = convert_chrom(index_to_segment_dict[int(seg2[:-1])].chr_name.replace('Chr', ''))  # should be strcitly intra-chr
            # SV-missed
            smap_status12 = edge_in_smap(smap_df, chrom, chrom, bp1, bp2, [], 0.0, 50000) if bp1 != -99 else True
            smap_status34 = edge_in_smap(smap_df, chrom, chrom, bp3, bp4, [], 0.0, 50000)
            if (not smap_status12) or (not smap_status34):
                partial_event.append(e)
                sv_missed_event.append(e)
                sv_missed_svedge.append(get_sv_edge(e[1], index_to_segment_dict, duplicated_segments))
                continue
            if contains_seg_with_cnv(duplicated_segments, index_to_segment_dict, segment_multiplicity, cnv_df, 'gain'):
                partial_event.append(e)
                cnv_missed_event.append(e)
                cnv_missed_svedge.append(get_sv_edge(e[1], index_to_segment_dict, duplicated_segments))
                continue
        elif e[1] == 'right_duplication_inversion':
            seg1 = e[2][0].split('.')[3]
            seg2 = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')[0]
            seg3 = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')[-1]
            seg4 = e[2][0].split('.')[4]
            duplicated_segments = e[2][0].split('.')[2].replace('mt(', '').replace(')', '').split(',')
            bp1 = index_to_segment_dict[int(seg1[:-1])].end if seg1[-1] == '+' else index_to_segment_dict[int(seg1[:-1])].start
            bp2 = bp1  # self-edge
            bp3 = index_to_segment_dict[int(seg3[:-1])].end if seg3[-1] == '+' else index_to_segment_dict[int(seg3[:-1])].start
            if seg4 == 'q-ter':
                bp4 = -99
            else:
                bp4 = convert_chrom(index_to_segment_dict[int(seg4[:-1])].chr_name.replace('Chr', ''))
            chrom = convert_chrom(index_to_segment_dict[int(seg2[:-1])].chr_name.replace('Chr', ''))  # should be strictly intra-chr
            # SV-missed
            smap_status12 = edge_in_smap(smap_df, chrom, chrom, bp1, bp2, [], 0.0, 50000)
            smap_status34 = edge_in_smap(smap_df, chrom, chrom, bp3, bp4, [], 0.0, 50000) if bp4 != -99 else True
            if (not smap_status12) or (not smap_status34):
                partial_event.append(e)
                sv_missed_event.append(e)
                sv_missed_svedge.append(get_sv_edge(e[1], index_to_segment_dict, duplicated_segments))
                continue
            if contains_seg_with_cnv(duplicated_segments, index_to_segment_dict, segment_multiplicity, cnv_df, 'gain'):
                partial_event.append(e)
                cnv_missed_event.append(e)
                cnv_missed_svedge.append(get_sv_edge(e[1], index_to_segment_dict, duplicated_segments))
                continue
        full_event.append(e)
    if debug:
        print('1221 test')
        if full_event:
            out_dict = {'count': len(full_event), 'events': full_event}
            print(f"full_event: {out_dict}")
        if sv_missed_event:
            out_dict = {'count': len(sv_missed_event), 'events': sv_missed_event}
            print(f"sv_missed_event: {out_dict}")
            print(f"sv_missed_svedge: {sv_missed_svedge}")
        if cnv_missed_event:
            out_dict = {'count': len(cnv_missed_event), 'events': cnv_missed_event}
            print(f"cnv_missed_event: {out_dict}")
            print(f"cnv_missed_svedge: {cnv_missed_svedge}")
    return full_event, partial_event


def parse_event_multiplicities(dict_list):
    combined_dict = {}

    for d in dict_list:
        for key, value in d.items():
            key = key.replace('left_duplication_inversion', 'duplication_inversion')
            key = key.replace('right_duplication_inversion', 'duplication_inversion')
            if key not in combined_dict:
                combined_dict[key] = value
            else:
                combined_dict[key] += value

    # sort dictionary in correct event order, and return complexity
    event_order = ['reciprocal_translocation',
                   'nonreciprocal_translocation',
                   'duplication_inversion',
                   'inversion',
                   'duplicated_insertion',
                   'tandem_duplication',
                   'deletion']
    complexity_mapping = {'reciprocal_translocation': 2,
                          'nonreciprocal_translocation': 3,
                          'duplication_inversion': 2,
                          'inversion': 2,
                          'duplicated_insertion': 2,
                          'tandem_duplication': 1,
                          'deletion': 1}

    ## we are displaying insertion as duplicated insertion
    combined_dict['duplicated_insertion'] = combined_dict['insertion']
    combined_dict.pop('insertion')

    ## tally complexity
    total_complexity = 0
    for event in event_order:
        if event in combined_dict and combined_dict[event] > 0:
            total_complexity += complexity_mapping[event] * combined_dict[event]

    return combined_dict, total_complexity


def get_ucsc_url(chrom, start_pos, end_pos, db='hg38'):
    prefix = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db={}' \
             '&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position='.format(db)
    suffix = '{}%3A{}%2D{}'.format(chrom.lower(), start_pos, end_pos)
    return prefix + suffix


def test_ucsc_url():
    x = get_ucsc_url('Chr2', 130271298, 131386307)
    print(x)
