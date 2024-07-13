import os

from Report_Genes import *
from KT_visualizer import *
from KarInterpreter import *
from KarUtils import *


def format_report(interpreted_events, aligned_haplotypes, index_to_segment_dict, debug=False):
    iscn_events = []
    gene_reports = []
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
                if bp2_band != bp3_band:
                    main_str = 'inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'inversion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
            elif event_type == 'tandem_duplication':
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
                # different report format if insertion is from different chr
                if 'Chr' + path_chr == bp2_chr:
                    # TODO: check ISCN syntax if bp2_band == bp3_band
                    main_str = 'ins({})({}{}{})'.format(path_chr, bp1_band, bp2_band, bp3_band)
                else:
                    main_str = 'ins({};{})({};{}{})'.format(path_chr, bp2_chr, bp1_band, bp2_band, bp3_band)
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
                raise RuntimeError('more than 2-breaks detected')
            # get breakpoints and determine the swaps by number of qter/pter
            c_event_info = event[2]
            o_event_info = o_event[2]
            event_bps = c_event_info[0].split('.')[3:5] + c_event_info[1].split('.')[3:5] + o_event_info[0].split('.')[3:5] + o_event_info[1].split('.')[3:5]
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
            iscn_interpretation = 'balanced non-reciprocal translocation of Chr{}: {} into Chr{}: {}({})' \
                .format(event_seg_chr, chr_range, ins_chr, ins_site_left_bp, ins_site_left_band)
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
    return iscn_events, gene_reports


def chr_range_tostr(bpa, bpb, bpa_band, bpb_band):
    return "{}-{} ({} - {})".format(format(bpa, ',d'), format(bpb, ',d'), bpa_band, bpb_band)


def batch_populate_contents(omkar_output_dir, image_dir, file_of_interest=None, compile_image=False, debug=False, skip=None, forbidden_region_file='KarUtils/Metadata/acrocentric_telo_cen.bed'):
    headers = []
    cases_with_events = []
    image_paths = []
    iscn_reports = []
    genes_reports = []
    debug_outputs = []  # list of dicts [{'segs': [], 'mt_haps': [], 'wt_haps': []}]
    files = [file for file in os.listdir(omkar_output_dir)]
    # files = sorted(files, key=int_file_keys)
    for file in files:
        if file_of_interest is not None:
            if file.split('.')[0] not in file_of_interest:
                continue
        if skip is not None:
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
        if len(events) == 0:
            continue
        else:
            cases_with_events.append(filename)
        dependent_clusters, cluster_events = form_dependent_clusters(events, aligned_haplotypes, index_to_segment_dict)
        print(dependent_clusters)
        ## iterate over all clusters
        n_clusters = len(dependent_clusters)
        for image_cluster_idx, (c_cluster, c_events) in enumerate(zip(dependent_clusters, cluster_events)):
            # to remove all later file names, check cluster_idx != 0
            headers.append('{}: cluster {} (out of {})'.format(filename, image_cluster_idx + 1, n_clusters))
            ## include all homologues
            event_chr = set()
            for cluster_idx in c_cluster:
                event_chr.add(aligned_haplotypes[cluster_idx].chrom)
            hap_idx_to_plot = []
            for hap_idx, hap in enumerate(aligned_haplotypes):
                if hap.chrom in event_chr:
                    hap_idx_to_plot.append(hap_idx)

            c_aligned_haplotypes = [aligned_haplotypes[i] for i in hap_idx_to_plot]

            ## generate report text
            c_events = sort_events(c_events)
            iscn_events, genes_report = format_report(c_events, aligned_haplotypes, index_to_segment_dict, debug=debug)
            ## generate image
            c_vis_input = generate_visualizer_input(c_events, c_aligned_haplotypes, segment_to_index_dict)

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
            relative_image_path = image_dir.replace('latex_reports/', '') + image_path.split('/')[-1]
            if compile_image:
                if len(c_vis_input) <= 4:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_VERTICAL_SPLIT)
                else:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_HORIZONTAL_SPLIT)

            image_paths.append(relative_image_path)
            iscn_reports.append(iscn_events)
            genes_reports.append(genes_report)

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
                        debug_mt_haps.append(aligned_haplotype.mt_hap)
                        debug_wt_haps.append(aligned_haplotype.wt_hap)
                        debug_mt_aligned.append(aligned_haplotype.mt_aligned)
                        debug_wt_aligned.append(aligned_haplotype.wt_aligned)
                        hap_found = True
                        break
                if not hap_found:
                    raise RuntimeError('hap not found')
            debug_outputs.append({'segs': debug_segs, 'mt_haps': debug_mt_haps, 'wt_haps': debug_wt_haps, 'IDs': debug_hap_ids,
                                  'mt_aligned': debug_mt_aligned, 'wt_aligned': debug_wt_aligned})

    return headers, cases_with_events, image_paths, iscn_reports, genes_reports, debug_outputs



def get_ucsc_url(chrom, start_pos, end_pos, db='hg38'):
    prefix = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db={}' \
             '&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position='.format(db)
    suffix = '{}%3A{}%2D{}'.format(chrom.lower(), start_pos, end_pos)
    return prefix + suffix


def test_ucsc_url():
    x = get_ucsc_url('Chr2', 130271298, 131386307)
    print(x)
