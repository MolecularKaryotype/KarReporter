from KT_interpreter_html_report import *
from KarUtils import *


# batch_post_process_OMKar_output('example_input/subset_test/', '0903_test_files/post_processed_omkar/')
# batch_post_process_OMKar_output('0903_test_files/input/', '0903_test_files/post_processed_omkar/')



# i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'Example Cases', \
#                                                                 'example_input/manuscript_files/', \
#                                                                 'example_output/tmp/', \
#                                                                 'example_output/example_cases/'
# i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'Test Cases', \
#                                                                 '0903_test_files/post_processed_omkar/', \
#                                                                 '0903_test_files/images/', \
#                                                                 '0903_test_files/html_output/'
# i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'Infertility Screening', \
#                                                                 'Dremsek/omkar_output/', \
#                                                                 'Dremsek/images/', \
#                                                                 'Dremsek/output/'

# generate_html_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir, debug=False)

### Keyhole
# i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'Keyhole', \
#                                                                 'v2_reports/keyhole_post_processed/', \
#                                                                 'v2_reports/keyhole_report_full/images/', \
#                                                                 'v2_reports/keyhole_report_full/report/'
# batch_post_process_OMKar_output('D:\\keyhole_0717_paths\\', i_omkar_output_dir)
# generate_html_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir, debug=False)

### Sunnyside
# i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'Sunnyside', \
#                                                                 'v2_reports/sunnyside_post_processed/', \
#                                                                 'v2_reports/sunnyside_report_full/images/', \
#                                                                 'v2_reports/sunnyside_report_full/report/'
# batch_post_process_OMKar_output('D:\\sunnyside_0717_paths\\', i_omkar_output_dir)
# generate_html_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir, debug=False)

### Dremsek
# i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'Dremsek', \
#                                                                 'v2_reports/dremsek_post_processed/', \
#                                                                 'v2_reports/dremsek_report_full/images/', \
#                                                                 'v2_reports/dremsek_report_full/report/'
# batch_post_process_OMKar_output('D:\\dremsek_OMKar_output_paths\\', i_omkar_output_dir)
# generate_html_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir, debug=False)

### Simulation
# i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'KarSim', \
#                                                                 'example_input/b17_karsim_postprocessed/', \
#                                                                 'b17_karsim_report/images/', \
#                                                                 'b17_karsim_report/report/'
# batch_post_process_OMKar_output('example_input/b17_karsim/', i_omkar_output_dir)
# generate_html_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir, debug=True)


### Manuscript Draft
i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir = 'Manuscript', \
                                                                '0903_test_files/manuscript/input/', \
                                                                '0903_test_files/manuscript/images/', \
                                                                '0903_test_files/manuscript/report/'
generate_html_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_dir, debug=True)
