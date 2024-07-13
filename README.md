# KarReporter

Given a Molecular Karyotype file, generate an HTML report
containing 
1) ISCN Nomenclature events
2) Genes of interest
3) Visualization of the event

## Development

**Dependencies**: please do not modify. If there are bugs, please let me know, thank you!
1) [KarInterpreter](https://github.com/MolecularKaryotype/KarInterpreter)
2) [KarUtils](https://github.com/MolecularKaryotype/KarUtils)

**KT_interpreter_html_report.py** and **template.html**
- for HTML formatting
- most work needed are within these two files

**generate_content.py**
- upstream for HTML formatting
- this code creates all the information needed in the HTML report
- for more complicated tasks, you may need to include additional input for the KT_interpreter_html_report.py.
However, the use of this file requires additional knowledge such as how the Interpreter algorithm works, so we
should talk before you dive into this file.

**KT_visualizer.py**
- code used to create current visualization image
- eventually, we are looking to retiring the image generation. However, we can retain the content generation to
use them directly in the HTML formatting.
