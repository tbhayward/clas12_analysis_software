# Review and full analysis notes

Build either reporting profile from the v2 Overleaf archive and a completed production run. The extraction CSV and final canvases are inputs. Event processing is not repeated.

```sh
make note PROFILE=review
make note PROFILE=full
```

The defaults use `~/Desktop/dvcs_cross_section` and `~/Desktop/DVCS_RGA_pass2_analysis_note_v2.zip`. Override them for another snapshot:

```sh
make note PROFILE=review RUN="/path/to/completed/dvcs_cross_section" NOTE="/path/to/note.zip"
```

`NOTE_DEST` defaults to `output/note_build`. Each profile produces a complete clean Overleaf directory (`overleaf_review` or `overleaf_full`), its matching ZIP and a compiled note with stable names. Temporary compilation files are removed automatically. Python requires NumPy, pandas, Matplotlib, SciPy and Gepard. The review build evaluates dense KM15/BH curves for the world-data overlays using the existing analysis helper. Compilation requires `pdflatex`, `bibtex` and the packages used by the note.

The equivalent direct command is:

```sh
python3 external_scripts/build_note.py review --run "/path/to/completed/dvcs_cross_section" --note "/path/to/note.zip" --dest output/note_build
```

## One existing Overleaf project

```sh
make note PROFILE=both
```

This builds both profiles and produces `DVCS-single-project-additions.zip` plus the `single_project_additions/note_profiles` folder. Upload that folder into your existing Overleaf project. Follow its `README.txt` for five small edits to your existing `main.tex`. Keep `main.tex` as the main document. No new Overleaf project is required.

Choose the version in `note_profiles/selection.tex` with `\ReviewNotetrue` or `\ReviewNotefalse`, then recompile. Switching in Overleaf does not require Python or a new analysis run. Shared prose is kept in one source. `\NoteChoice{full text}{review text}` holds the differing passages, numerical tables and figure blocks. Only differing images are added. The original project files remain available for comparison.

The overlay is based on the supplied v2 archive and the approved changes. Incorporate any later edits to the affected chapters into their shared replacements. Other chapters and the bibliography remain in their existing locations. `main_single_project_example.tex` is supplied separately as a reference for the five edits, not as a file that must replace your current main document.

## Reporting coverage

| Profile | Cross section at both energies | BSA at both energies |
| --- | --- | --- |
| review | Nominal bins 03–06, 0.25–1.00 GeV² | All six nominal bins |
| full | All six nominal bins | All six nominal bins |

Membership follows nominal CSV edges, not measured means or imported reporting coordinates. Both profiles copy final-result PNGs from the input run without changing their values or uncertainties. The two lowest cross-section canvases are absent from the review ZIP. All BSA canvases are retained. The inactive proton-energy-loss chapter and its figures remain archived, with its input commented out.

The review build preserves the v2 section structure and discussion. It makes the approved numerical substitutions in the abstract, systematic-uncertainty summary and global-comparison chapter. All seven global-comparison subsections, all original labels and all eleven graphic inclusions remain. The same physical low, intermediate and higher xB world-data canvases are retained, with the reported t intervals. Pass-1, world-data and model comparisons are recomputed on the selected pass-2 sample. The Jo–Defurne cross-checks and the native pass-1 benchmark are also restricted to the relevant momentum-transfer range.

Final cross-section coverage and the associated opening reference change in review. BSA retains all six intervals and uses the input-run figures. Its result text identifies statistical error bars. The lower-t study status appears briefly in the review abstract and results. The original calibration and extraction studies retain their production populations and corrections. Their source chapters are not rewritten.

Edit the baseline Overleaf archive to revise prose or stable figures. The builder targets the approved v2 passages and fails if required source passages are absent or ambiguous. When editing those passages, update their replacements deliberately. The compiled source directories are generated outputs. Preserve edits in the baseline archive before rebuilding. The full profile retains the baseline full-sample abstract and comparison discussion, with the input-run final figures and the common BSA result clarification.

Both builds load `booktabs` and remove an unused repeated subfigure label so that the baseline compiles without missing table commands or duplicate labels. Projection figures retain the production layout, colors and markers. Distinct phi intervals remain separate even when their centers coincide. The full-sample projection medians are checked against the production plotting CSVs before the review selection is applied.

## Validation and provenance

Before packaging, the builder checks nominal bin IDs for duplicates, reproduces both uncertainty summary tables, verifies final figures against the input run and checks the expected figure counts. Review comparison inputs are checked against nominal cross sections before refitting. The four-period consistency diagnostic is recomputed over the selected sample without changing the production uncertainty allocations. LaTeX compilation must finish without unresolved references, citations or multiply defined labels. Outputs replace the previous files only after successful compilation.

Each source directory and ZIP includes `build-info.json` with the source paths and checksums, selected cross-section IDs, uncertainty support counts, final-figure checksums, review comparison IDs and fit diagnostics. Undefined fit diagnostics are JSON null. The full profile preserves the baseline full-sample comparison chapter. It does not regenerate every historical diagnostic from the run. Rebuilding either profile requires the complete input run, including its comparison tables and plotting helpers.

For the restored November 17 run, regression values are:

| Quantity | Full | Review |
| --- | ---: | ---: |
| Positive 10.6 GeV cross sections | 1,445 | 1,141 |
| Mean statistical uncertainty | 9.60% | 10.40% |
| Mean point-to-point uncertainty | 9.05% | 9.50% |
| Positive 10.2 GeV cross sections | 1,276 | 1,020 |
| Mean statistical uncertainty | 10.97% | 11.79% |
| Mean point-to-point uncertainty | 8.77% | 9.17% |
| Matched 10.6 GeV pass-1 comparisons | 1,289 | 997 |
| Median pass-2/pass-1 ratio | 0.877 | 0.927 |

The twelve cross-section canvases and six 10.2 GeV BSA canvases in this run match the original v2 archive exactly. The six 10.6 GeV BSA canvases differ in roughly 0.10–0.22% of pixels. Over 1,400 isolated markers match in position, with Desktop vertical bars approximately 11–12% shorter at the median. This is an image comparison, not a numerical reconstruction of the archived BSA extraction. The input-run figures are used in both profiles.
