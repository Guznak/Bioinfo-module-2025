# AVA / Recovery training

Open AVA_training.html for the student exercises. It contains the full teaching
explanations, figures, boxes, examples, and tasks. It can be read offline.

## Work in R

1. Fully extract the archive into your persistent course folder, for example
   Documents/alien-lecture. Windows "Extract all" suggests a folder named after
   the ZIP; change it to the course folder so AVA_training sits beside Rlib.
2. Open AVA_training.Rproj.
3. Open student_workbook.R and save your work there.
4. Run source("check_setup.R"), then follow the HTML in order.

The check finds Rlib inside the project, beside it, one folder further up, or in
Documents/alien-lecture.
It checks packages and inputs without installing or updating anything.
Core packages are tidyverse, limma, and knitr.

Input files stay in data/bulk, data/ml, and data/annotation. Save results in
output. The main HTML explains working directories and paths under
“Where are we?”.

If your code is stuck, catch_up/after_exercise_3.R and
catch_up/after_exercise_9.R recreate the objects needed for the next step.
extra_practice holds two older R/tidyverse courses for optional practice after
Exercise 0; they load some files from a university server and need internet.

data/annotation/gene_annotation.csv (gene symbol, description, biotype) was
built from org.Hs.eg.db 3.20.0 and EnsDb.Hsapiens.v86 by
instructor/build_gene_annotation.R; IDs without a match are left empty.

## Separate solutions version

The instructor bundle additionally contains AVA_training_solutions.Rmd,
AVA_training_solutions.html, and reference_analysis.R. The solution HTML displays
all worked code, plots, tables, and interpretation notes without collapsed answers.

Both Rmd sources are ordinary editable documents. There is no visibility switch.
Knit the source you want to rebuild in RStudio. Shared explanations need editing
in both files. Rendering requires rmarkdown, knitr, the analysis packages, and
Pandoc, normally supplied with RStudio.

## Single-cell continuation

Extract single_cell_workshop beside this chapter's folder. The link from the
main HTML then opens single_cell_workshop/single_cell_workshop.html.
Open that chapter's own Rproj before running its scripts.
It has one route of about 60 minutes, starting from its student_workbook.R.

The course uses real transformed cancer data and an explicitly simulated
logistic-regression archive (300 people, instructor/generate_archive.R). The separate single-cell workshop uses real brain data.
The package library is supplied separately and must match the classroom setup.
