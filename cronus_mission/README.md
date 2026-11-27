# Cronus mission (Day 2)

This project holds the Day 2 rescue-mission blocks: the **bulk mission**
(mission_cronus.R, the cleaned March 2026 script, nine acts) and **single cells
from the crew** (about 35 minutes, guided). How the two are combined on Day 2
is not decided yet.

## Bulk mission

Open cronus_mission.Rproj, then open mission_cronus.R and run it act by act.
It finds Rlib like the single-cell block, needs tidyverse, patchwork, ggrepel,
limma and RColorBrewer, and installs nothing. Data are in data/bulk; the story
texts and two tables (Cooper's status records, the W-Y gene key) are encrypted
in transmissions/ and open with clearance codes released by Mission Control.
Students write their T3 groups in R; no CSV editing.

Instructor material is in instructor/mission: INSTRUCTOR_GUIDE.md (codes, release
rules, expected findings, what changed against March), clearance_codes.csv,
plaintexts/ and tables/ (editable sources), make_transmissions.R (re-encrypts
after edits), and the March originals in original_scripts/.

## Single cells from the crew

Open mission_single_cells.html for the student version. In R:

1. Fully extract the archive into your persistent university account folder,
   beside AVA_training and single_cell_workshop.
2. Open cronus_mission.Rproj.
3. Create a script and follow the HTML step by step.
   `source("start_single_cells.R")` loads the packages and builds the object `cells`.

start_single_cells.R finds Rlib inside the project, beside it, one folder
further up, or in Documents/alien-lecture. It needs Seurat, tidyverse, and patchwork, and installs nothing.

If a step stalls, two checkpoints let students continue:
checkpoints/cronus_cells_clustered.rds (after quality control, clustering, and
UMAP) and checkpoints/cronus_cells_annotated.rds (with cell-type names).

The data in data/single_cell are **simulated for the course**: tissue biopsies
at T3 from seven crew members whose bulk samples exist in the mission archive
(two controls, two antidote recipients, one creature, Cooper, and the survivor).
Marker genes have real human names; HVAR genes use dashes (HVAR-041) because
Seurat replaces underscores, and hvar_modules.csv gives the archive IDs
(HVAR_041).

## Instructor material

- mission_single_cells_solutions.Rmd/.html: all steps run, with worked
  interpretations and a note on what is built into the data.
- instructor/simulate_cronus_cells.R: the simulator (expression means and
  negative binomial counts only; no clusters or coordinates are supplied).
- instructor/make_course_data.R: writes data/single_cell (seed 2183) and both
  checkpoints by running the student Rmd's own code chunks.
- instructor/validate_cronus_cells.R: runs a blind Seurat workflow on three
  seeds and checks every designed finding against the hidden truth;
  results and figures are in instructor/feasibility.
- instructor/simulation_truth.csv: hidden cell types, states, damaged cells,
  and doublets. Never give this to students.

Story text in the two AVA transmissions is marked "edit freely". Shared text
must be edited in both Rmds. After changing the simulator, run
make_course_data.R, then knit both Rmds, then course_review/package_course_versions.py.

Validated locally with R 4.4.3 and Seurat 5.4.0; the university setup
(R 4.5.3, Seurat 6.0.0) has not been tested with this block yet.
