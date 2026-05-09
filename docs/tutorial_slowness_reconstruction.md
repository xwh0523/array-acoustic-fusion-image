# Tutorial: lithology-constrained slowness-envelope reconstruction

This tutorial demonstrates how to use the slowness-envelope reconstruction software.

## 1. Start the software

From the repository root directory, run:

python src/slowness_envelope_reconstruction_app.py

The graphical user interface will open.

## 2. Select the input file

Click `Select File`.

Select the dataset:

Data-Interpolation Reconstruction Algorithm for Slowness Curves

The input file should contain at least the following columns:

Depth, DTP, DTS, DTST, Lithology

Additional profile columns can be included after the lithology column.

## 3. Read the data

Click `Read Data`.

The software will load the input file and display the data in the preview table.

The software will also report:

- number of rows
- number of columns
- depth range
- detected profile columns

## 4. Check the default segment

After the data are loaded, a default segment covering the full depth range is created automatically.

The user can modify:

- top depth
- bottom depth
- sliding window length
- lithology constraint option
- negative-value clipping option
- display thresholds

## 5. Run the reconstruction

Click `Run`.

The software will reconstruct floating baselines for:

- P-wave slowness curve
- S-wave slowness curve
- Stoneley-wave slowness curve

The software then calculates:

- S_DTP
- S_DTS
- S_DTST

## 6. Interpret the output plot

The output plot contains:

1. P-wave slowness curve and floating baseline
2. S-wave slowness curve and floating baseline
3. Stoneley-wave slowness curve and floating baseline
4. P-wave slowness-envelope attribute
5. S-wave slowness-envelope attribute
6. Stoneley-wave slowness-envelope attribute
7. FDA lithology track
8. Optional profile curves

The depth axis increases downward.

## 7. Export results

Click `Export Results`.

The exported file contains:

Depth  
Segment_ID  
DTP  
DTP_BL  
S_DTP  
DTS  
DTS_BL  
S_DTS  
DTST  
DTST_BL  
S_DTST  
Lithology  
Optional profile columns  

## 8. Export figure

Click `Export Figure`.

The figure can be saved as:

- PNG
- PDF
- TIFF

## 9. Notes

The target lithology code can be modified in the software interface. The default target lithology is sandstone, but users may define their own lithology codes and lithology classes according to their datasets.