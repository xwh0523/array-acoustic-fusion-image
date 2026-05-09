# Example and test data

This folder contains example or desensitized datasets used to test the software modules provided in this repository.

The original field logging, core, electrical imaging, and production-test datasets used in the manuscript are not publicly available because of confidentiality agreements with the data provider. The files provided here are intended to demonstrate the input format, computational workflow, and expected software behavior.

## Files

### Data-Interpolation Reconstruction Algorithm for Slowness Curves

This dataset is used as the input dataset for the lithology-constrained slowness-envelope reconstruction module.

The corresponding software module is:

python src/slowness_envelope_reconstruction_app.py

The minimum required input format is:

Depth, DTP, DTS, DTST, Lithology

Additional profile columns can be included after the lithology column. These profile columns are optional and are used only for visualization in the profile-track panel.

The recommended input format is:

Depth, DTP, DTS, DTST, Lithology, Profile_1, Profile_2, Profile_3, ...

where:

Depth      = logging depth  
DTP        = P-wave slowness curve  
DTS        = S-wave slowness curve  
DTST       = Stoneley-wave slowness curve  
Lithology  = FDA-based lithology classification result  
Profile_i  = optional profile curve used for visualization  

The default lithology code is defined as:

1 = Sandstone  
2 = Limestone  
3 = Coal  
4 = Mudstone  

These codes are only the default setting used in the example dataset. Users may define their own lithology codes and lithology classes according to their datasets. The target lithology code used for baseline reconstruction can be modified in the software interface, and the number of lithology classes is not restricted to four.

In the reconstruction algorithm, sandstone intervals are treated as the target reservoir lithology. Local minima are detected in sandstone intervals and used as interpolation anchors. Non-sandstone samples are forced to act as interpolation anchors when the lithology constraint is enabled. This reduces false slowness-envelope anomalies caused by lithological variations.

The software outputs the reconstructed floating baselines and slowness-envelope attributes:

DTP_BL  
S_DTP  
DTS_BL  
S_DTS  
DTST_BL  
S_DTST  

where:

DTP_BL   = reconstructed floating baseline of the P-wave slowness curve  
S_DTP    = P-wave slowness-envelope attribute  
DTS_BL   = reconstructed floating baseline of the S-wave slowness curve  
S_DTS    = S-wave slowness-envelope attribute  
DTST_BL  = reconstructed floating baseline of the Stoneley-wave slowness curve  
S_DTST   = Stoneley-wave slowness-envelope attribute  

---

### Attribute dataset

This dataset is used as the input dataset for the multi-attribute fusion-image framework module.

The corresponding software module is:

python src/fusion_image_framework_app.py

The minimum required input format is:

Depth, Attribute_1

The default manuscript-related input format is:

Depth, S_DTP, 1-AMP, S_DTS, CAI, S_DTST

where:

Depth   = logging depth  
S_DTP   = P-wave slowness-envelope attribute  
1-AMP   = transformed average amplitude response  
S_DTS   = S-wave slowness-envelope attribute  
CAI     = coupled attenuation index  
S_DTST  = Stoneley-wave slowness-envelope attribute  

Additional attribute columns can be added after S_DTST. If the input file contains a header row, the software automatically reads the attribute names from the first line. If no header row is provided, default attribute names are assigned automatically.

The software converts the input attribute matrix into a continuous fusion image using inverse-distance weighting interpolation and Gaussian filtering. It then calculates a composite fracture indicator curve through sliding-window image-volume integration.

The exported result contains:

Depth  
LF  
Img_01  
Img_02  
...  
Img_N  

where:

LF      = composite fracture indicator curve  
Img_N   = interpolated and smoothed fusion-image column  

---

## Notes on data availability

The datasets in this folder are provided only for software testing, demonstration, and reproducibility of the computational workflow. They should not be interpreted as the complete original field dataset used in the manuscript.

The original field logging, core, electrical imaging, and production-test data cannot be publicly released due to confidentiality agreements with the data provider. Users may replace the example datasets with their own data as long as the input column order follows the formats described above.