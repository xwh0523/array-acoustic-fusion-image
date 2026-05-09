# User guide

This guide describes the input data, parameters, outputs, and expected behavior of the two software modules included in this repository.

The repository contains two graphical software modules:

1. `slowness_envelope_reconstruction_app.py`
2. `fusion_image_framework_app.py`

Both programs are designed for testing and reproducing the computational workflow described in the associated manuscript.

---

## 1. Lithology-constrained slowness-envelope reconstruction module

### Run the software

From the repository root directory, run:

python src/slowness_envelope_reconstruction_app.py

### Purpose

This module reconstructs floating baselines for three acoustic slowness curves:

- DTP: P-wave slowness curve
- DTS: S-wave slowness curve
- DTST: Stoneley-wave slowness curve

It then calculates the corresponding slowness-envelope attributes:

- S_DTP
- S_DTS
- S_DTST

The method is designed to reduce lithology-related false anomalies by incorporating lithology constraints during the baseline reconstruction process.

### Input format

The minimum required input format is:

Depth, DTP, DTS, DTST, Lithology

Additional profile columns can be included after the lithology column:

Depth, DTP, DTS, DTST, Lithology, Profile_1, Profile_2, Profile_3, ...

The profile columns are optional and are only used for visualization in the profile-track panel.

### Lithology code

The default lithology code used in the example dataset is:

1 = Sandstone  
2 = Limestone  
3 = Coal  
4 = Mudstone  

These codes are only the default setting used in the example dataset. Users may define their own lithology codes and lithology classes according to their datasets. The target lithology code used for baseline reconstruction can be modified in the software interface, and the number of lithology classes is not restricted to four.

### Main options

#### Sliding window length

This option controls the number of depth samples used for local-minimum detection.

A value of 3 means that a three-point moving window is used. The center point is selected as a local minimum if it is smaller than the neighboring points.

#### Enable lithology constraint

If this option is enabled, samples outside the target lithology are forced to act as interpolation anchors. This allows the floating baseline to follow lithology-related background changes and reduces false slowness-envelope anomalies.

If this option is disabled, local minima are detected without lithology constraints.

#### Clip negative values to zero

If this option is enabled, negative envelope values are set to zero.

#### Show interpolation anchors

If this option is enabled, the selected interpolation anchors are displayed on the slowness tracks.

#### Display thresholds

The threshold controls only the visualization of the slowness-envelope tracks. It does not change the calculated envelope values.

### Outputs

The exported result file contains:

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

where:

DTP_BL = reconstructed floating baseline of the P-wave slowness curve  
DTS_BL = reconstructed floating baseline of the S-wave slowness curve  
DTST_BL = reconstructed floating baseline of the Stoneley-wave slowness curve  
S_DTP = P-wave slowness-envelope attribute  
S_DTS = S-wave slowness-envelope attribute  
S_DTST = Stoneley-wave slowness-envelope attribute  

### Expected behavior

After loading a valid input file and clicking `Run`, the software displays:

- P-wave slowness curve and floating baseline
- S-wave slowness curve and floating baseline
- Stoneley-wave slowness curve and floating baseline
- P-wave slowness-envelope attribute
- S-wave slowness-envelope attribute
- Stoneley-wave slowness-envelope attribute
- FDA lithology track
- Optional profile curves

The depth axis increases downward, following the standard well-log display convention.

---

## 2. Multi-attribute fusion-image framework module

### Run the software

From the repository root directory, run:

python src/fusion_image_framework_app.py

### Purpose

This module converts multiple acoustic attributes into a continuous multi-attribute fusion image. The fusion image is generated using inverse-distance weighting interpolation and Gaussian filtering. A sliding-window integration method is then used to extract a composite fracture indicator curve.

### Input format

The minimum required input format is:

Depth, Attribute_1

The default manuscript-related input format is:

Depth, S_DTP, 1-AMP, S_DTS, CAI, S_DTST

where:

Depth = logging depth  
S_DTP = P-wave slowness-envelope attribute  
1-AMP = transformed average amplitude response  
S_DTS = S-wave slowness-envelope attribute  
CAI = coupled attenuation index  
S_DTST = Stoneley-wave slowness-envelope attribute  

Additional attribute columns can be added. If a header row is included, the software automatically reads the attribute names from the first line.

### Main options

#### Interpolated columns

This option controls the number of interpolated columns along the attribute axis.

The value must be greater than or equal to the number of input attributes.

#### IDW power

This option controls the power parameter used in inverse-distance weighting interpolation.

A larger value gives greater weight to nearby attributes along the attribute axis.

#### Gaussian sigma

This option controls the standard deviation of the Gaussian filter used for image smoothing.

A larger value produces a smoother fusion image.

#### Window size

This option controls the sliding-window length used for image-volume integration.

#### Step

This option controls the sliding step in samples.

#### Colormap

This option controls the colormap used for visualization.

#### Clip negative values to zero

If this option is enabled, negative values in the fusion image and composite indicator curve are set to zero.

### Outputs

The exported result file contains:

Depth  
LF  
Img_01  
Img_02  
...  
Img_N  

where:

LF = composite fracture indicator curve  
Img_01 to Img_N = interpolated and smoothed fusion-image columns  

### Expected behavior

After loading a valid input file and clicking `Run`, the software displays:

- A multi-attribute fusion image
- A composite fracture indicator curve

The depth axis increases downward, following the standard well-log display convention.

---

## Notes

The original field data used in the manuscript cannot be publicly released due to confidentiality agreements with the data provider. The datasets included in this repository are provided for software testing, demonstration, and reproducibility of the computational workflow.