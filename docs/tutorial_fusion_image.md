# Tutorial: multi-attribute fusion-image framework

This tutorial demonstrates how to use the multi-attribute fusion-image software.

## 1. Start the software

From the repository root directory, run:

python src/fusion_image_framework_app.py

The graphical user interface will open.

## 2. Select the input file

Click `Select Input File`.

Select the dataset:

Attribute dataset

The minimum required input format is:

Depth, Attribute_1

The default manuscript-related input format is:

Depth, S_DTP, 1-AMP, S_DTS, CAI, S_DTST

Additional attribute columns can be added. If the first row contains column names, the software reads the attribute names automatically.

## 3. Check detected attributes

After loading the file, the software reports:

- number of samples
- number of attributes
- attribute order

The number of input attributes is not fixed. The software automatically uses all columns after the depth column as input attributes.

## 4. Set parameters

The main parameters are:

### Interpolated columns

This controls the number of columns generated along the attribute axis.

The value must be greater than or equal to the number of input attributes.

### IDW power

This controls the power parameter of inverse-distance weighting interpolation.

A larger value gives more influence to nearby attributes.

### Gaussian sigma

This controls the smoothing strength of the Gaussian filter.

A larger value produces a smoother fusion image.

### Window size

This controls the sliding-window length used for image-volume integration.

### Step

This controls the sliding step in samples.

### Colormap

This controls the colormap used for visualization.

### Clip negative values to zero

If this option is enabled, negative values are set to zero.

## 5. Run the fusion-image workflow

Click `Run`.

The software will:

1. read the attribute matrix
2. interpolate the attribute axis using inverse-distance weighting
3. smooth the image using Gaussian filtering
4. generate a multi-attribute fusion image
5. calculate the composite fracture indicator curve using sliding-window integration

## 6. Interpret the output plot

The output plot contains:

1. Multi-attribute fusion image
2. Composite fracture indicator curve

The depth axis increases downward.

## 7. Export results

Click `Export Results`.

The exported file contains:

Depth  
LF  
Img_01  
Img_02  
...  
Img_N  

where:

LF = composite fracture indicator curve  
Img_01 to Img_N = interpolated and smoothed fusion-image columns  

## 8. Save figure

Click `Save Figure`.

The figure can be saved as:

- PNG
- PDF
- TIFF

## 9. Notes

The default manuscript-related attributes are:

S_DTP  
1-AMP  
S_DTS  
CAI  
S_DTST  

However, users can add more attributes to the input file. The software automatically reads all attribute columns after the depth column and uses them for fusion-image generation.