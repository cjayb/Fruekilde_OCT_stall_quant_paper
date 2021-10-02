# Capillary stall quantification from optical coherence tomography angiogram maximum intensity projections

By: Signe K. Fruekilde, Eugenio G. Jiménez, Kim R. Drasbek,
Christopher J. Bailey

This repository accompanies
([our bioRxiv preprint](https://doi.org/10.1101/2021.10.01.461840))
documenting the need and implementation of a method to obtain measures of
capillary blood flow stall events in pre-clinical optical coherence
tomography (OCT) scans.

Check out our OSF project page at
[https://osf.io/me6hf/](https://osf.io/me6hf/) for more information about
the method, as well as *access to data used to demonstrate it*.

## Basic workflow

1. Clone the repository and add it to your Matlab path
2. Run `a_create_first_angio_files.m`. Select a single `nii`-format file
containing a frame of data, and identify the depth at which the 3D images
are most in-focus.
3. Provide the output from the previous step to
`b_create_other_angios_destripe_frangi.m`, which will apply image
cleaning and vessel enhancements operations, as well as flattening the 3D
slabs into 2D maximum intensity-projections (MIPs). Note that this is
applied to all frames in the same folder as the one chosen above.
NB! The function `clean_angio_stripes_and_emph_vessels.m` contains fixed
parameter values that should be optimised for new imaging setups.
4. Select the folder with cleaned 2D MIPs, and let
`c_identify_cap_segments.m` automatically extract capillary segments from
from all the frames, applying correction of in-plane translation.
5. Run `d_set_stall_thresholds.m`.
Based on the image intensity of each capillary segment over time,
manually define intensity thresholds indicating the occurrence of _stalls_,
i.e., vessel portions in which blood flood is temporarily blocked.
6. Run `e_viz_stallogram_overlays.m` to extract 'stall statistics', such as
portion of extracted capillaries with at least one stall, the longest stall
of each capillary, etc. 

The automatic segmentation of capillary segments is based on three
important resources:

* Münch, B., Trtik, P., Marone, F., and Stampanoni, M. (2009).
Stripe and ring artifact removal with combined wavelet—Fourier filtering,
Opt. Express 17, 8567-8591 (2009);
[click here for more information.](https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-17-10-8567&id=179485)
* Kroon, D.-J. (2021). Hessian based Frangi Vesselness filter,
MATLAB Central File Exchange;
[click here for more information.](https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter)
* The `edgelink`-algorithm by Peter Kovesi; 
[click here for more information.](https://www.peterkovesi.com/matlabfns/)

## Detailed description of core scripts

### `NiiReader.m`

Calls `NII2RR.m` for conversion of reflexivity profile to image; see
parameters therein. Note that the easiest way to use our methods on data
from other systems is to implement a new reader class. The reader should
expose methods/attributes for

* `readFrame`: returns a single 2D (XZ) slice through the data.
* `zrange`: specify the start and stop indices in the Z-direction for slab-
extraction

Basic usage:

```matlab
    % >>> nii = NiiReader(fullfile(pathname, filename));
    % >>> frame = nii.readFrame(150);  % read 150th XZ-plane
    % >>> imagesc(frame)  % already converted to reflexivity
```

### `calculate_angio.m`

Creates 3D angiograms as the subtraction of two B-scans. Calls
`CorrSlicePhase.m`; see parameter choices therein and

> Srinivasan, V. J., Jiang, J. Y., Yaseen, M. A., Radhakrishnan, H.,
Wu, W., Barry, S., ... Boas, D. A. (2010). Rapid volumetric angiography of
cortical microvasculature with optical coherence tomography.
Optics Letters, 35(1), 43. https://doi.org/10.1364/ol.35.000043

### `RemoveStripesVertical.m`

Applies wavelet-filtering to remove stripe-artefacts, as described in
[Münch et al. (2009)](https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-17-10-8567&id=179485)

### `a_create_first_angio_files.m`

This is the beginning of the analysis pipeline described in Fruekilde et
al. (in preparation). The objective is first to reduce the amount of data
to be read by identifying the focus depth of the acquired data. Depending
on the lense used, our experience is that a stack of only around 20-30
slices is in focus, rendering most of the 1,024 slices obsolete. To
reduce the memory and disk space consumption, we here define the focus
'slab', and save the information to disk. The next stage of processing
will read the saved 'angiogram' and apply the same ROI selection to all
frames. Select a 'representative frame'.

### `b_create_other_angios_destripe_frangi.m`

After creating the reduced-depth angiogram slab in step a), we here apply
the selection to all 3D data frames, and write out the reduced angiograms
for all. This step is time-consuming, but can be performed with as little
as 8 GB of RAM on a laptop in a matter of minutes (60 frames).

The reduced-depth slabs are flattened using a maximum-intensity
projection (MIP): `log(max( ... , 3))`, where the dimension of the
max-operation is along the third (depth) dimension. Each frame is
projected independently, cleaned, and saved to disk as a 2D TIFF image
(both cleaned and original versions are saved).

Note that the script calls `clean_angio_stripes_and_emph_vessels.m`
This is where the parameters for 'destriping' and vessel emphasis using
the Frangi filter are defined (and should be adjusted manually).

Calls `clean_angio_stripes_and_emph_vessels.m`.

### `c_identify_cap_segments.m`

Previous analysis stages have resulted in n_frames TIFF images of cleaned
Ps. This script reads in the folder of TIFFs, and begins by calculating
an average across frames. It then identifies the frame which 'resembles'
the average the most (similarity is quantified as the cosine of the
'angle' between images after flattening each into a 1-dimensional vector)
Subsequently, each frame is coregistered to the 'most representative'
single frame using only translations (no rotations or shears), and a new
average is computed. This average is then skeletonised, from which edges
are extracted using the edgelink-algorithm. Each one of these edges is
dilated using a 3x3 image kernel to create n_edges 'capillary masks'.
The image intensity within each masks across time frames thus depends on
whether blood flows or is stalled. The 'stallogram' is a matrix of
dimensions (n_frames, n_edges), i.e., a stacking of the extracted mask
intensity time courses.

This script does not require any user interaction, and runs in less than
30 seconds on a standard modern laptop. The output is a figure showing
the final edge map of the image, and the stallogram. We recommend that
each image is inspected visually to ensure that a reasonable-looking
capillary network is extracted, with on the order of 120-180 edges for an
image dimension of 400 x 400 pixels. Adjusting the image margin to avoid
filter artefacts at the image edges should be performed initially, but
this is expected to stay constant for a given imaging setup.

### `d_set_stall_thresholds.m`

The 'stallograms' generated in the previous step reflect image intensity
profiles of the edges that approximate biological capillary segments. The
purpose of the present script is to identify those intensity fluctuations
that correspond to 'true' stall events, where capillary flow is blocked.
In an ideal situation, without movement or other artefacts, a simple
threshold could be pre-specified. With real data, however, the threshold
of each segment must be determined individually: only a fraction of an
edge may actually be blocked (see manuscript for discussion), movement
may render portions of the image plane blank and be confused with an
actual stall, etc. Depending on the noise level and number of stalls,
experienced users can perform thresholding on a 60-frame acquisition in
5-15 min (ca. 150 edges/capillaries).

The GUI in the present scripts accepts the following interactions:
- up/down arrow: select previous/next edge
- left/right arrow: select previous/next frame
- 2,w,s,x: adjust the stall-threshold of the current edge in steps of
           +0.1, +0.01, -0.01, -0.1, respectively
- 1: toggle showing the current edge overlay (red)
- g: jump to edge number (type in command window)
- f: jump to frame number (type in command window)
- space: mark current frame bad

*NB! Remember to save the results (last cell).*

### `e_viz_stallogram_overlays.m`

Once a stallogram has been created (c) and binarised (d), we are able to
visualise and quantify the results. This scripts demonstrates ways to do
this, though other statistics may be relevant in particular applications.

Choose single roi folder (TIFF, cleaned), reads stalls from step d),
e.g., WT07/ROI1/2D_MIP_clean
