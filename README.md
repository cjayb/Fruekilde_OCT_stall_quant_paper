# Capillary stall quantification from optical coherence tomography
angiogram maximum intensity projections

By: Signe K. Fruekilde, Eugenio G. Jiménez, Kim R. Drasbek,
Christopher J. Bailey

This repository accompanies
([our bioRxiv preprint FIX LINK](https://osf.io/me6hf/)) documenting the
need and implementation of a method to obtain measures of capillary blood
flow stall events in pre-clinical optical coherence tomography (OCT) scans.

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
4. Select the folder with cleaned 2D MIPs, and let
`c_identify_cap_segments.m` automatically extract capillary segments from
from all the frames, applying correction of in-plane translation.
5. Based on the image intensity of each capillary segment over time,
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

### `a_create_first_angio_files.m`

This is the beginning of the analysis pipeline described in Fruekilde et
al. (in preparation). The objective is first to reduce the amount of data
to be read by identifying the focus depth of the acquired data. Depending
on the lense used, our experience is that a stack of only around 20-30
slices is in focus, rendering most of the 1,024 slices obsolete. To
reduce the memory and disk space consumption, we here define the focus
"slab', and save the information to disk. The next stage of processing
will read the saved 'angiogram' and apply the same ROI selection to all
frames. Select a 'representative frame'.

### `b_create_other_angios_destripe_frangi.m`

After creating the reduced-depth angiogram slab in step a), we here apply
the selection to all 3D data frames, and write out the reduced angiograms
for all. This step is time-consuming, but can be performed with as little
as 8 GB of RAM on a laptop in a matter of minutes (60 frames).

The reduced-depth slabs are flattened using a maximum-intensity
projection (MIP): log(max( ... , 3)), where the dimension of the
max-operation is along the third (depth) dimension. Each frame is
projected independently, cleaned, and saved to disk as a 2D TIFF image
(both cleaned and original versions are saved).

Note that the script calls clean_angio_stripes_and_emph_vessels.m
This is where the parameters for 'destriping' and vessel emphasis using
the Frangi filter are defined (and should be adjusted manually).
