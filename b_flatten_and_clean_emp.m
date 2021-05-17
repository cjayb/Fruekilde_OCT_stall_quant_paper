%% b) Fruekilde et al. (in preparation)
% NB this script runs on the mat-files created by P-J. M. on the EPM
% datasets: 4D angios with n_y x n_y x 20 x n_frames voxels.
%
% After creating the reduced-depth angiogram slab in step a), we here apply
% the selection to all 3D data frames, and write out the reduced angiograms
% for all. This step is time-consuming, but can be performed with as little
% as 8 GB of RAM on a laptop in a matter of minutes (60 frames).
%
% The reduced-depth slabs are flattened using a maximum-intensity
% projection (MIP): log(max( ... , 3)), where the dimension of the
% max-operation is along the third (depth) dimension. Each frame is
% projected independently, cleaned, and saved to disk as a 2D TIFF image
% (both cleaned and original versions are saved).
%
% Note that the script calls clean_angio_strips_and_emph_vessels.m
% This is where the parameters for 'destriping' and vessel emphasis using
% the Frangi filter are defined (and should be adjusted manually).
%
% The implementation here is dependent on the following resources:
% https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-17-10-8567&id=179485
% https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
%%
clear all

OCT_base_folder = '/Users/au210321/data/Signe/OCT/EPM';
[filename, pathname] = uigetfile(fullfile(OCT_base_folder, '*.mat'));

file_base = split(filename, '.');
file_base = file_base{1};
out_folder_tif = fullfile(OCT_base_folder, file_base);
out_folder_tif_cleaned = fullfile(OCT_base_folder, file_base, 'cleaned');

mkdir(out_folder_tif)
mkdir(out_folder_tif_cleaned)

%%
load(fullfile(pathname, filename))

for dd = 1:size(angio, 4)
    frame_name = sprintf('frame_%03d.tif', dd);

    mip = double(log(squeeze(max(angio(:, :, :, dd), [], 1))));

    mip = (mip - min(mip(:))) / (max(mip(:)) - min(mip(:)));  % Normalization
    mip(isnan(mip(:))) = min(mip(:));  % NB!

    imwrite(mip, fullfile(out_folder_tif, frame_name));
    
    clean_mip = clean_angio_stripes_and_emph_vessels(mip);
    imwrite(clean_mip, fullfile(out_folder_tif_cleaned, frame_name));
end
