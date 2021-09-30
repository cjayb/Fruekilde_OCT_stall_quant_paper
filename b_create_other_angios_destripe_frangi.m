%% b) Fruekilde et al. (in preparation)
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
% Note that the script calls clean_angio_stripes_and_emph_vessels.m
% This is where the parameters for 'destriping' and vessel emphasis using
% the Frangi filter are defined (and should be adjusted manually).
%
% The implementation here is dependent on the following resources:
% https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-17-10-8567&id=179485
% https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
%%
clear all

% location of 3D_angio-folder
fprintf(1, 'Choose 3D angio file of first frame\n')
[first_angio_name, input_folder] = uigetfile('*.mat', 'Select first input');
load(fullfile(input_folder, first_angio_name), 'angio', 'pl_zrange', ...
                                               'nii_zrange', 'zstack', ...
                                               'nii_folder')

slab_thickness = 20;  % MIP thickness in pixels; 20 used in paper
show_cleaning_process = true;  % plot effect of filtering
                                           
%%

all_raws = dir(fullfile(nii_folder, '*.nii'));
pl_zrange_man = pl_zrange;

tic;
% NB don't skip first file, as we need to do the extra slabs
for fidx = 1:length(all_raws)
    cur_fname = fullfile(all_raws(fidx).folder, all_raws(fidx).name);
    nii = NiiReader(cur_fname);
    nii.zrange = nii_zrange;
    fprintf(1, 'Calculating angio %d of %d ... ', fidx, length(all_raws))
    raw_angio = calculate_angio(nii);  % G
    nii.close();
        
    pl_zrange_adjust = floor((slab_thickness - length(pl_zrange_man)) / 2);
    pl_zrange = pl_zrange_man(1) - pl_zrange_adjust:pl_zrange_man(end) + pl_zrange_adjust;
        
    out_folder_mip = fullfile(input_folder, '..', '2D_MIP');
    out_folder_mip_clean = fullfile(input_folder, '..', '2D_MIP_clean');
    
    if fidx == 1
        mkdir(out_folder_mip)
        mkdir(out_folder_mip_clean)
    end
    
    % save all angio
    new_angio_base = fullfile(input_folder, all_raws(fidx).name(1:end-4));
    angio = raw_angio(pl_zrange, :, :);  % just save the slab!
    
    % Don't overwrite the manual angio!
    new_angio_file = [new_angio_base, '_3D.mat'];
    if ~(exist(new_angio_file, 'file') == 2)
        fprintf(1, 'saving chosen slab (thickness %d) ... ', slab_thickness)
        save(new_angio_file, ...
            'angio','pl_zrange','nii_zrange', 'zstack', 'nii_folder', '-v7.3');
    else
        fprintf(1, '(skipping: existing manual angio) ... ')
    end
    
    % save to MIP folders and return log-scaled mip
    fprintf(1, 'along with MIPs in tif format\n')
    orig_mip = save_angio_as_MIP(angio, out_folder_mip,...
        [all_raws(fidx).name(1:end-4), '_2D_MIP.tif']);
    
    mip = orig_mip;
    mip(isnan(mip(:))) = min(orig_mip(:));  % NB!
    if show_cleaning_process
        [destripe, frangi, clean_mip] = ...
            clean_angio_stripes_and_emph_vessels(mip);
        plot_cleaning_progress_diagnostics(mip, destripe, frangi, clean_mip)
    else
        clean_mip = clean_angio_stripes_and_emph_vessels(mip);
    end
 
    outname_clean_mip_tiff = fullfile(out_folder_mip_clean,...
        [all_raws(fidx).name(1:end-4), '_2D_MIP.tif']);
    imwrite(clean_mip, outname_clean_mip_tiff);
    
    
end
t_60 = toc;
fprintf(1, '\n\n*** Processing took %.1f minutes ***\n\n', t_60 / 60)
