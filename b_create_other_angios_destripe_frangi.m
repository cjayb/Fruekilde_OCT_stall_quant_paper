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
% Note that the script calls clean_angio_strips_and_emph_vessels.m
% This is where the parameters for 'destriping' and vessel emphasis using
% the Frangi filter are defined (and should be adjusted manually).
%
% The implementation here is dependent on the following resources:
% https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-17-10-8567&id=179485
% https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
%%
clear all

scratch_folder = '/tmp/OCT';
OCT_base_folder = '/Users/au210321/data/Signe/OCT/nii';

% while we're at it, we might as well create 2D projections for various 
do_slabs = [20, 10, 30];

%%
do_rois = {'ROI1_angio', 'ROI2_angio'};

% This is OPTIONAL: mark the  animals and conditions you want to
% include. Otherwise the code will add everything it finds!
% do_animals = {'Mouse1', 'Mouse4', 'Mouse5'};
% do_conditions = {'baseline', 'day2'};
do_animals = {};
do_conditions = {};

all_ROI_folders = dir(fullfile(manual_folder, '**', 'ROI*_angio'));
for ii = 1:length(all_ROI_folders)

    bar = split(all_ROI_folders(ii).folder, filesep);  % file separator for this platform!
    animal = bar{end};
    condition = bar{end - 1};

    first_angio_name = fullfile(all_ROI_folders(ii).folder,...
                                all_ROI_folders(ii).name,...
                                'test_spectral_001_angio.mat');
    if ~exist(first_angio_name, 'file')
        error('Cannot find %s', first_angio_name)
    end
    
    if ~any(strcmp(do_animals, animal))
        do_animals{end + 1} = animal;
        fprintf(1, 'Added %s\n', animal)
    end
    if ~any(strcmp(do_conditions, condition))
        do_conditions{end + 1} = condition;
        fprintf(1, 'Added %s\n', condition)
    end
end
%%
clear animal condition
roi_folders_completed = {};
status_file = sprintf('batch_status-%s.mat', date);
if exist(status_file, 'file') == 2
    load(status_file, 'roi_folders_completed')
end

for this_condition = do_conditions
    condition = this_condition{1};
    cond_folder = fullfile(manual_folder, condition);
    for this_animal = do_animals
        animal = this_animal{1};
        for this_roi = do_rois
            roi = this_roi{1};

            % get the manually defined z-range
            in_folder = fullfile(cond_folder, animal, roi);
            
            if any(strcmp(roi_folders_completed, in_folder))
                fprintf(1, '%s, %s, %s already done, skipping\n', condition, animal, roi)
                continue
            end
            
            first_angio_name = fullfile(in_folder, 'test_spectral_001_angio.mat');
            if ~exist(first_angio_name, 'file')
                error('Cannot find %s', first_angio_name)
            end
            load(first_angio_name, 'angio', 'pl_zrange', 'nii_zrange', 'zstack')
            pl_zrange_man = pl_zrange;

            all_raws = dir(fullfile(OCT_base_folder, condition, ...
                animal, roi, 'test_spectral_*.nii'));
            if length(all_raws) ~= 60
                error('Found %d raw NII files (not 60)', length(all_raws))
            end

            tic;
            % NB don't skip first file, as we need to do the extra slabs
            for fidx = 1:length(all_raws)
                cur_fname = fullfile(all_raws(fidx).folder, all_raws(fidx).name);
                nii = NiiReader(cur_fname);
                nii.zrange = nii_zrange;
                fprintf(1, '%s->%s->%s->Calculating angio %d of %d ... ', condition, animal, roi, fidx, length(all_raws))            
                raw_angio = calculate_angio(nii);  % G
                nii.close();                
                
                % just redo the 20 as well...
                for slab_thickness = do_slabs
                    
                    pl_zrange_adjust = floor((slab_thickness - length(pl_zrange_man)) / 2);                    
                    pl_zrange = pl_zrange_man(1) - pl_zrange_adjust:pl_zrange_man(end) + pl_zrange_adjust;                    
                                        
                    zrange_str = sprintf('_%d_%d', pl_zrange(1) + (nii_zrange(1) - 1), ...
                        pl_zrange(end) + (nii_zrange(1) - 1));  % refer to the slices in RAW (nii) data!

                    OCT_folder = sprintf('OCT_%dpixels', slab_thickness);
                    out_folder = fullfile(scratch_folder, OCT_folder, this_condition{1}, animal, roi);
                    
                    out_folder_mip = fullfile(out_folder, ['MIPData', zrange_str]);
                    out_folder_mip_cleaned = fullfile(out_folder, ['MIPData', zrange_str], 'cleaned');

                    out_folder_tif = fullfile(out_folder, ['MIPDataTiff', zrange_str]);
                    out_folder_tif_cleaned = fullfile(out_folder, ['MIPDataTiff', zrange_str], 'cleaned');

                    if fidx == 1
                        mkdir(out_folder)
                        mkdir(out_folder_mip)
                        mkdir(out_folder_mip_cleaned)
                        mkdir(out_folder_tif)
                        mkdir(out_folder_tif_cleaned)
                    end
                        
                    % save all angio 
                    new_angio_base = fullfile(out_folder, all_raws(fidx).name(1:end-4));        
                    angio = raw_angio(pl_zrange, :, :);  % just save the slab!

                    % Don't overwrite the manual angio!
                    if ~(exist([new_angio_base, '_angio.mat'], 'file') == 2)
                        fprintf(1, 'and saving (only) chosen focus slab thickness %d ...\n', slab_thickness)
                        save([new_angio_base, '_angio.mat'], ...
                            'angio','pl_zrange','nii_zrange', 'zstack', '-v7.3');
                    else
                        fprintf(1, '(skipping existing manual angio ...)\n')
                    end

                    % save to MIPData folders!
                    fprintf(1, 'along with MIPData in mat and tif formats ...\n')
                    outname_mip = save_angios_as_MIPData(...
                        angio, out_folder_mip, out_folder_tif,...
                        [all_raws(fidx).name(1:end-4), '_angio.mat']);

                    % now 'angio' is actually the MIP (logarithmic)
                    load(outname_mip, 'angio');
                    orig_mip = angio;

                    mip = orig_mip;
                    mip(isnan(mip(:))) = min(mip(:));  % NB!
                    clean_mip = clean_angio_stripes_and_emph_vessels(mip);

                    fprintf(1, 'and cleaned MIPdata in mat and tiff format ...\n')                
                    outname_clean_mip = fullfile(out_folder_mip_cleaned,...
                        [all_raws(fidx).name(1:end-4), '_angio.mat']);
                    angio = clean_mip;
                    save(outname_clean_mip, 'angio');
                    outname_clean_mip_tiff = fullfile(out_folder_tif_cleaned,...
                        [all_raws(fidx).name(1:end-4), '_MIPz.tif']);
                    imwrite(clean_mip, outname_clean_mip_tiff);
                                        
                    if slab_thickness == 20
                        figure(601); clf; colormap bone; 
                        subplot(1,2,1)
                        imagesc(orig_mip); axis equal; axis tight
%                         set(gca,'CLim', [0, 0.75]);
                        xlabel('XY-plane')
                        set(gcf, 'Position', [15   759   624   278])
                        title('First frame MIP before...')
                        subplot(1,2,2)
                        imagesc(clean_mip); axis equal; axis tight
%                         set(gca,'CLim', [0, 0.75]);
                        xlabel('XY-plane')
                        set(gcf, 'Position', [15   759   624   278])
                        title('...and after cleaning')
                        drawnow
                    end                                
                    
                end % slab!
            end  % raw nii
            t_60 = toc;
            fprintf(1, '\n\n*** Processing %s, %s, %s took %.1f hours ***\n\n', condition, animal, roi, t_60 / 3600)

            roi_folders_completed{end + 1} = in_folder;
            save(status_file, 'roi_folders_completed')  % overwrite
        end % roi
    end % animal
end