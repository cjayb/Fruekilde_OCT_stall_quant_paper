%% c) Fruekilde et al. (in preparation)
% Previous analysis stages have resulted in n_frames TIFF images of cleaned
% MIPs. This script reads in the folder of TIFFs, and begins by calculating
% an average across frames. It then identifies the frame which 'resembles'
% the average the most (similarity is quantified as the cosine of the
% 'angle' between images after flattening each into a 1-dimensional vector)
% Subsequently, each frame is coregistered to the 'most representative'
% single frame using only translations (no rotations or shears), and a new
% average is computed. This average is then skeletonised, from which edges
% are extracted using the edgelink-algorithm. Each one of these edges is
% dilated using a 3x3 image kernel to create n_edges 'capillary masks'.
% The image intensity within each masks across time frames thus depends on
% whether blood flows or is stalled. The 'stallogram' is a matrix of
% dimensions (n_frames, n_edges), i.e., a stacking of the extracted mask
% intensity time courses.
%
% This script does not require any user interaction, and runs in less than
% 30 seconds on a standard modern laptop. The output is a figure showing
% the final edge map of the image, and the stallogram. We recommend that
% each image is inspected visually to ensure that a reasonable-looking
% capillary network is extracted, with on the order of 120-180 edges for an
% image dimension of 400 x 400 pixels. Adjusting the image margin to avoid
% filter artefacts at the image edges should be performed initially, but
% this is expected to stay constant for a given imaging setup.

clear all

rim = 5;  % N pixel margin during vessel extraction
rim_xfactor = 4;  % top (small x) is more affected, increase rim there
min_conn_voxels = 100;  % in binarized vessel-images; area in pixels
min_branch_len = 10;  % in skeleton-images; length in pixels
min_edge_len = 16;  % disregard edges shorter than this (after edgelink)

do_correct_motion = true;
xfm_type = 'translation';

% scratch_folder = tempdir;
scratch_folder = '/Volumes/LPSDATA2/scratch/OCT_wildtypes_20pixels_cjb';
caps_viz_matfile = 'capmap.mat';
prev_ROI_dir = [];

fprintf(1, 'Select TIFF folder\n')
if isempty(prev_ROI_dir)
    mip_folder = uigetdir(scratch_folder);
else
    mip_folder = uigetdir(fullfile(prev_ROI_dir, '../'));
end

fname_caps = fullfile(mip_folder, caps_viz_matfile);

tiffs = dir(fullfile(mip_folder, '*.tif'));
if length(tiffs) < 1
    error(['No TIF images found in ' mip_folder])
end

while tiffs(1).name(1) == '.'  % delete macOS zombie files (QuickLook)
    delete(fullfile(tiffs(1).folder, tiffs(1).name))
    tiffs = tiffs(2:end);
end
info = imfinfo(fullfile(tiffs(1).folder, tiffs(1).name));
%
eq_vessels = zeros(info.Height, info.Width, length(tiffs), 'double');

fprintf(1, 'Reading in %d TIFFs and averaging from:\n', length(tiffs))
fprintf(1, '%s\n', mip_folder)
for ii = 1:length(tiffs)   
    
    mip = imread(fullfile(tiffs(ii).folder, tiffs(ii).name));
%     vessels = clean_angio_stripes_and_emph_vessels(mip);
    vessels = double(mip);  % assume already cleaned and equalized
    
    vessels(1:rim * rim_xfactor, :) = 0;
    vessels(info.Height-rim:end, :) = 0;
    vessels(:, 1:rim) = 0;
    vessels(:, info.Width-rim:end) = 0;    

    clean_img = histeq(vessels / max(vessels(:)));
    eq = imadjust(clean_img, [0.7 1]);
    
    eq_vessels(:, :, ii) = eq;
end
avg_vessel_mask = imbinarize(mean(eq_vessels, 3), 'adaptive');

% THIS is only used when starting from unclean tiffs 
% for ii = 1:length(tiffs)   
%     vessels = clean_angio_stripes_and_emph_vessels(eq_vessels(:, :, ii),...
%                                                    avg_vessel_mask);
%     
%     vessels(1:rim, :) = 0;
%     vessels(info.Height-rim:end, :) = 0;
%     vessels(:, 1:rim) = 0;
%     vessels(:, info.Width-rim:end) = 0;    
%     clean_img = histeq(vessels / max(vessels(:)));
%     
%     eq_vessels(:, :, ii) = imadjust(clean_img, [0.7 1]); 
%     if ~mod(ii, 10)
%         fprintf(1, 'Second pass, filtered %d of %d\n', ii, length(tiffs))
%     end
% end
% avg_vessel_mask = imbinarize(mean(eq_vessels, 3), 'adaptive');

% REGISTER all frames to the one most similar to the average
if do_correct_motion    

    [optimizer, metric] = imregconfig('monomodal');
    
    % calculate "similarity" between each frame and the mean (cosine of angle)
    similarity = zeros(length(tiffs), 1);
    n_pixels = info.Width * info.Height;
    avg_vessel_vector = reshape(double(mean(eq_vessels, 3)), n_pixels, 1);
    for ii = 1:length(tiffs)
        cur_frame = reshape(eq_vessels(:, :, ii), n_pixels, 1);
        similarity(ii) = dot(cur_frame, avg_vessel_vector) /...
                            (norm(cur_frame, 2) * norm(avg_vessel_vector, 2));
    end

    % select the single frame most similar to the average (small angle => large cosine)
    [~, target_frame_idx] = max(similarity);
    if length(target_frame_idx) > 1
        target_frame_idx = target_frame_idx(1);
    end

    targetIm = eq_vessels(:, :, target_frame_idx);
    rsl_vessels = eq_vessels;

    for itiff = 1:length(tiffs)
        if itiff == 1
            steps = ceil(length(tiffs) / 10);
            steps_str = ''; for ii = 1:steps, steps_str = [steps_str '=']; end
            fprintf(1, '\nRegistering to frame %d | %s', target_frame_idx, steps_str)
        elseif ~mod(itiff, 10)
            fprintf(1, '\b')
        end
        if itiff == target_frame_idx
            continue
        end
        xfm = imregtform(eq_vessels(:, :, itiff), targetIm, xfm_type, ...
                         optimizer, metric);    
        rsl_vessels(:, :, itiff) = imwarp(eq_vessels(:, :, itiff), xfm,...
                                          'OutputView', imref2d(size(targetIm)));
    end

    % replace with motion corrected version!
    eq_vessels = rsl_vessels;
end

% recalculate average!
avg_vessel_mask = imbinarize(mean(eq_vessels, 3), 'adaptive');

%
% FIND CONNECTED COMPONENTS
CC = bwconncomp(avg_vessel_mask);
for uuu = 1:length(CC.PixelIdxList)
    if length(CC.PixelIdxList{uuu}) < min_conn_voxels  % at least this many voxels in 2D
        avg_vessel_mask(CC.PixelIdxList{uuu}) = 0;
    end
end

% these are actually quite different!
skel = bwskel(avg_vessel_mask, 'MinBranchLength', min_branch_len); 
% skel = bwmorph(avg_vessel_mask, 'skel', Inf);

[edgelist, edgeim, etype] = edgelink(skel);  % the magic happens here

figure(301); clf; 
subplot(2,2,1); colormap bone; imagesc(skel); hold on
ylabel('X'); xlabel('Y')

filt_edgelist = {};
iiel = 0;
filt_edgeim = edgeim;

cm = get(gca, 'ColorOrder');
% filter shortest (<min_edge_len) edges away!
for ii = 1:length(edgelist)
    ci = mod(ii, size(cm, 1)) + 1;
    el = edgelist{ii};
    if size(el, 1) > min_edge_len
        plot(el(:,2), el(:,1), '.', 'Color', cm(ci,:))
        iiel = iiel + 1;
        filt_edgelist{iiel} = el;
    else
        filt_edgeim(el(:,1), el(:,2)) = 0;
        adj_inds = filt_edgeim > iiel;
        filt_edgeim(adj_inds) = filt_edgeim(adj_inds) - 1;
    end
end

% figure(212); clf
subplot(2,2,3)
cmap = jet(length(filt_edgelist));
cmap = cmap(randperm(length(filt_edgelist)), :);
cmap = [0.3 0.3 0.3; cmap];
filt_edgeim_rgb = imdilate(ind2rgb(filt_edgeim,cmap), ones(3,3));
image(filt_edgeim_rgb);
ylabel('X'); xlabel('Y')

fprintf(1, '\n%d edges (capillaries) detected after filtering.\n', ...
        length(filt_edgelist))

%
stallogram = zeros(length(filt_edgelist), length(tiffs), 'double');
for itiff = 1:length(tiffs)
    
    BW = imbinarize(eq_vessels(:, :, itiff), 'adaptive');
    
    for iedge = 1:length(filt_edgelist)
        el = filt_edgelist{iedge};
        inds = sub2ind(size(BW), el(:, 1), el(:, 2));

        path_mask = false(size(BW)); path_mask(inds) = true;
        thick_solution_path = imdilate(path_mask, ones(3,3));
        stallogram(iedge, itiff) = 1 - mean(BW(thick_solution_path));
        
    end
    if itiff == 1
        steps = ceil(length(tiffs) / 10);
        steps_str = ''; for ii = 1:steps, steps_str = [steps_str '=']; end
        fprintf(1, '\nCreating stallogram | %s', steps_str)
    elseif ~mod(itiff, 10)
        fprintf(1, '\b')
    end
end
fprintf(1, '\n')
subplot(1, 2, 2)
imagesc(stallogram); colormap hot
xlabel('Time (frame number)')
ylabel('Capillary segment (edge)')

% save cleaned MIP maps and skeletonized edges (after filtering)
save(fname_caps, 'eq_vessels', 'filt_edgeim', 'filt_edgelist', ...
                 'target_frame_idx', 'skel', 'stallogram');
fprintf(1, 'Capillary network and stallogram saved to: %s\n', fname_caps)

prev_ROI_dir = mip_folder;
