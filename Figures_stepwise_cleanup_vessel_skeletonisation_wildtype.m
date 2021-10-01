animal = 'WT07';
roi = 'ROI1';
slice = '001';  % thesis; sharp stripe

scratch_base = '/Users/au210321/tmp';
image_out_folder = fullfile(scratch_base, animal, roi, 'Figures');
dpi = 300;
inc_titles = true;

mip_base = fullfile(scratch_base, animal, roi, '2D_MIP');
fname = fullfile(mip_base, ['test_spectral_' slice '_2D_MIP.tif']);
img = imread(fname);
%%
clear global
global cap_id filt_edgelist bin_stalls stallogram eq_vessels frame_id
global unique_thresholds show_cap_overlay bad_frames

nr = 2;
nc = 3;

figure(201); clf; colormap bone
set(gcf, 'Position', [ 39         490        1263         667])
subplot(nr, nc, 1); imagesc(img); axis equal; axis tight
tit = sprintf('Original, frame %s', slice);
% colorbar
set(gca, 'XTick', []); set(gca, 'YTick', []);
fname = sprintf('Fig1Ca_Orig_frame_%s.tif', slice);
exportgraphics(gca, fullfile(image_out_folder, fname), 'Resolution', dpi)
if inc_titles, title(tit), end
%% plot sequence of cleanup operations (Figure 1)
destripe_wavelet = 'db2';
destripe_keep_comps = 20;
destripe_sigma = 1;

% https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
cfg_opts = struct();
cfg_opts.FrangiScaleRange = [1, 3];  % The range of sigmas used, default [1 8]
cfg_opts.FrangiScaleRatio = 2;  % Step size between sigmas, default 2
cfg_opts.FrangiBetaOne = 0.5;  % Frangi correction constant, default 0.5
cfg_opts.FrangiBetaTwo = 10;  % Frangi correction constant, default 15
cfg_opts.BlackWhite = false;% Detect black ridges (default) set to true, for
                      % white ridges set to false.
cfg_opts.verbose = false;  % Show debug information, default true
destripe = RemoveStripesVertical(img, destripe_keep_comps, ...
                                 destripe_wavelet, destripe_sigma);
subplot(nr, nc, 2); imagesc(destripe); axis equal; axis tight
% colorbar
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, 'Fig1Cb_Destripe.tif'), 'Resolution', dpi)
if inc_titles, title('destripe'), end

med_destripe = medfilt2(destripe, [3 3]);
subplot(nr, nc, 3); imagesc(med_destripe); axis equal; axis tight
% colorbar
set(gca, 'XTick', []); set(gca, 'YTick', []);
if inc_titles, title('+ medfilt'), end

frangi = FrangiFilter2D(med_destripe, cfg_opts);
subplot(nr, nc, 4); imagesc(frangi); axis equal; axis tight
% colorbar
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, 'Fig1Cc_VesselFilter.tif'), 'Resolution', dpi)
if inc_titles, title('+ frangi'), end

clean_img = histeq(frangi / max(frangi(:)));
subplot(nr, nc, 5); imagesc(clean_img); axis equal; axis tight
% colorbar
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, 'Fig1Cd_IntensityEQ.tif'), 'Resolution', dpi)
if inc_titles, title('+ histeq (maxnorm)'), end


%% Vessel extraction
rim = 5;  % N pixel margin during vessel extraction
rim_xfactor = 4;  % top (small x) is more affected, increase rim there
min_conn_voxels = 100;  % in binarized vessel-images; area in pixels
min_branch_len = 10;  % in skeleton-images; length in pixels
min_edge_len = 16;  % disregard edges shorter than this (after edgelink)
init_stall_threshold = 0.5;  % initial cutoff, adjusted manually for each capillary

do_correct_motion = true;
xfm_type = 'translation';

mip_folder = [mip_base '_clean'];

tiffs = dir(fullfile(mip_folder, '*.tif'));
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
avg_vessel_mask = imbinarize(mean(eq_vessels, 3), 'adaptive');

%%
figure(202); clf; colormap bone
set(gcf, 'Position', [ 39         490        1263         667])
subplot(nr, nc, 1); imagesc(mean(eq_vessels, 3)); axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, 'Fig2A_AvgCleaned.tif'), 'Resolution', dpi)
title('Average of cleaned frames')

avg_vessel_mask = imbinarize(mean(eq_vessels, 3), 'adaptive');
subplot(nr, nc, 2); imagesc(avg_vessel_mask); axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, 'Fig2B_AvgBinarized.tif'), 'Resolution', dpi)
title('Binarized vessel mask')

%% Skeletonization and edge linkage

% FIND CONNECTED COMPONENTS!
CC = bwconncomp(avg_vessel_mask);
for uuu = 1:length(CC.PixelIdxList)
    if length(CC.PixelIdxList{uuu}) < min_conn_voxels  % at least this many voxels in 2D
        avg_vessel_mask(CC.PixelIdxList{uuu}) = 0;
    end
end

subplot(nr, nc, 3); imagesc(avg_vessel_mask); axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, 'Fig2C_ConnComps.tif'), 'Resolution', dpi)
title('Vessel mask after removal of smallest elements')

% these are actually quite different!
skel = bwskel(avg_vessel_mask, 'MinBranchLength', min_branch_len); 
% skel = bwmorph(avg_vessel_mask, 'skel', Inf);

subplot(nr, nc, 4); imagesc(skel); axis equal; axis tight
tit = sprintf('Skeletonized vessel mask; min branch len=%d vox', min_branch_len);
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, 'Fig2D_Skeletonization.tif'), 'Resolution', dpi)
title(tit)

[edgelist, edgeim, etype] = edgelink(skel);  % the magic happens here

filt_edgelist = {};
iiel = 0;
filt_edgeim = edgeim;

cm = get(gca, 'ColorOrder');
% filter shortest (<min_edge_len) edges away!
for ii = 1:length(edgelist)
    ci = mod(ii, size(cm, 1)) + 1;
    el = edgelist{ii};
    if size(el, 1) > min_edge_len
        iiel = iiel + 1;
        filt_edgelist{iiel} = el;
    else
        filt_edgeim(el(:,1), el(:,2)) = 0;
        adj_inds = filt_edgeim > iiel;
        filt_edgeim(adj_inds) = filt_edgeim(adj_inds) - 1;
    end
end

subplot(nr, nc, 5);
cmap = jet(length(filt_edgelist));
cmap = cmap(randperm(length(filt_edgelist)), :);
cmap = [0.2 0.2 0.2; cmap];
filt_edgeim_rgb = imdilate(ind2rgb(filt_edgeim,cmap), ones(3,3));
image(filt_edgeim_rgb); axis equal; axis tight
tit = sprintf('Capillary network; min edge len=%d vox', min_edge_len);
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, 'Fig2E_CapNetwork.tif'), 'Resolution', dpi)
title(tit)


subplot(nr, nc, 6);
underlay = max(eq_vessels, [], 3);
tit = sprintf('Capillaries on MAX over frames');
fname = sprintf('Fig2F_CapsOverlay_MaxFrames.tif', slice);

% underlay = mean(eq_vessels, 3);
% tit = sprintf('Capillaries on MEAN over frames');
% fname = sprintf('Fig2F_CapsOverlay_MeanFrames.tif', slice);

% underlay = eq_vessels(:, :, str2num(slice));
% tit = sprintf('Capillaries over frame %s', slice);
% fname = sprintf('Fig2F_CapsOverlay_Frame_%s.tif', slice);

imagesc(underlay); hold on; axis equal; axis tight

overlay = zeros(size(underlay));
for cap_id = 1:length(filt_edgelist)
    coords = filt_edgelist{cap_id};
    for jj = 1:size(coords, 1)
        overlay(coords(jj, 1), coords(jj, 2)) = cap_id;
    end
end
dil_overlay = imdilate(overlay, ones(3, 3));
cmap = jet(length(filt_edgelist));
cmap = cmap(randperm(length(filt_edgelist)), :);

overlay_h = image(ind2rgb(dil_overlay, cmap));
% inv_underlay = -underlay + 1;
% set(overlay_h, 'AlphaData', inv_underlay .* double(dil_overlay > 0))
set(overlay_h, 'AlphaData', 0.5 * double(dil_overlay > 0))
hold off
set(gca, 'XTick', []); set(gca, 'YTick', []);
exportgraphics(gca, fullfile(image_out_folder, fname), 'Resolution', dpi)
title(tit)

%%

figure(2003); clf
set(gcf, 'Position', [100   300   848   848])
subplot(1,2,2)
set(gca, 'FontSize', 12)

stallogram = zeros(length(filt_edgelist), length(tiffs), 'double');
for itiff = 1:length(tiffs)
    % binarize or not??
    BW = imbinarize(eq_vessels(:, :, itiff), 'adaptive');
    
    for iedge = 1:length(filt_edgelist)
        el = filt_edgelist{iedge};
        inds = sub2ind(size(BW), el(:, 1), el(:, 2));

        path_mask = false(size(BW)); path_mask(inds) = true;
        thick_solution_path = imdilate(path_mask, ones(3, 3));
        stallogram(iedge, itiff) = 1 - mean(BW(thick_solution_path));
        
    end
end
imagesc(stallogram); colormap hot
title('Capillary intensity plot')
xlabel('Time (frame no.)')
ylabel('Capillary ID')

fname = 'Fig3A_CapillaryIntensityPlot.tif';
exportgraphics(gca, fullfile(image_out_folder, fname), 'Resolution', dpi)

%% Check prev figure to make sure vessel extraction makes sense, then...
% examine one

iact_fig = figure(203); clf;
set(iact_fig, 'Position', [100   300   848   848])

cap_id = 1;
n_caps = size(stallogram, 1);
n_frames = size(stallogram, 2);
frame_id = 0;
show_cap_overlay = true;
bad_frames = false(n_frames, 1);

el = filt_edgelist{cap_id};
subplot(2,2,1); colormap gray; imagesc(imdilate(skel, ones(3,3))); hold on
plot(el(:,2), el(:,1), '.', 'Color', [1 0 0])
axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', [])
title('Capillary network')

if ~exist('unique_thresholds', 'var') || isempty(unique_thresholds)
    unique_thresholds = init_stall_threshold * ones(n_caps, 1);
    bin_stalls = imbinarize(stallogram, init_stall_threshold);
else
    cutoff = repmat(unique_thresholds, 1, n_frames);
    bin_stalls = stallogram > cutoff;
end

row_mask = zeros(size(bin_stalls), 'double'); row_mask(cap_id, :) = 1;

highlight_row = imoverlay(bin_stalls, row_mask.* bin_stalls, [1 0 0]);
subplot(1,2,2); imagesc(highlight_row)
title('Stallogram')
xlabel('Time (frame no.)')
ylabel('Capillary ID')

subplot(2,2,3); plot(stallogram(cap_id, :))
set(gca, 'Ylim', [0 1]); set(gca, 'YTick', [0 1])
set(gca, 'YTickLabel', {'Flow', 'Stall'})
ytickangle(90)
title('Stalling dynamics for selected capillary')
xlabel('Time (frame no.)')

axes = get(iact_fig, 'Children');
for an = 1:length(axes)
    set(axes(an), 'FontSize', 12)
end

set(iact_fig, 'KeyPressFcn', @setThresholds_callback);

fname = 'Fig3B_ThresholdingGUI_InitialView.tif';
exportgraphics(gcf, fullfile(image_out_folder, fname), 'Resolution', dpi)

%%
fname = 'Fig3B_ThresholdingGUI_NoisyExample.tif';
exportgraphics(gcf, fullfile(image_out_folder, fname), 'Resolution', dpi)
%%
fname = 'Fig3B_ThresholdingGUI_ClearExampleFlow.tif';
exportgraphics(gcf, fullfile(image_out_folder, fname), 'Resolution', dpi)
%%
fname = 'Fig3B_ThresholdingGUI_ClearExampleStall.tif';
exportgraphics(gcf, fullfile(image_out_folder, fname), 'Resolution', dpi)
