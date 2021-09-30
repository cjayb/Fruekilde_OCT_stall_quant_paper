%% d) Fruekilde et al. (in preparation)
% The 'stallograms' generated in the previous step reflect image intensity
% profiles of the edges that approximate biological capillary segments. The
% purpose of the present script is to identify those intensity fluctuations
% that correspond to 'true' stall events, where capillary flow is blocked.
% In an ideal situation, without movement or other artefacts, a simple
% threshold could be pre-specified. With real data, however, the threshold
% of each segment must be determined individually: only a fraction of an
% edge may actually be blocked (see manuscript for discussion), movement
% may render portions of the image plane blank and be confused with an
% actual stall, etc. Depending on the noise level and number of stalls,
% experienced users can perform thresholding on a 60-frame acquisition in
% 5-15 min (ca. 150 edges/capillaries).
%
% The GUI in the present scripts accepts the following interactions:
% - up/down arrow: select previous/next edge
% - left/right arrow: select previous/next frame
% - 2,w,s,x: adjust the stall-threshold of the current edge in steps of
%            +0.1, +0.01, -0.01, -0.1, respectively
% - 1: toggle showing the current edge overlay (red)
% - g: jump to edge number (type in command window)
% - f: jump to frame number (type in command window)
%
% NB! Remember to save the results (last cell)

clear all

init_stall_threshold = 0.67;  % initial cutoff, adjusted manually for each capillary

thresh_matfile = 'unique_stall_thresholds.mat';
stalls_matfile = 'stalls.mat';
caps_viz_matfile = 'capmap.mat';
prev_ROI_dir = [];
%%%

clear global
global cap_id filt_edgelist bin_stalls stallogram eq_vessels frame_id
global unique_thresholds show_cap_overlay bad_frames

fprintf(1, 'Select TIFF folder\n')
if isempty(prev_ROI_dir)
    mip_folder = uigetdir();
else
    mip_folder = uigetdir(fullfile(prev_ROI_dir, '../'));
end

fname_thresh = fullfile(mip_folder, thresh_matfile);
fname_caps = fullfile(mip_folder, caps_viz_matfile);

if exist(fname_thresh, 'file') == 2
    fprintf(1, 'Found previously defined thresholds, using them!\n')
    load(fname_thresh, 'unique_thresholds');
end

% get stallogram and cap segments
fprintf(1, 'Loading previously identified capillary segments...')
load(fname_caps, 'eq_vessels', 'filt_edgeim', 'filt_edgelist', ...
                 'target_frame_idx', 'skel', 'stallogram')
fprintf(1, 'done.\n')

tiffs = dir(fullfile(mip_folder, '*.tif'));
while tiffs(1).name(1) == '.'  % delete macOS zombie files (QuickLook)
    delete(fullfile(tiffs(1).folder, tiffs(1).name))
    tiffs = tiffs(2:end);
end
info = imfinfo(fullfile(tiffs(1).folder, tiffs(1).name));

iact_fig = figure(302); clf;
set(iact_fig, 'Position', [100, 100, 800, 800])

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
set(gca, 'XLim', [1 n_frames])
set(gca, 'YTickLabel', {'Flow', 'Stall'})
ytickangle(90)
title('Stalling dynamics for selected capillary')
xlabel('Time (frame no.)')

axes = get(iact_fig, 'Children');
for an = 1:length(axes)
    set(axes(an), 'FontSize', 12)
end

set(iact_fig, 'KeyPressFcn', @setThresholds_callback);

%% done? REMEMBER TO SAVE! (overwrites without warning)
save(fname_thresh, 'unique_thresholds', 'bad_frames');
fprintf(1, 'Thresholds saved to: %s\n', fname_thresh)

fname_stalls = fullfile(mip_folder, stalls_matfile);

% save binary stallograms, as well as all the variables needed to reproduce
save(fname_stalls, 'bin_stalls', 'stallogram', 'unique_thresholds',...
                   'bad_frames', 'filt_edgelist');
fprintf(1, 'Stalls saved to: %s\n', fname_stalls)

prev_ROI_dir = mip_folder;