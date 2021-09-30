% Choose single roi folder (TIFF, cleaned), reads stalls from step d)

scratch_base = '/Users/au210321/tmp';
image_out_folder = fullfile(scratch_base, 'wt07_roi1_figure4');
mkdir(image_out_folder)

thresh_matfile = 'unique_stall_thresholds.mat';
stalls_matfile = 'stalls.mat';
caps_matfile = 'capmap.mat';
min_stall_len = 2;

fprintf(1, 'Select TIFF folder\n')
mip_folder = uigetdir(scratch_base);

fname_stalls = fullfile(mip_folder, stalls_matfile);
fname_caps = fullfile(mip_folder, caps_matfile);

dpi = 300;
%%
clear global
global filt_edgelist sgram plot_mip


% save binary stallograms, as well as all the variables needed to reproduce
load(fname_stalls, 'bin_stalls', 'filt_edgelist');
load(fname_caps, 'eq_vessels');
             
sgram = filter_stallogram(bin_stalls, min_stall_len);
rowsum = sum(single(sgram), 2);
rowsum = rowsum(rowsum > 0);
fprintf(1, 'Mean stall rate: %.02f\n', 100 * mean(rowsum) / size(sgram, 2))
fprintf(1, 'Min stall rate: %.02f\n', 100 * min(rowsum) / size(sgram, 2))
fprintf(1, 'Max stall rate: %.02f\n', 100 * max(rowsum) / size(sgram, 2))
%
iact_fig = figure(303); clf;
% set(iact_fig, 'Position', [100, 100, 800, 800])

ax1 = axes;

[min_stalls, min_idx] = min(sum(bin_stalls, 1));
plot_mip = max(eq_vessels, [], 3);
% plot_mip = median(eq_vessels, 3);
% plot_mip = eq_vessels(:, :, min_idx);
viz_all_caps(plot_mip, sgram, filt_edgelist)

set(iact_fig, 'KeyPressFcn', @setOverlays_callback);

%% only run if want to export

fname = 'Fig4A_CapillariesStalling.tif';
% fname = 'Fig4B_StallRatePercent.tif';
% fname = 'Fig4C_LongestStall.tif';
exportgraphics(gcf, fullfile(image_out_folder, fname), 'Resolution', dpi)
%% NEW vs OLD comparison?
figure(601); clf
set(gcf, 'Position', [100   140   400  700 ])
new_sgram = bin_stalls + sgram;
cmap = [1 1 1; 0 0 1; 1 0 0];
imagesc(new_sgram); colormap(cmap)
axis equal; axis tight
title('Binary stallogram')
xlabel('Time (frame no.)')
ylabel('Capillary (index no.)')


fname = 'Fig4A_binary_stall_matrix.tif';
exportgraphics(gcf, fullfile(image_out_folder, fname), 'Resolution', dpi)



%%%%%%%%%%%%%%%
function viz_all_caps_orig(plot_mip, filt_edgeim, filt_edgelist)

clf
imagesc(plot_mip); colormap gray
hold on

cmap = jet(length(filt_edgelist));
cmap = cmap(randperm(length(filt_edgelist)), :);
dil_edgeim = imdilate(filt_edgeim, ones(3, 3));
filt_edgeim_rgb = ind2rgb(dil_edgeim, cmap);
overlay_h = image(filt_edgeim_rgb);
ylabel('X'); xlabel('Y')

set(overlay_h, 'AlphaData', 0.5 * double(dil_edgeim > 0))
th = title(sprintf('All capillaries (%d)', length(filt_edgelist)));
set(th, 'FontSize', 14)
end

function viz_all_caps(plot_mip, sgram, filt_edgelist, clear_image)

clear_image = false;
if nargin < 4
    clear_image = true;
end
if clear_image, clf; end

imagesc(plot_mip); colormap gray
hold on
overlay = zeros(size(plot_mip));
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
set(overlay_h, 'AlphaData', 0.5 * double(dil_overlay > 0))
axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);


th = title(sprintf('All capillaries (%d)', length(filt_edgelist)));
set(th, 'FontSize', 14)
end


function viz_any_stalls(plot_mip, sgram, filt_edgelist, clear_image)

clear_image = false;
if nargin < 4
    clear_image = true;
end
if clear_image, clf; end

overlay = zeros(size(plot_mip));
n_caps_stalling = 0;
for cap_id = 1:length(filt_edgelist)
    coords = filt_edgelist{cap_id};
    if sum(sgram(cap_id, :)) > 0
        for jj = 1:size(coords, 1)
            overlay(coords(jj, 1), coords(jj, 2)) = cap_id;
        end
        n_caps_stalling = n_caps_stalling + 1;
    end
end
imagesc(plot_mip); colormap gray; axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
hold on
dil_overlay = imdilate(overlay, ones(3, 3));
cmap = jet(length(filt_edgelist));
cmap = cmap(randperm(length(filt_edgelist)), :);

overlay_h = image(ind2rgb(dil_overlay, cmap));
set(overlay_h, 'AlphaData', 0.5 * double(dil_overlay > 0))

th = title(sprintf('Stalling capillaries (%d of %d; %.2f%%)', ...
    n_caps_stalling, length(filt_edgelist), ...
    100 * n_caps_stalling / length(filt_edgelist)));
set(th, 'FontSize', 14)

end

%%%%%%
function viz_stall_rate(plot_mip, sgram, filt_edgelist, clear_image)

clear_image = false;
if nargin < 4
    clear_image = true;
end
if clear_image, clf; end

n_frames = size(sgram, 2);

ax1 = axes;

overlay = zeros(size(plot_mip));
n_caps_stalling = 0;
for cap_id = 1:length(filt_edgelist)
    coords = filt_edgelist{cap_id};
    if sum(sgram(cap_id, :)) > 0
        for jj = 1:size(coords, 1)
            overlay(coords(jj, 1), coords(jj, 2)) = 100 * sum(sgram(cap_id, :)) / n_frames;
        end
        n_caps_stalling = n_caps_stalling + 1;
    end
end
imagesc(plot_mip); axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
colormap(ax1,'gray');
hold on
th = title(sprintf('Stall rate (percent of %d frames)', n_frames));
set(th, 'FontSize', 14)

ax2 = axes;

dil_overlay = imdilate(overlay, ones(3, 3));

imagesc(ax2, dil_overlay, 'alphadata', 0.75 * double(dil_overlay > 0));
axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
colormap(ax2, 'parula');
caxis(ax2, [0 max(nonzeros(dil_overlay))]);
linkprop([ax1 ax2], 'Position');
ax2.Visible = 'off';
colorbar;
set(ax1, 'Position', get(ax2, 'Position'))


end

%%%%%%
function viz_longest_stall(plot_mip, sgram, filt_edgelist, clear_image)

clear_image = false;
if nargin < 4
    clear_image = true;
end
if clear_image, clf; end
ax1 = axes;

overlay = zeros(size(plot_mip));
[lens, vals] = runLengthEncodeRows(sgram);
n_caps_stalling = 0;
for cap_id = 1:length(filt_edgelist)
    if sum(sgram(cap_id, :)) > 0
        coords = filt_edgelist{cap_id};
        stalls = find(vals{cap_id});
        runlen = lens{cap_id};
        stall_durations = runlen(stalls);  % length in frames of each stall
        for jj = 1:size(coords, 1)
            overlay(coords(jj, 1), coords(jj, 2)) = max(stall_durations);
        end
        n_caps_stalling = n_caps_stalling + 1;        
    end
end
imagesc(plot_mip); axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
colormap(ax1,'gray');
hold on
th = title(sprintf('Longest stall (# frames)'));
set(th, 'FontSize', 14)
ax2 = axes;

dil_overlay = imdilate(overlay, ones(3, 3));

imagesc(ax2, dil_overlay, 'alphadata', 0.75 * double(dil_overlay > 0));
axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
colormap(ax2, 'parula');
caxis(ax2, [0 max(nonzeros(dil_overlay))]);
linkprop([ax1 ax2], 'Position');
ax2.Visible = 'off';
colorbar;
set(ax1, 'Position', get(ax2, 'Position'))

end


function setOverlays_callback(fig_obj, eventData)

global filt_edgelist sgram plot_mip

ck = get(fig_obj, 'CurrentKey');        

if strcmp(ck, '1')
    viz_any_stalls(plot_mip, sgram, filt_edgelist)
elseif strcmp(ck, '2')
    viz_stall_rate(plot_mip, sgram, filt_edgelist)
elseif strcmp(ck, '3')
    viz_longest_stall(plot_mip, sgram, filt_edgelist)   
elseif strcmp(ck, '0')
    viz_all_caps(plot_mip, sgram, filt_edgelist)     
end
end