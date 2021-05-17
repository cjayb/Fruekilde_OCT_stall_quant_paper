%% a) Fruekilde et al. (in preparation)
% This is the beginning of the analysis pipeline described in Fruekilde et
% al. (in preparation). The objective is first to reduce the amount of data
% to be read by identifying the focus depth of the acquired data. Depending
% on the lense used, our experience is that a stack of only around 20-30
% slices is in focus, rendering most of the 1,024 slices obsolete. To
% reduce the memory and disk space consumption, we here define the focus
% "slab', and save the information to disk. The next stage of processing
% will read the saved 'angiogram' and apply the same ROI selection to all
% frames. Select a 'representative frame'.

%% locate the effective data range (here in Nifti format)

clear all

% all output will be written here
scratch_folder = '/tmp/OCT';

OCT_base_folder = '/Users/au210321/data/Signe/OCT/nii';
[filename, pathname] = uigetfile(fullfile(OCT_base_folder, '*.nii'));

export_QC_images = false;
dpi = 300;
image_out_folder = fullfile(scratch_folder, 'Figures');

avg_frames_y = [190, 210];  % visualise around the mid-y in a 400x400 stack

nii = NiiReader(fullfile(pathname, filename));
avg_ref = nii.readFrame(avg_frames_y(1));  % z-by-x
plot_nx = floor(size(avg_ref, 2) / 2);

for ns = avg_frames_y(1) + 1:avg_frames_y(2)
    avg_ref = avg_ref + nii.readFrame(ns);
end
avg_ref = avg_ref / (diff(avg_frames_y) + 1);
%% visualise XZ-plane
figure(501); clf
Bscan = abs(avg_ref(:, 1:plot_nx)) + abs(avg_ref(:, plot_nx+1:end));
Bscan = log(Bscan / max(Bscan(:)));
imagesc(Ascan); axis equal; axis tight
colormap parula
% set(gca, 'CLim', [0 0.5])
set(gcf, 'Position', [40   180   480   890])
set(gca, 'FontSize', 12)
tit = sprintf('Log mean magnitude image (y=[%d, %d])', avg_frames_y);
ylabel('Z (pixels)'); xlabel('X (pixels)')
title(tit)
if export_QC_images
    fname = 'Fig1A_MeanA-scan-logarithmic.tif';
    exportgraphics(gcf, fullfile(image_out_folder, fname), 'Resolution', dpi)    
end
% colorbar

%% Select range to reconstruct in 3D

% choose a generous range (not too tight)
nii.zrange = [50, 200];

nii_zrange = nii.zrange;

% phase correction takes a minute or two
angio = calculate_angio(nii);
nx = nii.dims(2);  % original x-width, twice actual image size

nii.close();
%% define and center the 3D slab to analyse
zstack = 20;  % thickness of slab to look at in z-direction

xstack = zstack; ystack = zstack;  % just use the same in x-y
angio_zsize = diff(nii.zrange);

% where to start looking?
global pl_zrange  % for slider to return the chosen range

% WRONG! pl_zstart = nii.zrange(1) + floor(diff(nii.zrange) / 2) - floor(zstack / 2);
% pl_zstart is relative to angio, where index 1 = nii.zrange(1) in raw nii
pl_zstart = floor(angio_zsize / 2) - floor(zstack / 2);

pl_zrange = pl_zstart:pl_zstart + zstack - 1;

pl_xstart = (nx/2) / 2 - floor(xstack / 2);  % halve the x-dim
pl_xrange = pl_xstart:pl_xstart + xstack - 1;  % select midslice in x
pl_yrange = pl_xstart:pl_xstart + xstack - 1;  % select midslice in y

figure(503); clf; colormap bone
set(gcf, 'Position', [100 200 800 320])
clim = [0, 0.5];

% subplot(1,3,1);
subplot(1,2,1);
imagesc(squeeze(mean(angio(pl_zrange, :, :), 1))); axis equal; axis tight
set(gca,'CLim', clim); xlabel('XY-plane (pixels)')
% subplot(1,3,2)
subplot(2,2,2)
imagesc(squeeze(mean(angio(:, pl_xrange, :), 2))); axis equal; axis tight
set(gca,'CLim', clim); xlabel('XZ-plane (pixels)'); ylabel('Z (depth)')
% subplot(1,3,3)
subplot(2,2,4)
imagesc(squeeze(mean(angio(:, :, pl_yrange), 3))); axis equal; axis tight
set(gca,'CLim', clim); xlabel('YZ-plane (pixels)'); ylabel('Z (depth)')
drawnow

add_slider_to_angio_plot(angio, zstack, pl_zstart, angio_zsize, clim)

axes = get(gcf, 'Children');
for an = 1:length(axes)
    set(axes(an), 'FontSize', 12)
end
%% export QC images showing selected slab depth (OPTIONAL)
if export_QC_images
    fname = 'Fig1B_SlabSelection-noGUI.tif';
    exportgraphics(gcf, fullfile(image_out_folder, fname), 'Resolution', dpi)
    fname = 'Fig1B_SlabSelection-GUI.tif';
    print(gcf, fullfile(image_out_folder, fname), '-dtiff', ['-r' num2str(dpi)])
end
%% write out the 3D angio stack for the single frame, ready for next stage

[Gnz,Gnx,Gny] = size(angio);
Izc = mean(mean(angio(:,round(Gnx/2)-5:round(Gnx/2)+4,round(Gny/2)-5:round(Gny/2)+4),2),3);
save(fullfile(scratch_folder,  [filename(1:end-4),'_Izc.mat']),'Izc')

fprintf(1, 'Saving angio...')

% just take the selected zrange
angio = angio(pl_zrange, :, :);

save(fullfile(scratch_folder,  [filename(1:end-4),'_angio.mat']),...
    'angio','pl_zrange','nii_zrange', 'zstack', '-v7.3');

fprintf(1, 'done\n')

%%
function varargout = add_slider_to_angio_plot(G, zstack, zinit, angio_zsize, clim)
global pl_zrange

slctr=uicontrol('style','slider','position',[50 60 20 200],...
    'min',1,'max',angio_zsize,'value', round(zinit), ...
    'sliderstep', [1, 10] / angio_zsize, 'callback',@callbackfn);

% hsttext=uicontrol('style','text',...
%     'position',[20 260 120 20], 'value', round(zinit), 'visible','on');

    function callbackfn(source,eventdata)
        center=floor(get(slctr,'value'));
        
        center = max(center, center + floor(zstack / 2));
        center = min(center, size(G, 1) - ceil(zstack / 2));
        
%         width=get(slwdt,'value');
        uistr = sprintf('z-center: %d', center);
%         set(hsttext,'visible','on','fontsize', 14, 'string',uistr)

%         subplot(1,3,1)
        pl_zstart = round(center - floor(zstack / 2));
        pl_zrange = pl_zstart:pl_zstart + zstack - 1;

        subplot(1,2,1);
        imagesc(squeeze(mean(G(pl_zrange, :, :), 1))); axis equal; axis tight
        set(gca,'CLim', clim)
        th = title(sprintf('z-center: %d', center));
        xh = xlabel('XY-plane (pixels)');
        set(th, 'FontSize', 14)
        set(xh, 'FontSize', 12)
        
%         subplot(1,3,2); hold on
        subplot(2,2,2); hold on
        children = get(gca, 'Children');
        if length(children) > 1
            delete(children(1))
            delete(children(2))
            
        end
        xlim = get(gca, 'XLim');
        lh(1) = plot(xlim, [pl_zrange(1), pl_zrange(1)], 'r', 'LineWidth', 2);
        lh(2) = plot(xlim, [pl_zrange(end), pl_zrange(end)], 'r', 'LineWidth', 2);
        hold off
%         imagesc(squeeze(mean(G(:, xrange, :), 2)))
%         set(gca,'CLim', [0, 1])
%         subplot(1,3,3); hold on
        subplot(2,2,4); hold on
        children = get(gca, 'Children');
        if length(children) > 1
            delete(children(1))
            delete(children(2))
        end
        xlim = get(gca, 'XLim');
        lh(3) = plot(xlim, [pl_zrange(1), pl_zrange(1)], 'r', 'LineWidth', 2);
        lh(4) = plot(xlim, [pl_zrange(end), pl_zrange(end)], 'r', 'LineWidth', 2);
        hold off
%         imagesc(squeeze(mean(G(:, :, yrange), 3)))
%         set(gca,'CLim', [0, 1])
        drawnow
        varargout{1} = pl_zrange;
    end
end