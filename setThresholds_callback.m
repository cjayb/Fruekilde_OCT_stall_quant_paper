function setThresholds_callback(fig_obj, eventData)

global cap_id filt_edgelist bin_stalls stallogram eq_vessels frame_id
global unique_thresholds show_cap_overlay bad_frames

n_caps = size(stallogram, 1);
n_frames = size(stallogram, 2);

ck = get(fig_obj, 'CurrentKey');        

if strcmp(ck, 'downarrow')
    if cap_id < n_caps
        subplot(2,2,1)
        el = filt_edgelist{cap_id};
        plot(el(:,2), el(:,1), '.', 'Color', [1 1 1])
        cap_id = cap_id + 1;
        el = filt_edgelist{cap_id};
        plot(el(:,2), el(:,1), '.', 'Color', [1 0 0])
        
    end
elseif strcmp(ck, 'uparrow')
    if cap_id > 1
        subplot(2,2,1)
        el = filt_edgelist{cap_id};
        plot(el(:,2), el(:,1), '.', 'Color', [1 1 1])
        cap_id = cap_id - 1;
        el = filt_edgelist{cap_id};
        plot(el(:,2), el(:,1), '.', 'Color', [1 0 0])
    end
    
elseif strcmp(ck, 's')
    cur_threshold = unique_thresholds(cap_id);
    cur_threshold = cur_threshold - 0.01;
    unique_thresholds(cap_id) = min(max(cur_threshold, 0), 1);
elseif strcmp(ck, 'w')
    cur_threshold = unique_thresholds(cap_id);
    cur_threshold = cur_threshold + 0.01;
    unique_thresholds(cap_id) = min(max(cur_threshold, 0), 1);
elseif strcmp(ck, 'x')
    cur_threshold = unique_thresholds(cap_id);
    cur_threshold = cur_threshold - 0.1;
    unique_thresholds(cap_id) = min(max(cur_threshold, 0), 1);
elseif strcmp(ck, '2')
    cur_threshold = unique_thresholds(cap_id);
    cur_threshold = cur_threshold + 0.1;
    unique_thresholds(cap_id) = min(max(cur_threshold, 0), 1);
    
elseif strcmp(ck, 'g')
    gotoline = str2num(input('Go to capillary # ', 's'));
    if gotoline > 0 && gotoline <= n_caps
        subplot(2,2,1)
        el = filt_edgelist{cap_id};
        plot(el(:,2), el(:,1), '.', 'Color', [1 1 1])
        cap_id = gotoline;
        el = filt_edgelist{cap_id};
        plot(el(:,2), el(:,1), '.', 'Color', [1 0 0])
    end
    figure(fig_obj)    
elseif strcmp(ck, 'f')
    gotocol = str2num(input('Go to frame # ', 's'));
    if gotocol > 0 && gotocol <= n_frames
        frame_id = gotocol;
    end
    figure(fig_obj)    
elseif strcmp(ck, '1')
    show_cap_overlay = ~show_cap_overlay;
elseif strcmp(ck, 'space')    
    bad_frames(frame_id) = ~bad_frames(frame_id);
    
end

%%%%%
    
if strcmp(ck, 'downarrow') || strcmp(ck, 'uparrow') || strcmp(ck, 'g') || ...
        strcmp(ck, '2') || strcmp(ck, 'w') || strcmp(ck, 's') || strcmp(ck, 'x')

    bin_stalls(cap_id, :) = stallogram(cap_id, :) > unique_thresholds(cap_id);

    row_mask = zeros(size(bin_stalls), 'double'); row_mask(cap_id, :) = 1; 
    
    highlight_row = imoverlay(bin_stalls, row_mask.*bin_stalls, [1 0 0]);
    subplot(1,2,2); imagesc(highlight_row);
    hold on; plot([1 n_frames], [cap_id cap_id], 'g--', 'LineWidth', 0.5); hold off
    hold on; plot([frame_id frame_id], [1 n_caps], 'g--', 'LineWidth', 0.5); hold off
    title('Stallogram')
    xlabel('Time (frame no.)')
    ylabel('Capillary ID')

end

if strcmp(ck, 'rightarrow')
    if frame_id < size(stallogram, 2)
        frame_id = frame_id + 1;
    end
elseif strcmp(ck, 'leftarrow')
    if frame_id > 1
        frame_id = frame_id - 1;
    end
end
if strcmp(ck, 'rightarrow') || strcmp(ck, 'leftarrow') || strcmp(ck, 'f')
    subplot(2,2,1); hold off
    cur_frame = eq_vessels(:, :, frame_id);
    imagesc(cur_frame); hold on
    el = filt_edgelist{cap_id};
    inds = sub2ind(size(cur_frame), el(:, 1), el(:, 2));
    path_mask = zeros(size(cur_frame)); path_mask(inds) = 1;
    thick_solution_path = imdilate(path_mask, ones(3,3));
%     highlight_path = imoverlay(cur_frame, thick_solution_path, [1 0 0]);
    red = cat(3, ones(size(thick_solution_path)), ...
                 zeros(size(thick_solution_path)), ...
                 zeros(size(thick_solution_path)));
    h = image(red);
    if show_cap_overlay
        set(h, 'AlphaData', 0.5*thick_solution_path);
    else
        set(h, 'AlphaData', 0);
    end
    axis equal; axis tight
    set(gca, 'XTick', []); set(gca, 'YTick', [])
    title('Capillary network')

    subplot(1,2,2); delete(findobj(gca, 'type', 'line'))
    hold on; plot([1 n_frames], [cap_id cap_id], 'g--', 'LineWidth', 0.5); hold off
    hold on; plot([frame_id frame_id], [1 n_caps], 'g--', 'LineWidth', 0.5); hold off

end

subplot(2,2,3); plot(stallogram(cap_id, :));
set(gca, 'Ylim', [0 1]);
set(gca, 'XLim', [1 n_frames])
for ii = 1:length(bad_frames)
    if bad_frames(ii)
        h = xline(ii);
        set(h, 'LineWidth', 4, 'Color', [1 0 1], 'Alpha', 0.2)
    end
end
xline(frame_id);
yline(unique_thresholds(cap_id), 'Color', 'r');

set(gca, 'Ylim', [0 1]); set(gca, 'YTick', [0 1])
set(gca, 'YTickLabel', {'Flow', 'Stall'})
ytickangle(90)
title('Stalling dynamics for selected capillary')
xlabel('Time (frame no.)')

axes = get(fig_obj, 'Children');
for an = 1:length(axes)
    set(axes(an), 'FontSize', 12)
end
