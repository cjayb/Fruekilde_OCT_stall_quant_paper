function plot_cleaning_progress_diagnostics(orig_img, destripe, frangi, ...
                clean_mip)

inc_titles = true;

nr = 2;
nc = 2;

figure(201); clf; colormap bone
set(gcf, 'Position', [ 39   524   645   633])
subplot(nr, nc, 1); imagesc(orig_img); axis equal; axis tight
tit = 'Original frame';
% colorbar
set(gca, 'XTick', []); set(gca, 'YTick', []);
if inc_titles, title(tit), end

% plot sequence of cleanup operations

subplot(nr, nc, 2); imagesc(destripe); axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
if inc_titles, title('destripe'), end

subplot(nr, nc, 3); imagesc(frangi); axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
if inc_titles, title('+ frangi'), end

subplot(nr, nc, 4); imagesc(clean_mip); axis equal; axis tight
set(gca, 'XTick', []); set(gca, 'YTick', []);
if inc_titles, title('+ histeq (maxnorm)'), end
