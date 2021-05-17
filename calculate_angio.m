function G = calculate_angio(nii)

nz = diff(nii.zrange) + 1;  % limited z-range
nx = nii.dims(2);  % original x-width
ny = nii.dims(3);  % original y-width

RR1 = ones(nz, nx/2, ny, 2) * (1+1i);  % 128-bit complex (2xdouble)
wh = waitbar(0, 'Correcting slice phase...');
for iy=1:ny
    curFrame = nii.readFrame(iy);
    RR1(:,:,iy,1) = curFrame(:, 1:nx / 2);
    RR1(:,:,iy,2) = curFrame(:, nx / 2 + 1:end);
    for ix = 1:nx/2
        RR1(:,ix,iy,:) = CorrSlicePhase(RR1(:,ix,iy,:),1,3);
    end
    waitbar(iy/ny, wh)
end
close(wh)

% difference along last dimension
G = abs(diff(RR1,1,4));
