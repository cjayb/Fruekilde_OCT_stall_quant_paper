% Convert spectrum to reflectivity profile
% RAM = FRG raw file X 4
% Provided by BOAS Lab, https://github.com/BUNPC

function RR = NII2RR(NII, intpDk, bavgfrm)

% choose parameter for lamda-k interpolation
if nargin < 2
	intpDk = -.23; %-.37 was the default value %-.22 Conrad
end
if nargin < 3
	bavgfrm = 0;
end

% substract the reference signal, Subtract mean
[nk,nx,nf] = size(NII);
if bavgfrm == 1
    NII = NII - repmat(mean(NII(:,:),2),[1 nx nf]);
else
    for ifr=1:nf
        NII(:,:,ifr) = NII(:,:,ifr) - repmat(mean(NII(:,:,ifr),2),[1 nx]);
    end
end

% lamda-k interpolation
if intpDk ~= 0
    k = linspace(1-intpDk/2, 1+intpDk/2, nk);
    lam = 1./fliplr(k);
    for ifr=1:nf
        NII(:,:,ifr) = interp1(lam, NII(:,:,ifr), linspace(min(lam),max(lam),length(lam)), 'spline');
        if (nf > 100) && (mod(ifr,ceil(nf/5)) == 0)  
            disp(['... NIItoRR ' num2str(ifr) '/' num2str(nf) '  ' datestr(now,'HH:MM')]);  
        end
    end			
end

% ifft
nz = round(nk/2);
RR = zeros(nz,nx,nf,'single');
for ifr=1:nf
    RRy = ifft(NII(:,:,ifr));
    RR(:,:,ifr) = RRy(1:nz,:);
end
