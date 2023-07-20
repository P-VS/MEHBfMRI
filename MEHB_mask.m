function mask = MEHB_mask(indata)

voldim = size(indata);
if numel(voldim)>1
    meandata = mean(indata,4);
else
    meandata = indata;
end

nmeandata = meandata(isfinite(meandata));

thr = opt_thr_corr(nmeandata);
   
mask = zeros(voldim(1),voldim(2),voldim(3));

tmp1 = find(isfinite(meandata));
tmp2 = meandata(tmp1)>thr;

mask(tmp1(tmp2)) = 1;

se = strel('sphere',1);
mask = imopen(mask,se);
mask = imclose(mask,se);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thr = opt_thr_corr(img)
costfunc = @(thr) -correlation(img, img > thr);
[thr ncc] = fminbnd(costfunc, min(img), max(img));
fprintf('Maximal correlation of %g found with threshold of %g\n', ...
    -ncc, thr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = correlation(x, y)
cs = corrcoef(x, double(y));
c = cs(1, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thr = opt_thr_antimode(img) %#ok, optfunc can eval as handle to it
% See appendix B of Luo and Nichols (2003), Neuroimage 19, 1014-1032
% http://www.sph.umich.edu/ni-stat/SPMd/SPMd.pdf
% http://dx.doi.org/10.1016/S1053-8119(03)00149-6
% Implemented from the description, any mistakes are my fault - Ged

% Hartigan method
srt = sort(img);
lo = srt(floor(numel(srt) * 0.1));
hi = srt(ceil(numel(srt) * 0.9));
srt(srt <= lo) = [];
srt(srt >= hi) = [];
dfs = diff(srt);
mx = max(dfs);
k = floor(mean(find(dfs == mx))); % in case non-unique
thr = mean([srt(k) srt(k+1)]);
fprintf('Anti-mode threshold of %g by Hartigan method\n', thr);

% Histogram method
iqrange = diff(srt(round(numel(srt)*[0.25 0.75])));
binwidth = 1.595 * iqrange * numel(srt)^(-1/5);
numbins = ceil((srt(end) - srt(1)) / binwidth);
[counts bins] = hist(srt, numbins);
mn = min(counts);
k = round(mean(find(counts == mn))); % in case non-unique
thr2 = bins(k);
fprintf('Anti-mode threshold of %g by histogram method\n', thr2);

% Assume Hartigan more accurate if similar, but less reliable if different
if abs(thr - thr2) > 3 * binwidth
    thr = thr2;
end
fprintf('Anti-mode threshold of %g chosen.\n', thr);