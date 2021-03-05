function Mv = pl2mv(A, B, varargin)
% Mv = pl2mv(A, ...)
%
% Extract pillar tracking data from skel2n... output A into Matlab movie struct Mv. By default
% Mv will consist of the original grayscale source data, with tracked pillars superimposed in
% color (blue when in base position, red when deflected, no outline for frames where a measurement
% problem was indicated by the pillar field's "isprob" subfield). Indices to the element of the
% pillars field will also be displayed (in gray). Other arguments:
%  - 'filename' : by default program looks for original movie in location given in A,
%      use this if the file has been moved or renamed.
%  - 'saveto' : by default movie is returned to memory, if this is given it is saved as an
%      AVI file to disk, and the saveto name is returned.
%  - 'framerange': use this to limit movie to specified frame range.
%  - 'pillars' : indices of pillars in B to highlight (by default all are).
%  - 'crop', [i1 i2 j1 j2]: restrict movie to subframe within rows i1 to i2 and columns j1 to j2. If
%      given as the string 'auto' then the crop box will be automatically determined by the pillars
%      that are to be highlighted.
%  - 'rescale', 0 < x < ...: rescale frames by specified amount. If a duple [m,n] >= 1 is given then
%      frames will be shrunk if it has more than m rows or n columns (preserving aspect ratio).
%  - 'colormap': use this to specify alternate colors for deflected and base modes.
%  - 'textcolor': specify alternate grayscale color (between 0 and 255) for text; use when
%      background luminance is too close to gray 127 (which will make numbers hard to read)
%  - 'showframe', {T|F}: show the current frame number in the top left corner; can alternately
%      be given as a specific position in the (original) frame in [x,y] coords.
%  - 'thicken', {T|F}: if set to true, pillar outlines are thickened (more visible, but relation to
%      image underneath may be obscured).

if nargin < 2
    B = A.pillars; S = getargs({});
elseif ~isstruct(B)
    S = getargs(cat(2, B, varargin)); B = A.pillars;
else   % B is a pillar array, use that
    S = getargs(varargin);
end
if isfield(A, 'source') % don't need any other parts of the output, just source info and pillar measurements
    A = A.source;
end

% collect the pillar base and measurement data from B (after which we don't need B anymore)
frng = A.framerange; nA = diff(frng) + 1; kp = length(B); plist = (1:kp)';
ctr_base = reshape([B.center], 2, kp)'; rad_base = [B.radius]';
hasdefl = false(nA, kp); defl_circ = cell(nA, kp);
for k = plist'
    if ~isempty(B(k).defl)
        ii = [B(k).defl.frame]; hasdefl(ii, k) = true;
        defl_circ(ii, k) = {B(k).defl.center};
    end
end
usebase = ~hasdefl;
if isfield(B, 'isprob')   % if there were in the run problems unmark the affected frames (so pillar outline will disappear in movie at these points)
    for k = find(~cellfun('isempty',{B.isprob}))
        ii = [B(k).isprob.frame]; usebase(ii,k) = false;
    end
end
if ~isempty(S.pillars)    % at this point, check if user put restriction on pillars they want shown
    plist = intersect(S.pillars(:), plist); kp = length(plist);
    ctr_base = ctr_base(plist, :); rad_base = rad_base(plist); usebase = usebase(:,plist);
    defl_circ = defl_circ(:,plist); hasdefl = hasdefl(:,plist);
end
clear B;

% now apply any other restrictions (framerange, cropping, etc) to the pillar data
ofs_fr = 0;
if ~isempty(S.framerange)
    frtmp = [max(S.framerange(1), frng(1)), min(S.framerange(2), frng(2))];
    ii = (frtmp(1) - frng(1) + 1):(frtmp(2) - frng(1) + 1);
    defl_circ = defl_circ(ii,:); usebase = usebase(ii,:); hasdefl = hasdefl(ii,:);
    frng = frtmp; ofs_fr = frng(1) - 1; nA = diff(frng) + 1;
end

% apply cropping + init crop vars.
bb = A.size; docrop = false;
if ~isempty(S.crop)
    if isequal(S.crop, 'auto')
        rad_tmp = 2*max(rad_base);
        b1 = max(1, round(fliplr(min(ctr_base - rad_tmp, [], 1))));
        b1 = [b1; min(bb, round(fliplr(max(ctr_base + rad_tmp, [], 1))))];
    else
        b1 = [max(1, S.crop([1,3])); min(S.crop([2,4]), bb)];
    end
    if ~isequal(b1, [ones(1,2); bb])
        docrop = true; bb = diff(b1, 1, 1) + 1;
        icr = b1(1,1):b1(2,1); jcr = b1(1,2):b1(2,2);
        xycr = b1(1,[2,1]) - 1;
        for i = 1:kp
            ctr_base(i,:) = ctr_base(i,:) - xycr;
            for j = find(hasdefl(:,i)')
                defl_circ{j,i} = defl_circ{j,i} - xycr;
            end
        end
    end
end

% apply rescale + init scale var. Since we need to predict exactly where features in each image will end up we use the
% more predictable [numrows, numcols] interface for ML's imresize function (even though we are given a single scale factor)
doresc = false;
if ~isempty(S.rescale)
    if numel(S.rescale) > 1
        r0 = min(min(S.rescale(:)./bb(:)), 1);
    else
        r0 = S.rescale;
    end
    if r0 ~= 1
        doresc = true; b0 = round(r0*bb);
        r2fact = fliplr(b0./bb); r1fact = mean(r2fact); bb = b0;
        for i = 1:kp
            ctr_base(i,:) = r2fact.*(ctr_base(i,:) - 1/2) + 1/2;   % in ML images everything is offsect from pt {1/2,1/2}
            rad_base(i) = r1fact*rad_base(i);    % proportions in r2fact should be close enough that we can still draw
            for j = find(hasdefl(:,i)')          % circles rather than ellipses
                defl_circ{j,i} = r2fact.*(defl_circ{j,i} - 1/2) + 1/2;
            end
        end
    end
end

% check other user args, if given arg vals by user override the defaults
fnam = A.filename;
if ~isempty(S.filename)
    fnam = S.filename;
end
cmap = gray(256); cmap(1:2,:) = [1 0 0; 0 0 1];
if ~isempty(S.colormap)
    cmap(1:2,:) = S.colormap;
end
txcol = [127 127 127]; txsz = min(72, max(8, ceil(min(rad_base))));
if ~isempty(S.textcolor)
    txcol = repmat(max(S.textcolor(1), 2), 1, 3);
end
doframeno = false;
if ~isequal(S.showframe, false)
    doframeno = true; framestr = 'FRAME: ';
    if ~isequal(S.showframe, true)
        framepos = S.showframe;
        if docrop
            framepos = framepos - xycr;
        end
        if doresc
            framepos = r2fact.*(framepos - 1/2) + 1/2;
        end
    else
        framepos = txsz*[1/2, 1];
    end
end
dothick = false;
if S.thicken
    dothick = true; str1 = strel('disk', 1);
end

% calculate indices for the base and deflected circles (note: won't need centers of deflected circles after this, so
% defl_circ contents can just be replaced)
base_circ = cell(1, kp);
for i = 1:kp
    base_circ{i} = m3circ(ctr_base(i,:), rad_base(i), bb)';
    for j = find(hasdefl(:,i)')
        defl_circ{j,i} = m3circ(defl_circ{j,i}, rad_base(i), bb)';
    end
end

% make the movie, (if saving avi file to disk do setup for that now too)
Av = VideoReader(fnam); dosave = false;
if ~isempty(S.saveto)
    dosave = true; Mv = S.saveto; Au = VideoWriter(Mv);
    Au.FrameRate = A.avinfo.FrameRate; open(Au);
else
    Mv = struct('cdata',cell(1,nA), 'colormap',cmap);
end
for i = 1:nA
    ifr = i + ofs_fr; bi = read(Av, [ifr, ifr]);
    if docrop
        bi = bi(icr,jcr);
    end
    if doresc
        bi = imresize(bi, bb);
    end

    bi = insertText(bi, ctr_base, plist, 'AnchorPoint','Center', 'BoxOpacity',0, 'TextColor',txcol, 'FontSize',txsz);
    bi = bi(:,:,1); bi(bi < 2) = 2;
    bi([defl_circ{i,hasdefl(i,:)}]) = 0; bi([base_circ{usebase(i,:)}]) = 1;

    if dothick
        bi(imdilate(bi == 0, str1)) = 0; bi(imdilate(bi == 1, str1)) = 1;
    end
    if doframeno   % need to do this after circles are drawn in case it lands on top of any
        bi = insertText(bi, framepos, [framestr, num2str(ifr)], 'AnchorPoint','LeftCenter', 'BoxOpacity',1, ...
            'BoxColor',[255 255 255], 'TextColor',[2 2 2], 'FontSize',txsz); bi = bi(:,:,1);
    end
    if dosave
        writeVideo(Au, ind2rgb(bi, cmap));
    else
        Mv(i).cdata = bi;
    end
end
if dosave
    close(Au);
end

return;

% ---------------------------------------------------------------------------------------%
function ij2 = m3circ(xy, rd, imext)
% ij2 = m3circ(xy, rd, imext)
% Calculates points on perimeter of a circle so that when they are rounded there are no
% gaps, then converts to indices into matrix with size imext.

x = xy(1); y = xy(2);
ix = ceil(x - rd):floor(x + rd); jy = ceil(y - rd):floor(y + rd);
xy1 = [ix', y + sqrt(rd^2 - (ix' - x).^2); ...
       x + sqrt(rd^2 - (fliplr(jy)' - y).^2), fliplr(jy)'; ...
    fliplr(ix)', y - sqrt(rd^2 - (fliplr(ix)' - x).^2); ...
       x - sqrt(rd^2 - (jy' - y).^2),  jy'];
ij1 = unique(round(xy1), 'rows');

ij1 = ij1((ij1(:,2) >= 1 & ij1(:,2) <= imext(1)) & (ij1(:,1) >= 1 & ij1(:,1) <= imext(2)),:);
ij2 = sub2ind(imext, ij1(:,2), ij1(:,1));

return;

% -------------------------------------------------------------------------------------------------------------------- %
function S1 = getargs(C)
% Parse input arguments.

% user args. [] is used to flag that default is to be applied.
sf =   {'thicken','showframe','filename','saveto','framerange','pillars','colormap','textcolor','crop','rescale'};
indic = 2;
S1 = cell2struct(cell(length(sf),1), sf);

nrgs = length(C); i = 1;
while i <= nrgs
    j = find(strncmpi(C{i}, sf, length(C{i})));
    if length(j) < 1
        error(['Unknown argument string "' C{i} '"']);
    elseif length(j) > 1
        error(['Cannot resolve argument string "' C{i} '" between "' sf{j(1)} '", "' sf{j(2)} '" ...']);
    end
    if j <= indic && (i == nrgs || ischar(C{i+1}))
        S1.(sf{j}) = true; i = i + 1;
    else
        S1.(sf{j}) = C{i+1}; i = i + 2;
    end
end
for i = 1:indic   % check indicators, if not set to any value then give the default.
    if isempty(S1.(sf{i}))
        S1.(sf{i}) = false;
    end
end

return;

