function Mv = wm2mv(A, B, varargin)
% Mv = wm2mv(skel_struct, worm_track_array, ...)
%
% Extract worm images from intermediate data worm_track_array (object tracking structure)
% into color movie Mv; 1st arg skel_struct (or just its source field) is required to get
% the movie frame size and number of frames. Other args work the same as for pl2mv:
%  - 'saveto' : by default movie is returned to memory, if this is given it is saved as an
%      AVI file to disk, and the saveto name is returned.
%  - 'framerange': use this to limit movie to specified frame range.
%  - 'worms' : indices of tracks in worm_track_array to include (by default all are).
%  - 'crop', [i1 i2 j1 j2]: restrict movie to subframe within rows i1 to i2 and columns j1 to j2. If
%      given as the string 'auto' then the crop box will be automatically calculated from the data.
%  - 'rescale', 0 < x < ...: rescale frames by specified amount. If a duple [m,n] >= 1 is given then
%      frames will be shrunk if they have more than m rows or n columns (preserving aspect ratio).
%  - 'colormap': use this to specify alternate colors for tracks. By default the standard line colors
%      are cycled.
%  - 'number', {T|F}: superimpose the track number on each object at its centroid location
%  - 'showframe', {T|F}: show the current frame number in the top left corner.

S = getargs(varargin);
if isfield(A, 'source') % don't need any other parts of the output, just source info and pillar measurements
    A = A.source;
end

% set up crossref between frames and data fields
frng = A.framerange; nA = diff(frng) + 1; kw = length(B); wlist = (1:kw)';
haswrm = zeros(nA, kw);
for k = wlist'
    ii = [B(k).data.frame]; haswrm(ii, k) = 1:length(B(k).data);
end

% apply restrictions across dimensions (frame range, worm selection, frame extents etc)
ofs_fr = 0;   % frame range ...
if ~isempty(S.framerange)
    frtmp = [max(S.framerange(1), frng(1)), min(S.framerange(2), frng(2))];
    ii = (frtmp(1) - frng(1) + 1):(frtmp(2) - frng(1) + 1); haswrm = haswrm(ii,:);
    frng = frtmp; ofs_fr = frng(1) - 1; nA = diff(frng) + 1;
end

% worm list
wlist = wlist(any(haswrm > 0, 1));    % remove tracks not present in requested frame range
if ~isempty(S.worms)    % at this point, check if user put restriction on pillars they want shown
    wlist = intersect(S.worms(:), wlist);
end
kw = length(wlist); haswrm = haswrm(:,wlist);

% % cropping...
bb0 = A.size; bb = bb0; docrop = false;
if ~isempty(S.crop)
    if isequal(S.crop, 'auto')   % autocrop extents are derived from bounding boxes of data elements
        b1 = [Inf(1,2); -Inf(1,2)];
        for i = 1:kw
            ii = haswrm(haswrm(:,i) > 0, i);
            btmp = reshape([B(wlist(i)).data(ii).ij], 4, length(ii))';
            b1(1,:) = min(b1(1,:), min(btmp(:,[1,3]),[],1));
            b1(2,:) = max(b1(2,:), max(btmp(:,[2,4]),[],1));
        end
    else
        b1 = [max(1, S.crop([1,3])); min(S.crop([2,4]), bb)];
    end
    if ~isequal(b1, [ones(1,2); bb])
        docrop = true; bb = diff(b1, 1, 1) + 1; xycr = b1(1,[2,1]) - 1;
        icr = b1(1,1):b1(2,1); jcr = b1(1,2):b1(2,2);
    end
end

% rescaling...
doresc = false;
if ~isempty(S.rescale)
    if numel(S.rescale) > 1
        r0 = min(min(S.rescale(:)./bb(:)), 1);
    else
        r0 = S.rescale;
    end
    if r0 ~= 1
        doresc = true; b0 = round(r0*bb);
        r2fact = fliplr(b0./bb); bb = b0;
    end
end

% set colormap and color assignments. We put a limit on the number of colors since we want to reserve some for other purposes
cmap = [1 1 1];
if ~isempty(S.colormap)
    nc = min(127, size(S.colormap, 1)); cmap = [cmap; S.colormap(1:nc,:)];
else
    nc = 7; cmap = [cmap; lines(7)];
end
wcols = uint8(mod((1:kw) - 1, nc) + 1);

doanynum = false; donumbers = false; doframeno = false;
if S.number
    donumbers = true; doanynum = true; numloc = cell(nA, kw);
    for i = 1:kw
        jj = find(haswrm(:,i) > 0); ii = haswrm(jj,i);
        ctmp = reshape([B(wlist(i)).data(ii).Centroid], 2, length(jj))';
        if size(ctmp, 1) > 2   % smooth trajectory so that numbers don't bump around too much
            ctmp = [csaps(ii, ctmp(:,1), 1/10, ii), csaps(ii, ctmp(:,2), 1/10, ii)];
        end
        if docrop
            ctmp = bsxfun(@minus, ctmp, xycr);
        end
        if doresc
            ctmp = bsxfun(@times, r2fact, ctmp - 1/2) + 1/2;
        end
        numloc(jj,i) = num2cell(ctmp, 2);
    end
end
if S.showframe
    doframeno = true; framestr = 'FRAME: '; doanynum = true;
end
if doanynum
    cmap = [cmap; 1 - cmap];  % use complementary colors of whatever is in background for numbering
end

dosave = false;
if ~isempty(S.saveto)
    dosave = true; Mv = S.saveto; Au = VideoWriter(Mv);
    Au.FrameRate = 10; open(Au);
else
    Mv = struct('cdata',cell(1,nA), 'colormap',cmap);
end
for i = 1:nA
    ifr = i + ofs_fr; bi = zeros(bb, 'uint8');
    for j = 1:kw
        k = haswrm(i,j);
        if k > 0
            wtmp = B(wlist(j)).data(k); btmp = false(bb0);
            ii = wtmp.ij(1):wtmp.ij(2); jj = wtmp.ij(3):wtmp.ij(4);
            btmp(ii,jj) = bwunpack(wtmp.Image, wtmp.BoundingBox(4));
            if docrop
                btmp = btmp(icr,jcr);
            end
            if doresc
                btmp = imresize(btmp, bb);
            end
            bi(btmp) = wcols(j);
        end
    end

    if doanynum
        btmp = zeros(bb, 'uint8'); ii = find(haswrm(i,:));
        if donumbers && ~isempty(ii)
            btmp = insertText(btmp, cell2mat(numloc(i,ii)'), wlist(ii), 'AnchorPoint','Center', 'BoxOpacity',0, ...
                'TextColor',[255 255 255], 'FontSize',15);
        end
        if doframeno
            btmp = insertText(btmp, [7.5 15], [framestr, num2str(ifr)], 'AnchorPoint','LeftCenter', 'BoxOpacity',0, ...
                'TextColor',[255 255 255], 'FontSize',15);
        end
        jj = find(btmp(:,:,1) > 127);
        bi(jj) = bi(jj) + nc + 1;
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

% -------------------------------------------------------------------------------------------------------------------- %
function S1 = getargs(C)
% Parse input arguments.

% user args. [] is used to flag that default is to be applied.
sf =   {'number','showframe','saveto','framerange','worms','colormap','crop','rescale'};
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

