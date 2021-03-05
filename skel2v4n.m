function a2 = skel2v4n(Ax, varargin)
% output_data = skel2v4n(source_data, ...)
%
% Streamlined skel program with pillar management. "source_data" will generally be a
% filename for an avi file to be processed, though other primary inputs are supported for
% debugging etc. purposes. Remaining arguments are string-value pairs (full descriptions
% to be added when program complete):
%   'verbose'  (user interface)
%   'framerange', 'strip', 'threshold', 'medianfilter'  (file and basic image processing)
%   'magnification', 'pixelsize', 'wormstage' (experimental parameters)
%   'diameter', 'spacing', 'tolspace'  (pillar grid parameters)

%   'verbose', (optional logical value): if true then report current stage of computation,
%      with runtime info and progress dots when .
%   'framerange', [1x2 int]: by default the entire movie is processed if possible, use
%      this to process only frames in the specified range.
%   'strip', [kx2 subscript list]: strip specified pixels (by assigning a neutral value)
%      from all frames before further processing. Use this to prevent embedded time-codes
%      or hot pixels from affecting the results.
%   'magnification', 'pixelsize': experimental parameters. Needed if pillar coords are
%      going to be given. Default (for now) is 1 for both.
%   'diameter','spacing': pillar parameters (in microns). If a diameter is given then
%      will be checked for the presence of pillars, and if a spacing is given then the
%      pillar placement will be checked for a grid structure. Detected pillars will be
%      removed from tracked objects after the fact, and their deflections measured.
%  More to come ...

% First check if an ordinary run (filename to final output) or starting in debug mode with output from previous run.
% Throughout the var holding the main intermediate output will be Ax; at points there will be a 2nd major intermediate
% output, Px, holding results related to pillar management, plus occasionally a 3rd, Qx; in debug mode these are
% included in the output struct in the fields 'Ax', 'Px', and 'Qx'
uargs = getargs(varargin); dbstart = 0;
if isstruct(Ax) && isfield(Ax, 'debug')
    a2 = rmfield(Ax, {'debug', 'Ax'}); dbstart = Ax.debug; Ax = Ax.Ax;
    if isfield(a2, 'Px')
        Px = a2.Px; a2 = rmfield(a2, 'Px');
    end
    if isfield(a2, 'Qx')
        Qx = a2.Qx; a2 = rmfield(a2, 'Qx');
    end
elseif ~ischar(Ax)
    error('Incorrect input; see program help.');
else
    a2 = struct('date',date);
end

% if debug level is set less than 1 then 
if uargs.debug < 1
    a2 = uargs; a2.auxiliary = Ax; return;
end
if uargs.verbose
    t1 = cputime; t2 = t1; 
end

% if given filename, read the file into memory
if dbstart < 1
    doprn = uargs.verbose && ischar(Ax);
    if doprn
        fprintf(1, '1. Reading file %s ...', Ax);
    end
    [Ax, a2.source] = getim2(Ax, uargs.framerange, uargs.strip);
    if isempty(Ax)
        fprintf(1, ' Problem encountered; exiting.\n'); return;
    elseif doprn
        t3 = t2; t2 = cputime; fprintf(1, ' done (%f sec)\n', t2 - t3);
    end
end
msz = a2.source.size; nfr = diff(a2.source.framerange) + 1;

% verify the experimental params + check against user input, then convert to pixel distances
a2.params = chkexparms(a2, uargs); resc = a2.params.mic2pix;
rad = a2.params.inpix.radius; sp = a2.params.inpix.grdspc; wmln = a2.params.inpix.lensc;
if uargs.debug == 1
    a2.debug = 1; a2.Ax = Ax; return;
end

% apply a median filter to the data, plus do threshold + background etc. calculations
if dbstart < 2
    if uargs.verbose
        fprintf(1, '2. Applying median filter and doing threshold calculations ...');
    end
    [Ax, bbg, bfg, thr, md] = procim2(Ax, uargs.medianfilter);
    a2.improc = struct('medfilt',md, 'Bg',bbg, 'Bf',bfg, 'thr',thr);
    if uargs.verbose
        t3 = t2; t2 = cputime; fprintf(1, ' done (%f sec)\n', t2 - t3);
    end
else
    bbg = a2.improc.Bg; bfg = a2.improc.Bf; thr = a2.improc.thr;
end
if ~isempty(uargs.threshold)
    thr = uargs.threshold; a2.improc.thr = thr;
end
if uargs.debug == 2
    a2.debug = 2; a2.Ax = Ax; return;
end

% next stage: logical mask for trackable objects + get pillar array and background objects that need to be monitored
tsp = uargs.tolspace;
if isempty(sp)
    tsp = resc*tsp;
end
if dbstart < 3
    if uargs.verbose
        fprintf(1, '3. Obtaining object mask and pillar array ...');
    end
    [Afg, obg, pchk, pgd] = initg2(bfg, bbg, thr, wmln, rad, sp, tsp);
    a2.array = struct('mask',Afg, 'ptrk',pchk, 'bgobj',obg, 'grid',pgd);
    if isempty(Afg)    % if no trackable objects found return info struct so user can check foreground, threshold, etc
        fprintf(1, ' Problem encountered; exiting.\n'); uargs.debug = 3;
    elseif uargs.verbose
        t3 = t2; t2 = cputime; fprintf(1, ' done (%f sec)\n', t2 - t3);
    end
else
    Afg = a2.array.mask; pchk = a2.array.ptrk; pgd = a2.array.grid; % obg = a2.array.bgobj;
end
if uargs.debug == 3
    a2.debug = 3; a2.Ax = Ax; return;
end
clear bfg bbg

% now scan the movie for trackable objects. There is no auxiliary field for this op, but after this point full movie no
% longer needed; Ax now becomes (temporarily) a cell array of regionprops structs
if dbstart < 4
    if uargs.verbose
        fprintf(1, '4. Scanning movie for trackable objects ...');
    end
    Ax = getobj2(Ax, thr, Afg, 1);    % last arg says to pack large subimages if encountered
    if uargs.verbose
        t3 = t2; t2 = cputime; fprintf(1, ' done (%f sec)\n', t2 - t3);
    end
end
if uargs.debug == 4
    a2.debug = 4; a2.Ax = Ax; return;
end
clear Afg

% tracking stage. Since we potentially have objects of interest a couple of orders of magnitude smaller than worms
% the tracker returns 2 arrays of tracked objects: Ax (worms and other large, presumably mobile, objects) and
% Px (pillars and other relatively small but non-trivial, presumably static, objects). For this sorting the tracker
% requires worm length and pillar diameter for initial categorization.
if dbstart < 5
    if uargs.verbose
        fprintf(1, '5. Sorting found objects into tracks ...');
    end
    [Ax, Px, a2.track] = trkobj2(Ax, wmln, 2*rad);
    if iscell(Ax)   % output type flags tracking prob. This won't be exactly the same cell array that was given as input
        fprintf(1, ' Problem encountered; exiting.\n'); uargs.debug = 5;
    elseif uargs.verbose
        t3 = t2; t2 = cputime; fprintf(1, ' done (%f sec)\n', t2 - t3);
    end
end
if uargs.debug == 5
    a2.debug = 5; a2.Px = Px; clear Px; a2.Ax = Ax; return;
end

% check the medium scale static objects for consistency with previous pillar location + radius estimates; when objects
% are verified to be pillars, refine the pillar center and radius estimate using the best subset of tracked image objects.
% 2nd return from obj2plr are leftover objects (either non-pillars, or pillar images containing extraneous material).
if dbstart < 6
    if uargs.verbose
        fprintf(1, '6. Checking tracked image objects against pillar array ...');
    end
    [Px, Qx] = obj2plr(msz, nfr, rad, Px, pchk, pgd);
    if uargs.verbose
        t3 = t2; t2 = cputime; fprintf(1, ' done (%f sec)\n', t2 - t3);
    end
end
if uargs.debug == 6
    a2.debug = 6; a2.Px = Px; a2.Qx = Qx; clear Px Qx; a2.Ax = Ax; return;
end
clear pchk pgd

% Do the deflection measurements, + remove pillar outlines from worm images. As well, incorporate leftover images in Qx
% with the worm images when they are close enough to the worm body (some of these may include tail or head sections
% disconnected from main worm by thresholding etc).
if dbstart < 7
    if uargs.verbose
        fprintf(1, '7. Doing deflection measurements + 1st worm image revision ...');
    end
    [Ax, Px, Qx] = trkplr(msz, nfr, rad, Ax, Px, Qx);
    if uargs.verbose
        t3 = t2; t2 = cputime; fprintf(1, ' done (%f sec)\n', t2 - t3);
    end
end
if uargs.debug == 7
    a2.debug = 7; a2.Px = Px; a2.Qx = Qx; a2.Ax = Ax; return;
end

a2.pillars = Px; clear Px Qx
for i = 1:length(Ax)
    Ax(i).data = rmfield(Ax(i).data, 'Image');
end
a2.wrmtmp = Ax; clear Ax

% 8. Revise tracks (if multiple worm tracks), revise pillar / worm crossref if needed, get 1st est pillar / worm touch pts
%  need obg at this point?

% 9. clean up worm outlines, get perimeters, revise pillar / worm touch pts

% 10. get worm skeletons

if uargs.verbose
    t2 = cputime; fprintf(1, 'Computations completed (%f seconds total)\n', t2 - t1);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, iminf] = getim2(A, frrange, stripx)
% [image_data, image_info] = getimdata(image_source, frame_range, strip_pixel_list)
%
% Simplified image retrieval to usable greyscale Matlab data.
% "image_source" is a string that designates an avi movie. Additional arguments:
%   - frame_range [1x2 int]: Process only frames within the given range (inclusive)
%   - strip_pixel_list [kx2]: list of matrix indices to pixels that we may want to set to
%    a neutral value (these may be a time-code embedded in each frame or known hot pixels
%    for the camera).
% Image data is returned as a negative image (darker pixels => higher values) 3D uint8
% array. 2nd return image_info is a struct containing basic information about the images.
%
% Future versions may (if needed) have support for large movies and various color etc
% movie formats. For now movies must be 8-bit grayscale requiring <= 400MB of storage.

% First, initialize the VideoReader and get some info about the avi
if ischar(A)
    d2 = VideoReader(A); d3 = get(d2);
    numf = d3.Duration * d3.FrameRate; ht = d3.Height; wd = d3.Width;
    iminf = struct('filename',fullfile(d3.Path, d3.Name), 'framerange',[], 'numframes',numf, 'size',[ht, wd], ...
        'stripixels',[], 'avinfo',d3);

    if nargin < 2 || (isempty(frrange) || (frrange(1) <= 1 && frrange(2) >= numf))
        ij = {}; iminf.framerange = [1 numf];
    else
        ij = {[frrange(1), min(numf, frrange(2))]}; iminf.framerange = ij{1};
    end

% At this point we are restricting ourselves to 8-bit grayscale movies requiring <= 400MB
% memory (this will probably change in the future).
    if ~(isequal(d3.VideoFormat, 'Mono8') || (isequal(d3.VideoFormat, 'Grayscale') && d3.BitsPerPixel == 8))
        disp(' '); warning('Only 8-bit grayscale avi data format currently supported. See output for more details.');
        A = []; return;   % alt (see ML docs) = RGB24, Indexed, Mono16, RGB48, and signed versions of Mono and RGB
    end
    A = squeeze(255 - read(d2, ij{:}));
else   % A is given as a 3D array of uint8's in memory
    A = 255 - A; [n1, n2, n3] = size(A);
    iminf = struct('filename',[], 'framerange',[1 n3], 'numframes',n3, 'size',[n1, n2], 'stripixels',[], 'avinfo',[]);
end

% if requested, strip specified pixels
if nargin > 2 && ~isempty(stripx)
    iminf.stripixels = stripx;
    na = size(A); a_mn = min(A, [], 3);    % get replacement value = min of unaffected area of each frame across entire movie
    mina = min(a_mn(setdiff((1:prod(na(1:2)))', sub2ind(na(1:2), stripx(:,1), stripx(:,2)))));
    ii = repmat(stripx(:,1), 1, na(3)); jj = repmat(stripx(:,2), 1, na(3));    % get indices to affected pixels across a and assign
    kk = repmat(1:na(3), size(stripx, 1), 1);
    A(sub2ind(na, ii(:), jj(:), kk(:))) = mina;
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [A, amn, amx, thr, mfsz] = procim2(A, mfsz, prndots)
% [preproc1_data, foreground, background, threshold, medfilt] = procim2(image_data, threshold)
%
% Apply median filter to data and then calculate background, foreground, and threshold.

if nargin < 2 || isempty(mfsz)
    mfsz = [3, 3];
elseif any(mfsz == 0)
    mfsz = [];
end
if nargin < 3
    prndots = 0;
end

[n1, n2, n3] = size(A);
if ~isempty(mfsz)
    for i = 1:n3
        A(:,:,i) = medfilt2(A(:,:,i), mfsz);
        if prndots > 0 && ~mod(i,prndots)
            fprintf(1,'.');
        end
    end
end

amn = min(A, [], 3); amx = max(A, [], 3);
thr = 255 * graythresh(A);
% if over 20% of the background thresholds to true then presumably there is something up (ie. fuzz, crud, or a large
% light or dark region). If this is due to a large dark region we don't want the threshold to change (much), but
% otherwise we want to bring it up, and using the foreground to revise the value will have the desired effect
if nnz(amn > thr) > 0.2*n1*n2
    thr = 255 * graythresh(amx);
end
if prndots > 0
    fprintf(1,'.');
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [Amsk, bgobj, plist, pgrid] = initg2(Afg, Abg, thr, trklen, prad, pspace, tolg)
% [A_mask, bkg_objects, chk_pillars, grid_info] = initg2(A_foreground, A_background, threshold,
%                           char_length, pillar_radius, pillar_spacing, grid_tolerance)
%
% Create a bitmask using the foreground and background images which can be applied to the
% movie to reduce the number of found objects, + check the background for a pillar array
% plus any large objects that need to be taken into account. Arguments:
%   A_foreground, A_background, threshold: foreground and background images with the
%     threshold that will be applied to the movie as a whole.
%   char_length: characteristic length scale (in pixels), should roughly equal length of
%     youngest worm in the movie. Default = 1000 (assuming adult with 1:1 rescaling).
%   pillar_radius: pillar radius (in pixels). Can be a range. If given then a pillar
%     list will be constructed by performing a Hough transform on the background image.
%     If given as a singleton then a 10% tolerance will be applied to the Hough transform;
%     to override give a range with identical values.
%   pillar_spacing: distance between pillar centers (in pixels). If this is given then
%     circles will be uniquely assigned to grid positions; conflicts etc will be discarded
%     from the pillar list.
%   grid_tolerance: amount of play (in pixels) allowed in grid construction. Default = 10%
%     of spacing. If no spacing is given (ie random pillar positions) then use this to set
%     a minimum inter-pillar spacing.
% Returns:
%   A_mask: a logical array the same size of a movie frame indicating where trackable
%     activity is taking place.
%   bkg_objects: a list containing large background objects (such as shadows) that the
%     foreground object(s) may interact with. Format: regionprops struct array.
%   chk_pillars: a list of pillars which impinge on the trajectory of the trackable
%     object. Contains fields giving initial estimates for center, radius, and other info.
%   grid_info: structure giving info about the grid (eg. rotation wrt view, extents) +
%     the full list of found pillars, for reference. Note: numeric id's for tracked pillars
%     will be based on the full list.

if nargin < 3 || isempty(thr)
    thr = 255*graythresh(Afg);
end
if nargin < 4 || isempty(trklen)
    trklen = 1000;    % this really should be provided, but if not, assume 1 micron per pixel
end
bgsz = trklen^2/16;   % objects big enough to consider tracking (approx area of adult worm)
if nargin < 5
    prad = [];
end
if nargin < 6
    pspace = [];
end
if nargin < 7 || isempty(tolg)
    tolg = 0.1*pspace;
end

% create primary version of mask. function imkeepobj keeps connected regions at least as large as specified area
Amsk = imkeepobj(Afg > thr, bgsz); bgobj = []; plist = []; pgrid = [];
if nnz(Amsk) == 0;
    Amsk = []; disp(' '); warning('No large enough trackable foreground objects found.'); return;
elseif isempty(Abg)
    return;
end

% check for large background objects (don't need as large as foreground, just large enough to worry about). Once these
% are removed any small flotsam resulting in the mask is removed as well.
tdiff = median(Afg(:) - Abg(:)); atmp = imkeepobj(Abg > (thr - tdiff), bgsz/4);
if nnz(atmp & Amsk) > 0
    bgobj = regionprops(imkeepobj(atmp, Amsk), {'Area','Centroid','BoundingBox','Image'});   % only return bg obj's that
    [bgobj.ij] = bnd2idx({bgobj.BoundingBox});                                               % affect the fg obj's trajectory
    Amsk = imkeepobj(Amsk & ~atmp, bgsz/16);    % save prev Amsk for checks below ...
% % depricated below:
%     aold = Amsk;   
% % want to add a bit of insurance to prevent cutting objects in half, losing pillars etc ...
%     ad1 = Amsk & imdilate(atmp, strel('disk', 1));   % boundary between the retained mask and removed background objects
%     dval1 = ceil(trklen/30); ad2 = imheritconn(aold, Amsk, dval1);    % reconnect object(s) that may have been split etc.
%     Amsk = Amsk | imdilate(ad1 | ad2, strel('disk', dval1));          % (dval1 = about 1/2 target worm width)
end

% If a pillar radius is provided, do the Hough transform on the thresholded background to find pillar candidates.
% The thresholded image is skeletonized to reduce the effect of rim thickness on the procedure, then dilated a bit
% to make the circles more prominent. If large background objects exist these are removed first, since they would
% skeletonize into weird arrays of lines, but the objects are first eroded a bit in order to allow circles on the
% periphery to come out.
if isempty(prad)
    return;
elseif length(prad) < 2
    prad = prad + 0.1*prad*[-1,1];
elseif prad(1) == prad(2)    % if caller really wants it, no tol version is also available
    prad = prad(1);
end
a2mp = Abg > thr & ~atmp; % ~imerode(atmp, strel('disk', 2));
plist1 = chuff(imdilate(bwmorph(a2mp, 'skel', Inf), strel('disk', 1)), prad);   % chuff = std hough transform ...
if isempty(plist1)
    disp(' '); warning('No pillars with specified radius found.'); return;
end

% grid check. If no spacing given this will still check whether pillars are too close to each other
pgrid = pgrid2(plist1, pspace, size(Abg), tolg);
if isempty(pgrid) || isempty(pgrid.plist)
    disp(' '); warning('Could not detect pillar array with specified radius and spacing.'); return;
elseif pgrid.isprob > 0
    disp(' '); warning('Issues encountered during array detection; check output.'); return;
end

% if we have a valid grid, check pillar locations against mask and select out those that interact with the trackable
% object(s). Function mkmask creates a logical mask for a circle with a given center and radius which can be safely
% checked against a subregion of the main image. If the pillar does touch the object we augment the mask so that we can
% detect everything in the immediate area.
plist = [];
for i = 1:length(pgrid.plist)
    tmpl = pgrid.plist(i); [a3mp, ij] = mkmask(tmpl.center, tmpl.radius, size(Amsk), -1);
    ii = ij(1):ij(2); jj = ij(3):ij(4);   % ij gives extents of pillar mask wrt Amsk
    if nnz(Amsk(ii,jj) & a3mp) > 0
        Amsk(ii,jj) = Amsk(ii,jj) | imdilate(a3mp, strel('disk', round(tmpl.radius)));
        tmpl.imij = ij; tmpl.pnum = i; plist = [plist; tmpl]; % tmpl.msk = a3mp;
    end
end
Amsk = imfill(imdilate(Amsk & ~atmp, strel('disk',1)), 'holes');

return;

% -------------------------------------------------------------------------------------- %
function bw2 = imkeepobj(bw1, szlim)
% Works like bwmorph 'clean' except that it is possible to set an arbitrary object size;
% returns the image with only objects larger than szlim retained. If szlim is -ve works in
% the opposite direction (removing objects bigger than -szlim).
%   Alternately, if the 2nd argument is a binary image, then objects that intersect with
% any non-zero region of this reference are retained.

if nargin < 2 || (isempty(szlim) || isequal(szlim, 0))
    bw2 = bw1; return;
elseif isequal(szlim, 1)
    bw2 = bwmorph(bw1, 'clean'); return;   % for the actual clean op this is more efficient
end
bw2 = false(size(bw1)); ob1 = regionprops(bw1, {'Area','PixelIdxList'});
if isequal(size(szlim), size(bw1))
    nob = length(ob1); kk = false(1, nob);
    for k = 1:nob
        kk(k) = any(szlim(ob1(k).PixelIdxList));
    end
    kk = find(kk);
elseif szlim > 0
    kk = find([ob1.Area] > szlim);
else
    kk = find([ob1.Area] < -szlim);
end
if isempty(kk)
    return;
end
for k = kk
    bw2(ob1(k).PixelIdxList) = true;
end

return;

% -------------------------------------------------------------------------------------- %
function A_conn = imheritconn(A1, A2, dstol)
% Determine whether objects in A2, which are subobjects of those in A1, may be connected
% based on whether they are close enough to each other and whether they belonged to the
% same A1 object. As well, check if any of the objects in A1 extend to the image borders
% based on whether the containing object in A1 extended to the border. Returns: an image
% of the same size with the minimal objects needed to re-establish any probable
% connections. 3rd arg dstol is the distance tolerance.

mn = size(A1); 
B1 = bwlabel(A1); C1 = bwlabel(A2); n2 = max(C1(:)); 

% get the parent relationship
r2par = zeros(1, n2);
for j = 1:n2
    ii = B1(C1 == j); r2par(j) = ii(1);
end

% check the border. function subconn checks pairs of regions for points within given tolerance
klist = []; bflag = false(1, 2); B1(2:(mn(1)-1),2:(mn(1)-1)) = 0;
for i = 1:n2
    ktmp = subconn(C1 == i, B1 == r2par(i), 3*dstol, dstol, true);
    bflag(1) = bflag(1) | ~isempty(ktmp); klist = [klist; ktmp];
end

% check pairs of regions
for i = 1:n2
    for j = (i+1):n2
        if r2par(i) == r2par(j)
            ktmp = subconn(C1 == i, C1 == j, 5*dstol, dstol, false);
            bflag(2) = bflag(2) | ~isempty(ktmp); klist = [klist; ktmp];
        end
    end
end

nk = size(klist, 1); A_conn = false(mn);
if any(bflag)
    [i1,j1] = ind2sub(mn, klist(:,1)); [i2,j2] = ind2sub(mn, klist(:,2));
    mx = max(abs(i1 - i2), abs(j1 - j2)) + 1;
    for i = 1:nk
        ii = round(linspace(i1(i), i2(i), mx(i))); jj = round(linspace(j1(i), j2(i), mx(i)));
        A_conn(sub2ind(mn, ii, jj)) = true;
    end
end

return;
    
% ---------------------------------------------------------------------------------------%
function kk = subconn(c1, c2, dtol, ntol, onbdry)
% subfunction for above, check not only whether 2 areas are within dtol of one another, but
% also whether there are possible multiple connections. c1 and c2 are bitmasks for the two
% regions, dtol is the distance tolerance, ntol is the population tolerance.

% initial one-sided check, if nothing turns up we don't need to do more
kk = []; [d1, id1] = bwdist(c1); ok2 = c2 & d1 <= dtol;
if nnz(ok2) < 1
    return;
end

% check other side and preproc. If c2 is boundary of parent region it is flagged to get special treatment
[d2, id2] = bwdist(c2); ok1 = c1 & d2 <= dtol;
L1 = bwlabel(ok1); r1 = regionprops(L1, {'Area', 'PixelIdxList'}); q1 = [r1.Area] >= ntol;
if onbdry   % if checking boundary do not want to a) split by connectivity and b) constrain by population
    L2 = double(ok2); r2 = regionprops(L2, {'Area', 'PixelIdxList'}); q2 = true;
else
    L2 = bwlabel(ok2); r2 = regionprops(L2, {'Area', 'PixelIdxList'}); q2 = [r2.Area] >= ntol;
end

% check for connections from adequately populated components of first set to 2nd; when any
% valid ones found mark the appropriate component of the 2nd set for ref below
done2 = false(size(q2));
for i = find(q1);
    ij = r1(i).PixelIdxList;
    while ~isempty(ij)
        [~, k] = min(d2(ij)); m = double(id2(ij(k))); l = L2(m);
        if q2(l)
            done2(l) = true; kk = [kk; ij(k), m]; break;
        else
            ij(L2(id2(ij)) == l) = [];
        end
    end
end

% check components of 2nd set that have not yet been matched to anything
for i = find(q2 & ~done2);
    ij = r2(i).PixelIdxList;
    while ~isempty(ij)
        [~, k] = min(d1(ij)); m = double(id1(ij(k))); l = L1(m);
        if q1(l)
            kk = [kk; m, ij(k)]; break;
        else
            ij(L1(id1(ij)) == l) = [];
        end
    end
end

return;
    
% -------------------------------------------------------------------------------------------------------------------- %
function [clist, Ac] = chuff(A, rad, tol, res, inbnds, gstd, mbox) %, npt, ptol)
% [circle_struct_list, hough_matrix] = chuff(BWimage, radius, tolerance, resolution,
%                                             inbounds, gaussfilt_std, medfilt_box)
%
% Perform circular hough transform on logical (ie thresholded) image BWimage. This differs
% from ML's imfindcircles in that it is designed to find rims as opposed to filled (bright
% or dark) circles, by acting directly on salient features as defined by the caller.
% Remaining arguments (use an empty array to accept the default):
%   radius: target radius. Does not need to be integer. Can be a range (as with imfindcircles).
%   tolerance: cutoff wrt Hough metric for retention of circles. If given as a value
%      between 0 and 1 then circles with Hough values below the cutoff are discarded. If
%      an integer value k >= 1 is given then the best k circles are returned. Default
%      (roughly 0.15) corresponds to acceptance of any arc segment decently outlining at
%      least one radian.
%   resolution: by default, if a range is given then values at 1/2 pixel increments are checked.
%      (resolution = 1/2). Use this to set a different division.
%   inbounds: if true, return only circles with centers within the image boundaries (ie.
%      what imfindcircles does). Default = false.
%   gaussfilt_std: standard deviation for gaussian filter. Default = 1.5 pixels. Set to 0
%      to omit gaussian filter step.
%   medfilt_box: box size for median filter. Default = [5 5]. Set to 0 to omit median
%      filter step.
%
% Returns a struct array with fields 'center', 'radius', and 'houghval' (ie. the Hough
% metric). Unlike imfindcircles, order of elements is not according to the Hough metric
% but rather by position (standard ML left-to-right). Second return is the actual Hough
% accumulator matrix.

% Note: code here is implementation of standard circular Hough transform, but order of ops
% etc. is derived from study of ML's imfindcircles (which uses a slightly different method,
% optimized for finding filled disks but not so good at finding circle outlines).

% REM'ED: going with vector at whatever increment caller desires
%   refine_radii: by default, if a radius range is given then a finite set of radius values
%      are checked and the best of these returned for each circle (step size is 1/2 pixel).
%      If refine_radii is set to true then additional calculations are done to determine with
%      greater precision what the actual radii are.

% process args first, empty means use default ...
if nargin < 7 || isempty(mbox)    % median filter, default is same as used by imfindcircles
    mbox = [5, 5];
elseif isequal(mbox, 0)
    mbox = [];
elseif length(mbox) < 2
    mbox = [mbox, mbox];
end
if nargin < 6 || isempty(gstd)    % gaussian filter, default is just large enough to have an effect
    gstd = 1.5;
elseif gstd <= 0
    gstd = [];
end
if nargin < 5 || isempty(inbnds)  % flag to throw away out-of-bounds centers (imfindcircles' behavior)
    inbnds = false;
end
nres = 2;                         % division for range-of-radii checking. Default is same as imfindcircles
if nargin >= 4 && ~isempty(res)
    nres = 1/res;
end
tol0 = 1/2/pi;                    % tolerance used for suppression of trivial local maxima in Hough matrix prior to finding
if nargin < 3 || isempty(tol)     % circle centers. This should be generally the same as the explicit tol. If instead we
    tol = tol0;                   % are getting the best circle (or k best circles), we usually don't need the suppression
elseif tol >= 1                   % so by default we skip it, but we can activate it by giving an non-integer tol > 1
    tol0 = mod(tol, 1); tol = floor(tol);
else
    tol0 = tol;
end

% get the list of radii to check.
mn = size(A);
if length(rad) > 1
    rlist = (floor(nres*rad(1)):ceil(nres*rad(end)))/nres; rlen = length(rlist);
else
    rlist = rad; rlen = 1;
end
rcircf = round(2*pi*rlist);
impad = ceil(rlist(end));   % padding on array to check for arc segments with centers outside frame

% create and fill the accumulator. accum is weighted by reciprocal circumference, otherwise
% when there are multiple radii the larger ones would dominate. For m
[iy, jx] = find(A);
if isempty(iy)
    clist = []; Ac = zeros(size(A)); return;
end
mnc = mn + 2*impad;
if rlen > 1   % multidimensional HT, we keep the accumulator 2D by tracking where best vals came from.
    Ac = zeros(mnc); kc = Ac;
    for i = 1:rlen
        kz = linspace(0, 2*pi*(1 - 1/rcircf(i)), rcircf(i));
        xx = bsxfun(@plus, jx, rlist(i)*cos(kz)); yy = bsxfun(@plus, iy, rlist(i)*sin(kz));
        ij = round([yy(:), xx(:)]) + impad; clear xx yy;
        actmp = accumarray(ij, 1/rcircf(i), mn + 2*impad); clear ij;
        kk = actmp > Ac; Ac(kk) = actmp(kk); kc(kk) = i; clear actmp kk;
    end
else
    kz = linspace(0, 2*pi*(1 - 1/rcircf), rcircf);
    xx = bsxfun(@plus, jx, rlist*cos(kz)); yy = bsxfun(@plus, iy, rlist*sin(kz));
    ij = round([yy(:), xx(:)]) + impad; clear xx yy;
    Ac = accumarray(ij, 1/rcircf, mnc); kc = 1; clear ij;
end

Ac2 = Ac;    % apply smoothing filters.
if ~isempty(gstd)
    gbox = ceil(3*gstd); gbox = gbox + 1 - mod(gbox, 2); gbox = [gbox, gbox];
    Ac2 = imfilter(Ac2, fspecial('gaussian', gbox, gstd), 'same');
end
if ~isempty(mbox)
    Ac2 = medfilt2(Ac2, mbox);
end
% suppress local maxima that are too low, then find the remainder. Suppression threshold
% is rescaled by the degree to which the smoothing has reduced the height of the peaks.
if tol0 > 0
    Ac2 = imhmax(Ac2, tol0*max(Ac2(:))/max(Ac(:)));
end

% find maximal regions in the (smoothed) hough matrix
mxlst = regionprops(imregionalmax(Ac2), Ac, {'WeightedCentroid','PixelValues','PixelIdxList'});

% collect centroids from mx list. If the region in the original Ac matrix is uniformly 0 (which we don't want)
% centroid will be nan, so remove these
xy = reshape([mxlst.WeightedCentroid], 2, length(mxlst))'; jj = any(isnan(xy), 2);
if any(jj)    % uniformly zero within original Ac matrix, so maxima is artifact of imhmax procedure, remove
    if all(jj)
        clist = []; return;
    end
    jj = find(~jj); xy = xy(jj,:); mxlst = mxlst(jj);
end

% Convert xy to indices for extraction of best radius, then put xy back into original space
if rlen > 1
    kk = kc(sub2ind(size(kc), round(xy(:,2)), round(xy(:,1)))); jj = kk < 1;
    if any(jj)        % if xy location in kc matrix = 0 then we have a horseshoe or other non-concave region, want to reject
        if all(jj)    % these since they are imhmax artifacts as well
            clist = []; return;
        end
        jj = find(~jj); xy = xy(jj,:); kk = kk(jj); mxlst = mxlst(jj);
    end
else
    kk = ones(length(mxlst), 1);
end
xy = xy - impad;

% if requested, check bounds now
if inbnds
    jj = ~any([xy < 1/2, xy(:,1) >= mn(2) + 1/2, xy(:,2) >= mn(1) + 1/2], 2);
    if ~any(jj)
        clist = []; return;
    end
    xy = xy(jj,:); kk = kk(jj); mxlst = mxlst(jj);
end

% at this point, apply tolerance using the hough values in the flagged maximal regions. Issue: we have 3 non-equivalent
% versions of the hough metric to choose from: the smoothed (Ac2) value, the mean Ac val in the region, and the max Ac val
% in the region. This is particularly problematic when choosing the best circle centerpoint (tol = 1). When the metrics
% agree we are certainly ok, but when they give different choices there is no a priori way to resolve them (time consuming
% visual inspection has turned up cases where each gives the best answer against the others). For now: we will use the
% Nate Silver approach (to political polling) and aggregate.
ii = cellfun(@(x) x(1), {mxlst.PixelIdxList}');        % smoothed matrix has uniform val in max region, can use 1st index
hall = [Ac2(ii), cellfun(@mean, {mxlst.PixelValues})', cellfun(@max, {mxlst.PixelValues})'];
hwt = mean(hall, 2);
if tol == 1
    ii = find(hwt == max(hwt));
    if length(ii) > 1
%         disp(' '); warning('(*) non-unique max houghvals');   % strategy: find region where max individual metric is most dominant
        hchk = bsxfun(@rdivide, hall(ii,:), max(hall(ii,:),[],1)); h2 = sum(hchk, 2);
        jj = find(h2 == max(h2));
        if length(jj) > 1
%             warning('(*) really non-unique!');   % can't imagine this happening, but ...
            clist = []; return;
        end
        ii = ii(jj);
    end
elseif tol > 1
    [~, ii] = sort(hwt, 'descend'); ii = sort(ii(1:min(tol, length(ii))));   % not going to worry about ties for 4th place ...
else
    ii = find(hwt >= tol);   % apply tol (< 1) as cutoff, 0 keeps everything
end
if isempty(ii)
    clist = []; return;
end
xy = xy(ii,:); kk = kk(ii); hwt = hwt(ii); ni = length(ii);

% get radii via the kc (log of which iteration gave the best houghval) matrix
if rlen > 1
    rr = reshape(rlist(kk), ni, 1);
elseif ni > 1
    rr = repmat(rlist, ni, 1);
else
    rr = rlist;
end

% construct the output
if ni == 1
    clist = struct('center',xy, 'radius',rr, 'houghval',hwt);
else
    clist = struct('center',num2cell(xy, 2), 'radius',num2cell(rr), 'houghval',num2cell(hwt));
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function grinf = pgrid2(clist0, spc, imext, stol, htol)
% grid_struct = pgrid2(circle_structs, spacing, spacing_tol, hough_tol)
%
% Checks chuff.m output given in circle_structs against grid positions implied by given
% spacing. Returns a grid struct with a list of likely pillars plus other grid info. The
% likely pillar list will contain circles from the input as well as pillars inferred from
% gaps in the detected grid (the latter can be recognized by NaN houghval fields).
%
% Additional arguments give full image extents and set tolerances for grid construction.
% spacing_tol is used to determine when circles should be retained as pillars or discarded.
% Default is computed from the spacing and pillar radii so as to ensure that pillars do not
% touch (ie. it is not very strict). hough_tol is used to select good circles on which to
% base the initial grid estimate; default = 0.65.
%
% If spacing is not given then circles are checked against each other, and those which
% overlap any circles which are better according to the hough metric are discarded. In this
% case the any tolerances given are ignored. 

grinf = struct('plist',[], 'size',[], 'id',[], 'spacing',[], 'angle',[], 'offset',[], ...
    'origid',[], 'discard',[], 'hghtol',[], 'spctol',[], 'isprob',0);
if isempty(clist0)
    return;
end
if nargin < 2
    spc = [];
end
if nargin < 3   % image extents. If spacing is non-null then this should really exist as well, but allow nothing here for now
    imext = [];
end
if nargin < 4
    stol = [];
end
if nargin < 5 || isempty(htol)
    htol = 0.65;
end
clist = clist0; n = length(clist);
pts = cell2mat({clist.center}'); xy = complex(pts(:,1),pts(:,2));
rds = [clist.radius]'; wts = [clist.houghval]';

% if not creating grid just check for too-close circles and return
if isempty(spc)
    [~, ii] = sort(wts, 'descend');
    for i = 1:n
        if ~isnan(xy(ii(i)))
            dxy = abs(xy(ii(i)) - xy); dxy(ii(i)) = Inf;
            xy(~isnan(xy) & dxy < (rds(ii(i)) + rds(:) + 1)) = NaN;
        end
    end
    jj = isnan(xy); grinf.discard = clist(jj);
    clist(1).grij = []; grinf.plist = clist(~jj);
    grinf.origid = find(~jj); grinf.spctol = stol; return;
end
if length(spc) > 1
    disp(' '); error('Unequal spacing in x and y directions not currently supported.');
end
if isempty(stol)   % default spacing tolerance is very loose, maximum dist without pillars being allowed to touch
    stol = max(spc/2 - (max(rds) + 1), 2.5);
end
grinf.hghtol = htol; grinf.spctol = stol;

% ow we have a grid spec, so find the actual grid. This requires finding the offset of the
% top right corner point and determining if the view is rotated at all (it is quite difficult
% to get a perfectly aligned picture of the micro-environment). To do both we need to use
% the pillar centers we have been given, but only the more trustworthy ones.
% Start with the angle estimate ...
iok = wts >= htol; dxy = bsxfun(@minus, xy(iok).', xy(iok));     % use likely neighbours among ok points, triu mask removes reverse angles
ii = find(triu(ones(size(dxy)),1) & (abs(dxy) > 0.8*spc & abs(dxy) < 1.2*spc));
if isempty(ii)
    grinf.discard = clist; grinf.isprob = 1; return;   % probably error, caller should deal with this
end
axy = mod(angle(dxy(ii)) + pi/4, pi/2) - pi/4;    % this has the effect of converting left, up, and down angles to associated rightward angle
ang1 = median(axy); spc1 = median(abs(dxy(ii)));  % if actual spacing is not exactly what user gives it could cause problems, so get better estimate now

% get estimate for the offset of the top right grid point and size of grid
xtmp = xy(iok)*exp(complex(0,-ang1)); minx = complex(min(real(xtmp)), min(imag(xtmp)));
xtmp = round((xtmp - minx)/spc1); mn1 = [max(imag(xtmp)), max(real(xtmp))] + 1;
ofs1 = minx*exp(complex(0,ang1)); gref1 = bsxfun(@complex, (1:mn1(2))-1, (1:mn1(1))'-1);

% refine the values by using estimates as starting point for simultaneous sol'n using fminsearch
prms2 = fminsearch(@(a) mean(min(abs(bsxfun(@minus, xy(iok), a(1)*gref1(:).'*exp(complex(0,a(2))) + ...
    complex(a(3),a(4)))),[],2)), [spc1, ang1, real(ofs1), imag(ofs1)]);
ofs2 = prms2(3:4); ang2 = prms2(2); spc2 = prms2(1); grinf.angle = ang2; grinf.spacing = spc2;

% check if extents need to be modified to incorporate all points, and then obtain indices into the full grid
xtmp = xy;
if ~isempty(imext)   % dummy vals to push grid to extent of entire frame
    xtmp = [xtmp; complex([1 1 imext(2) imext(2)]', [1 imext(1) 1 imext(1)]')];
end
xtmp = round((xtmp - complex(ofs2(1),ofs2(2)))*exp(complex(0,-ang2))/spc2); minx = complex(min(real(xtmp)),min(imag(xtmp)));
if minx ~= 0
    cmnx = spc2*minx*exp(complex(0,ang2)); ofs2 = ofs2 + [real(cmnx), imag(cmnx)]; xtmp = xtmp - minx;
end
xtmp = [imag(xtmp),real(xtmp)] + 1; mn2 = max(xtmp, [], 1);
ij = sub2ind(mn2, xtmp(1:n,1), xtmp(1:n,2));    % using 1:n here will prevent inclusion of dummies
grinf.offset = ofs2; grinf.size = mn2;

% use a new reference grid to check whether points are close enough to grid locations, and check whether
% multiple elements are binned to the same location; save the closet/best to the grid point.
gref2 = spc2*bsxfun(@complex, (1:mn2(2))-1, (1:mn2(1))'-1)*exp(complex(0,ang2)) + complex(ofs2(1),ofs2(2));
gid = zeros(mn2); jok = false(size(xy));
for i = 1:n
    k = ij(i); j = gid(k);
    if abs(xy(i) - gref2(k)) > stol
        continue;
    elseif j == 0
        gid(k) = i; jok(i) = true;
    elseif abs(xy(i) - gref2(k))/wts(i) < abs(xy(j) - gref2(k))/wts(j)    % this takes into account both dist and houghval
        jok(j) = false; gid(k) = i; jok(i) = true;
    end
end
grinf.discard = clist(~jok);
if nnz(~jok) > nnz(jok)
    grinf.isprob = 2;
end

% check grid for missing intermediate elements; in return val clist these will be filled in with placeholders
% under assumption that pillars are in these locations, only obscured or deformed for whatever reason.
tmpid = gid == 0;   % following calculation of extents will have same effect if done by columns
tmp1 = [sum(cumprod(tmpid, 2), 2) + 1, mn2(2) - sum(cumprod(fliplr(tmpid), 2), 2)];
tmp2 = tmp1; cur1 = tmp1(1,:); cur2 = tmp2(end,:);
for i = 1:mn2(1)
    cur1 = [min(cur1(1),tmp1(i,1)), max(cur1(2),tmp1(i,2))];
    cur2 = [min(cur2(1),tmp2(mn2(1)-i+1,1)), max(cur2(2),tmp2(mn2(1)-i+1,2))];
    tmp1(i,:) = cur1; tmp2(mn2(1)-i+1,:) = cur2;
end
tid = [max(tmp1(:,1),tmp2(:,1)), min(tmp1(:,2),tmp2(:,2))];
for i = 1:mn2(1)
    ii = tid(i,1):tid(i,2); ii = ii(gid(i,ii) == 0); gid(i,ii) = -1;
end
grinf.origid = gid;
if ~grinf.isprob && nnz(gid < 0) > nnz(gid > 0)/2
    grinf.isprob = 3;
end

% put together the return values. Start with tmp cell array so we can insert placeholders
jj = find(gid > 0); ctmp = cell(mn2); ctmp(jj) = num2cell(clist(gid(jj)));
jj = find(gid < 0); ctmp(jj) = num2cell(struct('center',num2cell([real(gref2(jj)),imag(gref2(jj))], 2), 'radius',mean(rds), 'houghval',NaN));
for i = 1:mn2(1)    % set the grid position field
    for j = 1:mn2(2)
        if ~isempty(ctmp{i,j})
            ctmp{i,j}.grij = [i, j];
        end
    end
end
ctmp = ctmp'; clist = [ctmp{:}]';

% save the indices for later back ref from the grid positions
g2id = zeros(mn2); i2j = cell2mat({clist.grij}');
g2id(sub2ind(mn2, i2j(:,1), i2j(:,2))) = 1:length(clist);
grinf.id = g2id; grinf.plist = clist;

return;

% -------------------------------------------------------------------------------------- %
function mx = ksmode(x, d)
% mx = kmode(x, d)
%
% Uses ksdensity to get an estimate for the mode wrt a continuous distribution (ML's mode function
% will only work when values come from a discrete distribution). 2nd arg d gives dimension to
% work along (as with other ML stats functions); caller can also give NaN here, in which case
% a single mode will be returned for all values. For complex numbers the real and imaginary parts
% are evaluated separately and then returned in complex form. Finally, any NaN's in data are stripped,
% so in this sense works like ML's nanmean etc.

sz = size(x); isz = find(sz > 1);
if isempty(sz)
    mx = x; return;
elseif nargin < 2 || isempty(d)
    d = isz(1);
end
if ~isreal(x)
    mx = complex(ksmode(real(x), d), ksmode(imag(x), d));
    return;
end
if ~isnan(d) && (length(isz) > 1 || isz ~= d)
    cx = num2cell(x, d); mx = zeros(size(cx));
    for i = 1:numel(cx)
        mx(i) = ksmode(cx{i});
    end
    return;
end

x = x(:); x = x(~isnan(x)); n = length(x);
if n <= 2
    mx = mean(x); return;
end
bw = (median(abs(x - median(x)))/0.6745)*(4/(3*n))^(1/5);
[ypdf, xpdf] = ksdensity(x, 'width', bw/2);    % use half ML's default bandwidth, don't want too PDF smooth
fpdf = spline(xpdf, -ypdf);                    % use spline utility fnmin to get the peak
[~, mx] = fnmin(fpdf);

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [amsk, ijrng] = mkmask(cxy, rd, imsz, dnum, fullsz)
% [amsk, ijrng] = mkmask(cxy, rd, imsz, dnum, fullsz)
%
% Given center, radius, and subimage bounds, create a mask for the circle (a dark circle
% in a subimage with extents 2x that of the circle itself. Optional params: imsz gives outer
% bounds, if any part of generated image is outside this range it is truncated to fit; dnum
% specifies that the circle be dilated (+ve) or filled (-ve); and if fullsz is set then the
% size of the mask will be set exactly to the bounds given by imsz (expanding if necessary).
% Returns the logical matrix amsk, plus its bounds wrt the given center.

if nargin < 5 || isempty(fullsz)
    fullsz = false;
end
if nargin < 4 || isempty(dnum)
    dnum = 0;
end
if nargin < 3
    imsz = [];
elseif length(imsz) == 2
    imsz = [1 imsz(1) 1 imsz(2)];    % put in ij bounds format
end

% start by constructing a circle outline so as to have no gaps or unnecessary pixels
x = cxy(1); y = cxy(2);
ix = ceil(x - rd):floor(x + rd); jy = ceil(y - rd):floor(y + rd);
a1 = [ix', y + sqrt(rd^2 - (ix' - x).^2); ...
       x + sqrt(rd^2 - (fliplr(jy)' - y).^2), fliplr(jy)'; ...
    fliplr(ix)', y - sqrt(rd^2 - (fliplr(ix)' - x).^2); ...
       x - sqrt(rd^2 - (jy' - y).^2),  jy'];

% round and remove redundancies (points that are in the same position)
a1 = unique(round(a1), 'rows');
csz = fliplr([min(a1,[],1); max(a1,[],1)]);
asz = diff(csz) + 1;

% create the mask, with padding
amsk = false(asz); amsk(sub2ind(asz, a1(:,2)-csz(1,1)+1, a1(:,1)-csz(1,2)+1)) = true;
pval = round(rd); amsk = padarray(amsk, [pval, pval], false);
csz(1,:) = csz(1,:) - pval; csz(2,:) = csz(2,:) + pval;

% apply any desired ops
if dnum < 0
    amsk = imfill(amsk, round(size(amsk)/2));
elseif dnum > 0
    amsk = cdilate(amsk, cxy - fliplr(csz(1,:)) + 1, rd, dnum);
end

% truncate the mask and adjust the bounds if need be
ijrng = reshape(csz, 1, 4);
if ~isempty(imsz)
    [i2, j2] = size(amsk); ijfix = max((imsz - ijrng).*[1 -1 1 -1], 0);
    if any(ijfix)
        amsk = amsk((1+ijfix(1)):(i2-ijfix(2)),(1+ijfix(3)):(j2-ijfix(4)));
        ijrng = ijrng + ijfix.*[1 -1 1 -1];
    end
    if fullsz    % padding specified as well
        pdsz = (ijrng - imsz).*[1 -1 1 -1];
        if any(pdsz > 0)
            amsk = padarray(padarray(amsk, pdsz(1:2:3), 'pre'), pdsz(2:2:4), 'post');
            ijrng = imsz;
        end
    end
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function Ac2 = cdilate(Ac1, cpt, rad, val)
% Ac2 = cdilate(Ac1, cpt, val)
% Perform "circular dilation" on the image. Since image is always going to be a circle, we can get a better (as in more
% circular) effect than straightforward dilation by retaining points if and only if they are within exactly val + 1/2 of
% the distance rad from the given center cpt. An added benefit is that fractional values can be used, allowing a
% pseudo-continuous increase in rim width as desired.

% use ordinary dilation on slightly bigger val to start, then we'll peal away the excess
v2 = val + 1/2; dval = ceil(v2);
Ac2 = imdilate(Ac1, strel('disk', dval));
ij = find(Ac2); [iy, jx] = ind2sub(size(Ac1), ij);
dxy = abs(complex(jx, iy) - complex(cpt(1), cpt(2))) - rad;
Ac2(ij(dxy < -v2 | dxy > v2)) = false;

return;

% -------------------------------------------------------------------------------------------------------------------- %
function C = getobj2(A, th, amsk, pk, prndots)
% objlist = getobj2(movie_data, threshold, obj_mask, do_pack)
%
% Scan movie data for objects. Arguments threshold and obj_mask are the threshold value to use
% and a logical mask for foreground objects respectively. Argument do_pack controls object
% packing: 0 (default) = no packing, 1 = straight packing of large objects, and 2 = agressive
% packing of all objects. Return value is a cell array (1 cell per frame) containing the object
% list (regionprops structs).

[n1, n2, nn] = size(A); C = cell(nn,1);
if nargin < 3 || isempty(amsk)
    amsk = true(n1, n2);
end
if nargin < 4 || isempty(pk)
    pk = 0;
end
if nargin < 5 || isempty(prndots)
    prndots = 0;
end

% get the objects ...
for i = 1:nn
    ci = regionprops(A(:,:,i) > th & amsk, {'Area','Centroid','BoundingBox','Image'});
    if ~isempty(ci)
        [ci.ij] = bnd2idx(ci); [ci.packed] = deal(0);
        if pk ~= 0
            for j = 1:length(ci)
                if ci(j).BoundingBox(4) >= 32
                    ci(j).Image = bwpack(ci(j).Image);
                    ci(j).packed = 1;
                elseif ci(j).BoundingBox(3) >= 32
                    ci(j).Image = bwpack(ci(j).Image')';
                    ci(j).packed = -1;
                end
            end
        end
    end
    C{i} = ci;
    if prndots > 0 && ~mod(i, prndots)
        fprintf(1, '.');
    end
end
    
return;

% -------------------------------------------------------------------------------------------------------------------- %
function varargout = bnd2idx(xywh)
% Convert bounding box bounds given by regionprops [x1, y1, w, h] into subscript ranges [i1, i2, j1, j2]

if isstruct(xywh)
    xywh = cell2mat({xywh.BoundingBox}');
elseif iscell(xywh)
    xywh = cell2mat(xywh(:));
end
rxy12 = [ceil(xywh(:,1:2)), floor(xywh(:,1:2) + xywh(:,3:4))];
ij12 = rxy12(:,[2 4 1 3]);
if nargout < size(ij12, 1)
    varargout{1} = ij12;
else
    varargout = num2cell(ij12, 2)';
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [Tr, Ts, tdata] = trkobj2(C, lgsc, mdsc) %, obg)
% [Tr, Ts, tdata] = trkobj2(C, lgsc, mdsc)
%
% Convert per-frame object lists in cell array C to per-track object lists in struct arrays Tr and Ts. Additional
% arguments:
%   lgsc - length scale for trackable objects (eg. worms)
%   mdsc - length scale for medium scale objects (eg. pillar diameter)
% The scale arguments are used to sort objects into categories requiring different actions depending on whether we
% expect the objects to be moving and possibly entering or leaving the frame, or standing still and occasionally being
% touched or occluded.
%
% Return values:
%   Tr - struct array of tracks for mobile large objects of interest
%   Ts - struct array of tracks for static medium scale objects of potential interest
%   tdata - structure containing general information which can be used for debugging etc.
% Note: per-frame tracking info for Tr and Ts will be contained in a sub-array in the field 'data'. The information at
% the primary level will be statistical in nature (ie. median position or area and range of movement etc.)

if nargin < 2 || isempty(lgsc)    % maximum length scale = average length of long thin trackable object
    lgsc = sqrt(15*10^4.5);
end
if nargin < 3 || isempty(mdsc)
    mdsc = 100/pi;                % medium length scale = min height and width of non-negligible hollow object
end
nc = length(C);
tdata = struct('numobj',[], 'byframe',[], 'largescale',lgsc, 'medscale',mdsc, 'areadiv',[], 'byarea',[], ...
    'trklg',[], 'trkmd',[]);

% apply scales to object categorization by area. Large area scale =~ worm size, medium area scale =~ open circle with
% diameter given by mdsc and rim at least 1-pixel wide. Defaults set above reduce this to orders of magnitude in
% usual sense eg 1, 10, 100, 1000, 10000, ...
s1 = log10(pi*mdsc^2/min(mdsc,25)); s2 = log10(lgsc^2/15);
if s2 - s1 < 1e-2
    pm = [0, 5];    % if bad vals for our length scales then have no choice but to track everything. Want actual diff between 2 and 3 or so
else
    pm = 2.5/(s2 - s1); pm = [pm, 3 - pm*s1];
end
tdata.areadiv = round(10.^(((1:6) - pm(2))/pm(1)));

% go through the object list and sort by area
cc = cell(nc, 6);
for i = 1:nc
    ci = C{i}; ni = length(ci);
    if ni == 0
        continue;
    end
% categories: these are 3 = medium, static; 5 = large, trackable; + 1 (flotsam), 2, 4 (intermediate) & 6 = too large
    ri = polyval(pm, log10([ci.Area])); ii = max(1, min(floor(ri), 6));

% since category 3 is essentially meant for pillars (ie thin-rimmed annuli with diameter approximately equal to the medium
% scale), demote cat 3 objects if they are too compact and promote then if too extended
    kk = find(ii == 3); bb = reshape([ci(kk).BoundingBox], 4, length(kk))';
    ii(kk(all(bb(:,3:4) < mdsc/2, 2))) = 2; ii(kk(all(bb(:,3:4) > 2*mdsc, 2))) = 4;
    for j = 1:6
        cc{i,j} = ci(ii == j);
    end
end
numc = cellfun('length', cc); tdata.byframe = sum(numc, 2);
tdata.byarea = sum(numc, 1); tdata.numobj = sum(tdata.byarea);
isc = numc > 0; cc(~isc) = {[]};

% check packing, make consistent per category (want 4 & 5 packed, others not), and remove packing indicator
ptmp = cell(size(cc)); ptmp(isc) = cellfun(@(x) [x.packed], cc(isc), 'UniformOutput', false);
for j = 1:3    % unpack anything that got packed in small-to-medium object lists
    if any([ptmp{:,j}])
        for i = find(cellfun(@any, ptmp(:,j)))'
            for k = find(ptmp{i,j} < 0)
                cc{i,j}(k).Image = bwunpack(cc{i,j}(k).Image', cc{i,j}(k).BoundingBox(3))';
            end
            for k = find(ptmp{i,j} > 0)
                cc{i,j}(k).Image = bwunpack(cc{i,j}(k).Image, cc{i,j}(k).BoundingBox(4));
            end
        end
    end
end
for j = 4:6    % pack anything that didn't get packed in larger object lists (unlikely to find anything here, but ...)
    if any([ptmp{:,j}] < 1)
        for i = find(cellfun(@(x) any(x < 1), ptmp(:,j)))'
            for k = find(ptmp{i,j} < 0)
                cc{i,j}(k).Image = bwunpack(cc{i,j}(k).Image', cc{i,j}(k).BoundingBox(3))';
            end
            for k = find(ptmp{i,j} < 1)
                cc{i,j}(k).Image = bwpack(cc{i,j}(k).Image);
            end
        end
    end
end
cc(isc) = cellfun(@(x) rmfield(x, 'packed'), cc(isc), 'UniformOutput', false);

% lets start with sanity check for situations possibly arising from input problems (either bad movie or user params) or
% bad behavior (worms clustering). Caller should return with this output so user can inspect it as soon as possible.
if any(isc(:,6))
    disp(' '); warning('Extremely large objects encountered, check output.');
    Tr = cc; Ts = []; return;
end

% first pass for trackable objects (we'll return to these later when checking against other categories); again we start
% with the sanity check. We then check a few unambiguous cases which should not require explicit tracking before ...
tdata.trklg = struct('tracks',[], 'persistent',false, 'area_tol',0.1, 'sub_box',mdsc);
if nnz(isc(:,4:6)) < 1
    disp(' '); warning('No objects large enough for tracking found, check output.');
    Tr = cc; Ts = []; return;
end
jj = find(any(isc(:,4:5), 2));
if all(sum(numc(:,4:5), 2) <= 1) && ~any(diff(jj) > 1)   % at most 1 primary object at any one time, with no gaps (but can
    Tr = cc(:,5); Tr(~isc(:,5)) = cc(~isc(:,5),4);       % be leaving / entering frame), so don't need to do explicit track
else
    Tr = mtrack(cc(:,4:5), 0.1);    % 0.1 here means area is allowed to change by at most 10% between frames
end
cc(:,4:5) = {[]}; clear isc numc

% track the medium sized (presumably static) objects, then check smaller objects to see if they can be associated (as eg
% pieces disconnected from the main object because of close-to-threshold pixel values in original movie). Note: caller
% will verify tracks and associations at a later stage, using pillar array info.
tdata.trkmd = struct('tracks',[], 'persistent',true, 'dist_tol',mdsc/2, 'sub_dist',mdsc);
Ts = ptrack(cc(:,3), mdsc/2); cc(:,3) = {[]};
[ts2, cc(:,1:2)] = rtrack(cc(:,1:2), Ts, mdsc, false);    % last arg flags that we want tol to apply to centroid distances
Ts = mkcompact(Ts, ts2); clear ts2;                                 % convert cells to struct arrays for easier handling

% now check leftovers of previous check against tracking objects. In this case we'll collect anything that is within or
% close enough to the worm's bounding box (last arg to rtrack says to use the bounding box tolerance
[tr2, cc(:,1:2)] = rtrack(cc(:,1:2), Tr, mdsc, true);
Tr = mkcompact(Tr, tr2); clear tr2;
tdata.trklg.tracks = rmfield(Tr, 'data');    % save record of tracks found (though not images or trajectories etc)

% Any leftovers in category 1 will be specks etc, so we don't care about them, but objects in category 2 could be
% relevant, so we compile these into track form and append to the medium scale array
if any(~cellfun('isempty',cc(:,2)))
    tmp1 = ptrack(cc(:,2), mdsc/2); tmp2 = rtrack(cc(:,1), tmp1, mdsc, false);
    Ts = cat(2, Ts, mkcompact(tmp1, tmp2)); clear tmp1 tmp2;
end
tdata.trkmd.tracks = rmfield(Ts, 'data');

return;

% -------------------------------------------------------------------------------------------------------------------- %
function tc = mtrack(c, atol)
% tc = mtrack(c, atol)
% Perform standard nearest neighbor tracking on assumption of mobile objects. c is a cell array of image objects,
% and atol is the tolerance with respect to area change. Returns length(c) X num_tracks cell array tc, where the j-th
% column tracks the j-th distinct object, and the elements of the i-th row come from the i-th cell of c.

n = size(c, 1); tc = cell(n, 0);
isc = ~cellfun('isempty', c); kl = 0; kk = []; km = 0;
for j = 1:n
% cj = next set of candidates. If the entering set is empty then we clear the current set of tracks; otherwise, if the
% current track set (cl) is empty we promote all the entering objects into new tracks (note: this also works to
% initialize the track list in the 1st iteration). If both are non-empty then we proceed with the comparison as usual.
    cj = cell2mat(c(j,isc(j,:))'); kj = length(cj);
    if kj == 0
        kl = 0; kk = [];
        continue;
    elseif kl == 0;
        tc = cat(2, tc, cell(n, kj)); kk = km + (1:kj); km = kk(end);
        tc(j,kk) = num2cell(cj'); kl = kj;
        continue;
    else
        cl = cell2mat(tc(j-1, kk));
    end
    
% Criteria: all-to-all Euclidean distance comparison with verification by all-to-all area comparison (this is the
% simplest way to pre-apply tolerances). We are counting on the current set of tracks to not be too large (for test
% cases what we do here is still faster than eg. knnsearch, probably because of the minimal overhead required)
    xl = reshape([cl.Centroid], 2, kl)'; xj = reshape([cj.Centroid], 2, kj)';
    d12 = abs(bsxfun(@minus, complex(xj(:,1), xj(:,2)), complex(xl(:,1), xl(:,2)).'));
    da = abs(bsxfun(@rdivide, [cj.Area]', [cl.Area]) - 1); d12(da > atol) = NaN;
    [d2, j2] = unqmin(d12);   % version of row-min that only gives non-NaN vals for unique columns

% assign to tracks, + end tracks without any current addition, + create new if necessary
    ii = ~isnan(d2); jj = j2(ii);
    if ~isempty(jj)
        tc(j,kk(jj)) = num2cell(cj(ii)); kk = sort(kk(jj)); kl = length(kk);
    else
        kk = []; kl = 0;
    end
    ii = ~ii; ki = nnz(ii);
    if ki > 0
        tc = cat(2, tc, cell(n, ki)); ktmp = km + (1:ki); kk = [kk, ktmp]; km = kk(end);
        tc(j, ktmp) = num2cell(cj(ii)'); kl = kl + ki;
    end
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function tc = ptrack(c, dtol)
% tc = ptrack(c, dtol)
% Perform persistent version of nearest neighbor tracking operation on arrays of regionprops structs in cell array c,
% where each cell corresponds to a frame of the data source. Additional argument dtol gives the distance tolerance to
% be applied.
%   As opposed to the standard mobile object tracking implemented in mtrack (above), this tracking method assumes that
% objects are relatively static in location and persistent across time; which means that they are allowed to disappear
% temporarily and then reappear later. Entering objects are associated with the track object having the most similar
% size within the allowed distance range. For track objects both centroid and area are aggregate measures.
%   Like mtrack this returns a length(c) X num_tracks cell array tc, where the j-th column tracks the j-th distinct
% object, and the elements of the i-th row come from the i-th cell of c.

% since tracks are persistent (so the current set can never go back to zero) we start by promoting the 1st non-empty
% cell in c to the current set of representatives. As opposed to above, the current track reps are maintained in a
% struct array, with relevant values updated as needed.
n = size(c, 1); isc = ~cellfun('isempty', c);
j1 = find(any(isc, 2), 1, 'first');
if isempty(j1)
    tc = cell(n, 0); return;
end
cl = cell2mat(c(j1,find(isc(j1,:)))); kl = length(cl); kcnt = ones(kl, 1);
tc = cell(n, kl); tc(j1, :) = num2cell(cl);

for j = (j1+1):n
% next set of candidates. Unlike above, if nothing in this frame we just go to the next
    cj = cell2mat(c(j,isc(j,:))'); kj = length(cj);
    if kj == 0
        continue;
    end

% criteria are essentially the same as above, though applied in the opposite order.
    xl = reshape([cl.Centroid], 2, kl)'; xj = reshape([cj.Centroid], 2, kj)';
    d12 = abs(bsxfun(@minus, complex(xj(:,1), xj(:,2)), complex(xl(:,1), xl(:,2)).'));
    al = [cl.Area]'; aj = [cj.Area]'; da = abs(bsxfun(@rdivide, aj, al') - 1);
    da(d12 > dtol) = NaN; [d2, j2] = unqmin(da);    % as above we want unique assignments to tracks

% assign to tracks, + create new if necessary (this is still possible event though tracks can't end)
    ii = ~isnan(d2); jj = j2(ii);
    if ~isempty(jj)
        tc(j,jj) = num2cell(cj(ii));   % below: update aggregate metrics for current track list
        xtmp = num2cell((repmat(kcnt(jj),1,2) .* xl(jj,:) + xj(ii,:))./repmat(kcnt(jj)+1,1,2), 2);
        atmp = num2cell((kcnt(jj) .* al(jj) + aj(ii))./(kcnt(jj) + 1)); kcnt(jj) = kcnt(jj) + 1;
        [cl(jj).Centroid] = deal(xtmp{:}); [cl(jj).Area] = deal(atmp{:});
    end
    ii = ~ii; ki = nnz(ii);
    if ki > 0    % add the new tracks. 
        tc = cat(2, tc, cell(n, ki)); tc(j, kl + (1:ki)) = num2cell(cj(ii));
        cl = [cl; cj(ii)]; kcnt = [kcnt; ones(ki,1)]; kl = kl + ki;
    end
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [tc, sc] = rtrack(c, rc, dtol, dob)
% [tc, sc] = rtrack(c, rc, dtol, dob)
% Track regionprops objects in cell array c against existing track objects given in reference cell array rc. Remaining
% args dtol and dob give distance tolerance and flag to apply this to bounding box (versus centroid) respectively.
%   In this case object locations are compared to locations of reference objects, and if within tol, assigned to the
% nearest as a subobject. Unique assignments are not necessary.
%   Returns length(c) X t cell array, where t is the number of tracks in rc, in which the j,i-th element are all
% subobjects associated with the i-th track in frame i. The second return value sc is a cell array the same size as c
% giving leftover unassigned objects in c.

[n, m] = size(c); tc = cell(size(rc)); sc = cell(n, m);
isc = ~cellfun('isempty', rc); numc = cellfun('length', c);
if nargin < 4 || isempty(dob)
    dob = false;
end

for j = 1:n
    cj = cell2mat(c(j,numc(j,:)>0)'); kj = length(cj);
    if kj == 0
        continue;
    end
    kk = find(isc(j,:)); kl = length(kk);
    if kl == 0    % no active tracks in current frame, send entire current set to leftovers
        sc(j,:) = c(j,:); continue;
    end
    
    cl = cell2mat(rc(j,kk));
    xl = reshape([cl.Centroid], 2, kl)'; xj = reshape([cj.Centroid], 2, kj)';
    d12 = abs(bsxfun(@minus, complex(xj(:,1), xj(:,2)), complex(xl(:,1), xl(:,2)).'));
    if dob   % compare bounding boxes. 
        blwtol = bsxfun(@plus, reshape([cl.BoundingBox], 4, kl)', [-[1,1]*dtol, [2,2]*dtol]);
        bj = reshape([cj.BoundingBox], 4, kj)'; db = rectint(bj, blwtol);
        d12(bsxfun(@lt, db, prod(bj(:,3:4), 2))) = NaN;    % cj box needs to be completely contained within dtol of cl's box
%         d12(db == 0) = NaN;     % any edge of cj's box needs to be within dtol of cl's box
    else
        d12(d12 > dtol) = NaN;
    end
    [d2, j2] = min(d12, [], 2);   % assign each object to closest track object without constraints

    ii = ~isnan(d2); jj = find(ii);   % assign (non-uniquely) good matches to relevant tracks
    for i = 1:length(jj)
        tc{j,kk(j2(jj(i)))} = [tc{j,kk(j2(jj(i)))}; cj(jj(i))];
    end
    ii = mat2cell(~ii, numc(j,:), 1);   % copy the non-assigned objects to the leftover output
    for i = 1:m
        sc{j,i} = c{j,i}(ii{i});
    end
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [mval, mi] = unqmin(A)
% [mval, mi] = unqmin(A, dim)
% This is wrapper for min function as applied specifically to rows of matrices, that also checks for uniqueness wrt
% columns. Conflicts are resolved by choosing the row which provides the most minimal value for the column; losers are
% assigned NaN values and 0 indices in the output.

[mval, mi] = min(A, [], 2);
[mchk, ii] = sortrows([mi, mval], [1 -2]);   % second col of mchk insures that last in set of consecutive indices
jj = find(diff(mchk(:,1))' == 0);            % corresponds to min of mins
mval(ii(jj)) = NaN; mi(ii(jj)) = 0;

return;

% -------------------------------------------------------------------------------------------------------------------- %
function aobj = mkcompact(cobj, bobj)
% aobj = mkcompact(cobj, bobj)
% Collect data in n x m cell array of sorted objects cobj into a 1 x m track structure array aobj. Individual per-frame
% data structs are held in subfield data, other fields contain median values + other stats. 2nd arg bobj is same-sized
% cell array as cobj containing associated smaller objects; these are added as sub-fields to corresponding elements of
% the data arrays (caller can decide whether to incorporate into Images or not).

[n, m] = size(cobj); ic = ~cellfun('isempty', cobj);
if nargin < 2 || isempty(bobj)
    bobj = cell(n,m);
else
    bobj = cellfun(@transpose, bobj, 'UniformOutput',false);    % orient everything in subs field horizontally
end                                                             % (more convenient orientation down the road)
ib = cellfun('length', bobj);
% if any(ib(:) & ~ic(:))
%     disp(' '); warning('(*) something went wrong during subobject allocation');
%     disp(['Frames: ' num2str(find(any(ib & ~ic, 2))')]);
% end

aobj = struct('Area',cell(1,m), 'Centroid',[], 'Range',[], 'ij',[], 'data',[], 'frames',[], 'subs',[]);
for i = 1:m
    ii = find(ic(:,i))'; ni = length(ii); ai = cell2mat(cobj(ii,i));
    [ai.sub] = deal(bobj{ii,i}); aobj(i).subs = sum(ib(:,i));

    ar = [ai.Area]; ac = reshape([ai.Centroid], 2, ni)'; aj = reshape([ai.ij], 4, ni)';
    aobj(i).Area = median(ar); aobj(i).Centroid = median(ac, 1);
    acm = min(ac, [], 1); aobj(i).Range = [acm, max(ac, [], 1) - acm];
    aobj(i).ij = reshape([min(aj(:,[1 3]), [], 1); max(aj(:,[2 4]), [], 1)], 1, 4);

    jj = find(diff(ii) > 1); aobj(i).frames = [ii(1) ii(jj+1); ii(jj) ii(ni)]';
    ii = num2cell(ii); [ai.frame] = deal(ii{:}); aobj(i).data = ai;
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [nplst, reobj] = obj2plr(imext, numf, rad, objs, plst, grdp)
% [nplst, reobj] = obj2plr(imext, numf, rad, objs, plst, grdp)
%
% Compare the array of medium-scale static found objects (objs) with previously obtained estimates of trackable pillars
% (plst) and grid parameters (grdp); 1st 3 args imext, numf, and rad give image extents, number of frames, and
% pillar radius respectively. Function checks objects against pillar locations / dimensions and grid coords, and if
% there is a satisfactory match uses the obj info to refine the original estimate of the location and radius. Returns
% nplist (revised pillar list) and reobj (cell array of object lists for each frame which either did not match a pillar
% or had extraneous non-pillar material; this will be readded to tracking object images during the next stage.
%
% If there were no pillars in the experiment, or they were not arrayed, plst and/or gridp can be empty; in the latter
% case the 1st output nplist will also be empty.

nplst = []; reobj = cell(numf, 1);
if nargin < 6
    grdp = [];
end
if nargin < 5
    plst = [];
end
if nargin < 4 || isempty(objs)    % without objects we can't do anything
    return;
end

% combine obj tracks if centroids are closer then a pillar diameter. Original tracking phase used tighter tolerance
% since we didn't want track centroids (which were dynamically updated) to wander too much; since we now have final
% centroid vals we can afford to be a bit more forgiving. Note: anything in stobj which does not match a pillar rim
% fairly well will be folded back into reobj by end.
stobj = objs; nt = length(stobj);
ic = num2cell(1:nt); [stobj.track] = deal(ic{:});   % log orig track #'s for back-ref
xyb = reshape([stobj.Centroid], 2, nt)';

% clusterdata is ML's hierarchical clustering func; called with these params it will simply associate objects closer
% than the cutoff (but still make sure all pairs in each cluster are close enough, so some possible assoc's are dropped)
if nt > 1
    chkdst = clusterdata(xyb, 'cutoff',2*rad, 'criterion','distance', 'linkage','complete');
    if max(chkdst) < nt            % if everyone further than 1 diam from each other all clusters are trivial, nothing to do
        kk = true(1, nt); chkdst = sortrows([chkdst, (1:nt)']);
        chkdst = mat2cell(chkdst(:,2), diff(find([1; diff(chkdst(:,1)) > 0; 1])), 1);  % convert chkdst from cluster index
        for i = find(cellfun('length', chkdst)' > 1)                                   % vector to cells of indices to cluster
            ii = chkdst{i}; stobj(ii(1)) = mrgtrks(stobj(ii)); kk(ii(2:end)) = false;  % members; kk marks the stobj's to rem afterwards
        end
        stobj = stobj(kk); nt = nnz(kk);
    end
end

% Next: incorporate subobjects into primaries if they are close enough to have any portion within the general object's
% bounding box. subs that are too far out are put back into the reobj
for i = 1:nt
    if stobj(i).subs > 0
        bbi = [stobj(i).ij([3 1]) - 1/2, stobj(i).ij([4 2]) - stobj(i).ij([3 1]) + 1];
        jj = find(~cellfun('isempty', {stobj(i).data.sub}));
        for j = jj
            tmpd = stobj(i).data(j); tmpb = tmpd.sub; tmpd.sub = [];
            for k = 1:length(tmpb)
                rv = rectint(bbi, tmpb(k).BoundingBox);
                if rv > 0
                    tmpd = addobj(tmpd, tmpb(k));    % addobj combines disjoint images and updates relevant fields
                else
                    reobj{tmpd.frame} = [reobj{tmpd.frame}, tmpb(k)];
                end
            end
            stobj(i).data(j) = tmpd;
        end
    end
end
stobj = rmfield(stobj, 'subs');

% Next: for each track obtain a composite real-valued image based on superimposing the images for each frame. In the
% process we also partition the data field into equivalence sets wrt bounding boxes (a quick stand in for morphological
% variation). These will be used differently whether the object is determined to be a pillar or not.
[stobj.composite] = deal([]);
for i = 1:nt
    [stobj(i).composite, stobj(i).data] = collectbb(stobj(i).data);
end

% Next: generate pillar/grid data and do primary (distance) check of objects against pillar / grid locations (if provided).
% Matching objects will be subjected to more scrutiny below, followed by revision of the pillar coords. Var xyp2 contains
% center, radius, pillar index (from general list), and (lastly) matching object index and subindex
xyb = reshape([stobj.Centroid], 2, nt)'; np = length(plst);
nplst = cell(np, 1); chkc = cell(nt, 3); kchk = [];
if ~isempty(grdp)   % have grid, use grid position to match objects to pillars via reverse application of grid params
    R_a = [cos(grdp.angle) sin(grdp.angle); -sin(grdp.angle) cos(grdp.angle)];      % grid's rotation (usually slight but can't be ignored)
    xyp2 = [reshape([plst.center],2,np); repmat(rad,1,np); plst.pnum; NaN(2,np)]';  % this contains centers + other relevant info for each pillar

    gijb = fliplr(round((bsxfun(@minus, xyb, grdp.offset)*(R_a'))/grdp.spacing)) + 1;   % pt to grid conversion
    gp2 = reshape([plst.grij], 2, np)'; [ism, kchk] = ismember(gijb, gp2, 'rows');
    if ~all(ism)   % if we find anything outside of known pillar positions generate temporary grid elements
        [gp3, ~, kk] = unique(gijb(~ism,:), 'rows'); ng = size(gp3, 1);
        xytmp = bsxfun(@plus, grdp.spacing * (fliplr(gp3) - 1) * R_a, grdp.offset);
        xyp2 = [xyp2; xytmp, repmat(rad, ng, 1), NaN(ng, 3)];
        kchk(~ism) = np + kk; np = np + ng;
    end
elseif ~isempty(plst)   % have pillars but no grid, so we can only check the pillars we have.
    xyp2 = [reshape([plst.center], 2, np); plst.radius; plst.pnum; NaN(2, np)]';
    kchk = knnsearch(xyp2(:,1:2), xyb); % kchk = [jj, (1:nt)'];
end

% take a closer look at object to pillar matches. We do this check in a separate pass before generating revised pillars
% in order to a) remove from consideration non-pillar objects that might have wandered close to a pillar location, and
% b) select from the data objects the one that best matches the target circle specs.
kt = true(nt, 1);
for k2 = 1:length(kchk)
    k = kchk(k2);
    if sqrt(sum((xyb(k2,1:2) - xyp2(k,1:2)).^2)) > 2*rad    % object too distant from pillar/grid pt, entire track goes to reobj
        continue;
    end
    tmpd = stobj(k2).data; nd = length(tmpd); chk1 = [];
    for j = 1:nd
        chk1 = [chk1, circhk(tmpd(j), xyp2(k,1:2), xyp2(k,3), imext)]; % disp([k2, j]); % for debugging
    end
    [~, ibst] = sortrows([chk1.ismatch; chk1.rval]'); ibst = ibst(end);
    if chk1(ibst).ismatch == 0                                % close enough, but non-match, so ditch track like above
        continue;
    end

    ktr = stobj(k2).track; jfr = [];  % collect track and frame numbers for back-ref. If ismatch value is very low non-zero
    for j = 1:nd                      % pillar possibly touched by worm, so give affected data elements back to reobj (and
        jj = tmpd(j).frame;           % don't save frames for that data element)
        if chk1(j).ismatch < 1/6
            reobj(jj) = cellfun(@(x1, x2) [x1, x2], reobj(jj), num2cell(tmpd(j).sub), 'UniformOutput',false);
        else
            jfr = [jfr, jj];
        end
    end
    chk1 = chk1(ibst); kt(k2) = false;
    if ~isnan(xyp2(k,5))   % check if not only match to this pillar (this can actually happen, if touching long tail section)
        k3 = xyp2(k,5); chktmp = chkc{k3,3};
        if chktmp.ismatch > chk1.ismatch || (chktmp.ismatch == chk1.ismatch && chktmp.rval > chk1.ismatch)
            chkc{k3,1} = [chkc{k3,1}, ktr]; chkc{k3,2} = [chkc{k3,2}, jfr]; % prev better match, log these frames and tracks and go to next object
            continue;
        else
            ktr = [chkc{k3,1}, ktr]; jfr = [chkc{k3,2}, jfr];               % else collect the prev frames and tracks to this one
        end
    end
    xyp2(k,5) = k2; xyp2(k,6) = ibst;
    chkc{k2,1} = ktr; chkc{k2,2} = jfr; chkc{k2,3} = chk1;  % need frames for pillar measurements, chk struct for pillar revision below
end

% create the revised pillar list; func getplr2 reruns chuff with higher resolution after preproc on obj's image. Note:
% we reject pillars at this point if they are too far off the screen (these tend to be more trouble then they are worth)
for i = 1:np
    k1 = xyp2(i,5); k2 = xyp2(i,6);
    if any(xyp2(i,1:2) - rad/2 < 1/2 | xyp2(i,1:2) + rad/2 > fliplr(imext) + 1/2)
        continue;
    elseif isnan(k1)  % no matching tracked object. If at an earlier stage we found one here then use getplr2 to create a stub
        if ~isnan(xyp2(i,4)) && ~isnan(plst(i).houghval)
            nplst{i} = getplr2([], xyp2(i,1:4), [], imext);
        end
    else
        tmpd = rmfield(stobj(k1).data(k2), 'sub');
        tmpd.composite = struct('Image',stobj(k1).composite, 'ij',stobj(k1).ij, 'frames',stobj(k1).frames);
        tmpd.track = sort(chkc{k1,1}); tmpd.frame = sort(chkc{k1,2});
        nplst{i} = getplr2(tmpd, xyp2(i,1:4), chkc{k1,3}, imext);
    end
end
nplst = [nplst{:}]';

% for now if anything survives to this point just add back to reobjs. Alternatively: we could go with original idea
% of a category of static non-pillar objects which get removed from worm images too, but need to see a case or two
% before proceeding in that direction ...
stobj = stobj(kt); nt = nnz(kt);
for i = 1:nt
    for j = 1:length(stobj(i).data)
        jj = stobj(i).data(j).frame;
        reobj(jj) = cellfun(@(x1, x2) [x1, x2], reobj(jj), num2cell(stobj(i).data(j).sub), 'UniformOutput',false);
    end
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function ps = getplr2(obj, crad, chkref, imext)
% ps = getplr2(obj, radobj, crad, chkref, imext)
% Refine pillar location with image object. Since we know (sort of) that this really has a tight image of a single
% pillar which is presumably not large, we can be a bit more aggressive with pre-proc and use a higher resolution.
% Output ps is an augmented pillar object (chuff output). Arguments: obj = image object; crad = 4-tuple giving relevant
% pillar info (ctr (1st 2 places), target radius, and pillar identifier); chkref = circhk output, contains info
% necessary for finding the interior and exterior boundaries; imext = full image extents, needed for pillar mask.

% Note: this is a critical subroutine
% Main issue: because of shadows, thresholding etc pillar rim will a) almost always be thicker than 1 pixel and b)
% through most but not all of its circumference thicker than 2 pixels. Because shadows in particular are directional the
% midpoint of this set of pixels is not necessarily the actual physical rim location. For compatibility with our solution
% to the "yin-yang" problem of tracking pillar locations when they are being displaced by the worm(s), we want a pretty
% accurate radius measurement here for something close to but not exactly the interior boundary of the thresholded pillar
% rim. Thus: before applying hough transform we extract the interior perimeter (where rim <= 3 pixels thick) and 2nd
% (ie. iterated) interior perimeter where the rim > 3 pixels thick.

ijbnds = reshape([max(1, round(crad([2 1]) - 2*crad(3))); min(imext, ceil(crad([2 1]) + 2*crad(3)))], 1, 4);
if isempty(obj)    % alt call: no matching tracked object found, so create a stub
    ps = struct('center',crad(1:2), 'radius',crad(3), 'houghval',NaN, 'pnum',crad(4), 'msk',[], 'ij',[], 'trkobj',[]);
    [ps.msk, ps.ij] = mkmask(ps.center, ps.radius + 2, ijbnds, -1, true);
    return;
end

imb = obj.Image; sz = size(imb); ijb = obj.ij;
pkpt = min(sz, max(1, round(crad([2 1]) - ijb([1 3]) + 1)));
rad = crad(3); pnum = crad(4);

% in order to get interior + iterated perimeters we need to isolate the hole of the obj (and remove any dangling
% projections ect). During circle - object check (circhk) gaps and projections have been identified, and the circhk
% output (chkref) contains needed info to clean up obj before proceeding.
if chkref.ismatch < 1
    imb(chkref.int.rmpts) = false; imb(chkref.ext.rmpts) = false;   % remove any projections etc from imb profile
    imb2 = imb; imb2(chkref.circ.adpts) = true;   % this fills gaps in imb's circumference so pick point fill will work
    imt = imfill(imb2, pkpt); imho = imfill(imt, 'holes') & ~imt;
    imt = imt & ~imb2; imb = imb | imho;    % get the main hole, + fill any other extraneous ones
else
    imt = imfill(imb, pkpt) & ~imb;
end

% % get the augmented interior boundary. Sections that are wider than 3 pixels still have material after we have removed
% % one outer boundary and two inner boundaries (equiv, remove the general perimeter + the twice dilated hole). We can
% % then use dilations to find where we want to use the 2nd boundary as opposed to the 1st boundary.
% strl1 = strel('disk', 1); imt1 = imdilate(imt, strl1); imt2 = imdilate(imt1, strl1);
% imp2 = bwmorph((imb & ~bwperim(imb)) & ~imt2, 'clean');
% impt = imt2 & imdilate(imp2, strl1); impt = ((imt1 & ~imdilate(impt, strl1)) | impt) & imb;

strl1 = strel('disk', 1); impt = imdilate(imb & imdilate(imt, strl1), strl1);

% run chuff with high resolution, retaining only the circle with the highest hough val. The mask field will be used to
% remove pillars from worm images once their deflection has been measured, so the circle it represents is a couple of
% pixels bigger in radius (ijbnds insures that the bounding box remains what it would be with no expansion).
ps = chuff(impt, rad*[8,10]/9, 1, 1/100);
if isempty(ps)
%     disp(' '); warning('(*) Unexpected empty chuff output during pillar location refinement.');
    return;
end
ps.center = ps.center + ijb([3 1]) - 1; ps.pnum = pnum;
[ps.msk, ps.ij] = mkmask(ps.center, ps.radius + 2, ijbnds, -1, true);
ps.trkobj = obj; ps.trkobj.circhk = chkref;

return;

% -------------------------------------------------------------------------------------------------------------------- %
function t2 = mrgtrks(t1s)
% t2 = mrgtracks(t1s)
% Merge tracks in array t1s into single track t2. For any particular frame, if more than one object is present than the
% lesser ones are made into subobjects of the major one. Once the data fields are merged the track metrics are
% recalculated as needed.

% start by ordering tracks so the most persistent object is primary; we cell out contents of data subfield by frame #
[~, ii] = sortrows([cellfun('length',{t1s.data}); t1s.Area], [-1 -2]);
t2 = t1s(ii); m = length(t2); newd = cell(0, m); n = 0;
for i = 1:m
    ii = [t2(i).data.frame];
    if ii(end) > n
        newd = cat(1, newd, cell(ii(end) - n, m)); n = ii(end);
    end
    newd(ii,i) = num2cell(t2(i).data);
end
isd = ~cellfun('isempty', newd);

% idea: shift everything leftward into col 1, putting further left objects into subobject fields as needed (also making
% sure to keep subobjects that are already present). btot keeps track of the # of sub's we've generated
jj = find(any(isd, 2))'; mj = length(jj); btot = 0;
for j = jj
    ii = find(isd(j,:)); tmpd = [newd{j,ii}]; tmpb = [tmpd.sub];
    if length(ii) > 1
        tmpb = [rmfield(tmpd(2:end), {'sub', 'frame'}), tmpb]; tmpd = tmpd(1);
    end
    tmpd.sub = tmpb; newd{j,1} = tmpd; btot = btot + length(tmpb);
end

% create the new track, and update fields if any of the lesser tracks had elements promoted
t2(1).track = sort([t2.track]); t2(1).subs = btot;
t2 = t2(1); t2.data = cell2mat(newd(jj, 1));
m1 = nnz(isd(:,1)); m2 = mj - m1;
if m2 > 0
    j2 = find(~isd(jj,1));
    t2.Area = (m1*t2.Area + sum([t2.data(j2).Area]))/mj;
    c2 = cell2mat({t2.data(j2).Centroid}'); ij2 = cell2mat({t2.data(j2).ij}');

    t2.Centroid = (m1*t2.Centroid + sum(c2, 1))/mj;
    r1 = t2.Range(1:2) + t2.Range(3:4); t2.Range(1:2) = min(t2.Range(1:2), min(c2,[],1));
    t2.Range(3:4) = max(r1, max(c2,[],1)) - t2.Range(1:2);

    t2.ij([1 3]) = min(t2.ij([1 3]), min(ij2(:,[1 3]), [], 1));
    t2.ij([2 4]) = max(t2.ij([2 4]), max(ij2(:,[2 4]), [], 1));

    ii = find(diff(jj) > 1); t2.frames = [jj(1) jj(ii+1); jj(ii) jj(end)]';
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function objnw = addobj(objld, objad)
% objnw = addobj(objld, objad)
% Combine the images in objld and objad under assumption that they are disjoint. objnw inherits its fields from objld;
% any that are affect by the op (Area, Centroid, etc) are updated as needed;

objnw = objld; imld = objnw.Image;
% check bounds, in most cases should just be able to OR objad's image onto a sub-section of imld
ijl = objld.ij; ija = objad.ij; ijdf = (ijl - ija) .* [1 -1 1 -1];
if any(ijdf > 0)    % need to expand the target image, requires padding + updating bounds etc.
    ijd2 = max(ijdf, 0); ijdf = min(ijdf, 0);
    imld = padarray(padarray(imld, ijd2([1 3]), 'pre'), ijd2([2 4]), 'post');
    ijl = ijl + ijd2 .* [-1 1 -1 1]; objnw.ij = ijl;
    objnw.BoundingBox = [ijl([3 1]) - 1/2, ijl([4 2]) - ijl([3 1]) + 1];
end
szd = size(imld); ii = (1-ijdf(1)):(szd(1)+ijdf(2)); jj = (1-ijdf(3)):(szd(2)+ijdf(4));
imld(ii,jj) = imld(ii,jj) | objad.Image; objnw.Image = imld;

objnw.Area = objld.Area + objad.Area;
objnw.Centroid = (objld.Area * objld.Centroid + objad.Area * objad.Centroid) / objnw.Area;

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [imd1, d2] = collectbb(d1)
% [imd1, d2] = collectbb(d1)
% Aggregate image data across frames by superimposition of pixels. This is done in two stages: collect data objects with
% the same bounding box and extract the count matrix for that range, then (after adjustment for different extents) take
% the sum of the count matrices. Returns the normalized count matrix imd1 as well as a modified data array d2, where
% each element contains the representative BW image for the bb-set, and the sub field contains the original data elements.

% Start by sorting the set by ij coords (equiv to bb's, but easier to manage)
nd = length(d1); [ijs, ~, kk] = unique(reshape([d1.ij], 4, nd)', 'rows', 'stable');
ij1 = [min(ijs(:,1)), max(ijs(:,2)), min(ijs(:,3)), max(ijs(:,4))];
sz1 = ij1([2 4]) - ij1([1 3]) + 1; imd1 = zeros(sz1);
nj = size(ijs, 1); d2 = repmat(d1(1), nj, 1);
for j = 1:nj
% get the count matrix (tmpim)
    jj = find(kk == j); nj = length(jj); tmpd = d1(jj(1));
    tmpd.frame = [d1(jj).frame]; tmpd.sub = rmfield(d1(jj), {'sub', 'frame'});
    tmpim = sum(cell2mat(reshape({d1(jj).Image}, [1 1 nj])), 3);
% update the overall track image
    dij = ij1 - ijs(j,:); ii = (1-dij(1)):(sz1(1)-dij(2));
    jj = (1-dij(3)):(sz1(2)-dij(4)); imd1(ii,jj) = imd1(ii,jj) + tmpim;

% get centroid for this bb-set (for d2) using the normalized count matrix
    tmpim = tmpim/nj; a1 = sum(tmpim(:));
    sx = sum(sum(bsxfun(@times, (ijs(j,3):ijs(j,4)), tmpim), 2), 1);
    sy = sum(sum(bsxfun(@times, (ijs(j,1):ijs(j,2))', tmpim), 1), 2);
    tmpd.Centroid = [sx, sy]/a1;
% convert the normalized count to BW image
    tmpim = tmpim > graythresh(tmpim);
    tmpd.Image = tmpim; tmpd.Area = nnz(tmpim);
    d2(j) = tmpd;
end    

% normalize the overall count matrix but do not convert
imd1 = imd1/nd;

return;

% -------------------------------------------------------------------------------------------------------------------- %
function sval = circhk(obj1, ctr, rad, imext)  % , rwid)
% sval  = circhk(obj1, ctr, rad, imext)
%
% Check how well image object obj1 matches to a circle defined by center ctr and radius rad. 4th argument imext
% gives full frame extents (to check if the circle would truncated within the frame, so that border hugging objects are
% not penalized). Returns a struct sval giving information about the quality of the match.
%  -- Field "ismatch": qualitative match, 0 = not a circle, 1 = full match (including location), values between 0 and 1
%   indicate formal deviations such as gaps in the circumference. If < 1/4 then object may be in contact with a worm
%   or other external object.
%  -- Field "rval" gives proportion of intersection area to union area wrt obj1 and an area-equivalent circle image,
%   which can be used to distinguish between objects with the same "ismatch" value.
% Remaining fields give specific information for recovering a proper closed circle from obj1 if 0 < ismatch <= 1/2.

% ismatch values: if ismatch < 1, then the set bits of 1/ismatch indicate the following cases:
%   5: significant material collected into a relatively small region outside of obj2 ( => being touched by external obj)
%   4: signif. gap in circumference ( => thresholding issue, missing section possibly being touched)
%   3: signif. material in interior, collected ( => possibly tail tip getting under bottom of pillar, or garbage/eggs)
%   2: insignificant gaps or extraneous material, or extra material close to obj2 ( => most likely dirt or image quality)
%   1: match indicator. Note: any number of these can be set, so minimum value for ismatch is 1/31.

sval = struct('ismatch',0, 'rval',0, 'iminfo',[], 'circ',[], 'int',[], 'ext',[]);  % , 'aratio',[0, 0]
if nargin < 1 || isempty(obj1)
    return;
elseif nargin < 4
    imext = [];
end
ismtyp = false(1,4);    % keep track of applicable categories, used later

% get relevant fields from image object
im1 = obj1.Image; bb1 = obj1.BoundingBox;
ij1 = obj1.ij; sz1 = size(im1); sval.iminfo.size = sz1;
a1 = nnz(im1); sval.iminfo.area = a1;
if rectint(bb1, [ctr - rad, 2*rad*ones(1,2)]) == 0    % no possible overlap, return (with non-match)
    return;
end

% to start: get angular range wrt circle params (need to know if there is significant truncation due to frame bounds),
% and create a reference circle image, then do the basic intersection area checks
c2ref = chkexts(ctr, rad, imext); sval.iminfo.circref = c2ref;
im2 = mkmask(ctr, rad, ij1, 1, true); a2 = nnz(im2); sval.iminfo.aref = a2;  % ref circle has rim width 3 (ie once dilated)
a3 = nnz(im1 & im2); sval.rval = a3/(a1 + a2 - a3);                          % compute the rval
if a3 <= (a1 * a2)/prod(sz1)    % ie match is no better than random over given area, reject
    return;
end
sval.circ = struct('npt',a3, 'ratio',a3/a2, 'gaps',[], 'inimg',[], 'adpts',[]);
if sval.circ.ratio < 0.40
    return;
end

% check im1 against the interior (imh) and exterior (imx) of our ref circle image as well
pkpt = min(max(1, round(fliplr(ctr) - ij1([1 3]) + 1)), sz1);   % pick point will always work on im2 regardless of truncation
imh = imfill(im2, pkpt); imx = ~imh; imh = imh & ~im2;

ah = nnz(imh); sval.iminfo.intref = ah; a1h = nnz(im1 & imh);
sval.int = struct('npt',a1h, 'ratio',a1h/ah, 'dmin',[], 'subobj',[], 'rmpts',[]);
if ~(sval.int.ratio <= 0.60)
    return;
end

ax = nnz(imx); sval.iminfo.extref = ax; a1x = nnz(im1 & imx);
sval.ext = struct('npt',a1x, 'ratio', a1x/ax, 'dmax',[], 'subobj',[], 'rmpts',[]);
if ~(sval.ext.ratio <= 0.60)   % this and above are based on empirical observation (a very thick rim can easily get above 1/2)
    return;
end

% preproc: since extraction of the interior boundary involves converting a hole boundary (returned by bwboundaries) into
% the inner object boundary, life will be a lot easier if we don't have other holes to worry about, so remove them now.
imho = imfill(im1, 'holes') & ~im1; nho = nnz(imho);
if nho > 0
    if nnz(imh & imho) > (ah - a1h)/2                            % large hole at center, this has to be removed from imho
        hotmp = regionprops(imho, {'Area','PixelIdxList'});
        if length(hotmp) > 1                       % multiple holes, remove the one with the most material in the center
            [~, ih] = max(cellfun(@(x) nnz(imh(x)), {hotmp.PixelIdxList}));
            imho(hotmp(ih).PixelIdxList) = false;
            im1 = im1 | imho; nho = nho - hotmp(ih).Area;
        else
            nho = 0; clear imho;                            % only one (big) hole, so no subsidiary holes to worry about
        end
        clear hotmp;
    else               % most of the empty part of center of im1 is not in imho, so all holes in imho are all removeable
        im1 = im1 | imho;
    end
else             % no holes whatsoever (gaps in circumference prevented imfill from working, but this can be dealt with)
    clear imho;
end

% idea: to check for projections we need points in order around interior and exterior boundaries. So while bwboundaries
% output format is kind of annoying, this is the place to start.
[idc, bmat, ncmp, cgraf] = bwboundaries(im1); chkgaps = false;
if any(full(cgraf(:)))
    [pnt, pxt] = getinterior(idc, full(cgraf), ij1, ctr);  % an interior exists, just need convert it from an outer hole boundary to an inner object boundary
else
    [pnt, pxt, trgaps] = findgaps(idc, ij1, ctr, rad);     % no interior, need to find the gaps and then separate the boundary ourselves.
    chkgaps = true;                                        % 3rd return val is gap due to truncation at image edge (no filler points added in this range)
end
if isempty(pnt) || isempty(pxt)                            % empty boundaries flag definitively non-circular input
    return;
end
mknt = ismember(pnt(1).ij, pxt(1).ij) & ~pnt(1).mrk; pnt(1).mrk(mknt) = 1;   % if points had to be added above they and some neighbors are marked,
mkxt = ismember(pxt(1).ij, pnt(1).ij) & ~pxt(1).mrk; pxt(1).mrk(mkxt) = 1;   % here we also mark any other parts where width goes down to 1 pixel (we want to protect these points too)

% there were gaps, need to record this to the circ subfields. gapinimg is length of gaps within the image (total gaps include
% those due to truncation to fit within the image, minus any gaps that result from pillar being on overall frame boundary)
if chkgaps
    imtmp = false(sz1); imtmp(pnt(1).ij(pnt(1).mrk == 2)) = true;
    objtmp = getsubobjs(imtmp, ctr, ij1); gtmp = runion(objtmp.Angle);
    sval.circ.inimg = gtmp; gapinimg = rsum(gtmp); angrange = 2*pi - rsum(trgaps);
    if ~isempty(trgaps)  % truncated image, have to check against ref truncation and add the gap ranges for full gap set
        sval.circ.gaps = rtol(1/2/rad, runion(gtmp, rdiff(trgaps, rcomp(c2ref))));
    else
        sval.circ.gaps = gtmp;
    end
    sval.circ.adpts = cell2mat({objtmp.PixelIdxList}');

%     if gapinimg/angrange > 1/3
%         disp(' '); warning('(*) test case for gap-angle ratio cutoff (%f), should check', gapinimg/angrange);
%     end
    if angrange - gapinimg < 1       % less than 1 radian of circum to work with, can't trust results
        return;
    elseif gapinimg/angrange > 0.5   % portion within the image is way too gappy, probably not a circle
        return;
    elseif rmax(sval.circ.gaps) > pi/3      % at least 1 significant gap, missing material may have been touching worm
        ismtyp(3) = true;
    elseif ~isempty(sval.circ.gaps)
        ismtyp(1) = true;   % small gaps, due to thresholding + difference in thickness of shadow around rim
    end
end

% check for interior and exterior projections. At this point the image essentially been accepted as a circle, we just
% need to know how much it deviates from the ideal ...
hprj = findproj(pnt(1), rad, -1);
if ~isempty(hprj) || length(pnt) > 1     % have interior projections or floating subobjects
    im1h = false(size(im1));
    for i = 1:size(hprj, 1)
        if isnan(hprj(i,1))    % flags extreme crappiness on boundary, so return with no match
            return;
        elseif hprj(i,1) == 0  % flags crappiness on boundary which doesn't amount to a real projection; effectively puts
            continue;          % into category 1 (messy but not pathological) if no other projections
        end
        [hid, hmrk] = proj2ind(hprj(i,:), pnt(1), ij1, ctr);
        im1h(hid) = true; pnt(1).mrk(hmrk) = 3;
    end
    for i = 2:length(pnt)
        im1h(pnt(i).ij) = true;
    end
    im1h = imfill(im1h, 'holes') & im1;   % need to AND here since sub-curve closure or fill might add some new points
    im1 = im1 & ~im1h;
    if nho > 0
        imtmp = im1h & imho; im1h = im1h & ~imtmp;   % if small holes filled earlier re-add now
        imho = imho & ~imtmp; nho = nnz(imho);
    end

    nlst = getsubobjs(bwlabel(im1h, 4), ctr, ij1);    % get 4-connected sub-objects, only want to log substantial items
    sval.int.rmpts = cell2mat({nlst.PixelIdxList}');  % as full subobjects (all pixel locations are saved though)
    sval.int.dmin = min([nlst.Distance]);
    ii = find([nlst.Area] > 3);
    if ~isempty(ii)
        sval.int.subobj = rmfield(nlst(ii), 'PixelIdxList'); ismtyp(2) = true;   % non-trivial interior projections
    else
        ismtyp(1) = true;                                                        % messy stuff sticking out
    end
end
if im1(pkpt(1), pkpt(2))   % if clearing interior projections leaves pick point occupied then shape is too weird to take
    return;                % a chance on using
end

xprj = findproj(pxt(1), rad, 1);
if ~isempty(xprj) || length(pxt) > 1   % pretty well the same basic procedure as above for interior objects
    im1x = false(size(im1));
    for i = 1:size(xprj, 1)
        if isnan(xprj(i,1))
            return;
        elseif xprj(i,1) == 0
            continue;
        end
        [xid, xmrk] = proj2ind(xprj(i,:), pxt(1), ij1, ctr);
        im1x(xid) = true; pxt(1).mrk(xmrk) = 3;
    end
    for i = 2:length(pxt)
        im1x(pxt(i).ij) = true;
    end
    im1x = imfill(im1x, 'holes') & im1; im1 = im1 & ~im1x;
    if nho > 0
        imtmp = im1x & imho; im1x = im1x & ~imtmp;
        imho = imho & ~imtmp; nho = nnz(imho);
    end

    xlst = getsubobjs(bwlabel(im1x, 4), ctr, ij1);
    sval.ext.rmpts = cell2mat({xlst.PixelIdxList}');
    sval.ext.dmax = max([xlst.Distance]);
    ii = find([xlst.Area] > 3);
    if ~isempty(ii)
        sval.ext.subobj = rmfield(xlst(ii), 'PixelIdxList'); ismtyp(4) = true;
    else
        ismtyp(1) = true;
    end
end
if all(pxt(1).mrk(im1(pxt(1).ij))) || all(pnt(1).mrk(im1(pnt(1).ij)))   % sanity check
%     disp(' '); warning('(*) only marked points left on boundary');
    return;
end

if nho > 0   % any remaining extra holes are on the circumference, update circ field + mark for messiness
    sval.circ.adpts = [sval.circ.adpts; find(imho)]; ismtyp(1) = true;
end
% now set the category
if any(ismtyp)
    sval.ismatch = 1/sum(2.^find(ismtyp));
else
    sval.ismatch = 1;
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [idx, inpp] = proj2ind(ids12, pp, ij0, ct0)
% [idx, inpp] = proj2ind(ids12, pp, ij0, ct0)
% Given indices ids12 to points on curve pp, along with subscript bounds ij0 and centerpoint ct0, calculate indices idx
% for a curve within ij0 containing points in pp between ids12(1) and ids12(2), plus whatever extra points are needed
% to close the curve. Added points are calculated to be consistent with curvature implied by values of pp's axy and dxy
% fields (so that they are convex for outward projections and concave for inward projections). 2nd return inpp is
% logical vector flagging points in pp which are also in the returned curve.

sz0 = diff(reshape(ij0, 2, 2)) + 1; np = length(pp.ij);
i1 = ids12(1); i2 = ids12(2); ii = mod(i1 + (0:mod(i2 - i1, np)) - 1, np) + 1;

% section belonging to pp (simple)
idx = pp.ij(ii); inpp = false(np, 1); inpp(ii) = true;

% points needed to close curve. Closing section is reversed since it will be appended to end of current idx.
pij = p1top2(pp.iy([i1; i2]), pp.jx([i1; i2]), ct0 - ij0([3 1]) + 1, sz0);
idnew = sub2ind(sz0, pij(2:end-1,1), pij(2:end-1,2));   % don't need 1st and last points
idx = [idx; flipud(idnew)];

return;

% -------------------------------------------------------------------------------------------------------------------- %
function intrpts = p1top2(pk, pj, ct1, sz1)
% intrpts = p1top2(pi, pj, ct1, sz1)
% Interpolate connecting curve between subscript points {pi(1),pj(1)} and {pi(2),pj(2)}, bowing out wrt center point
% ct1 (in xy coords) within image bounds sz1. pi and pj should be ints, but ct1 can be real and does not have to lie
% within the image. Returns subscripts (ij coords) for points in order between pij(1) and pij(2).

% Idea: we add midpoints in remaining gaps until the closing curve is connected. Points are interpolated wrt angle and
% distance from ct1.
cxy = complex(ct1(1), ct1(2)); pxy = complex(pj, pk); dj = find(abs(diff(pxy)) > 1.5);
while ~isempty(dj)
    intrp = repmat([1; 0], 1, length(pxy)); intrp(2, dj) = 2; intrp = intrp(intrp > 0);  % index locations of old,new pts
    ptmp = NaN(length(intrp), 1); ptmp(intrp == 1) = pxy;

    dtmp = [pxy(dj) - cxy, pxy(dj+1) - cxy];                            % interpolate points themselves
    dm = mean(abs(dtmp), 2); am = angle(dtmp); am = am(:,1) + mod(diff(am,1,2), 2*pi)/2; 
    ptmp(intrp == 2) = round(cxy + dm .* exp(complex(0, am)));
    pxy = ptmp; dj = find(abs(diff(pxy)) > 1.5);
end
intrpts = [max(1, min(sz1(1), imag(pxy))), max(1, min(sz1(2), real(pxy)))];
dchk = find(all(diff(intrpts, 1, 1) == 0, 2));
if ~isempty(dchk)
    intrpts(dchk, :) = [];
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function impts = findproj(pp, rd, dr)
% impts = findproj(pp, rd, dr)
% Check curve given by point set pts for projections in the direction perpendicular to pts' axy field. Returns a list of
% indices to projections (locations where projection starts and ends); if no projections are found then the list is
% empty. Remaining args rd and dr give standard height and direction to check to check (positive for exterior, negative
% for interior).

% idea: estimate a (normalized) derivative of height vs. angle; where the deriv spikes a projection is likely pushing
% out. Normalization involves mult by direction dr (so that projections are always preceeded by a +ve deriv), div by
% height rd (rescales for different magnifications etc), and suppression of -ve angles (to preserve sign of height
% diff, inf vals that will result in deriv are fine for our purposes)
impts = []; np = length(pp.ij); dd = dr*diff([pp.dxy; pp.dxy(1)]);
da = diff(unwrap([pp.axy; pp.axy(1)])); dv = (dd/rd) ./ max(da, 0);
dnorm = dr*(pp.dxy - rd); ismrk = pp.mrk > 0;
if any(dnorm(ismrk) > rd/2)      % dnorm = normalized heights used further down; check here whether worth pursuing in 1st case
    impts = NaN(1,2); return;
end
dnorm(ismrk) = NaN; jj = find(~ismrk);

% we start by dividing the point set into subsets of consecutive unmarked points with rising and failing heights.
jout = jj(dv(mod(jj - 2, np) + 1) > 0); jin = jj(dv(jj) <= 0);  % Note: since rising sets start after the diff we back the indices up by one position for the eval.
if isempty(jout) || isempty(jin)
%     disp(' '); warning('(*) no unmarked points in boundary profile');
    impts = NaN(1,2); return;
end
jout = mat2cell(jout, diff(find([1; diff(jout) > 1; 1])));
if length(jout) > 1 && (jout{1}(1) == 1 && jout{end}(end) == np)
    jout = cat(1, jout(2:end-1), {[jout{end}; jout{1}]});
end
jin = mat2cell(jin, diff(find([1; diff(jin) > 1; 1])));
if length(jin) > 1 && (jin{1}(1) == 1 && jin{end}(end) == np)
    jin = cat(1, jin(2:end-1), {[jin{end}; jin{1}]});
end

% sections of rising and failing dv vals should have matching last/first resp pt, in simplest case a single pair will make
% a small projection. Note: matchup in particular cases may be prevented by marked points, so may have to create null sets
kout = cellfun(@(x) x(end), jout); kin = cellfun(@(x) x(1), jin);
nogo = false(length(jout), 2);   % keep track of any null cells generated (which will be given nominal content)
if isequal(kout, kin)
    jall = cat(2, jout, jin); kall = kout;
elseif isequal(kout, circshift(kin, -1))
    jall = cat(2, jout, circshift(jin, -1)); kall = kout;
else
    jall = cat(2, jout, cell(size(jout))); [k1, i1] = ismember(kin, kout);
    jall(i1(k1), 2) = jin(k1); jall = cat(1, jall, cat(2, cell(nnz(~k1),1), jin(~k1)));
    [kall, ii] = sort([kout; kin(~k1)]); jall = jall(ii, :);
    nogo = cellfun('isempty', jall);    % give the null cells some content so that we don't have to worry about indexing etc
    jall(nogo(:,1),1) = num2cell(kall(nogo(:,1))); jall(nogo(:,2),2) = num2cell(kall(nogo(:,2)));
end
nj = size(jall, 1);

% collect relevant values for determining whether there are projections and where they are. dchk holds combined test
% results; these are for spikes in the height derivative (ie sudden change in boundary direction perpendicular to
% circumference) and sustained increases or decreases in height (boundary veering away from circumference).
% Note: tolerances are based on max neighbour and once-removed dists
tol1 = sqrt(2); tol2 = 2*tol1; tolhalf = tol1/2;
dvals = cell(nj,2); ddvls = cell(nj,2);
for i = 1:nj
    jtmp = mod(jall{i,1} - 2, np) + 1;   % for upward again need to back indices 1 pos
    dvals{i,1} = dv(jtmp); ddvls{i,1} = dd(jtmp);
    dvals{i,2} = -dv(jall{i,2}); ddvls{i,2} = -dd(jall{i,2});
end
dchk = cellfun(@(x1, x2) max(x1) > tol1 | sum(x2) > tol2, dvals, ddvls) & ~nogo;
imsct = find(any(dchk, 2));
if isempty(imsct)
    return;
end

% extract the best index to start / stop each projection from. This is where the first single pt rise of at least half
% basic tol (or first accumulated rise of 1/2 that tol) occurs, if such a rise does occur (it should for each
% flagged section in dchk, as well as some others). For falling sections this of course needs to be done backwards
ivals = zeros(nj, 2);
for i = imsct'
    tmpl = [dvals{i,1} > tolhalf | cumsum(ddvls{i,1}) > tol1; true];
    jtmp = [jall{i,1}; kall(i)]; ivals(i,1) = jtmp(find(tmpl, 1));
    tmpl = [flipud(dvals{i,2}) > tolhalf | cumsum(flipud(ddvls{i,2})) > tol1; true];
    jtmp = flipud([kall(i); jall{i,2}]); ivals(i,2) = jtmp(find(tmpl, 1));
end

% complex projections (which are not rare) may contain multiple rising / falling sections, so if there are multiple sections
% flagged we need to check if these belong to individual projections or should be combined. We do this by considering valley
% points (between falling and rising sections) which should be lower than points in projections. When actually greater
% we know we have a complex projection that needs to be constructed by merging smaller projections.
if length(imsct) > 1
    bvals = zeros(nj, 2); bvals(imsct,:) = dnorm(ivals(imsct,:));
    cvals = dnorm(mod(cellfun(@(x) x(end), jall(:,2)), np) + 1);         % cvals = points outside of either rises or falls
    
    isok = -double(bsxfun(@and, ~dchk(:,1), ~dchk(:,2)'));   % find matchups crossing cval pts (or not, if i = j)
    for i = find(dchk(:,1)')                                 % isok is real matrix so we can use -1 to mark "do not use
        for j = find(dchk(:,2)')                             % or check this match" (since better matchup overrides)
            if isok(i,j) < 0
                continue;
            elseif i == j
                isok(i,j) = 1; continue;
            end    % criterion: accept if the region where start > cval reachs the region where end > cvals (coming from the other side) 
            mc = mod(j - i, nj); ii = mod(i + (0:mc) - 1, nj) + 1;
            ic = find(~([cvals(ii(1:mc)); -Inf] >= bvals(ii(1),1)), 1) - 1;
            jc = find(~([cvals(ii(mc:-1:1)); -Inf] >= bvals(ii(end),2)), 1) - 1;
            if ic + jc > mc
                isok(ii,ii) = -1; isok(i,j) = 1;                  % if accept mark any combos within this pair as no-use
            end
        end
    end
    for i = find(~any(isok,2))'   % rising sections with no corresponding falling section (nor contained in any larger
        isok(i,i) = 1;            %  projection); promote the complementary section to make a full pair
    end
    for j = find(~any(isok,1))    % falling sections that still need corresponding rising sections, as above
        isok(j,j) = 1;
    end
    [i, j] = find(isok > 0); imsct = sortrows([i, j]);
else
    imsct = [imsct, imsct];
end
if any(any(diff(sort(imsct, 1), 1, 1) == 0))
%     disp(' '); warning('(*) problem encountered during circle verification (projection check)');  % shouldn't happen, so if this gets flagged there was a bug in the code somewhere
    impts = NaN(1,2); return;
end

% obtain the points for each section pair, and then check for certain unworkable (as is) projections. The 1st 2 situations
% arise from interactions with the edge of the image, and so we downgrade from projection to "messy" (the caller will
% make sure that the boundary gets flagged for non-optimality, but won't try to remove anything) unless it actually projects
% too far (then the caller just rejects the image). The 3rd can happen with thin meandering projections (which we want to
% remove), so we back up the pts a bit for a better aperture (if no backup works then tell caller to reject as well).
impts = [ivals(imsct(:,1),1), ivals(imsct(:,2),2)];
for i = 1:size(impts, 1)
    if nogo(imsct(i,1),1) || nogo(imsct(i,2),2)
        if max(dnorm(jall{imsct(i,1),1})) > rd/2 || max(dnorm(jall{imsct(i,2),2})) > rd/2
            impts = NaN(1,2); return;
        else
            impts(i,:) = zeros(1,2); continue;
        end
    end
    if mod(diff(impts(i,:)), np) > np/2
        ii = mod((impts(i,1):impts(i,1)+mod(diff(impts(i,:)), np)) - 1, np) + 1; achk = unwrap(pp.axy(ii));
        if diff(achk([1 end])) > pi
            if max(dnorm(ii)) > rd/2
                impts = NaN(1,2); return;
            else
                impts(i,:) = zeros(1,2); continue;
            end
        end
    end
    if diff(unwrap(pp.axy(impts(i,:)))) < 0
        ij1 = find(jall{imsct(i,1),1} == impts(i,1)); ij2 = find(jall{imsct(i,2),2} == impts(i,2)); 
        nj2 = length(jall{imsct(i,2),2}); mindy = -Inf; impts(i,:) = NaN(1,2);   % if nothing found can't save, so tell caller to reject
        for k1 = jall{imsct(i,1),1}(ij1:-1:1)'
            for k2 = jall{imsct(i,2),2}(ij2:nj2)'
                achk = diff(unwrap(pp.axy([k1,k2]))); hchk = min(dnorm([k2,k2]));  % 2nd criterion: want the new pair
                if achk > 0 && hchk > mindy                                        % to reduce the min height of the proj
                    mindy = hchk; impts(i,:) = [k1, k2];                           % as little as possible
                end
            end
        end
    end
end

return;

% ---------------------------------------------------------------------------------------------------------------------- 
function [pint, pext, truncd] = findgaps(blst, ij, ct, rd)
% [pint, pext] = getinterior(blst, ij, ct, rd)
% Given bwboundaries blst for an object with no hole (plus image extents ij, ref center ct, and ref radius rd), finds
% the boundary(s) belonging to the main object, divides the interior and exterior portions, fills the gaps, and returns
% arrays of interior and exterior boundary points (in the same manner as getinterior). The mrk field of the point
% structures will indicate points added to fill gaps. 3rd return truncd indicates truncation because the reference
% circle goes out of the image extents; in this range no filler points are added, but we will need the ranges for
% diagnostic purposes.

truncd = [];

% idea: we will first convert everything to point sets, since we will need distances and angles to do the calcs. We sort
% the boundaries into definitely interior subcomponents, definitely exterior subcomponents, and components that need to
% be split into interior and exterior sides; at the same time get minimum and maximum distances from ct point.
nb = length(blst); ptmp = []; pcat = zeros(1, nb);
drange = zeros(nb, 2); idr = zeros(nb, 2); dtol = rd/10;
for i = 1:nb
    pnxt = getpset(blst{i}, ij, ct);
    [drange(i,1), idr(i,1)] = min(pnxt.dxy);
    [drange(i,2), idr(i,2)] = max(pnxt.dxy);
    pcat(i) = (drange(i,1) - dtol > rd) - (drange(i,2) + dtol < rd);   % categories: -ve for interior, +ve for ext, 0 for needs split
    ptmp = [ptmp, pnxt];
end
if ~any(pcat == 0)   % don't think this can happen if program has made it this far, but just in case
%     disp(' '); warning('(*) No boundaries near target radius, should check');
    pint = []; pext = []; return;
end

% get angular extents and indices for category 0 boundaries. Initial point choice is conditioned on distance field (this
% will be revised later, when we have both angular endpoints together for comparison)
arange = NaN(nb, 2); ira = zeros(nb, 2); kk = find(~pcat);
for i = kk
    achk = unwrap(ptmp(i).axy); jj = find(abs(ptmp(i).dxy - rd) <= dtol);
    [a1, j1] = min(achk(jj)); [a2, j2] = max(achk(jj));
    arange(i,:) = mod([a1, a2] + pi, 2*pi) - pi; ira(i,:) = [jj(j1), jj(j2)];
end

% find order of the 0-cat boundaries around the center point, and shift so the largest gap is between
% the last and the first boundary
nk = length(kk);
if nk > 1
    [~, ii] = sort(arange(kk,1));
    [~, jj] = max(mod(arange(kk(circshift(ii, -1)),1) - arange(kk(ii),2) , 2*pi));
    kk = kk(circshift(ii, nk - jj));
end

% check for truncation, and if there is anything significant use alternate connecting sets at end sections. Note: it
% _is_ possible for truncation to happen at eg opposite sides at the same time, but in fact it is more likely that the
% circle will fall slightly out of the image due to off-centeredness etc. Hence we will only include truncation ranges
% that intersect with the maximum gap (set between 1st and last ptmp(kk) element above). Note: rfmt converts angle
% interval lists between different formats relating to how intervals wrapping around -pi/pi are treated (0 is standard).
istrunc = false; trtmp = rfmt(rcomp(chkexts(ct, rd, ij)), 1);
if ~isempty(trtmp)
    gp1 = rfmt([arange(kk(end),2), arange(kk(1),1)], 1);   % max angular gap in boundary sets
    tmprint = risintr(trtmp, [gp1 - 2*pi; gp1]); ii = find(tmprint > 1/2);
    if ~isempty(ii)
        if length(ii) > 1   % if multiple intervals have overlap combine them (so need to find where multi-interval effectively starts and ends)
            [~, itmp] = max(mod(trtmp(ii,1) - circshift(trtmp(ii,2), 1), 2*pi));
            trtmp = rfmt([trtmp(ii(itmp),1), trtmp(ii(mod(itmp-2,length(ii))+1),2)], 2);
        else
            trtmp = rfmt(trtmp(ii,:), 2);
        end
        istrunc = true;
    end
end

% find minimal point sets which connect consecutive boundary sections in their revised order
plink = cell(1, nk); idfix = cell(nk, 2);  % idfix index list is cell array in case multiple values need to be stored
for i = 1:nk-1
    k1 = kk(i); k2 = kk(i+1);
    [btmp, idfix{i,2}, idfix{i+1,1}] = connectingpset(ptmp(k1), ira(k1,2), ptmp(k2), ira(k2,1), ij, ct);
    plink{i} = getpset(btmp, ij, ct, 2);   % last arg: connecting sets are marked so they don't affect other boundary calcs
end
if ~istrunc   % if not truncated cycle goes full 360, need to connect last to first elements.
    k1 = kk(nk); k2 = kk(1);
    [btmp, idfix{nk,2}, idfix{1,1}] = connectingpset(ptmp(k1), ira(k1,2), ptmp(k2), ira(k2,1), ij, ct);
    plink{nk} = getpset(btmp, ij, ct, 2);
else    % have truncated circle, check boundary + estimated ends against image extents
    k1 = kk(1); k2 = kk(nk); plink = cat(2, {[]}, plink);
    [btmp, id1, id2] = connectingpset([], [], ptmp(k1), ira(k1,1), ij, ct, trtmp(2));
    plink{1} = getpset(btmp, ij, ct, 2); idfix{1,1} = [id1, id2];   % @ order of these pts ok?
    [btmp, id1, id2] = connectingpset(ptmp(k2), ira(k2,2), [], [], ij, ct, trtmp(1));
    plink{end} = getpset(btmp, ij, ct, 2); idfix{nk,2} = [id1, id2];
end

% now we have actual dividing points for each set, split them by interior / exterior division. Last arg of splitpset
% function = -1 ==> reverse order of points (so everything in consistent ccw (cw in image coords) direction)
pint = cell(1, nk); pext = cell(1,nk);
for i = 1:nk   % note: truncated circles might return 2 indices, in ccw order (so which is used where depends on if int or ext)
    k = kk(i);
    pext{i} = splitpset(ptmp(k), idfix{i,1}(end), idfix{i,2}(1), 1);
    pint{i} = splitpset(ptmp(k), idfix{i,2}(end), idfix{i,1}(1), -1);
end

% finally, combine the links with both pint and pext to convert into single connected point sets, then append any extra components, and return
if length(plink) > nk
    pint = cat(2, reshape(cat(1, plink(1:nk), pint), 1, 2*nk), plink(end));
    pext = cat(2, reshape(cat(1, plink(1:nk), pext), 1, 2*nk), plink(end));
else
    pint = reshape(cat(1, pint, plink), 1, 2*nk); pext = reshape(cat(1, pext, plink), 1, 2*nk);
end
pint = [combinepsets(pint), ptmp(pcat < 0)];
pext = [combinepsets(pext), ptmp(pcat > 0)];
if istrunc   % pint's extremal points will give most accurate range, after correction for pixel size
    truncd = rfmt([pint(1).axy(end), pint(1).axy(1)] + [1, -1]/2/rd, 0);
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [idset, ip1, ip2] = connectingpset(ps1, i1, ps2, i2, ij0, ct0, trang)
% [ptstruct, ip1, ip2] = connectingpset(ps1, i1, ps2, i2, ij0, ct0, trang)
% Find a line joining point sets ps1 (starting close to point i1) to ps2 (ending close to point i2). Exact connecting
% points are chosen to minimize point-to-point distance. Returns a subscript list (similar to bwboundaries output) which
% can be converted to a point structure with bounding box and center information. If either ps1 or ps2 are empty the
% line is drawn to the relevant edge of the image (image size being given in i1 or i2 as needed).

% Start by getting points in neighborhoods of i1 (far end of ps1) and i2 (near end of ps2) to consider for endpoints.
% tol is chosen to include points that are at most one step back and 1/10 deviation in distance from the ctr.
jj1 = []; jj2 = [];
if ~isempty(ps1)
    n1 = length(ps1.ij); ii = mod(i1 + (0:n1-1) - 1, n1) + 1;
    d1 = ps1.dxy(i1); achk = (ps1.axy(i1) - unwrap(ps1.axy(ii)) < 1.5/d1) & (abs(ps1.dxy(ii) - d1) <= d1/10);
    jj1 = ii([n1-find(flipud(~achk),1)+2:n1, 1:find(~achk, 1)-1]);
    if isempty(jj1)  % ie. all points satisfy tol; will only happen with small bits of cruft on circum
        jj1 = ii;
    end
end

% this is essentially ident to above, except that here we have the min-ang point instead of max
if ~isempty(ps2)
    n2 = length(ps2.ij); ii = mod(i2 + (0:n2-1) - 1, n2) + 1;
    d2 = ps2.dxy(i2); achk = (unwrap(ps2.axy(ii)) - ps2.axy(i2) < 1.5/d2) & (abs(ps2.dxy(ii) - d2) <= d2/10);
    jj2 = ii([n2-find(flipud(~achk),1)+2:n2, 1:find(~achk, 1)-1]);
    if isempty(jj2)
        jj2 = ii;
    end
end

sz0 = diff(reshape(ij0, 2, 2)) + 1; ct0b = ct0 - ij0([3 1]) + 1; 
if isempty(jj2)   % this and 2nd case involve truncated circle, will deal with below
    ps = ps1; jj0 = jj1; n0 = n1; d0 = d1; asgn = 1;
elseif isempty(jj1)
    ps = ps2; jj0 = jj2; n0 = n2; d0 = d2; asgn = -1;
else   % if have both ps1 and ps2 find the closest pair between the subsets and draw the line (using p1top2) between them
    dpairs = abs(bsxfun(@minus, ps1.xy(jj1), ps2.xy(jj2).'));
    [i, j] = find(dpairs == min(dpairs(:)));
    if length(i) > 1   % break ties by minimizing difference from original distances
        [~, k] = min(sum(abs([ps1.dxy(i), ps2.dxy(j)] - (d1+d2)/2), 2));
        i = i(k); j = j(k);
    end
    ip1 = jj1(i); ip2 = jj2(j);
    idset = p1top2([ps1.iy(ip1); ps2.iy(ip2)], [ps1.jx(ip1); ps2.jx(ip2)], ct0b, sz0);
    idset = idset(2:end-1,:); return;   % need to remove the 1st & last because already in ps
end

% If we're still here we're looking at one of the ends of the bdry for a truncated circle. Usually the boundary reaches
% the image edge, in which case we will just return the 2 indices (where the interior and exterior first encounter the
% edge). As well, since no filler pixels are required idset will be empty.
[~, i] = max(asgn * unwrap(ps.axy(jj0))); ip0 = jj0(i);   % since not doing compare, get the 
ang0 = ps.axy(ip0); id0 = [ps.iy(ip0), ps.jx(ip0)];

jjb = jj0(((ps.iy(jj0) == 1 | ps.iy(jj0) == sz0(1)) | (ps.jx(jj0) == 1 | ps.jx(jj0) == sz0(2))) & ...
    (ps.iy(jj0) == id0(1) | ps.jx(jj0) == id0(2)));
if ~isempty(jjb)  % ie pt ip0 or some near neighbor is on the image edge
    idset = []; [~, i] = max(asgn * unwrap(ps.axy(jjb)));
    ip0 = jjb(i); id0 = [ps.iy(ip0), ps.jx(ip0)];
    
    idchk = (id0 == 1) | (id0 == sz0); ii = mod(ip0 + (0:n0-1) - 1, n0) + 1;
    bchk = ~((ps.iy == id0(1) & idchk(1)) | (ps.jx == id0(2) & idchk(2)));
    kk = ii([(n0-find(bchk(fliplr(ii)),1)+2):n0, 1:(find(bchk(ii), 1)-1)]);
    if isempty(kk)   % all pts in bdry set on image edge, must short strip along edge
        if n0 > 2
%             disp(' '); warning('(*) test case, > 2 pixel section all on image edge');
            [~, itmp] = unique(ps.ij, 'stable');
            ip1 = itmp(1); ip2 = itmp(end); return;
        else
            ip1 = 1; ip2 = n0; return;  % need to put these at appropriate ends of segment, presumably with min(asgn * unwrap(ps.axy))
        end
    end
    tst1 = [ps.iy(kk(1)), ps.jx(kk(1))]; rechk1 = all(tst1 == 1 | tst1 == sz0);
    tst2 = [ps.iy(kk(end)), ps.jx(kk(end))]; rechk2 = all(tst2 == 1 | tst2 == sz0);
    if rechk1 || rechk2    % edge subset has reached corner, so extend around ...
        bchk = bchk & ~((ps.iy == tst1(1) | ps.jx == tst1(2)) & rechk1);
        bchk = bchk & ~((ps.iy == tst2(1) | ps.jx == tst2(2)) & rechk2);
        kk = ii([(n0-find(bchk(fliplr(ii)),1)+2):n0, 1:(find(bchk(ii), 1)-1)]);
    end
    ip1 = kk(1); ip2 = kk(end); return;
end

% we're connecting to the image edge, but we need to add some points to get there. So that downstream gap measurements
% are consistent we want to connect to a point on the edge close to where the reference circle is truncated; we then
% find the connecting points using the p1top2 function.
if asgn*diff(unwrap([ang0, trang])) <= 0   % precalculated truncation greater than actual, need to bring angle in a bit (don't want to go backwards)
    trang = ang0 + asgn/d0;
end
ebox = [1,1; sz0(2),1; sz0([2,1]); 1,sz0(1); 1,1];  % box around image edges
lintr = [ct0b; ct0b + sqrt(sum(sz0.^2))*[cos(trang), sin(trang)]];  % line that will intercept ebox at trunc angle
[idx2, idx1] = polyxpoly(lintr(:,1), lintr(:,2), ebox(:,1), ebox(:,2));  % use ML polyxpoly to find the intersection point
if asgn > 0
    idset = p1top2([id0(1); round(idx1)], [id0(2); round(idx2)], ct0b, sz0);
    idset = idset(2:end, :); ip1 = ip0; ip2 = [];
else  % in this case we flip the idset since we are going towards id0 (rather than coming from)
    idset = p1top2([round(idx1); id0(1)], [round(idx2); id0(2)], ct0b, sz0);
    idset = idset(1:end-1, :); ip1 = []; ip2 = ip0;
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [pint, pext] = getinterior(blst, amat, ij, ct)
% [pint, pext] = getinterior(blst, amat, ij, ct)
% If blst and amat are bwboundaries output for an image with a well-defined hole at the center, this separates the
% interior and exterior objects, converts the bwb hole boundary to an interior boundary of the object, and converts
% everything to point sets. Returns the interior and exterior boundaries as point set arrays, the 1st element of which
% will be the interior and exterior boundaries of the main object, with additional elements being boundaries of extra
% components.

% start by getting exterior boundaries via the adjacency matrix. Since we have already filled all holes but the circle
% interior this is rather simple
ii = find(~any(amat, 2))'; ll = any(amat(:, ii));    % main exterior boundary will be only blst element with children
ixt1 = ii(ll); pext = getpset(blst{ixt1}, ij, ct);   % but no parent, any other elements with no parent are extra components
for i = ii(~ll)
    pext = [pext, getpset(blst{i}, ij, ct)];
end

% main interior boundary consists of points in object immediately adjacent to the bwboundaries hole boundary; the hole
% is the main exterior boundary's child, and its children will be extra interior components. To get the object interior
% boundary we create the binary image with the boundary we want and then use bwtraceboundary on it.
in1 = find(amat(:,ixt1)); bj = blst{in1}; pt1 = bj(1,:) - [0,1];   % pt1 = starting point for obj inner boundary
sz = diff(reshape(ij, 2, 2)) + 1; imchk = false(sz);
imchk(sub2ind(sz, bj(:,1), bj(:,2))) = true;          % binary im of hole boundary
imchk = imdilate(imchk, strel('disk', 1)) & ~imchk;   % = 2 disconn curves, one containing pt1 is obj inner boundary
b2st = bwtraceboundary(imchk, pt1, 'E', 8, Inf, 'counterclockwise');   % 'E' => trace interior of curve

pint = getpset(b2st, ij, ct);
for i = find(amat(:,in1))'   % get any other components within interior
    pint = [pint, getpset(blst{i}, ij, ct)];
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function ptstruct = combinepsets(pcell)
% ptstruct = combinepsets(pcell)
% Combine all point sets given in cell array pcell into a single set ptstruct. Points + other info are all incorporated
% in the order given. As well, markers are set for elements now neighboring previously marked elements. If the set is
% closed then the points will be cycled such that the starting element is the leftmost then topmost (consistent with
% other boundaries returned by bwboundaries etc).

ii = find(~cellfun('isempty', pcell));    % for trunctated circles pointsets at ends could be empty
if isempty(ii)
    ptstruct = []; return;
end
pfc = fieldnames(pcell{ii(1)}); np = length(pfc);
nc = length(ii); C = cell(np, 1);
for i = 1:nc
    ctmp = struct2cell(pcell{ii(i)});
    for j = 1:np
        C{j} = cat(1, C{j}, ctmp{j});
    end
end

% check if closed curve, and if so cycle values so that left-top pixel is first. 4th field = "xy", complex points.
isclsd = false;
if abs(C{4}(end) - C{4}(1)) < 1.5
    isclsd = true;
    ii = find(C{3} == min(C{3}));                 % 3rd field = "jx", columns
    [~, jj] = min(C{2}(ii)); shft = 1 - ii(jj);   % 2nd field = "iy", row
    if shft ~= 0
        for j = 1:np
            C{j} = circshift(C{j}, shft);
        end
    end
end
ptstruct = cell2struct(C, pfc);

% adjust marker field; we want to add neighbors of mark 2 pts, plus any other pts adjacent to them that do not differ
% enough wrt angle (this is to prevent false positives during projection check)
nx = length(ptstruct.ij); mrkchk = (ptstruct.mrk == 2);
if isclsd
    mrkmo = any(bsxfun(@circshift, mrkchk, [-1 1]), 2) | mrkchk;
else
    mrkmo = ([mrkchk(2:end); 1] | [1; mrkchk(1:end-1)]) | mrkchk;   % if curve not closed mark endpoints as well
end
jj = find(mrkmo & ~mrkchk); achk = repmat(ptstruct.axy(jj),1,2);
atol = repmat(2/3./ptstruct.dxy(jj),1,2); ii = (1:numel(achk))';
for i = 1:nx
    jchk = mod([jj - i, jj + i] - 1, nx) + 1;
    if ~isclsd    % if not closed ends already marked, and can't go around them
        ii(jchk(ii) == 1 | jchk(ii) == nx) = [];
    end
    ii(mrkmo(jchk(ii))) = [];    % already marked ...
    ii(abs(diff(unwrap([ptstruct.axy(jchk(ii(:))), achk(ii(:))]'))') > atol(ii(:))) = [];   % outside of tol vals ...
    if isempty(ii)
        break;
    end
    mrkmo(jchk(ii)) = true;
end
ptstruct.mrk(mrkmo & ~mrkchk) = 1;

return;

% -------------------------------------------------------------------------------------------------------------------- %
function ptstruct = splitpset(ps0, i1, i2, direct)
% ptstruct = splitpset(ps0, i1, i2)
% Extract the point subset of ps0 between indices i1 and i2 (inclusive). If i1 > i2 then set wraps around the 1st index.
% Last arg direct, if negative, specifies that the set should be returned in reverse order (ie equiv to selecting i2:-1:i1)

np = length(ps0.ij); ptmp = struct2cell(ps0);
ii = mod(i1 + (0:mod(i2 - i1, np)) - 1, np) + 1;
if nargin > 3 && direct < 0
    ii = fliplr(ii);
end
for i = 1:length(ptmp)
    ptmp{i} = ptmp{i}(ii);
end
ptstruct = cell2struct(ptmp, fieldnames(ps0));

return;

% -------------------------------------------------------------------------------------------------------------------- %
function ptstruct = getpset(sbij, ijbnds, ct0, mval)
% ptstruct = getpset(sbij, ijbnds, ct0)
% Convert subscripts returned by bwboundaries into a point-distance-angle struct for an image contained in subscript
% range ijbnds wrt center ct0. Struct will contain indices, original subscripts, real points, angles, and distances in
% original order, as well as a marker vector (used to indicate pts added to eg close circumferential gaps; last arg
% allows caller to set a value for all marker elements on creation).

if isempty(sbij)
    ptstruct = []; return;
end
if nargin < 4 || isempty(mval)
    mval = 0;
end
szg = diff(reshape(ijbnds, 2, 2)) + 1; np = size(sbij,1);
if np > 1 && isequal(sbij(1,:), sbij(np,:))
    np = np - 1;
end
iy = sbij(1:np,1); jx = sbij(1:np,2); ij = sub2ind(szg, iy, jx);
xy = complex(jx + ijbnds(3) - 1, iy + ijbnds(1) - 1);
dxy = xy - complex(ct0(1), ct0(2)); axy = angle(dxy); dxy = abs(dxy);
ptstruct = struct('ij',ij, 'iy',iy, 'jx',jx, 'xy',xy, 'axy',axy, 'dxy',dxy, 'mrk',repmat(mval, np, 1));

return;

% -------------------------------------------------------------------------------------------------------------------- %
function objlst = getsubobjs(img, ct0, ij0)
% objlst = getsubobjs(img, ct0)
% Get connected subobjects from image img (can be bw, or label matrix) using regionprops. Each subobject struct is
% augmented with information about the angle and distance range wrt point ct0.

cxy = ct0 - ij0([3, 1]) + 1;
objlst = regionprops(img, {'Area', 'BoundingBox', 'Centroid', 'Extrema', 'PixelIdxList'});
[objlst(:).Angle] = deal([]); [objlst(:).Distance] = deal([]);
for i = 1:length(objlst)
    bb = objlst(i).BoundingBox; mind = Inf;
    if all(cxy >= bb(1:2)) && all(cxy <= bb(1:2) + bb(3:4))       % center point in bounding box, can't use extrema
        [iy, jx] = ind2sub(size(img), objlst(i).PixelIdxList);    % so we eval corners of all pixels
        xy = complex([jx-1/2, jx-1/2, jx+1/2, jx+1/2], [iy-1/2, iy+1/2, iy-1/2, iy+1/2]);
        xy = unique(xy(:));
        if img(round(cxy(2)), round(cxy(1)))  % if rounded location is occupied then min distance is 0
            mind = 0;
        end
    else
        xy = complex(objlst(i).Extrema(:,1), objlst(i).Extrema(:,2));
    end
    dxy = xy - complex(cxy(1), cxy(2));
    [axy, ii] = sort(angle(dxy)); dxy = abs(dxy(ii));
    mind = min(mind, min(dxy)); maxd = max(dxy); objlst(i).Distance = [mind, maxd];
    [mxang, j] = max(diff([axy; axy(1) + 2*pi]));
    if mxang < 1.5/maxd   % material in every direction from ct0
        objlst(i).Angle = [-pi, pi];
    elseif j == length(axy)
        objlst(i).Angle = [axy(1), axy(end)];
    else
        objlst(i).Angle = [-pi, axy(j); axy(j+1), pi];
    end
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function gps = chkexts(ct0, rd0, cij)
% gps = chkexts(ct0, rd0, imij)
% Check that the circle with center ct0 and radius rd0 is within the bounds given by imij (in image coords).
% Returns the angular range of the circle ([-pi pi] if circle is entirely within the subframe)

gps = [-pi, pi];
if nargin < 3 || isempty(cij)
    return;
elseif length(cij) < 4
    cij = [1 cij(1) 1 cij(2)];
end
tst1 = [round(ct0 - rd0) < cij([3 1]), round(ct0 + rd0) > cij([4 2])];
if ~any(tst1)
    return;
end

cxy = complex(ct0(1), ct0(2)); chgp = cell(1,4); sl = [Inf, 0, Inf, 0];
cbb = cij([3 1 4 2]); dang = pi*[2,-1,0,1]/2;
for i = find(tst1)
    [x1, y1] = linecirc(sl(i), cbb(i), ct0(1), ct0(2), rd0);
    a1 = sort(angle(complex(x1, y1) - cxy));
    if a1(1) < dang(i) && dang(i) < a1(2)
        chgp{i} = [a1(1) a1(2)];
    else
        chgp{i} = [-pi a1(1); a1(2) pi];
    end
end
gps = rcomp(runion(chgp{:}));

return;

% -------------------------------------------------------------------------------------------------------------------- %
%  Angle range interval set functions: Angular ranges are represented as [k x 2] ordered interval sets with all values
% between -pi and pi (an empty range is represented by an empty matrix). When a range straddles the pi,-pi value it will
% be represented as 2 sets (the 1st and the last). Operations on these sets are:
%     runion(a1, a2, ...): union of any number of angle range sets
%     rcomp(a1): complement of a1
%     rdiff(a1, a2): set theoretical a1 - a2
%     rtol: suppress intervals and gaps between intervals that are <= tol
%     rsum: scalar sum of interval lengths
%     rmax: scalar length of longest interval (intervals meeting at pi,-pi are treated as a single interval)
% ------------------------------
function urng = runion(varargin)    % union of ranges, redundancy check done by rtol is actually most involved aspect

urng = rtol(0, sortrows(cell2mat(varargin')));

return;

% ------------------------
function crng = rcomp(rng)    % complement of ranges, swap ends of intervals after appending -pi, pi

nr = size(rng, 1);
if nr == 0
    crng = [-pi, pi]; return;
elseif nr == 1 && (rng(1) == -pi && rng(2) == pi)
    crng = []; return;
end

tmp1 = rng(:,1); tmp2 = rng(:,2);
if tmp1(1) > -pi
    tmp2 = [-pi; tmp2];
else
    tmp1 = tmp1(2:end);
end
if tmp2(end) < pi
    tmp1 = [tmp1; pi];
else
    tmp2 = tmp2(1:end-1);
end
crng = [tmp2, tmp1];

return;

% -------------------------------
function drng = rdiff(rng1, rng2)    % set diff, = intersection of set 1 with complement of set 2

drng = rcomp(runion(rcomp(rng1), rng2));

return;

% ----------------------------
function trng = rtol(tol, rng)    % suppresses intervals and gaps <= tol. Intervals are assumed to already be in sorted
                                  % order, and in standard format (all within -pi,pi range)
if isempty(rng)
    trng = []; return;
end
trng = rng; ii = false(size(trng,1),1); k = 1;
for i = 2:size(trng, 1)    % start by removing the gaps (this has precedence since main application is runion overlap check)
    if trng(i,1) - trng(k,2) <= tol
        if trng(i,2) > trng(k,2)
            trng(k,2) = trng(i,2);
        end
        ii(i) = true;
    else
        k = i;
    end
end
trng = trng(~ii, :);
if trng(1,1) + 2*pi - trng(end,2) <= tol   % if end intervals have small gap effectively close it by bringing ends to pi val
    trng(1,1) = -pi; trng(end,2) = pi;
end

ii = diff(trng, 1, 2) <= tol;             % remove short intervals. The combined end intervals need to be separately
if trng(1,1) == -pi && trng(end,2) == pi  % check if there is wrapping
    if trng(1,2) + 2*pi - trng(end,1) > tol
        ii([1 end]) = false;
    end
end
trng = trng(~ii, :);

if isempty(trng)
    trng = [];
end

return;

% ---------------------
function trng = rfmt(rng, flg)    % convert between 3 alternate formats: 0 (standard) => all intervals fall between
                                  % -pi and pi (if there is an interval containing pi/-pi it is split); 1 => all
trng = rng;                       % intervals are represented on a single row with end >= beginning; 2 => all intervals
if isempty(rng)                   % on single rows, all angles <= pi, wrapping intervals have end < beginning
    return                        % Note: it is assumed that otherwise intervals are disjoint etc ...
end
% first put into standard format, then convert out if requested
if flg == 0
    if trng(end,2) > pi
        tmpr = [-pi, trng(end,2) - 2*pi]; trng(end,2) = pi; trng = [tmpr; trng];
    elseif trng(end,2) < trng(end,1)
        tmpr = [-pi, trng(end,2)]; trng(end,2) = pi; trng = [tmpr; trng];
    end
elseif flg == 1
    if trng(end,2) < trng(end,1)
        trng(end,2) = trng(end,2) + 2*pi;
    elseif size(trng,1) > 1 && (trng(1,1) == -pi && trng(end,2) == pi)
        trng(end,2) = trng(1,2) + 2*pi; trng = trng(2:end,:);
    end
else  % flg == 2
    if trng(end,2) > pi
        trng(end,2) = trng(end,2) - 2*pi;
    elseif size(trng,1) > 1 && (trng(1,1) == -pi && trng(end,2) == pi)
        trng(end,2) = trng(1,2); trng = trng(2:end,:);
    end
end

return;

% -------------------------------
function idrng = risintr(rng1, rng2)     % gives proportion of each interval of rng1 that intersect with rng2.

if isempty(rng1)
    idrng = logical([]); return;
elseif isempty(rng2)
    idrng = false(size(rng1,1), 1); return;
end
m1 = size(rng1); br1 = reshape([rng1(:,1), diff(rng1,1,2); ones(m1)], m1(1), 4);
m2 = size(rng2); br2 = reshape([rng2(:,1), diff(rng2,1,2); ones(m2)], m2(1), 4);
idrng = max(rectint(br1, br2), [], 2) ./ br1(:,3);

return;

% ---------------------
function sr = rsum(rng)    % scalar total length of intervals

if isempty(rng)
    sr = 0; return;
end
sr = sum(diff(rng, 1, 2));

return;

% ---------------------
function sr = rmax(rng)    % maximum interval length. Note possibility of wrapped intervals

if isempty(rng)
    sr = 0; return;
end
drng = diff(rng, 1, 2); sr = max(drng);
if size(rng,1) > 1 && (rng(1,1) == -pi && rng(end,2) == pi)
    sr = max(sr, drng(1) + drng(end));
end

return;

% -------------------------------------------------------------------------------------------------------------------- %
function [Tr, rplst, reobj] = trkplr(imext, numf, cdst, Tr, rplst, reobj)
% [Tr, rplst, reobj] = trkplr(imext, numf, cdst, Tr, rplst, reobj)
% Function does 2 things:
%  1. For worm tracks in track struct Tr, incorporates all data subobjects + all random small image objects in cell
%    array reobjs into the main data image for each frame if the subimage is close enough to the actual worm body in
%    the main image (this is done in order to reclaim eg tails etc that have been disconnected due to threshold effects).
%  2. for each revised pillar in rplst, for each frame where that pillar is being touched by a worm (so its image is
%    with the worm image in Tr) the deflection of the pillar is measured, and the pillar image is removed from the
%    worm image (masked out).
% Args imext, numf, and cdst give the image extents, number of frames, and distance tolerance for retaining sub- and
% re-objs with data images, respectively. Returns the modified worm track, the pillars augmented with deflection
% measurements, and any non-trivial reobjs that were not claimed by any of the worm images.

if nargin < 6 || isempty(reobj)
    reobj = cell(numf, 1);
end
dopillars = true;
if nargin < 5 || isempty(rplst)
    dopillars = false;
else
% setup for pillar deflection measurements: cpdefl <= valid deflection measurements; isprob <= low houghval measurements
% we don't want in official output; lpfr = per-frame indicator that measurement is needed; skipd <= log of skipped
% measurements (nothing found in image to measure even though indicated by lpfr, this could be benign)
    np = length(rplst); cpdefl = cell(numf, np); isprob = cell(numf, np);
    skipd = cell(1, np); lpfr = true(numf, np); 
    kk = ~isnan([rplst.houghval]);           % NaN houghval means no tracked object found for grid/pillar position ...
    if ~any(kk)
        disp(' '); warning('No base location could be established for any pillars.');
        rdrng = rplst(1).radius + [-1, 1];   % ... so we don't have base locations and radii vals for the pillars ...
    else
        rtmp = [rplst(kk).radius]; rdrng = [min(rtmp), max(rtmp)];
    end
    for i = find(kk)
        lpfr(rplst(i).trkobj.frame, i) = false;  % ... hence measurements need to be made for all frames ...
    end
    for i = find(~kk)
        rplst(i).radius = rdrng;                 % ... and multiple radii will need to be checked (by default chuff will
    end                                          %   use 1/2 pixel increments if given a radius range)
% tmp, log images belonging to pillars without base measurements
%     p4im = cell(1, np); 
end

numt = length(Tr); ltfr = false(numf, numt);           % grid out the worm data + sub objects, with indicator
dataobj = cell(numf, numt); subobj = cell(numf, numt);
for j = 1:numt
    tmpd = Tr(j).data; ii = [tmpd.frame]; ltfr(ii, j) = true;
    subobj(ii,j) = {tmpd.sub}; tmpd = rmfield(tmpd, 'sub');
    if dopillars
        [tmpd.pillar] = deal([]);                      % if tracking pillars add a cross-ref field
    end
    dataobj(ii,j) = num2cell(tmpd); Tr(j).data = [];   % temp clear orig data field to save memory
end

mc = ceil(cdst); strlx = strel('disk', mc); strl1 = strel('disk', 1);
for i = 1:numf
    tj = find(ltfr(i,:)); nt = length(tj);
    if ~nt
        continue;
    end
    tmpd = [dataobj{i,tj}]; im3 = false([imext, nt]);
    for j = 1:nt                              % collect the active track images for this frame
        ii = tmpd(j).ij(1):tmpd(j).ij(2); jj = tmpd(j).ij(3):tmpd(j).ij(4);
        im3(ii, jj, j) = bwunpack(tmpd(j).Image, tmpd(j).BoundingBox(4));
    end
    
% check the subobjects etc against the dilated worms to find those hovering just off the shoreline. If more than one
% worm gets a hit the sub is given to the worm whose dilated image has the most overlap with the sub. Check is
% cumulative (ie if a sub is added to a worm, the sub's dilation is added to the check images, thereby increasing the
% worms extent and the number of subs it can reach)
    tmpb = [subobj{i,:}, reobj{i}]; subobj(i,:) = {[]}; reobj{i} = [];
    mb = length(tmpb); done = false(1, mb); rechk = 1:nt; imchk = false([imext, nt]);
    while ~all(done) && ~isempty(rechk)
        imchk(:,:,rechk) = imdilate(im3(:,:,rechk), strlx);
        gotone = false(1, nt);
        for j = find(~done)
            ij = tmpb(j).ij; ii = ij(1):ij(2); jj = ij(3):ij(4);
            imtmp = tmpb(j).Image; imtst = bsxfun(@and, imchk(ii,jj,rechk), imtmp);
            [val, k] = max(sum(sum(imtst, 1), 2)); k = rechk(k);
            if val
                im3(ii,jj,k) = im3(ii,jj,k) | imtmp; gotone(k) = true; done(j) = true;
            end
        end
        rechk = find(gotone);
    end
    reobj{i} = tmpb(~done);

% pillar measurements. impchk is the perimeter image of all worms, on which chuff measurements will be done (padarray
% then reduction is done in order to suppress image bounds in perimeter); lpfr indicates which pillars need a measurement
% for this frame.
    if dopillars
        impchk = bwperim(padarray(any(im3, 3), [1 1], true, 'both'));
        impchk = imdilate(impchk(2:end-1,2:end-1), strl1);
        for j = find(lpfr(i,:))
            ij = rplst(j).ij; ii = ij(1):ij(2); jj = ij(3):ij(4);  % start by determining which worm pillar is touching
            [pval, jp] = max(sum(sum(bsxfun(@and, im3(ii,jj,:), rplst(j).msk), 1), 2));
            if pval < cdst                           % nothing at this location (could be that eg. image was returned
                skipd{j} = [skipd{j}, i]; continue;  % to reobj during prev stage but was not actually touching any worms)
            end
%             if isempty(rplst(j).trkobj)                    % if we don't have a base location (no associated tracked-object)
%                 p4im{j} = cat(3, p4im{j}, im3(ii,jj,jp));  % grab the sub-images for later examination
%             end

            ptmp = chuff(impchk(ii,jj), rplst(j).radius, 1);      % use circular hough transform to find pillar location
            if isempty(ptmp)               % divert untrustworthy results to isprob
                isprob{i,j} = struct('center',[], 'radius',[], 'houghval',[], 'frame',i); continue;
            elseif ptmp.houghval < 1/2/pi  % temporary measure, more stringent that prev 1/2/pi default
                isprob{i,j} = setfield(ptmp, 'frame', i); continue;
            end
            ptmp.center = ptmp.center + ij([3 1]) - 1; defl = ptmp.center - rplst(j).center;
            ptmp.wmid = tj(jp); ptmp.frame = i; ptmp.defl = defl;
            ptmp.dist = sqrt(sum(defl.^2)); ptmp.angle = atan2(defl(2), defl(1));
            cpdefl{i,j} = ptmp; tmpd(jp).pillar = [tmpd(jp).pillar, j];

            kd = round(fliplr(-defl)); kd = reshape([min(kd, 0); max(kd, 0)], 1, 4);   % + clear the deflected pillar from
            ii = ii(1-kd(1):end-kd(2)); jj = jj(1-kd(3):end-kd(4));                    % the worm image while we're here
            im3(ii,jj,jp) = im3(ii,jj,jp) & ~rplst(j).msk(1+kd(2):end+kd(1), 1+kd(4):end+kd(3));
        end
    end

    for j = 1:nt   % replace the original worm images with updated versions and place back into the dataobj array
        ii = find(any(im3(:,:,j), 2)); jj = find(any(im3(:,:,j), 1));
        tmpd(j).Image = bwpack(im3(ii(1):ii(end), jj(1):jj(end), j));
        ij = [ii(1) ii(end) jj(1) jj(end)]; tmpd(j).ij = ij;
        tmpd(j).BoundingBox = [ij([3 1]) - 1/2, ij([4 2]) - ij([3 1]) + 1];
        dataobj{i,tj(j)} = tmpd(j);
    end
end
for j = 1:numt
    Tr(j).data = [dataobj{:,j}]';
end
Tr = rmfield(Tr, 'subs'); clear dataobj subobj

% collect the measurements and collapse them into arrays.
if dopillars
    rplst = rmfield(rplst, 'msk');  % don't need this now we have the measurements and have removed the pillar from the worm images
    [rplst.defl] = deal([]);
    for j = 1:np
        rplst(j).defl = [cpdefl{:,j}];
        if ~isempty(skipd{j})             % skipped measurements: these are all likely due to images that were added back
            if ~isempty(rplst(j).trkobj)  % to reobj but not actually touching the worm, so revise frame lists
                rplst(j).trkobj.frame = union(rplst(j).trkobj.frame, skipd{j});
            else
                rplst(j).trkobj = struct('frame', skipd{j});
            end
        end
        if isnan(rplst(j).houghval)   % if no base location (for now) use measurements now to revise original estimate
            if ~isempty(rplst(j).defl)
                [~, k] = min([rplst(j).defl.dist]);
                ctk = rplst(j).defl(k).center; rplst(j).center = ctk;
                rplst(j).radius = mode([rplst(j).defl.radius]);

                ctmp = reshape([rplst(j).defl.center], 2, length(rplst(j).defl))';  % update deflection, dist, and angle ...
                dtmp = bsxfun(@minus, ctmp, ctk); atmp = num2cell(atan2(dtmp(:,2), dtmp(:,1)));
                rtmp = num2cell(sqrt(sum(dtmp.^2, 2))); dtmp = num2cell(dtmp, 2);
                [rplst(j).defl.defl] = deal(dtmp{:}); [rplst(j).defl.dist] = deal(rtmp{:});
                [rplst(j).defl.angle] = deal(atmp{:});
            else   % flag pathological case with NaN radius (to distinguish from NaN houghval pillars which gave legit measurements)
                rplst(j).radius = NaN;
            end
        end
        if ~isempty(rplst(j).defl)   % field is redundant once NaN-hougval elements have had radius field updated
            rplst(j).defl = rmfield(rplst(j).defl, 'radius');
        end
        jprob = [isprob{:,j}];   % don't add the isprob subfield unless actual problem measurements encountered
        if ~isempty(jprob)
            disp(' '); warning('Problems encountered with pillar %d; low houghval measurements have been diverted to "isprob" subfield', j);
            rplst(j).isprob = jprob;
        end
    end
    clear cpdefl isprob
%     [rplst.p4im] = deal(p4im{:});  % for now ...
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s2 = chkexparms(s1, u1)
% Check experimental parameters subfield of info struct. Basically, if it doesn't exist, create it with default values,
% and if non-empty args in user arguments struct u1 exist for any fields, update these to specified values. Put here to
% avoid distracting from flow of code in main body.

if ~isfield(s1, 'params') || isempty(s1.params)
    s2 = struct('mag',[], 'pixsz',[], 'mic2pix',[], 'wrmlen',[], 'diam',[], 'pdist',[], 'inpix',[], 'userargs',[]);
else
    s2 = s1.params;
end
m1 = u1.magnification; p1 = u1.pixelsize; w1 = u1.wormstage; d1 = u1.diameter; c1 = u1.spacing;

if ~isempty(m1)
    s2.mag = m1;
elseif isempty(s2.mag)
    s2.mag = 1;    % default magnification (for now) = 1:1
end
if ~isempty(p1)
    s2.pixsz = p1;
elseif isempty(s2.mag)
    s2.pixsz = 1;  % default pixel size (with above => 1:1 correspondence between pixels and microns)
end
s2.mic2pix = s2.mag / s2.pixsz;    % conversion from microns to pixels

% at this point need to convert from age to length scale (needed for tracking). Values derived from the Worm Atlas, also
% available in many std. texts eg the Oxford Handbook of Developmental Behavioral Neuroscience, Biology of Aging, etc...
if ~isempty(w1)
    stg2len = mean([250,250; 360,380; 490,510; 620,650; 900,940; 1110,1150], 2);
    if ischar(w1)
        j = find(strncmpi(w1, {'L1','L2','L3','L4','youngadult','adult'}, length(w1)));
        if length(j) < 1
            error(['Unknown worm stage "' w1 '"']);
        elseif length(j) > 1
            error(['Cannot resolve worm stage "' w1 '" between L1,L2,...']);
        end
        s2.wrmlen = stg2len(j);
    else
        s2.wrmlen = stg2len(w1);   % stage given as index, # => L#, with 5 and 6 indicating young adult and adult resp.
    end
elseif isempty(s2.wrmlen)
    s2.wrmlen = 1000;    % rough length of adult at 1:1 micron to pixel ...
end

if ~isempty(d1)   % default pillar diameter and spacing are empty (if they're not given we don't check for them)
    s2.diam = d1;
end
if ~isempty(c1)
    s2.pdist = c1;
end
s2.inpix = struct('radius',s2.mic2pix*s2.diam/2, 'grdspc',s2.mic2pix*s2.pdist, 'lensc',s2.mic2pix*s2.wrmlen);

% lastly, save the new user args (but keep any previous ones if nothing specified). Note: don't need the debugging fields.
u2 = rmfield(u1, {'verbose','debug','auxiliary'});
if ~isfield(s2, 'userargs') || isempty(s2.userargs)
    s2.userargs = u2;
else
    unam = fieldnames(u2);
    for i = 1:length(unam)
        if ~isempty(u2.(unam{i}))
            s2.userargs.(unam{i}) = u2.(unam{i});
        end
    end
end
  
return;

% ----------------------------------------------------------------------------------------------------------------------
function S1 = getargs(C)
% Parse input arguments.

% user args. [] is used to flag that default is to be applied.
sf =   {'verbose','debug','framerange','strip','magnification','pixelsize','wormstage','diameter','spacing','tolspace','threshold','medianfilter','auxiliary'};
indic = 1;
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
% set defaults relating to console output and debugging (these aren't saved or compared to anything)
if isempty(S1.verbose)
    S1.verbose = false;
end
if isempty(S1.debug)
    S1.debug = Inf;    % do all stages from beginning to end
end

return;

