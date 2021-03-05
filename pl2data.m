function D = pl2data(A, getdata)
% D = pl2data(A, getdata)
%
% Extract pillar tracking data from skel2n output A into a # frames by # pillars matrix D.
% By default deflections (pillars.defl.dist, converted to microns) are extracted. If any
% other quantity held in the defl subfield is desired this can be requested in the 2nd arg
% getdata (centers and deflections will be returned as complex points, centers in original
% pixel scale).

if nargin < 2 || isempty(getdata)
    getdata = 'dist';
else
    dnams = {'center','houghval','wmid','frame','defl','dist','angle'};
    j = find(strncmpi(getdata, dnams, length(getdata)));
    if length(j) ~= 1
        error(['Invalid field specification: ', getdata]);
    end
    getdata = dnams{j};
end

frng = A.source.framerange; nA = diff(frng) + 1;
np = length(A.pillars); D = NaN(nA, np);
for j = 1:np
    ii = []; jj = [];
    if ~isempty(A.pillars(j).defl)
        ii = [A.pillars(j).defl.frame]';
    else
        continue;
    end
    if isfield(A.pillars, 'isprob') && ~isempty(A.pillars(j).isprob)
        jj = [A.pillars(j).isprob.frame]';
    end
    
    if strcmp(getdata, 'center')
        D(:,j) = complex(A.pillars(j).center(1), A.pillars(j).center(2));
        xytmp = reshape([A.pillars(j).defl.center], 2, length(A.pillars(j).defl))';
        D(ii,j) = complex(xytmp(:,1), xytmp(:,2));
        D(jj,j) = complex(NaN,NaN);
    elseif strcmp(getdata, 'houghval')
        D(ii,j) = [A.pillars(j).defl.houghval]';
    elseif strcmp(getdata, 'wmid')
        D(ii,j) = [A.pillars(j).defl.wmid]';
    elseif strcmp(getdata, 'frame')
        D(ii,j) = ii + frng(1) - 1;
    elseif strcmp(getdata, 'defl')
        D(:,j) = complex(0, 0);
        xytmp = reshape([A.pillars(j).defl.defl], 2, length(A.pillars(j).defl))'/A.params.mic2pix;
        D(ii,j) = complex(xytmp(:,1), xytmp(:,2));
        D(jj,j) = complex(NaN,NaN);
    elseif strcmp(getdata, 'dist')
        D(:,j) = 0; D(ii,j) = [A.pillars(j).defl.dist]'; D(jj,j) = NaN;
    elseif strcmp(getdata, 'angle')
        D(ii,j) = [A.pillars(j).defl.angle]';
    end
end

return;
