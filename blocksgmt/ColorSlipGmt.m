function ColorSlipGmt(s, comp, file, varargin)
%
%  Write colored slip rates to a file to be plotted in GMT using:
%  psxy -M
%

if nargin == 4
   lims = varargin{:};
   sin = union(find(inpolygon(s.lon1, s.lat1, [min(lims(1:2)) max(lims(1:2)) max(lims(1:2)) min(lims(1:2))], [min(lims(3:4)) min(lims(3:4)) max(lims(3:4)) max(lims(3:4))])), ...
               find(inpolygon(s.lon2, s.lat2, [min(lims(1:2)) max(lims(1:2)) max(lims(1:2)) min(lims(1:2))], [min(lims(3:4)) min(lims(3:4)) max(lims(3:4)) max(lims(3:4))])));
   [s.lon1, s.lon2, s.lat1, s.lat2, s.ssRate, s.dsRate, s.tsRate] = deal(s.lon1(sin), s.lon2(sin), s.lat1(sin), s.lat2(sin), s.ssRate(sin), s.dsRate(sin), s.tsRate(sin));            
end

[path, name, ext, f] = fileparts(file);

% Make sure that the segments are in order


fid = fopen(file, 'w');
if comp == 1
   lims = [-1 1]*max(abs([min(s.ssRate) max(s.ssRate)]));
	cmap = redwhiteblue(256, lims);
	cidx = ceil(255*(s.ssRate + max(lims))./diff(lims) + 1);
	pnc2cpt(cmap, lims, [path filesep name '.cpt'])
	cvec = 255*cmap(cidx,:);
   for i = 1:numel(s.lon1)
      fprintf(fid, '> -W5p/%d/%d/%d\n%d %d\n%d %d\n', round(cvec(i, 1)), round(cvec(i, 2)), round(cvec(i, 3)), s.lon1(i), s.lat1(i), s.lon2(i), s.lat2(i));
   end
else
   nslips = s.dsRate - s.tsRate;
   lims = [-1 1]*max(abs([max(nslips) min(nslips)]));
	diffRate = diff(lims);
	nslips(find(nslips > max(nslips))) = max(nslips);
	nslips(find(nslips < min(nslips))) = min(nslips);
	cmap = bluewhitered(256, lims);
	cidx = ceil(255*(nslips + max(lims))./diffRate + 1);
	pnc2cpt(cmap, lims, [path filesep name '.cpt'])
	cvec = 255*cmap(cidx,:);
   for i = 1:numel(s.lon1)
      fprintf(fid, '> -W5p/%d/%d/%d\n%d %d\n%d %d\n', round(cvec(i, 1)), round(cvec(i, 2)), round(cvec(i, 3)), s.lon1(i), s.lat1(i), s.lon2(i), s.lat2(i));
   end
end

fclose(fid);
