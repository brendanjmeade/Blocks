function [Segment] = MoveIntersection(Segment)
%%  Graphically move an intersection

%%  Dynamically highlight fault endpoints
lonEnd                   = [Segment.lon1 ; Segment.lon2];
latEnd                   = [Segment.lat1 ; Segment.lat2];
set(gcf, 'WindowButtonDownFcn', 'ButtonDown');
done                     = 0;
setappdata(gcf, 'doneClick', done);
minDIdxOld               = 0;
while ~done
   done                  = getappdata(gcf, 'doneClick');
   [x, y]                = GetCurrentAxesPosition;
   set(findobj('Tag', 'Seg.pszCoords'), 'string', sprintf('(%7.3f)  %7.3f  ; %7.3f', npi2pi(x), x, y));

   %%  Find the closest intersection
   d                     = sqrt((lonEnd - x).^2 + (latEnd - y).^2);
   [minDVal minDIdx]     = min(d);
   [lonClose, latClose]  = deal(lonEnd(minDIdx), latEnd(minDIdx));
   
   %%  Draw circle if neccesary
   if minDIdxOld == 0
      plot(lonClose, latClose, 'ro', 'Tag', 'HighlightIntersection');
   elseif (minDIdxOld ~= minDIdx) & (minDIdxOld ~= 0)
      set(findobj('Tag', 'HighlightIntersection'), 'XData', lonClose, 'Ydata', latClose);
   end
   drawnow;
   minDIdxOld            = minDIdx;
end
set(gcf, 'WindowButtonDownFcn', '');

%%  Find all of the segments that touch this point
matchLon1                = find(Segment.lon1 == lonClose);
matchLon2                = find(Segment.lon2 == lonClose);
matchLat1                = find(Segment.lat1 == latClose);
matchLat2                = find(Segment.lat2 == latClose);
matchIdx                 = intersect([matchLon1 ; matchLon2], [matchLat1 ; matchLat2]);

%%  Find the other endpoint coordinates
nMatch                   = numel(matchIdx);
xOther                   = zeros(size(matchIdx));
yOther                   = zeros(size(matchIdx));
for iMatch = 1 : nMatch
   if (Segment.lon1(matchIdx(iMatch)) == lonClose) & (Segment.lat1(matchIdx(iMatch)) == latClose)
      xOther(iMatch)     = Segment.lon2(matchIdx(iMatch));
      yOther(iMatch)     = Segment.lat2(matchIdx(iMatch));
   else
      xOther(iMatch)     = Segment.lon1(matchIdx(iMatch));
      yOther(iMatch)     = Segment.lat1(matchIdx(iMatch));
   end
end

%%  Draw the initial dynamic lines and delete the intersection marker
for iMatch = 1 : nMatch
   plot([xOther(iMatch) lonClose], [yOther(iMatch) latClose], 'r-', 'Tag', strcat('lineMove', num2str(iMatch)), 'LineWidth', 1);
end
delete(findobj('Tag', 'HighlightIntersection'));

%%  Move the lines till the next click
set(gcf, 'WindowButtonDownFcn', 'ButtonDown');
done                     = 0;
setappdata(gcf, 'doneClick', done);
while ~done
   done                  = getappdata(gcf, 'doneClick');
   [x, y]                = GetCurrentAxesPosition;
   set(findobj('Tag', 'Seg.pszCoords'), 'string', sprintf('(%7.3f)  %7.3f  ; %7.3f', npi2pi(x), x, y));
   for iMatch = 1 : nMatch
      set(findobj('Tag', strcat('lineMove', num2str(iMatch))), 'xData', [xOther(iMatch) x], 'yData', [yOther(iMatch) y], 'erasemode', 'xor', 'linestyle', '-', 'visible', 'on', 'linewidth', 1.0, 'color', 'r');
   end
   drawnow;
end
set(gcf, 'WindowButtonDownFcn', '');

%%  Update and move positions of old segment lines
for iMatch = 1 : nMatch
   if (Segment.lon1(matchIdx(iMatch)) == lonClose) & (Segment.lat1(matchIdx(iMatch)) == latClose)
      Segment.lon1(matchIdx(iMatch)) = x;
      Segment.lat1(matchIdx(iMatch)) = y;
   else
      Segment.lon2(matchIdx(iMatch)) = x;
      Segment.lat2(matchIdx(iMatch)) = y;
   end
   set(findobj('Tag', strcat('Segment.', num2str(matchIdx(iMatch)))), 'xData', [Segment.lon1(matchIdx(iMatch)) Segment.lon2(matchIdx(iMatch))], 'yData', [Segment.lat1(matchIdx(iMatch)) Segment.lat2(matchIdx(iMatch))], 'erasemode', 'xor', 'linestyle', '-', 'visible', 'on', 'linewidth', 1.0, 'color', 'k');
end

%%  Delete dynamic lines
for iMatch = 1 : nMatch
   delete(findobj('Tag', strcat('lineMove', num2str(iMatch))));
end



function [x, y] = GetCurrentAxesPosition
%%  GetCurrentAxesPosition
%%  Returns pointer position on current axes in units of pixels
%%  Authors: David Liebowitz, Seeing Machines
%%           Tom Herring, MIT

%%  Get dimension information
scnsize             = get(0, 'ScreenSize');
figsize             = get(gcf, 'Position');
axesize             = get(gca, 'Position');  % Could get CurrentAxes from gcf
llsize              = [get(gca, 'Xlim') get(gca, 'Ylim')];
asprat              = get(gca, 'DataAspectRatio');

%%  Based on the aspect ratio, find the actual coordinates coordinates covered by the axesize.
ratio               = (llsize(2) - llsize(1)) * asprat(2) / (llsize(4) - llsize(3));
if ratio > 1,   % Longitude covers the full pixel range
    xoff            = figsize(1) + axesize(1);
    xscl            = (llsize(2) - llsize(1)) / axesize(3); 
    %%  For Latitude, compute height of axes
    dyht            = (axesize(4) - axesize(4) / ratio) / 2;
    yoff            = figsize(2) + axesize(2) + dyht;
    yscl            = (llsize(4) - llsize(3)) / (axesize(4) / ratio);
else
    dxwd            = (axesize(3) - axesize(3) * ratio) / 2;
    xoff            = figsize(1) + axesize(1) + dxwd;
    xscl            = (llsize(2) - llsize(1)) / (axesize(3) * ratio); 
    yoff            = figsize(2) + axesize(2);
    yscl            = (llsize(4) - llsize(3)) / axesize(4);
end
xin                 = llsize(1);
yin                 = llsize(3);

% Construct the mapping array
pix2ll              = [xoff xscl xin ; yoff yscl yin];

%%  Get the pointer's screen position
pix                 = get(0, 'PointerLocation');
x                   = (pix(1) - pix2ll(1, 1)) * pix2ll(1, 2) + pix2ll(1, 3);
y                   = (pix(2) - pix2ll(2, 1)) * pix2ll(2, 2) + pix2ll(2, 3);



function ButtonDown
%%  Set an application data variable to indicate that a button has been clicked
setappdata(gcf, 'doneClick', 1); 