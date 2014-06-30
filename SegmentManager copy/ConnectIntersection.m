function [Segment] = ConnectIntersection(Segment)
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

%%  Draw initial dynamic line and second intersection
plot([lonClose lonClose], [latClose latClose], '-r', 'lineWidth', 1, 'Tag', 'lineMove')
plot(lonClose, latClose, 'ro', 'Tag', 'HighlightIntersection02');

%%  Dynamically highlight possible connections
set(gcf, 'WindowButtonDownFcn', 'ButtonDown');
done                     = 0;
setappdata(gcf, 'doneClick', done);
while ~done
   done                  = getappdata(gcf, 'doneClick');
   [x, y]                = GetCurrentAxesPosition;
   set(findobj('Tag', 'Seg.pszCoords'), 'string', sprintf('(%7.3f)  %7.3f  ; %7.3f', npi2pi(x), x, y));
   
   %%  Find the closest intersection
   d                     = sqrt((lonEnd - x).^2 + (latEnd - y).^2);
   [minDVal minDIdx]     = min(d);
   [lonCon, latCon]      = deal(lonEnd(minDIdx), latEnd(minDIdx));
   
   set(findobj('Tag', 'lineMove'), 'xData', [lonClose lonCon], 'yData', [latClose latCon], 'erasemode', 'xor', 'linestyle', '-', 'visible', 'on', 'linewidth', 1.0, 'color', 'r');
   set(findobj('Tag', 'HighlightIntersection02'), 'xData', [lonCon], 'yData', [latCon]);
   drawnow;
end
set(gcf, 'WindowButtonDownFcn', '');

%%  Add new segment to structure Segment
newSegmentName           = char(inputdlg('New segment name:'));
Segment                  = AddGenericSegment(Segment, newSegmentName, lonClose, latClose, lonCon, latCon);

%%  Update the blue Segment lines
nSegment                 = numel(Segment.lon1);
plot([Segment.lon1(nSegment) Segment.lon2(nSegment)], [Segment.lat1(nSegment) Segment.lat2(nSegment)], '-b', 'Tag', strcat('Segment.', num2str(nSegment)), 'LineWidth', 2);

%%  Delete dynamic line and intersection markers
delete(findobj('Tag', 'lineMove'));
delete(findobj('Tag', 'HighlightIntersection'));
delete(findobj('Tag', 'HighlightIntersection02'));



function [x, y] = GetCurrentAxesPosition
%%  GetCurrentAxesPosition
%%  Returns pointer position on current axes in units of pixels
%%  Authors: David Liebowitz, Seeing Machines
%%           Tom Herring, MIT

%%  Get dimension information
scnsize             = get(0, 'ScreenSize');
figsize             = get(gcf, 'Position');
axesize             = get(gca, 'Position');
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