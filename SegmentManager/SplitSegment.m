function Segment = SplitSegment(Segment)
%%  Select a segment and split it in two

%%  Loop until button is pushed and redraw lines
set(gcf, 'WindowButtonDownFcn', 'ButtonDown');
done                     = 0;
setappdata(gcf, 'doneClick', done);
minDIdxOld               = 0;

while ~done
   done                  = getappdata(gcf, 'doneClick');
   [x, y]                = GetCurrentAxesPosition;

   %%  Find the closest line (midpoint)
   lonMid                = (Segment.lon1 + Segment.lon2) / 2;
   latMid                = (Segment.lat1 + Segment.lat2) / 2;
   d                     = sqrt((lonMid - x).^2 + (latMid - y).^2);
   [minDVal minDIdx]     = min(d);
   
   %%  Change color if neccesary
   if minDIdxOld == 0
      set(findobj('Tag', strcat('Segment.', num2str(minDIdx))), 'Color', 'r');
   elseif (minDIdxOld ~= minDIdx) & (minDIdxOld ~= 0)
      set(findobj('Tag', strcat('Segment.', num2str(minDIdxOld))), 'Color', 'b');
      set(findobj('Tag', strcat('Segment.', num2str(minDIdx))), 'Color', 'r');
   end
   set(findobj(gcf, 'Tag', 'Seg.modSegList'), 'Value', minDIdx + 2);
   minDIdxOld            = minDIdx;
   drawnow;
end
set(gcf, 'WindowButtonDownFcn', '');
segIdx                   = minDIdx;

%%  Get the midpoint coordinates
lonMid                   = lonMid(segIdx);
latMid                   = latMid(segIdx);

%%  Copy the segment with a new name
Segment                  = CopySegmentProp(Segment, segIdx, strcat(Segment.name(segIdx, :), 'b'), lonMid, latMid, Segment.lon2(segIdx), Segment.lat2(segIdx));
Segment                  = CopySegmentProp(Segment, segIdx, strcat(Segment.name(segIdx, :), 'a'), Segment.lon1(segIdx), Segment.lat1(segIdx), lonMid, latMid);
Segment                  = DeleteSegment(Segment, segIdx);

%%  Let RedrawSegments Handle the deletion



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