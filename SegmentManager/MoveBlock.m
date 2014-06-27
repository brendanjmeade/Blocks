function [Block] = MoveIntersection(Block)
%  Graphically move an intersection

%  Dynamically highlight fault endpoints
lon                      = Block.interiorLon;
lat                      = Block.interiorLat;
set(gcf, 'WindowButtonDownFcn', 'ButtonDown');
done                     = 0;
setappdata(gcf, 'doneClick', done);
minDIdxOld               = 0;
while ~done
   done                  = getappdata(gcf, 'doneClick');
   [x, y]                = GetCurrentAxesPosition;
   set(findobj('Tag', 'Seg.pszCoords'), 'string', sprintf('(%7.3f)  %7.3f  ; %7.3f', npi2pi(x), x, y));

   % Find the closest intersection
   d                     = sqrt((lon - x).^2 + (lat - y).^2);
   [minDVal minDIdx]     = min(d);
   [lonClose  latClose]  = deal(lon(minDIdx), lat(minDIdx));
   
   % Draw circle if neccesary
   if minDIdxOld == 0
      plot(lonClose, latClose, 'ro', 'Tag', 'HighlightIntersection');
   elseif (minDIdxOld ~= minDIdx) & (minDIdxOld ~= 0)
      set(findobj('Tag', 'HighlightIntersection'), 'XData', lonClose, 'Ydata', latClose);
   end
   drawnow;
   minDIdxOld            = minDIdx;
   idx                   = minDIdx;
end
set(gcf, 'WindowButtonDownFcn', '');

%  Move the lines till the next click
set(gcf, 'WindowButtonDownFcn', 'ButtonDown');
done                     = 0;
setappdata(gcf, 'doneClick', done);
while ~done
   done                  = getappdata(gcf, 'doneClick');
   [x, y]                = GetCurrentAxesPosition;
   set(findobj('Tag', 'Seg.pszCoords'), 'string', sprintf('(%7.3f)  %7.3f  ; %7.3f', npi2pi(x), x, y));
   set(findobj('Tag', 'HighlightIntersection'), 'xData', x, 'yData', y, 'erasemode', 'xor', 'visible', 'on');
   drawnow;
end
set(gcf, 'WindowButtonDownFcn', '');

%  Update and move positions of old interior point
Block.interiorLon(idx) = x;
Block.interiorLat(idx) = y;
set(findobj('Tag', strcat('Block.', num2str(idx))), 'xData', Block.interiorLon(idx), 'yData', Block.interiorLat(idx), 'erasemode', 'xor', 'visible', 'on');
delete(findobj('Tag', 'HighlightIntersection'));



function [x, y] = GetCurrentAxesPosition
%  GetCurrentAxesPosition
%  Returns pointer position on current axes in units of pixels
%  Authors: David Liebowitz, Seeing Machines
%           Tom Herring, MIT

%  Get dimension information
scnsize             = get(0, 'ScreenSize');
figsize             = get(gcf, 'Position');
axesize             = get(gca, 'Position');  % Could get CurrentAxes from gcf
llsize              = [get(gca, 'Xlim') get(gca, 'Ylim')];
asprat              = get(gca, 'DataAspectRatio');

%  Based on the aspect ratio, find the actual coordinates coordinates covered by the axesize.
ratio               = (llsize(2) - llsize(1)) * asprat(2) / (llsize(4) - llsize(3));
if ratio > 1,   % Longitude covers the full pixel range
    xoff            = figsize(1) + axesize(1);
    xscl            = (llsize(2) - llsize(1)) / axesize(3); 
    %  For Latitude, compute height of axes
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

%  Get the pointer's screen position
pix                 = get(0, 'PointerLocation');
x                   = (pix(1) - pix2ll(1, 1)) * pix2ll(1, 2) + pix2ll(1, 3);
y                   = (pix(2) - pix2ll(2, 1)) * pix2ll(2, 2) + pix2ll(2, 3);



function ButtonDown
%  Set an application data variable to indicate that a button has been clicked
setappdata(gcf, 'doneClick', 1); 