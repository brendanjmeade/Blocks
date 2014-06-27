function [Block] = NewBlock(Block)
% Draw a new segment

% Select the first point
set(gcf, 'WindowButtonDownFcn', 'ButtonDown');
done                     = 0;
setappdata(gcf, 'doneClick', done);
while ~done
   done                  = getappdata(gcf, 'doneClick');
   [x, y]                = GetCurrentAxesPosition;
   set(findobj('Tag', 'Seg.pszCoords'), 'string', sprintf('(%7.3f)  %7.3f  ; %7.3f', npi2pi(x), x, y));
   drawnow;
end
[lonClose latClose]      = deal(x, y);

% Add new segment to structure Segment
newBlockName             = char(inputdlg('New block name:'));
Block.name               = strvcat(Block.name, newBlockName);
Block.eulerLon           = [Block.eulerLon ; 0];
Block.eulerLat           = [Block.eulerLat ; 0];
Block.eulerLonSig        = [Block.eulerLonSig ; 0];
Block.eulerLatSig        = [Block.eulerLatSig ; 0];
Block.interiorLon        = [Block.interiorLon ; lonClose];
Block.interiorLat        = [Block.interiorLat ; latClose];
Block.rotationRate       = [Block.rotationRate ; 0];
Block.rotationRateSig    = [Block.rotationRateSig ; 0];
Block.rotationInfo       = [Block.rotationInfo ; 0];
Block.aprioriTog         = [Block.aprioriTog ; 0];
Block.other1             = [Block.other1 ; 0];
Block.other2             = [Block.other2 ; 0];
Block.other3             = [Block.other3 ; 0];
Block.other4             = [Block.other4 ; 0];
Block.other5             = [Block.other5 ; 0];
Block.other6             = [Block.other6 ; 0];


% Update the blue Segment lines
nBlock                 = numel(Block.interiorLon);
plot(Block.interiorLon(nBlock), Block.interiorLat(nBlock), 'go', 'Tag', strcat('Block.', num2str(nBlock)), 'LineWidth', 1, 'MarkerFaceColor', 'g');



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
    % For Latitude, compute height of axes
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