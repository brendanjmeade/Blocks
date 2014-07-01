function SegmentManager
warning off; % This seems to be neccesary for Matlab2014 to prevent massive CPU usage
             % Warning: The EraseMode property is no longer supported and will error in a future release. Use the ANIMATEDLINE function for animating
             % lines and points instead of EraseMode 'none'. Removing instances of EraseMode set to 'normal', 'xor', and 'background' has minimal impact. 
             % > In MoveIntersection at 83
             %   In SegmentManagerFunctions at 405 

% Color variables
white                            = [1 1 1];
lightGrey                        = 0.85 * [1 1 1];
fn                               = 'Lucida';
fs                               = 9;

% I/O options
global GLOBAL ul cul st;
GLOBAL.filestream                = 1;
ul                               = 10; % number of navigation undo levels
cul                              = ul - 1; % current undo level
st                               = 2; % where to start counting the undo levels

% Open figure
h                                = figure('Position', [0 0 1200 850], 'Color', lightGrey, 'menubar', 'figure', 'toolbar', 'figure');
set(gcf, 'MenuBar', 'none');
set(gcf, 'ToolBar', 'none');

% Load command file and all listed files

commandYOffset = 100;
Seg.loadCommandFrame             = uicontrol('style', 'frame',          'position', [5 785 290 54],        'visible', 'on', 'tag', 'Seg.navFrame', 'BackgroundColor', lightGrey);
Seg.loadTextCommand              = uicontrol('style', 'text',           'position', [10 731+commandYOffset 50 15],    'visible', 'on', 'tag', 'Seg.loadTextCommand', 'string', 'Command File', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.loadEditCommand              = uicontrol('style', 'edit',           'position', [10 710+commandYOffset 280 20],   'visible', 'on', 'tag', 'Seg.loadEditCommand', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'Fontsize', 8);
Seg.loadPushCommand              = uicontrol('style', 'pushbutton',     'position', [10 690+commandYOffset 70 20],    'visible', 'on', 'tag', 'Seg.loadPushCommand', 'callback', 'SegmentManagerFunctions(''Seg.loadPushCommand'')', 'string', 'Load', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);


% (Segment Manager) File I/O controls
segmentYOffset = 35;
Seg.loadSegFrame                 = uicontrol('style', 'frame',          'position', [5 545+segmentYOffset 290 194],        'visible', 'on', 'tag', 'Seg.navFrame', 'BackgroundColor', lightGrey);
Seg.loadText                     = uicontrol('style', 'text',           'position', [10 731+segmentYOffset 63 15],    'visible', 'on', 'tag', 'Seg.loadText', 'string', 'Segment File', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.loadEdit                     = uicontrol('style', 'edit',           'position', [10 710+segmentYOffset 260 20],   'visible', 'on', 'tag', 'Seg.loadEdit', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'Fontsize', 8);
Seg.dispCheck                    = uicontrol('style', 'checkbox',       'position', [270 710+segmentYOffset 20 20],    'visible', 'on', 'tag', 'Seg.dispCheck', 'callback', 'SegmentManagerFunctions(''Seg.dispCheck'')', 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs, 'enable', 'off');
Seg.loadPush                     = uicontrol('style', 'pushbutton',     'position', [10 690+segmentYOffset 70 20],    'visible', 'on', 'tag', 'Seg.loadPush', 'callback', 'SegmentManagerFunctions(''Seg.loadPush'')', 'string', 'Load', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.clearPush                    = uicontrol('style', 'pushbutton',     'position', [80 690+segmentYOffset 70 20],    'visible', 'on', 'tag', 'Seg.clearPush', 'callback', 'SegmentManagerFunctions(''Seg.clearPush'')', 'string', 'Clear', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.savePush                     = uicontrol('style', 'pushbutton',     'position', [150 690+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.savePush', 'callback', 'SegmentManagerFunctions(''Seg.savePush'')', 'string', 'Save', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);

% (Segment Manager) Modify segment controls
%Seg.modText                      = uicontrol('style', 'text',           'position', [10 666+segmentYOffset 78 15],    'visible', 'on', 'tag', 'Seg.modText', 'string', 'Modify Segment', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.modSegList                   = uicontrol('style', 'popupmenu',      'position', [10 645+20+segmentYOffset 205 20],   'visible', 'on', 'tag', 'Seg.modSegList', 'callback', 'SegmentManagerFunctions(''Seg.modSegList'')', 'string', {''}, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.modSegPush                   = uicontrol('style', 'pushbutton',     'position', [220 645+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modSegPush', 'callback', 'SegmentManagerFunctions(''Seg.modSegPush'')','string', 'Update', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
segModPropList                   = {'< none >', 'Longitude 1', 'Latitude 1', 'Longitude 2', 'Latitude 2', 'Dip', 'Dip sigma', 'Dip flag', 'Locking depth', 'Locking depth sigma', 'Locking depth flag', 'Strike slip rate', 'Strike slip rate sigma', 'Strike slip rate flag', 'Dip slip rate', 'Dip slip rate sigma', 'Dip slip rate flag', 'Tensile slip rate', 'Tensile slip rate sigma', 'Tensile slip rate flag', 'Resolution', 'Resolution flag', 'Resolution other', 'Patch file', 'Patch toggle', 'Other 3', 'Patch slip file', 'Patch slip toggle', 'Other 6', 'Other 7', 'Other 8', 'Other 9', 'Other 10', 'Other 11', 'Other 12'};
Seg.modPropList                  = uicontrol('style', 'popupmenu',      'position', [10 620+20+segmentYOffset 140 20],   'visible', 'on', 'tag', 'Seg.modPropList', 'callback', 'SegmentManagerFunctions(''Seg.modPropList'')', 'string', segModPropList, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.modPropEdit                  = uicontrol('style', 'edit',           'position', [155 620+20+segmentYOffset 135 20],  'visible', 'on', 'tag', 'Seg.modPropEdit', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.modGSelect                   = uicontrol('style', 'pushbutton',     'position', [10 595+20+segmentYOffset 70 20],    'visible', 'on', 'tag', 'Seg.modGSelect', 'callback', 'SegmentManagerFunctions(''Seg.modGSelect'')','string', 'GSelect', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modGSelectM                  = uicontrol('style', 'pushbutton',     'position', [10 575+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modGSelectBox', 'callback', 'SegmentManagerFunctions(''Seg.modGSelectBox'')','string', 'GSelectBox', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modClear                     = uicontrol('style', 'pushbutton',     'position', [220 595+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modClear', 'callback', 'SegmentManagerFunctions(''Seg.modClear'')','string', 'Clear', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);

% Graphical Modify Segment Controls
Seg.modNewPush                   = uicontrol('style', 'pushbutton',     'position', [150 595+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modNewPush', 'callback', 'SegmentManagerFunctions(''Seg.modNewPush'')','string', 'New', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modDeletePush                = uicontrol('style', 'pushbutton',     'position', [80 595+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modDeletePush', 'callback', 'SegmentManagerFunctions(''Seg.modDeletePush'')','string', 'GDelete', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modDeletePushBox             = uicontrol('style', 'pushbutton',     'position', [80 575+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modDeletePushBox', 'callback', 'SegmentManagerFunctions(''Seg.modDeletePushBox'')','string', 'GDeleteBox', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modExtendPush                = uicontrol('style', 'pushbutton',     'position', [150 575+20+segmentYOffset 70 20],  'visible', 'on', 'tag', 'Seg.modExtendPush', 'callback', 'SegmentManagerFunctions(''Seg.modExtendPush'')','string', 'Extend', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modMovePush                  = uicontrol('style', 'pushbutton',     'position', [220 575+20+segmentYOffset 70 20],  'visible', 'on', 'tag', 'Seg.modMovePush', 'callback', 'SegmentManagerFunctions(''Seg.modMovePush'')','string', 'Move', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modConnectPush               = uicontrol('style', 'pushbutton',     'position', [150 555+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modConnectPush', 'callback', 'SegmentManagerFunctions(''Seg.modConnectPush'')','string', 'Connect', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modSplitPush                 = uicontrol('style', 'pushbutton',     'position', [220 555+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modSplitPush', 'callback', 'SegmentManagerFunctions(''Seg.modSplitPush'')','string', 'Split', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);

% Pulldown to show values
segModShowList                   = {'< none >', 'Name', 'Longitude 1', 'Latitude 1', 'Longitude 2', 'Latitude 2', 'Dip', 'Dip sigma', 'Dip flag', 'Locking depth', 'Locking depth sigma', 'Locking depth flag', 'Strike slip rate and sigma', 'Strike slip rate', 'Strike slip rate sigma', 'Strike slip rate flag', 'Dip slip rate and sigma', 'Dip slip rate', 'Dip slip rate sigma', 'Dip slip rate flag', 'Tensile slip rate and sigma', 'Tensile slip rate', 'Tensile slip rate sigma', 'Tensile slip rate flag', 'Resolution', 'Resolution flag', 'Resolution other', 'Patch file', 'Patch toggle', 'Other 3', 'Patch slip file', 'Patch slip toggle', 'Other 6', 'Other 7', 'Other 8', 'Other 9', 'Other 10', 'Other 11', 'Other 12'};
Seg.modShowList                  = uicontrol('style', 'popupmenu',      'position', [10 530+20+segmentYOffset 140 20],   'visible', 'on', 'tag', 'Seg.modShowList', 'callback', 'SegmentManagerFunctions(''Seg.modShowList'')', 'string', segModShowList, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.modGSelectL                  = uicontrol('style', 'pushbutton',     'position', [10 555+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modGSelectLasso', 'callback', 'SegmentManagerFunctions(''Seg.modGSelectLasso'')','string', 'GSelectLasso', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modGDeleteL		             = uicontrol('style', 'pushbutton',     'position', [80 555+20+segmentYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modGDeleteLasso', 'callback', 'SegmentManagerFunctions(''Seg.modGDeleteLasso'')','string', 'GDeleteLasso', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);


% (Block Manager) File I/O controls
blockYOffset = 50;
Seg.loadBlockFrame               = uicontrol('style', 'frame',          'position', [5 345+blockYOffset 290 174],        'visible', 'on', 'tag', 'Seg.navFrame', 'BackgroundColor', lightGrey);
Seg.loadTextBlock                = uicontrol('style', 'text',           'position', [10 511+blockYOffset 46 15],    'visible', 'on', 'tag', 'Seg.loadTextBlock', 'string', 'Block File', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.loadEditBlock                = uicontrol('style', 'edit',           'position', [10 490+blockYOffset 260 20],   'visible', 'on', 'tag', 'Seg.loadEditBlock', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'Fontsize', 8, 'FontName', fn, 'FontSize', fs);
Seg.dispCheckBlock               = uicontrol('style', 'checkbox',       'position', [270 490+blockYOffset 20 20],    'visible', 'on', 'tag', 'Seg.dispCheckBlock', 'callback', 'SegmentManagerFunctions(''Seg.dispCheckBlock'')', 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs, 'enable', 'off');
Seg.loadPushBlock                = uicontrol('style', 'pushbutton',     'position', [10 470+blockYOffset 70 20],    'visible', 'on', 'tag', 'Seg.loadPushBlock', 'callback', 'SegmentManagerFunctions(''Seg.loadPushBlock'')', 'string', 'Load', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.clearPushBlock               = uicontrol('style', 'pushbutton',     'position', [80 470+blockYOffset 70 20],    'visible', 'on', 'tag', 'Seg.clearPushBlock', 'callback', 'SegmentManagerFunctions(''Seg.clearPushBlock'')', 'string', 'Clear', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.savePushBlock                = uicontrol('style', 'pushbutton',     'position', [150 470+blockYOffset 70 20],   'visible', 'on', 'tag', 'Seg.savePushBlock', 'callback', 'SegmentManagerFunctions(''Seg.savePushBlock'')', 'string', 'Save', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);

% (Block Manager) Modify block controls
% Seg.modTextBlock                      = uicontrol('style', 'text',           'position', [10 666-220 61 15],    'visible', 'on', 'tag', 'Seg.modTextBlock', 'string', 'Modify Block', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.modSegListBlock                   = uicontrol('style', 'popupmenu',      'position', [10 645-200+blockYOffset 205 20],   'visible', 'on', 'tag', 'Seg.modSegListBlock', 'callback', 'SegmentManagerFunctions(''Seg.modSegListBlock'')', 'string', {''}, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.modSegPushBlock                   = uicontrol('style', 'pushbutton',     'position', [220 645-200+blockYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modSegPushBlock', 'callback', 'SegmentManagerFunctions(''Seg.modSegPushBlock'')','string', 'Update', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
segModPropListBlock                   = {'< none >', 'Name', 'Euler longitude', 'Euler latitude', 'Rotation rate', 'Euler longitude sigma', 'Euler latitude sigma', 'Rotation rate sigma', 'Rotation rate info', 'Interior longitude', 'Interior latitude', 'a priori flag', 'Other 1', 'Other 2', 'Other 3', 'Other 4', 'Other 5', 'Other 6'};
Seg.modPropListBlock                  = uicontrol('style', 'popupmenu',      'position', [10 620-200+blockYOffset 140 20],   'visible', 'on', 'tag', 'Seg.modPropListBlock', 'callback', 'SegmentManagerFunctions(''Seg.modPropListBlock'')', 'string', segModPropListBlock, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.modPropEditBlock                  = uicontrol('style', 'edit',           'position', [155 620-200+blockYOffset 135 20],  'visible', 'on', 'tag', 'Seg.modPropEditBlock', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
% <--- Modified by JPL, 1/24/2008 
Seg.modAddBlock                   = uicontrol('style', 'pushbutton',     'position', [150 595-200+blockYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modAddBlock', 'callback', 'SegmentManagerFunctions(''Seg.modAddBlock'')','string', 'Add', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modDeleteBlock                = uicontrol('style', 'pushbutton',     'position', [80 595-200+blockYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modDeleteBlock', 'callback', 'SegmentManagerFunctions(''Seg.modDeleteBlock'')','string', 'GDelete', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modDeleteBlockBox             = uicontrol('style', 'pushbutton',     'position', [80 575-200+blockYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modDeleteBlockBox', 'callback', 'SegmentManagerFunctions(''Seg.modDeleteBlockBox'')','string', 'GDeleteBox', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modMoveBlock                  = uicontrol('style', 'pushbutton',     'position', [220 575-200+blockYOffset 70 20],  'visible', 'on', 'tag', 'Seg.modMoveBlock', 'callback', 'SegmentManagerFunctions(''Seg.modMoveBlock'')','string', 'Move', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modGSelectBlock               = uicontrol('style', 'pushbutton',     'position', [10 595-200+blockYOffset 70 20],    'visible', 'on', 'tag', 'Seg.modGSelectBlock', 'callback', 'SegmentManagerFunctions(''Seg.modGSelectBlock'')','string', 'GSelect', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modGSelectBlockBox            = uicontrol('style', 'pushbutton',     'position', [10 575-200+blockYOffset 70 20],  'visible', 'on', 'tag', 'Seg.modGSelectBlockBox', 'callback', 'SegmentManagerFunctions(''Seg.modGSelectBlockBox'')','string', 'GSelectBox', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modClearBlock                 = uicontrol('style', 'pushbutton',     'position', [220 595-200+blockYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modClearBlock', 'callback', 'SegmentManagerFunctions(''Seg.modClearBlock'')','string', 'Clear', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
% Pulldown to show values
segModShowListBlock               = {'< none >', 'Name', 'Euler longitude', 'Euler latitude', 'Rotation rate', 'Euler longitude sigma', 'Euler latitude sigma', 'Rotation rate sigma', 'Rotation rate info', 'Interior longitude', 'Interior latitude', 'a priori flag', 'Other 1', 'Other 2', 'Other 3', 'Other 4', 'Other 5', 'Other 6'};
Seg.modShowListBlock              = uicontrol('style', 'popupmenu',      'position', [10 570-220+blockYOffset 140 20],   'visible', 'on', 'tag', 'Seg.modShowListBlock', 'callback', 'SegmentManagerFunctions(''Seg.modShowListBlock'')', 'string', segModShowListBlock, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);



% Start interface for loading meshes
meshYOffset = -360;
Seg.loadMeshFrame             = uicontrol('style', 'frame',          'position', [5 685+meshYOffset 290 54],    'visible', 'on', 'tag', 'Seg.meshFrame', 'BackgroundColor', lightGrey);
Seg.loadTextMesh              = uicontrol('style', 'text',           'position', [10 731+meshYOffset 50 15],    'visible', 'on', 'tag', 'Seg.loadTextMesh', 'string', 'Mesh File', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.loadEditMesh              = uicontrol('style', 'edit',           'position', [10 710+meshYOffset 260 20],   'visible', 'on', 'tag', 'Seg.loadEditMesh', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'Fontsize', 8);
Seg.dispCheckMesh             = uicontrol('style', 'checkbox',       'position', [270 710+meshYOffset 20 20],    'visible', 'on', 'tag', 'Seg.dispCheckMesh', 'callback', 'SegmentManagerFunctions(''Seg.dispCheckMesh'')', 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs, 'enable', 'off');
Seg.loadPushMesh              = uicontrol('style', 'pushbutton',     'position', [10 690+meshYOffset 70 20],    'visible', 'on', 'tag', 'Seg.loadPushMesh', 'callback', 'SegmentManagerFunctions(''Seg.loadPushMesh'')', 'string', 'Load', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.clearPushMesh             = uicontrol('style', 'pushbutton',     'position', [80 690+meshYOffset 70 20],    'visible', 'on', 'tag', 'Seg.clearPushMesh', 'callback', 'SegmentManagerFunctions(''Seg.clearPushMesh'')', 'string', 'Clear', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.snapPushMesh              = uicontrol('style', 'pushbutton',     'position', [150 690+meshYOffset 70 20],    'visible', 'on', 'tag', 'Seg.snapPushMesh', 'callback', 'SegmentManagerFunctions(''Seg.snapPushMesh'')', 'string', 'Snap segs.', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);

% Block checking options
integrityYOffset = -65;
Seg.navFrame                      = uicontrol('style', 'frame',          'position', [5 274+integrityYOffset 290 40],        'visible', 'on', 'tag', 'Seg.navFrame', 'BackgroundColor', lightGrey);
Seg.modCheckText                  = uicontrol('style', 'text',           'position', [10 541-240+integrityYOffset 60 20],    'visible', 'on', 'tag', 'Seg.modCheck', 'string', 'File integrity', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.modCheckSegsBlock             = uicontrol('style', 'pushbutton',     'position', [10 520-240+integrityYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modCheckSegsBlock', 'callback', 'SegmentManagerFunctions(''Seg.modCheckSegsBlock'')','string', 'Check closure', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modCheckIpsBlock              = uicontrol('style', 'pushbutton',     'position', [80 520-240+integrityYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modCheckIpsBlock', 'callback', 'SegmentManagerFunctions(''Seg.modCheckIpsBlock'')','string', 'Check-Ips', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modSegmentChecker             = uicontrol('style', 'pushbutton',     'position', [150 520-240+integrityYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modSegmentChecker', 'callback', 'SegmentManagerFunctions(''Seg.modSegmentChecker'')','string', 'Seg. Checker', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.modClearSegmentChecks         = uicontrol('style', 'pushbutton',     'position', [220 520-240+integrityYOffset 70 20],   'visible', 'on', 'tag', 'Seg.modClearSegmentChecks', 'callback', 'SegmentManagerFunctions(''Seg.modClearSegmentChecks'')','string', 'Clear checks', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
% --->

% (Segment Manager) Additional features frame
dispYOffset = 35;
Seg.navFrame                     = uicontrol('style', 'frame',          'position', [5 74+dispYOffset-5 290 95],        'visible', 'on', 'tag', 'Seg.navFrame', 'BackgroundColor', lightGrey);
Seg.dispText                     = uicontrol('style', 'text',           'position', [10 156+dispYOffset 36 15],     'visible', 'on', 'tag', 'Seg.dispText', 'string', 'Display', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);

% (Segment Manager) Load line file
Seg.dispEditLine                 = uicontrol('style', 'edit',           'position', [10 135+dispYOffset 70 20],    'visible', 'on', 'tag', 'Seg.dispEditLine', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'Fontsize', 8, 'FontName', fn, 'FontSize', fs);
Seg.dispPushLine                 = uicontrol('style', 'pushbutton',     'position', [85 135+dispYOffset 30 20],    'visible', 'on', 'tag', 'Seg.dispPushLine', 'string', 'Load', 'callback', 'SegmentManagerFunctions(''Seg.dispPushLine'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.dispTextLine                 = uicontrol('style', 'text',           'position', [120 135+dispYOffset 100 18],   'visible', 'on', 'tag', 'Seg.dispTextLine', 'string', 'line file', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.dispCheckLine                = uicontrol('style', 'checkbox',       'position', [170 135+dispYOffset 20 20],    'visible', 'on', 'tag', 'Seg.dispCheckLine', 'callback', 'SegmentManagerFunctions(''Seg.dispCheckLine'')', 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);

% (Segment Manager) Load xy file
Seg.dispEditXy                   = uicontrol('style', 'edit',           'position', [10 115+dispYOffset 70 20],    'visible', 'on', 'tag', 'Seg.dispEditXy', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'Fontsize', 8, 'FontName', fn, 'FontSize', fs);
Seg.dispPushXy                   = uicontrol('style', 'pushbutton',     'position', [85 115+dispYOffset 30 20],    'visible', 'on', 'tag', 'Seg.dispPushXy', 'string', 'Load', 'callback', 'SegmentManagerFunctions(''Seg.dispPushXy'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.dispTextXy                   = uicontrol('style', 'text',           'position', [120 115+dispYOffset 100 18],   'visible', 'on', 'tag', 'Seg.dispTextXy', 'string', 'xy file', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.dispCheckXy                  = uicontrol('style', 'checkbox',       'position', [170 115+dispYOffset 20 20],    'visible', 'on', 'tag', 'Seg.dispCheckXy', 'callback', 'SegmentManagerFunctions(''Seg.dispCheckXy'')', 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);

% (Segment Manager) Load station file
Seg.dispEditSta                  = uicontrol('style', 'edit',           'position', [10 95+dispYOffset 70 20],     'visible', 'on', 'tag', 'Seg.dispEditSta', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'Fontsize', 8, 'FontName', fn, 'FontSize', fs);
Seg.dispPushSta                  = uicontrol('style', 'pushbutton',     'position', [85 95+dispYOffset 30 20],     'visible', 'on', 'tag', 'Seg.dispPushSta', 'string', 'Load', 'callback', 'SegmentManagerFunctions(''Seg.dispPushSta'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.dispTextSta                  = uicontrol('style', 'text',           'position', [120 95+dispYOffset 100 18],    'visible', 'on', 'tag', 'Seg.dispTextSta', 'string', 'station file', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.dispCheckSta                 = uicontrol('style', 'checkbox',       'position', [170 95+dispYOffset 20 20],     'visible', 'on', 'tag', 'Seg.dispCheckSta', 'callback', 'SegmentManagerFunctions(''Seg.dispCheckSta'')', 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
% *** new addition of station vector viewing ***
Seg.dispTextStaVec               = uicontrol('style', 'text',           'position', [120 75+dispYOffset 100 18],    'visible', 'on', 'tag', 'Seg.dispTextSta', 'string', 'vectors', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.dispCheckStaVec		         = uicontrol('style', 'checkbox',       'position', [170 75+dispYOffset 20 20],     'visible', 'on', 'tag', 'Seg.dispCheckStaVec', 'callback', 'SegmentManagerFunctions(''Seg.dispCheckStaVec'')', 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.dipsScaleStaVec              = 0;
Seg.dispTextStaName              = uicontrol('style', 'text',           'position', [10 75+dispYOffset 100 18],    'visible', 'on', 'tag', 'Seg.dispTextStaNames', 'string', 'names', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.dispCheckStaName             = uicontrol('style', 'checkbox',       'position', [60 75+dispYOffset 20 20],     'visible', 'on', 'tag', 'Seg.dispCheckStaNames', 'callback', 'SegmentManagerFunctions(''Seg.dispCheckStaNames'')', 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);

% Slider to control velocity scaling (operates on ALL plotted vectors identically)
%Seg.velSlider                    = uicontrol('style', 'slider',         'position', [200 75+dispYOffset 90 20], 'min', 1e-6, 'max', 2, 'value', 1.0, 'visible', 'on', 'tag', 'Seg.velSlider', 'callback', 'SegmentManagerFunctions(''Seg.velSlider'')', 'BackgroundColor', white, 'HorizontalAlignment', 'left', 'Fontsize', 8, 'enable', 'on');
Seg.velPushUp                    = uicontrol('style', 'push', 'position', [200 75+dispYOffset 20 20], 'String', '+',  'visible', 'on', 'tag', 'Seg.velPushUp', 'callback', 'SegmentManagerFunctions(''Seg.velPushUp'')', 'BackgroundColor', white, 'HorizontalAlignment', 'left', 'Fontsize', 12, 'enable', 'on');
Seg.velPushDown                  = uicontrol('style', 'push', 'position', [220 75+dispYOffset 20 20], 'String', '-', 'visible', 'on', 'tag', 'Seg.velPushDown', 'callback', 'SegmentManagerFunctions(''Seg.velPushDown'')', 'BackgroundColor', white, 'HorizontalAlignment', 'left', 'Fontsize', 12, 'enable', 'on');

% (Segment Manager) Topo, Grid lines, and meridian popups
%Seg.dispTopo                     = uicontrol('style', 'popupmenu', 'position', [195 135+dispYOffset 100 20],    'visible', 'on', 'tag', 'Seg.dispTopo', 'callback', 'SegmentManagerFunctions(''Seg.dispTopo'')', 'string', {'Topo off', 'Topo on'}, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
%Seg.dispGrid                     = uicontrol('style', 'popupmenu', 'position', [195 115+dispYOffset 100 20],    'visible', 'on', 'tag', 'Seg.dispGrid', 'callback', 'SegmentManagerFunctions(''Seg.dispGrid'')', 'string', {'Grid off', 'Grid on'}, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
%Seg.dispDips                     = uicontrol('style', 'popupmenu', 'position', [195 95+dispYOffset 100 20],     'visible', 'on', 'tag', 'Seg.dispDips', 'callback', 'SegmentManagerFunctions(''Seg.dispDips'')', 'string', {'Dips off', 'Dips on'}, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.dispTopo                     = uicontrol('style', 'checkbox', 'position', [195 135+dispYOffset 80 20],    'visible', 'on', 'tag', 'Seg.dispTopo', 'callback', 'SegmentManagerFunctions(''Seg.dispTopo'')', 'string', {'Topo'}, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.dispGrid                     = uicontrol('style', 'checkbox', 'position', [195 115+dispYOffset 80 20],    'visible', 'on', 'tag', 'Seg.dispGrid', 'callback', 'SegmentManagerFunctions(''Seg.dispGrid'')', 'string', {'Grid'}, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);
Seg.dispDips                     = uicontrol('style', 'checkbox', 'position', [195 95+dispYOffset 80 20],     'visible', 'on', 'tag', 'Seg.dispDips', 'callback', 'SegmentManagerFunctions(''Seg.dispDips'')', 'string', {'Dips'}, 'BackgroundColor', lightGrey, 'FontName', fn, 'FontSize', fs);

% (Segment Manager) Navigate frame and navigation rose buttons
Seg.navFrame                     = uicontrol('style', 'frame',          'position', [5 5 290 90],        'visible', 'on', 'tag', 'Seg.navFrame', 'BackgroundColor', lightGrey);
Seg.navText                      = uicontrol('style', 'text',           'position', [10 87 44 15],       'visible', 'on', 'tag', 'Seg.navText', 'string', 'Navigate', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.navSW                        = uicontrol('style', 'pushbutton',     'position', [230 10 20 20],      'visible', 'on', 'tag', 'Seg.navSW', 'string', 'SW', 'callback', 'SegmentManagerFunctions(''Seg.navSW'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navS                         = uicontrol('style', 'pushbutton',     'position', [250 10 20 20],      'visible', 'on', 'tag', 'Seg.navS', 'string', 'S', 'callback', 'SegmentManagerFunctions(''Seg.navS'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navSE                        = uicontrol('style', 'pushbutton',     'position', [270 10 20 20],      'visible', 'on', 'tag', 'Seg.navSE', 'string', 'SE', 'callback', 'SegmentManagerFunctions(''Seg.navSE'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navW                         = uicontrol('style', 'pushbutton',     'position', [230 30 20 20],      'visible', 'on', 'tag', 'Seg.navW', 'string', 'W', 'callback', 'SegmentManagerFunctions(''Seg.navW'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navC                         = uicontrol('style', 'pushbutton',     'position', [250 30 20 20],      'visible', 'on', 'tag', 'Seg.navC', 'string', 'C', 'callback', 'SegmentManagerFunctions(''Seg.navC'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navE                         = uicontrol('style', 'pushbutton',     'position', [270 30 20 20],      'visible', 'on', 'tag', 'Seg.navE', 'string', 'E', 'callback', 'SegmentManagerFunctions(''Seg.navE'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navNW                        = uicontrol('style', 'pushbutton',     'position', [230 50 20 20],      'visible', 'on', 'tag', 'Seg.navNW', 'string', 'NW', 'callback', 'SegmentManagerFunctions(''Seg.navNW'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navN                         = uicontrol('style', 'pushbutton',     'position', [250 50 20 20],      'visible', 'on', 'tag', 'Seg.navN', 'string', 'N', 'callback', 'SegmentManagerFunctions(''Seg.navN'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navNE                        = uicontrol('style', 'pushbutton',     'position', [270 50 20 20],      'visible', 'on', 'tag', 'Seg.navNE', 'string', 'NE', 'callback', 'SegmentManagerFunctions(''Seg.navNE'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);

% (Segment Manager) longitude and latitude ranges
Seg.navTextLonRange              = uicontrol('style', 'text',           'position', [10 82-10 70 15],       'visible', 'on', 'tag', 'Seg.navText', 'string', 'Lon Range', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.navTextLatRange              = uicontrol('style', 'text',           'position', [80 82-10 44 15],       'visible', 'on', 'tag', 'Seg.navText', 'string', 'Lat Range', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'left', 'FontName', fn, 'FontSize', fs);
Seg.navEditLonMax                = uicontrol('style', 'edit',           'position', [10 50 70 20],       'visible', 'on', 'tag', 'Seg.navEditLonMax', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'right', 'FontName', fn, 'FontSize', fs);
Seg.navEditLonMin                = uicontrol('style', 'edit',           'position', [10 30 70 20],       'visible', 'on', 'tag', 'Seg.navEditLonMin', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'right', 'FontName', fn, 'FontSize', fs);
Seg.navEditLatMax                = uicontrol('style', 'edit',           'position', [80 50 70 20],       'visible', 'on', 'tag', 'Seg.navEditLatMax', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'right', 'FontName', fn, 'FontSize', fs);
Seg.navEditLatMin                = uicontrol('style', 'edit',           'position', [80 30 70 20],       'visible', 'on', 'tag', 'Seg.navEditLatMin', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'right', 'FontName', fn, 'FontSize', fs);
Seg.navUpdate                    = uicontrol('style', 'pushbutton',     'position', [150 50 70 20],      'visible', 'on', 'tag', 'Seg.navUpdate', 'string', 'Update', 'callback', 'SegmentManagerFunctions(''Seg.navUpdate'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navBack                      = uicontrol('style', 'pushbutton',     'position', [150 30 70 20],      'visible', 'on', 'tag', 'Seg.navBack', 'string', 'Back', 'callback', 'SegmentManagerFunctions(''Seg.navBack'')', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);

% (Segment Manager) Zoom options
Seg.navZoomIn                    = uicontrol('style', 'pushbutton',     'position', [10 10 70 20],       'visible', 'on', 'tag', 'Seg.navZoomIn', 'callback', 'SegmentManagerFunctions(''Seg.navZoomIn'')', 'string', 'Zoom In', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navZoomOut                   = uicontrol('style', 'pushbutton',     'position', [80 10 70 20],       'visible', 'on', 'tag', 'Seg.navZoomOut',  'callback', 'SegmentManagerFunctions(''Seg.navZoomOut'')', 'string', 'Zoom Out', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);
Seg.navZoomRange                 = uicontrol('style', 'pushbutton',     'position', [150 10 70 20],      'visible', 'on', 'tag', 'Seg.navZoomRange', 'callback', 'SegmentManagerFunctions(''Seg.navZoomRange'')', 'string', 'Zoom Range', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'FontName', fn, 'FontSize', fs);

% (Segment Manager) Segment file figure axes
Seg.axHandle                     = axes('parent', gcf, 'units', 'pixels', 'position', [360 100 1000 700],     'visible', 'on', 'Tag', 'Seg.axHandle', 'Layer', 'top', 'xlim', [0 360], 'ylim', [-90 90], 'FontName', fn);
SegmentManagerFunctions('DrawClean');
set(gca, 'Fontname', fn, 'FontSize', fs)

% (Segment Manager) Print, Save, Zumax
Seg.pszCoords                    = uicontrol('style', 'edit',           'position', [980 10 190 20],     'visible', 'on', 'tag', 'Seg.pszCoords', 'BackgroundColor', lightGrey, 'HorizontalAlignment', 'center', 'string', ' ', 'FontName', fn, 'FontSize', fs);

% Create handles structure for easy use in the callback later
Handles.Seg                      = Seg;
set(h, 'userdata', Handles);

% Making the GUI visible and give it a name
set(h, 'visible', 'on', 'name', 'Segment Manager');
set(gcf, 'DoubleBuffer', 'on');

