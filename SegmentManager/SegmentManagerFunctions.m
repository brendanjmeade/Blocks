function SegmentManagerFunctions(option)
% SegmentManagerFunctions
%
% Functions called by Segment Manager GUI

% **** Big change: Switched from unlimited storage of the window bounds
% to only retaining the 10 most recent views.  Changed by JPL, 8 Jan 08.


% Declare variables
global GLOBAL ul cul st;
translateScale                     = 0.2;

% Parse callbacks
switch(option)
   
   % Start File I/O commands
   
   % Load .command file
   case 'Seg.loadPushCommand'
      % Delete all the childrean of the current axes
      fprintf(GLOBAL.filestream, '%s\n', option);
      
      % Get the name of the segment file
      ha                           = findobj(gcf, 'Tag', 'Seg.loadEditCommand');
      filename                     = get(ha, 'string');
      if exist(filename, 'file')
         filenameFull              = strcat(pwd, '\', filename);
      else
         [filename, pathname]      = uigetfile({'*.command'}, 'Load command file');
         if filename == 0
            set(ha, 'string', '');
            return;
         else
            set(ha, 'string', filename);
            filenameFull           = strcat(pathname, filename);
         end
      end
      
      % Read in the command file
      Command                      = ReadCommand(filenameFull)

      % Process .segment file
      ha                           = findobj(gcf, 'Tag', 'Seg.loadEdit');
      set(ha, 'string', Command.segFileName);
      Segment                      = ReadSegmentTri(Command.segFileName);
      setappdata(gcf, 'Segment', Segment);
      SegmentManagerFunctions('DrawSegments');
      set(findobj(gcf, 'Tag', 'Seg.modSegList'), 'string', cellstr(strvcat('< none >', 'Multiple', Segment.name)));
      set(findobj(gcf, 'tag', 'Seg.dispCheck'), 'enable', 'on', 'value', 1);

      % Process .block file
      ha                           = findobj(gcf, 'Tag', 'Seg.loadEditBlock');
      set(ha, 'string', 'Command.blockFileName');
      Block                        = ReadBlocksStruct(Command.blockFileName);
      Block                        = AlphaSortBlock(Block);
      setappdata(gcf, 'Block', Block);
      set(findobj(gcf, 'Tag', 'Seg.modSegListBlock'), 'string', cellstr(strvcat(' ', 'Multiple', Block.name)));
      nBlocks                      = numel(Block.interiorLon);
      blnames                      = cellstr([repmat('Block.', nBlocks, 1) strjust(num2str([1:nBlocks]'), 'left')]);
      plotbls                      = line([Block.interiorLon'; Block.interiorLon'], [Block.interiorLat'; Block.interiorLat'], 'MarkerFaceColor', 'm', 'MarkerSize', 5, 'marker', 'o', 'linestyle', 'none', 'MarkerEdgeColor', 'k');
      set(plotbls, {'tag'}, blnames);
      set(findobj(gcf, 'tag', 'Seg.dispCheckBlock'), 'enable', 'on', 'value', 1);

      % Process the .sta.data file
      ha                           = findobj(gcf, 'Tag', 'Seg.dispEditSta');
      set(ha, 'String', Command.staFileName);
      Station = PlotSta(Command.staFileName);
      PlotStaVec(Station)
      setappdata(gcf, 'Station', Station);
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckSta');
      set(hb, 'Value', 1);

      % Process the .msh files...only works for one right now
      P = ReadPatches(Command.patchFileNames);
      h = patch('Vertices', P.c, 'faces', P.v, 'facecolor', 'g', 'edgecolor', 'black');
      if ~isempty(Command.patchFileNames)
         set(findobj(gcf, 'tag', 'Seg.dispCheckMesh'), 'enable', 'on', 'value', 1);
      end


   
   % Load segment file
   case 'Seg.loadPush'
      % Delete all the childrean of the current axes
      fprintf(GLOBAL.filestream, '%s\n', option);
      
      % Get the name of the segment file
      ha                           = findobj(gcf, 'Tag', 'Seg.loadEdit');
      filename                     = get(ha, 'string');
      if exist(filename, 'file')
         filenameFull              = strcat(pwd, '\', filename);
      else
         [filename, pathname]      = uigetfile({'*.segment; *.segment.xml'}, 'Load segment file');
         if filename == 0
            set(ha, 'string', '');
            return;
         else
            set(ha, 'string', filename);
            filenameFull           = strcat(pathname, filename);
         end
      end
      
      % Read in the segment file
      Segment                      = ReadSegmentTri(filenameFull);

      setappdata(gcf, 'Segment', Segment);
      
      % Plot segments file
      SegmentManagerFunctions('DrawSegments');
      
      % Add the segment names to the segment pulldown list
      set(findobj(gcf, 'Tag', 'Seg.modSegList'), 'string', cellstr(strvcat('< none >', 'Multiple', Segment.name)));
      
      % Enable the display check box
      set(findobj(gcf, 'Tag', 'Seg.dispCheck'), 'enable', 'on', 'value', 1);
   
   	% Segment display toggle
   case 'Seg.dispCheck'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, '-regexp', 'Tag', '^Segment.');
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheck');
      if get(hb, 'Value') == 0;
      	set(ha, 'Visible', 'off');
      else
         set(ha, 'Visible', 'on');
      end   
      
   % Clear segment file
   case 'Seg.clearPush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                              = findobj(gcf, 'Tag', 'Seg.loadEdit');
      set(ha, 'string', '');
      setappdata(gcf, 'Segment', []);
      delete(findobj(gcf, '-regexp', 'tag', 'Segment.\d'));
      % Disable the display check box
      set(findobj(gcf, 'Tag', 'Seg.dispCheck'), 'enable', 'off', 'value', 0);
      
      
      
   % Save segment file
   case 'Seg.savePush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                     = getappdata(gcf, 'Segment');
      if size(Segment, 1) == 0
         return;
      else
         [filename, pathname]     = uiputfile({'*.segment; *.segment.xml'}, 'Save segment file');
         if filename == 0
            return;
         else
            filenameFull          = strcat(pathname, filename);
            WriteSegmentStruct(filenameFull, Segment);
         end
         set(findobj(gcf, 'tag', 'Seg.loadEdit'), 'string', ['  ' filename]);
      end
   
      
     
   %  Start Modify Commands
   case 'Seg.modSegPush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         valueStr                  = get(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'String');
         try % Convert string to numerical value
            value                  = str2double(valueStr);            
         catch
            fprintf(GLOBAL.filestream, 'Could not convert to a numerical value\n');
         end

         % Figure out which field has and which segments have been selected...then update
         segIdx                 = getappdata(gcf, 'segIdx');
         propIdx                = get(findobj(gcf, 'Tag', 'Seg.modPropList'), 'Value');
         if propIdx <= 2
            fprintf(GLOBAL.filestream, 'Select a property with a numerical value\n');
         end

         % Set and save properties
         if propIdx == 2; Segment.lon1(segIdx) = value; end;
         if propIdx == 3; Segment.lon2(segIdx) = value; end;
         if propIdx == 4; Segment.lat1(segIdx) = value; end;
         if propIdx == 5; Segment.lat2(segIdx) = value; end;
         if propIdx == 6; Segment.dip(segIdx) = value; end;
         if propIdx == 7; Segment.dipSig(segIdx) = value; end;
         if propIdx == 8; Segment.dipTog(segIdx) = value; end;
         if propIdx == 9; Segment.lDep(segIdx) = value; end
         if propIdx == 10; Segment.lDepSig(segIdx) = value; end;
         if propIdx == 11; Segment.lDepTog(segIdx) = value; end;
         if propIdx == 12; Segment.ssRate(segIdx) = value; end;
         if propIdx == 13; Segment.ssRateSig(segIdx) = value; end;
         if propIdx == 14; Segment.ssRateTog(segIdx) = value; end;
         if propIdx == 15; Segment.dsRate(segIdx) = value; end;
         if propIdx == 16; Segment.dsRateSig(segIdx) = value; end;
         if propIdx == 17; Segment.dsRateTog(segIdx) = value; end;
         if propIdx == 18; Segment.tsRate(segIdx) = value; end;
         if propIdx == 19; Segment.tsRateSig(segIdx) = value; end;
         if propIdx == 20; Segment.tsRateTog(segIdx) = value; end;
         if propIdx == 21; Segment.res(segIdx) = value; end;
         if propIdx == 22; Segment.resOver(segIdx) = value; end;
         if propIdx == 23; Segment.resOther(segIdx) = value; end;
         if propIdx == 24; Segment.patchFile(segIdx) = value; end;
         if propIdx == 25; Segment.patchTog(segIdx) = value; end;
         if propIdx == 26; Segment.other3(segIdx) = value; end;
         if propIdx == 27; Segment.patchSlipFile(segIdx) = value; end;
         if propIdx == 28; Segment.patchSlipTog(segIdx) = value; end;
         if propIdx == 29; Segment.other6(segIdx) = value; end;
         if propIdx == 30; Segment.other7(segIdx) = value; end;
         if propIdx == 31; Segment.other8(segIdx) = value; end;
         if propIdx == 32; Segment.other9(segIdx) = value; end;
         if propIdx == 33; Segment.other10(segIdx) = value; end;
         if propIdx == 34; Segment.other11(segIdx) = value; end;
         if propIdx == 35; Segment.other12(segIdx) = value; end;
         setappdata(gcf, 'Segment', Segment);
         SegmentManagerFunctions('RedrawSegments');
      end

   
   case 'Seg.modSegList'
      ha                           = findobj(gcf, 'Tag', 'Seg.modSegList');
      value                        = get(ha, 'Value');
      idxSeg                       = value - 2;
      Segment                      = getappdata(gcf, 'Segment');
      nSegments                    = length(Segment.lon1);
      if idxSeg < 1
         set(findobj(gcf, '-regexp', 'Tag', 'Segment.\d'), 'color', 'k');
      else
			set(findobj(gcf, '-regexp', 'Tag', 'Segment.\d'), 'color', 'k');
         set(findobj(gcf, 'Tag', strcat('Segment.', num2str(idxSeg))), 'Color', 'r');
         setappdata(gcf, 'segIdx', idxSeg);
      end

   case 'Seg.modPropList'
      disp('show segment properties')
      Segment                      = getappdata(gcf, 'Segment');
      ha                           = findobj(gcf, 'Tag', 'Seg.modSegList');
      value                        = get(ha, 'Value');
      segIdx                       = value - 2;
      propIdx                      = get(findobj(gcf, 'Tag', 'Seg.modPropList'), 'Value');   
      
      if (value > 1)
	        % Show selected property in edit box
			if propIdx == 2; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.lon1(segIdx))]); end;
			if propIdx == 3; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.lat1(segIdx))]); end;
			if propIdx == 4; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.lon2(segIdx))]); end;
			if propIdx == 5; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.lat2(segIdx))]); end;
			if propIdx == 6; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.dip(segIdx))]); end;
			if propIdx == 7; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.dipSig(segIdx))]); end;
			if propIdx == 8; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.dipTog(segIdx))]); end;
			if propIdx == 9; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.lDep(segIdx))]); end;
			if propIdx == 10; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.lDepSig(segIdx))]); end;
			if propIdx == 11; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.lDepTog(segIdx))]); end;
			if propIdx == 12; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.ssRate(segIdx))]); end;
			if propIdx == 13; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.ssRateSig(segIdx))]); end;
			if propIdx == 14; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.ssRateTog(segIdx))]); end;
			if propIdx == 15; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.dsRate(segIdx))]); end;
			if propIdx == 16; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.dsRateSig(segIdx))]); end;
			if propIdx == 17; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.dsRateTog(segIdx))]); end;
			if propIdx == 18; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.tsRate(segIdx))]); end;
			if propIdx == 19; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.tsRateSig(segIdx))]); end;
			if propIdx == 20; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.tsRateTog(segIdx))]); end;
			if propIdx == 21; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.res(segIdx))]); end;
			if propIdx == 22; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.resOver(segIdx))]); end;
			if propIdx == 23; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.resOther(segIdx))]); end;
			if propIdx == 24; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.patchFile(segIdx))]); end;
			if propIdx == 25; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.patchTog(segIdx))]); end;
			if propIdx == 26; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.other3(segIdx))]); end;
			if propIdx == 27; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.patchSlipFile(segIdx))]); end;
			if propIdx == 28; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.patchSlipTog(segIdx))]); end;
			if propIdx == 29; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.other6(segIdx))]); end;
			if propIdx == 30; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.other7(segIdx))]); end;
			if propIdx == 31; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.other8(segIdx))]); end;
			if propIdx == 32; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.other9(segIdx))]); end;
			if propIdx == 33; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.other10(segIdx))]); end;
			if propIdx == 34; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.other11(segIdx))]); end;
			if propIdx == 35; set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ['   ' num2str(Segment.other12(segIdx))]); end;
      else
			set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'string', ' ');
	  end	
      
   case 'Seg.modGSelect'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         set(findobj(gcf, '-regexp', 'Tag', 'Segment.\d'), 'color', 'k');
         
         % Select a segment graphically
         [minIdx]                  = GetSegmentSingle(Segment);

         % Set the listbox and save the segment index segment
         set(findobj(gcf, 'Tag', 'Seg.modSegList'), 'Value', minIdx + 2);
         setappdata(gcf, 'segIdx', minIdx);
         set(findobj(gcf, 'Tag', 'Seg.modPropList'), 'Value', 1);
         set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'String', ' ');
      end
      

   case 'Seg.modGSelectBox' 
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      % Calculate segment midpoints here, assuming they'll be needed in the future
      Segment.midlon               = (Segment.lon1 + Segment.lon2)/2;
      Segment.midlat               = (Segment.lat1 + Segment.lat2)/2;
      set(findobj(gcf, '-regexp', 'Tag', 'Segment.\d'), 'color', 'k');
      
      % Find the segments that are inside a clickable box
      fprintf(GLOBAL.filestream, 'Starting Box Select\n');
      segRange                     = GetRangeRbbox(getappdata(gcf, 'Range'));
      segPolyX                     = [min(segRange.lon) max(segRange.lon) max(segRange.lon) min(segRange.lon)];
      segPolyY                     = [min(segRange.lat) min(segRange.lat) max(segRange.lat) max(segRange.lat)];
      segIdx                       = find(inpolygon(Segment.midlon, Segment.midlat, segPolyX, segPolyY) == 1);
      for i = 1:numel(segIdx)
         fprintf(GLOBAL.filestream, '%s\n', Segment.name(segIdx(i), :));
         set(findobj('Tag', strcat('Segment.', num2str(segIdx(i)))), 'Color', 'r');
      end
      setappdata(gcf, 'segIdx', segIdx);
      
  case 'Seg.modGSelectLasso'
  	  fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');  		
      Range 					   = getappdata(gcf, 'Range');
      % Calculate segment midpoints here, assuming they'll be needed in the future
      Segment.midlon               = (Segment.lon1 + Segment.lon2)/2;
      Segment.midlat               = (Segment.lat1 + Segment.lat2)/2;
	   set(findobj(gcf, '-regexp', 'Tag', 'Segment.\d'), 'color', 'k');
	   mp = plot(Segment.midlon, Segment.midlat, '.', 'visible', 'off');
      %ignore = setdiff(get(gca, 'children'), mp);
      segIdx							  = myselectdata('sel', 'lasso', 'ignore', mp);
      for i = 1:numel(segIdx)
         fprintf(GLOBAL.filestream, '%s\n', Segment.name(segIdx(i), :));
         set(findobj(gcf, 'Tag', strcat('Segment.', num2str(segIdx(i)))), 'Color', 'r');
      end
      delete(mp); clear mp
      SetAxes(Range);
      setappdata(gcf, 'segIdx', segIdx);
      
  case 'Seg.modClear'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         nSegments                    = length(Segment.lon1);
         set(findobj(gcf, '-regexp', 'Tag', 'Segment.\d'), 'color', 'k');
         set(findobj(gcf, 'Tag', 'Seg.modSegList'), 'Value', 1);
         setappdata(gcf, 'segIdx', []);
         set(findobj(gcf, 'Tag', 'Seg.modPropList'), 'Value', 1);
         set(findobj(gcf, 'Tag', 'Seg.modPropEdit'), 'String', ' ');
      end
   
      
      
   case 'Seg.modNewPush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         Segment                   = NewSegment(Segment);
         setappdata(gcf, 'Segment', Segment);
         SegmentManagerFunctions('RedrawSegments');
      end

      
   case 'Seg.modDeletePush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         Segment                   = DeleteSegmentSingle(Segment);
         setappdata(gcf, 'Segment', Segment);
         SegmentManagerFunctions('RedrawSegments');
      end
      
      
   case 'Seg.modDeletePushBox'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      % Calculate segment midpoints here, assuming they'll be needed in the future
      Segment.midlon               = (Segment.lon1 + Segment.lon2)/2;
      Segment.midlat               = (Segment.lat1 + Segment.lat2)/2;
		set(findobj(gcf, '-regexp', 'Tag', 'Segment.\d'), 'color', 'k');
		
      if ~isempty(Segment)         
         fprintf(GLOBAL.filestream, 'Starting Box Select\n');
         segRange                     = GetRangeRbbox(getappdata(gcf, 'Range'));
         segPolyX                     = [min(segRange.lon) max(segRange.lon) max(segRange.lon) min(segRange.lon)];
         segPolyY                     = [min(segRange.lat) min(segRange.lat) max(segRange.lat) max(segRange.lat)];
         segIdx                       = find(inpolygon(Segment.midlon, Segment.midlat, segPolyX, segPolyY) == 1);
         for i = 1:numel(segIdx)
            fprintf(GLOBAL.filestream, 'Deleting %s\n', Segment.name(segIdx(i), :));
            set(findobj('Tag', strcat('Segment.', num2str(segIdx(i)))), 'Color', 'r');
         end
         Segment                      = DeleteSegment(Segment, segIdx);
         setappdata(gcf, 'Segment', Segment);
         SegmentManagerFunctions('RedrawSegments');
      end
      
      
   case 'Seg.modGDeleteLasso'
  	  fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');  		
      Range 					   = getappdata(gcf, 'Range');
      % Calculate segment midpoints here, assuming they'll be needed in the future
      Segment.midlon               = (Segment.lon1 + Segment.lon2)/2;
      Segment.midlat               = (Segment.lat1 + Segment.lat2)/2;
	   set(findobj(gcf, '-regexp', 'Tag', 'Segment.\d'), 'color', 'k');
	   mp = plot(Segment.midlon, Segment.midlat, '.', 'visible', 'off');
      %ignore = setdiff(get(gca, 'children'), mp);
      segIdx							  = myselectdata('sel', 'lasso', 'ignore', mp);
      for i = 1:numel(segIdx)
         fprintf(GLOBAL.filestream, '%s\n', Segment.name(segIdx(i), :));
         set(findobj(gcf, 'Tag', strcat('Segment.', num2str(segIdx(i)))), 'Color', 'r');
      end
      delete(mp); clear mp
      Segment                      = DeleteSegment(Segment, segIdx);
		setappdata(gcf, 'Segment', Segment);
      SetAxes(Range);
      setappdata(gcf, 'segIdx', segIdx);
      SegmentManagerFunctions('RedrawSegments');

      
   case 'Seg.modMovePush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         Segment                   = MoveIntersection(Segment);
         setappdata(gcf, 'Segment', Segment);
      end

      
   case 'Seg.modConnectPush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         Segment                   = ConnectIntersection(Segment);
         setappdata(gcf, 'Segment', Segment);
         SegmentManagerFunctions('RedrawSegments');
      end
      
      
   case 'Seg.modExtendPush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         Segment                   = ExtendIntersection(Segment);
         setappdata(gcf, 'Segment', Segment);
         SegmentManagerFunctions('RedrawSegments');
      end
   
      
   case 'Seg.modSplitPush'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         Segment                   = SplitSegment(Segment);
         setappdata(gcf, 'Segment', Segment);
         SegmentManagerFunctions('RedrawSegments');
      end
      
      
   %%%   Show segment properties   %%%
   % Show segment properties
   case 'Seg.modShowList'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         switch get(findobj(gcf, 'Tag', 'Seg.modShowList'), 'Value')
            case 1
               DeletePropertyLabels;
            case 2
               ShowPropertyLabels(strjust(Segment.name, 'center')) ; BigTitle('Segment names');
            case 3
               ShowPropertyLabels(num2str(Segment.lon1)) ; BigTitle('Longitude 1');
            case 4
               ShowPropertyLabels(num2str(Segment.lat1)) ; BigTitle('Latitude 1');
            case 5
               ShowPropertyLabels(num2str(Segment.lon2)) ; BigTitle('Longitude 2');
            case 6
               ShowPropertyLabels(num2str(Segment.lat2)) ; BigTitle('Latitude 2');
            case 7
               ShowPropertyLabels(num2str(Segment.dip)) ; BigTitle('Dip');
            case 8
               ShowPropertyLabels(num2str(Segment.dipSig)) ; BigTitle('Dip sigma');
            case 9
               ShowPropertyLabels(num2str(Segment.dipTog)) ; BigTitle('Dip flag');
            case 10
               ShowPropertyLabels(num2str(Segment.lDep)) ; BigTitle('Locking depth');
            case 11
               ShowPropertyLabels(num2str(Segment.lDepSig)) ; BigTitle('Locking depth sigma');
            case 12
               ShowPropertyLabels(num2str(Segment.lDepTog)) ; BigTitle('Locking depth flag');
            case 13
               labels = [num2str(Segment.ssRate) repmat('+/-', numel(Segment.ssRate), 1) strjust(num2str(Segment.ssRateSig), 'left')];
               ShowPropertyLabels(labels) ; BigTitle('Strike slip rate \pm sigma');
            case 14
               ShowPropertyLabels(num2str(Segment.ssRate)) ; BigTitle('Strike slip rate');
            case 15
               ShowPropertyLabels(num2str(Segment.ssRateSig)) ; BigTitle('Strike slip sigma');
            case 16
               ShowPropertyLabels(num2str(Segment.ssRateTog)) ; BigTitle('Strike slip flag');
            case 17
               labels = [num2str(Segment.dsRate) repmat('+/-', numel(Segment.dsRate), 1) strjust(num2str(Segment.dsRateSig), 'left')];
               ShowPropertyLabels(labels) ; BigTitle('Dip slip rate \pm sigma');
            case 18
               ShowPropertyLabels(num2str(Segment.dsRate)) ; BigTitle('Dip slip rate');
            case 19
               ShowPropertyLabels(num2str(Segment.dsRateSig)) ; BigTitle('Dip slip sigma');
            case 20
               ShowPropertyLabels(num2str(Segment.dsRateTog)) ; BigTitle('Dip slip flag');
            case 21
               labels = [num2str(Segment.tsRate) repmat('+/-', numel(Segment.tsRate), 1) strjust(num2str(Segment.tsRateSig), 'left')];
               ShowPropertyLabels(labels) ; BigTitle('Tensile slip rate \pm sigma');
            case 22
               ShowPropertyLabels(num2str(Segment.tsRate)) ; BigTitle('Tensile slip rate');
            case 23
               ShowPropertyLabels(num2str(Segment.tsRateSig)) ; BigTitle('Tensile slip sigma');
            case 24
               ShowPropertyLabels(num2str(Segment.tsRateTog)) ; BigTitle('Tensile slip flag');
            case 25
               ShowPropertyLabels(num2str(Segment.res)) ; BigTitle('Resolution');
            case 26
               ShowPropertyLabels(num2str(Segment.resOv)) ; BigTitle('Resolution override');
            case 27
               ShowPropertyLabels(num2str(Segment.resOther)) ; BigTitle('Resolution other');
            case 28
               ShowPropertyLabels(num2str(Segment.patchFile)) ; BigTitle('Patch name');
            case 29
               ShowPropertyLabels(num2str(Segment.patchTog)) ; BigTitle('Patch toggle');
            case 30
               ShowPropertyLabels(num2str(Segment.other3)) ; BigTitle('Patch Flag 2');
            case 31
               ShowPropertyLabels(num2str(Segment.patchSlipFile)) ; BigTitle('Other 1');
            case 32
               ShowPropertyLabels(num2str(Segment.patchSlipTog)) ; BigTitle('Other 2');
            case 33
               ShowPropertyLabels(num2str(Segment.other6)) ; BigTitle('Other 3');
            case 34
               ShowPropertyLabels(num2str(Segment.other7)) ; BigTitle('Other 4');
            case 35
               ShowPropertyLabels(num2str(Segment.other8)) ; BigTitle('Other 5');
            case 36
               ShowPropertyLabels(num2str(Segment.other9)) ; BigTitle('Other 6');
            case 37
               ShowPropertyLabels(num2str(Segment.other10)) ; BigTitle('Other 7');
            case 38
               ShowPropertyLabels(num2str(Segment.other11)) ; BigTitle('Other 8');
            case 39
               ShowPropertyLabels(num2str(Segment.other12)) ; BigTitle('Other 9');
               
         end
         
      end      
      
      
   % Start File I/O commands (block-file)
   % Load block file
   case 'Seg.loadPushBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      
      % Get the name of the block file
      ha                           = findobj(gcf, 'Tag', 'Seg.loadEditBlock');
      filename                     = get(ha, 'string');
      if exist(filename, 'file')
         filenameFull              = strcat(pwd, '\', filename);
      else
         [filename, pathname]      = uigetfile({'*.block; *.block.xml'}, 'Load block file');
         if filename == 0
            return;
            set(ha, 'string', '');
         else
            set(ha, 'string', filename);
            filenameFull           = strcat(pathname, filename);
         end
      end
      
      % Read in the block file
      Block                        = ReadBlocksStruct(filenameFull);
      setappdata(gcf, 'Block', Block);
      
      % Plot blocks file
	  SegmentManagerFunctions('RedrawBlocks')
	  set(findobj(gcf, 'Tag', 'Seg.dispCheckBlock'), 'enable', 'on', 'value', 1);
	
	% Block display toggle
   case 'Seg.dispCheckBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, '-regexp', 'Tag', '^Block.');
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckBlock');
      if get(hb, 'Value') == 0;
      	set(ha, 'Visible', 'off');
      else
         set(ha, 'Visible', 'on');
      end

      
   % Clear block file
   case 'Seg.clearPushBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                              = findobj(gcf, 'Tag', 'Seg.loadEditBlock');
      set(ha, 'string', '');
      setappdata(gcf, 'Block', []);
      delete(findobj(gcf, '-regexp', 'tag', '^Block.'));
      set(findobj(gcf, 'Tag', 'Seg.dispCheckBlock'), 'enable', 'off', 'value', 0);
      
   % Save block file
   case 'Seg.savePushBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Block                       = getappdata(gcf, 'Block');
      if size(Block, 1) == 0
         return;
      else
         [filename, pathname]     = uiputfile({'*.block; *.block.xml'}, 'Save block file');
         if filename == 0
            return;
         else
            filenameFull          = strcat(pathname, filename);
            WriteBlocksStruct(filenameFull, Block);
         end
         set(findobj(gcf, 'tag', 'Seg.loadEditBlock'), 'string', ['  ' filename]);
      end
      
      
   case 'Seg.modSegListBlock'
      ha                           = findobj(gcf, 'Tag', 'Seg.modSegListBlock');
      value                        = get(ha, 'Value');
      blockIdx                     = value - 2;
      Block                        = getappdata(gcf, 'Block');
      nBlocks                      = numel(Block.interiorLon);
      if blockIdx < 1
         set(findobj(gcf, '-regexp', 'Tag', '^Block.'), 'Color', 'g');
      else
         set(findobj(gcf, '-regexp', 'Tag', '^Block.'), 'Color', 'g');
         set(findobj(gcf, 'Tag', strcat('Block.', num2str(blockIdx))), 'Color', 'r');
         setappdata(gcf, 'blockIdx', blockIdx);
      end

      
   case 'Seg.modDeleteBlockBox'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Block                        = getappdata(gcf, 'Block');
		set(findobj(gcf, '-regexp', 'Tag', '^Block.'), 'color', 'g');
		
      if ~isempty(Block)         
         fprintf(GLOBAL.filestream, 'Starting Box Select\n');
         segRange                     = GetRangeRbbox(getappdata(gcf, 'Range'));
         segPolyX                     = [min(segRange.lon) max(segRange.lon) max(segRange.lon) min(segRange.lon)];
         segPolyY                     = [min(segRange.lat) min(segRange.lat) max(segRange.lat) max(segRange.lat)];
         blockIdx                     = find(inpolygon(Block.interiorLon, Block.interiorLat, segPolyX, segPolyY) == 1);
         for i = 1:numel(blockIdx)
            fprintf(GLOBAL.filestream, 'Deleting %s\n', Block.name(blockIdx(i), :));
            set(findobj('Tag', strcat('Block.', num2str(blockIdx(i)))), 'Color', 'r');
         end
			Block                     	  = DeleteBlock(Block, blockIdx);         
         setappdata(gcf, 'Block', Block);
         SegmentManagerFunctions('RedrawBlocks');
      end     

      
	case 'Seg.modGSelectBlock'
       fprintf(GLOBAL.filestream, '%s\n', option);
       Block                        = getappdata(gcf, 'Block');
       set(findobj(gcf, '-regexp', 'Tag', '^Block.'), 'color', 'g');
       blockIdx							  = SelectBlockSingle(Block);
       setappdata(gcf, 'blockIdx', blockIdx);
		
       % Set the listbox and save the segment index segment
       set(findobj(gcf, 'Tag', 'Seg.modPropListBlock'), 'Value', 1);
       set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'String', ' ');
	
        
	case 'Seg.modGSelectBlockBox'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Block                        = getappdata(gcf, 'Block');
      set(findobj(gcf, '-regexp', 'Tag', '^Block.'), 'color', 'g');
		
      if ~isempty(Block)         
         fprintf(GLOBAL.filestream, 'Starting Box Select\n');
         segRange                  = GetRangeRbbox(getappdata(gcf, 'Range'));
         segPolyX                  = [min(segRange.lon) max(segRange.lon) max(segRange.lon) min(segRange.lon)];
         segPolyY                  = [min(segRange.lat) min(segRange.lat) max(segRange.lat) max(segRange.lat)];
         blockIdx                  = find(inpolygon(Block.interiorLon, Block.interiorLat, segPolyX, segPolyY) == 1);
         for i = 1:numel(blockIdx)
            fprintf(GLOBAL.filestream, 'Selected %s\n', Block.name(blockIdx(i), :));
            set(findobj('Tag', strcat('Block.', num2str(blockIdx(i)))), 'Color', 'r');
         end
         setappdata(gcf, 'blockIdx', blockIdx);
      end
  		set(findobj(gcf, 'Tag', 'Seg.modPropListBlock'), 'Value', 1);
		set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'String', ' ');
      
        
	case 'Seg.modClearBlock'     
		fprintf(GLOBAL.filestream, '%s\n', option);
		SegmentManagerFunctions('RedrawBlocks');
		setappdata(gcf, 'blockIdx', []);
  
        
   case 'Seg.modAddBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Block                        = getappdata(gcf, 'Block');
      if ~isempty(Block)
         Block                     = NewBlock(Block);
         setappdata(gcf, 'Block', Block);
      end
		SegmentManagerFunctions('RedrawBlocks');
      
        
   case 'Seg.modDeleteBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Block                        = getappdata(gcf, 'Block');
      if ~isempty(Block)
         Block                     = DeleteBlockSingle(Block);
         setappdata(gcf, 'Block', Block);
      end
      SegmentManagerFunctions('RedrawBlocks');    

      
   case 'Seg.modMoveBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Block                        = getappdata(gcf, 'Block');
      if ~isempty(Block)
         Block                     = MoveBlock(Block);
         setappdata(gcf, 'Block', Block);
      end
      
   % Start Modify Commands
   case 'Seg.modSegPushBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Block                      = getappdata(gcf, 'Block');
      if ~isempty(Block)
         valueStr                  = get(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'String');
         try % Convert string to numerical value
            value                  = str2double(valueStr);            
         catch
            fprintf(GLOBAL.filestream, 'Could not convert to a numerical value\n');
         end

         % Figure out which field has and which blocks have been selected...then update
         blockIdx                  = getappdata(gcf, 'blockIdx');
         propIdx                   = get(findobj(gcf, 'Tag', 'Seg.modPropListBlock'), 'Value');
         if propIdx <= 2
            fprintf(GLOBAL.filestream, 'Select a property with a numerical value\n');
         end

         % Set and save properties
         if propIdx == 3; Block.eulerLon(blockIdx) = value; end;
         if propIdx == 4; Block.eulerLat(blockIdx) = value; end;
         if propIdx == 5; Block.rotationRate(blockIdx) = value; end;
         if propIdx == 6; Block.eulerLonSig(blockIdx) = value; end;
         if propIdx == 7; Block.eulerLatSig(blockIdx) = value; end;
         if propIdx == 8; Block.rotationRateSig(blockIdx) = value; end;
         if propIdx == 9; Block.rotationInfo(blockIdx) = value; end;
         if propIdx == 10; Block.interiorLon(blockIdx) = value; end;
         if propIdx == 11; Block.interiorLat(blockIdx) = value; end
         if propIdx == 12; Block.aprioriTog(blockIdx) = value; end;
         if propIdx == 13; Block.other1(blockIdx) = value; end;
         if propIdx == 14; Block.other2(blockIdx) = value; end;
         if propIdx == 15; Block.other3(blockIdx) = value; end;
         if propIdx == 16; Block.other4(blockIdx) = value; end;
         if propIdx == 17; Block.other5(blockIdx) = value; end;
         if propIdx == 18; Block.other6(blockIdx) = value; end;
         setappdata(gcf, 'Block', Block);
         SegmentManagerFunctions('RedrawBlocks');
      end
   
   case 'Seg.modPropListBlock'
      disp('show block properties')
      Block                     	  = getappdata(gcf, 'Block');
      ha                           = findobj(gcf, 'Tag', 'Seg.modSegListBlock');
      value                        = get(ha, 'Value');
      blockIdx                     = value - 2;
      propIdx                      = get(findobj(gcf, 'Tag', 'Seg.modPropListBlock'), 'Value');   
      
      if value > 1
			% Show selected property in edit box
			if propIdx == 2; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' Block.name(blockIdx, :)]); end;
			if propIdx == 3; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.eulerLon(blockIdx))]);end;
			if propIdx == 4; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.eulerLat(blockIdx))]);end;
			if propIdx == 5; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.rotationRate(blockIdx))]);end;
			if propIdx == 6; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.eulerLonSig(blockIdx))]);end;
			if propIdx == 7; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.eulerLatSig(blockIdx))]);end;
			if propIdx == 8; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.rotationRateSig(blockIdx))]);end;
			if propIdx == 9; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.rotationInfo(blockIdx))]);end;
			if propIdx == 10; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.interiorLon(blockIdx))]);end;
			if propIdx == 11; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.interiorLat(blockIdx))]);end;
			if propIdx == 12; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.aprioriTog(blockIdx))]);end;
			if propIdx == 13; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.other1(blockIdx))]);end;
			if propIdx == 14; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.other2(blockIdx))]);end;
			if propIdx == 15; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.other3(blockIdx))]);end;
			if propIdx == 16; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.other4(blockIdx))]);end;
			if propIdx == 17; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.other5(blockIdx))]);end;
			if propIdx == 18; set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ['   ' num2str(Block.other6(blockIdx))]);end;
		else
      	set(findobj(gcf, 'Tag', 'Seg.modPropEditBlock'), 'string', ' ');
      end	

      
   % Show block properties
   case 'Seg.modShowListBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Block                       = getappdata(gcf, 'Block');
      if ~isempty(Block)
         switch get(findobj('Tag', 'Seg.modShowListBlock'), 'Value');
            case 1
               DeletePropertyLabelsBlock;
            case 2
               ShowPropertyLabelsBlock(strjust(Block.name, 'center')) ; BigTitle('Block names');
            case 3
               ShowPropertyLabelsBlock(num2str(Block.eulerLon)) ; BigTitle('Euler longitude');
            case 4
               ShowPropertyLabelsBlock(num2str(Block.eulerLat)) ; BigTitle('Euler latitude');
            case 5
               ShowPropertyLabelsBlock(num2str(Block.rotationRate)) ; BigTitle('Rotation rate');
            case 6
               ShowPropertyLabelsBlock(num2str(Block.eulerLonSig)) ; BigTitle('Euler longitude sigma');
            case 7
               ShowPropertyLabelsBlock(num2str(Block.eulerLatSig)) ; BigTitle('Euler latitude sigma');
            case 8
               ShowPropertyLabelsBlock(num2str(Block.rotationRateSig)) ; BigTitle('Rotation rate sigma');
            case 9
               ShowPropertyLabelsBlock(num2str(Block.interiorLon)) ; BigTitle('Interior longitude');
            case 10
               ShowPropertyLabelsBlock(num2str(Block.interiorLat)) ; BigTitle('Interior latitude');
            case 11
               ShowPropertyLabelsBlock(num2str(Block.other1)) ; BigTitle('Other 1');
            case 12
               ShowPropertyLabelsBlock(num2str(Block.other2)) ; BigTitle('Other 2');
            case 13
               ShowPropertyLabelsBlock(num2str(Block.other3)) ; BigTitle('Other 3');
            case 14
               ShowPropertyLabelsBlock(num2str(Block.other4)) ; BigTitle('Other 4');
            case 15
               ShowPropertyLabelsBlock(num2str(Block.other5)) ; BigTitle('Other 5');
            case 16
               ShowPropertyLabelsBlock(num2str(Block.other6)) ; BigTitle('Other 6');
         end
      end
      
      
   % Load mesh file
   case 'Seg.loadPushMesh'
      % Delete all the childrean of the current axes
      fprintf(GLOBAL.filestream, '%s\n', option);
      
      % Get the name of the mesh file
      ha                           = findobj(gcf, 'Tag', 'Seg.loadEditMesh');
      filename                     = get(ha, 'string');
      if exist(filename, 'file')
         filenameFull              = strcat(pwd, '\', filename);
      else
         [filename, pathname]      = uigetfile({'*.msh'}, 'Load mesh file');
         if filename == 0
            return;
            set(ha, 'string', '');
         else
            set(ha, 'string', filename);
            filenameFull           = strcat(pathname, filename);
         end
      end
      
      % Read in the mesh file and plot
      P = ReadPatches(filenameFull);
      h = patch('Vertices', P.c, 'faces', P.v, 'facecolor', 'g', 'edgecolor', 'black', 'tag', 'Patch');
      % Enable the display check box
      set(findobj(gcf, 'Tag', 'Seg.dispCheckMesh'), 'enable', 'on', 'value', 1);
	
	% Mesh display toggle
   case 'Seg.dispCheckMesh'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'Patch');
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckMesh');
      if get(hb, 'Value') == 0;
      	set(ha, 'Visible', 'off');
      else
         set(ha, 'Visible', 'on');
      end

      
   % Clear mesh file
   case 'Seg.clearPushMesh'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                              = findobj(gcf, 'Tag', 'Seg.loadEditMesh');
      set(ha, 'string', '');
      delete(findobj(gcf, 'tag', 'Patch'));
      set(findobj(gcf, 'Tag', 'Seg.dispCheckMesh'), 'enable', 'off', 'value', 0);
      
      
% File integrity functions 

   % Check segment boundaries
   case 'Seg.modCheckSegsBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      if ~isempty(Segment)
         try
            [~, ~]  = BlockLabel(Segment.lon1, Segment.lat1, Segment.lon2, Segment.lat2, Segment.name);
            htemp                  = msgbox('All segments form closed blocks');
         catch
            htemp                  = msgbox('Segments do not form closed blocks!');
         end
      end
   
      
   % Check interior points
   case 'Seg.modCheckIpsBlock'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                      = getappdata(gcf, 'Segment');
      Block                        = getappdata(gcf, 'Block');

      if ~isempty(Segment) && ~isempty(Block)
         try
            [Segment.eastLabel Segment.westLabel]            = BlockLabel(Segment.lon1, Segment.lat1, Segment.lon2, Segment.lat2, Segment.name);
            [Block.associateLabel Block.exteriorBlockLabel]  = BlockAssociate(Segment.lon1, Segment.lat1, Segment.lon2, Segment.lat2, Segment.westLabel, Segment.eastLabel, Block.interiorLon, Block.interiorLat, Block.name);
            htemp                  = msgbox(sprintf('Interior points uniquely identify blocks %n blocks', numel(Block.associateLabel)));
         catch
            htemp                  = msgbox('Interior points do not uniquely identify blocks!');
         end
      end
      
      
   % Check for problematic segments
   case 'Seg.modSegmentChecker'
      fprintf(GLOBAL.filestream, '%s\n', option);
      SegmentManagerFunctions('Seg.modClearSegmentChecks');
		Segment                      = getappdata(gcf, 'Segment');
		h 									  = findobj(gcf, '-regexp', 'tag', 'Segment.\d');
      if ~isempty(Segment)
			SegmentCheckerForGui(Segment, h)
      end
    
    
   % Clear segment checks
   case 'Seg.modClearSegmentChecks'
      fprintf(GLOBAL.filestream, '%s\n', option);
      SegmentManagerFunctions('RedrawSegments');
		legend('deletelegend')
		hg = findobj(gcf, 'tag', 'hang');
		if ~isempty(hg)
			delete(hg)
        end

        
% Start Display Commands
   % Change topo settings
   case 'Seg.dispTopo'
      fprintf(GLOBAL.filestream, '%s\n', option);
      value                        = get(findobj(gcf, 'Tag', 'Seg.dispTopo'), 'Value');
      if value == 1
         set(findobj(gcf, 'Tag', 'topo'), 'visible', 'off');
         delete(findobj(gcf, 'Tag', 'topo'));
      elseif value == 2
          % plot_google_map('MapType', 'terrain', 'ShowLabels', 0, 'AutoAxis', 0)
         xlim = get(gca, 'xlim');
         ylim = get(gca, 'ylim');
         if xlim(1) > 180
            xlimtemp = xlim - 360;
         else
            xlimtemp = xlim;
         end

         set(gca, 'xlim', xlim);
         set(gca, 'ylim', ylim);
         [lonvec, latvec, map] = my_plot_google_map(xlimtemp, ylim);

         if xlim(1) > 180
            lonvec = lonvec + 360;
         end

         h = image(lonvec, latvec, map);    
         set(gca,'YDir','Normal');
         set(h, 'tag', 'topo');
         uistack(h, 'bottom');
         set(gca, 'xlim', xlim);
         set(gca, 'ylim', ylim);
         set(findobj(gcf, 'Tag', 'topo'), 'visible', 'on');
      else
         return;
      end
   
        
   % Change grid settings  
   case 'Seg.dispGrid'
      fprintf(GLOBAL.filestream, '%s\n', option);
      value                        = get(findobj(gcf, 'Tag', 'Seg.dispGrid'), 'Value');
      if value == 1
         set(gca, 'XGrid', 'off', 'YGrid', 'off');
      elseif value == 2
         set(gca, 'XGrid', 'on', 'YGrid', 'on');
      else
         return;
      end

   % Show, or don't show Dips  
   case 'Seg.dispDips'
      fprintf(GLOBAL.filestream, '%s\n', option);
      value = get(findobj(gcf, 'Tag', 'Seg.dispDips'), 'Value');
      Segment = getappdata(gcf, 'Segment');
      if value == 1
          % Remove surface projection of dipping structures
          h = (findobj(gcf, 'Tag', 'Dips'));
          if ~isempty(h)
              delete(h);
          end
      elseif value == 2
          % Plot surface projection of dipping structures
          PlotDips(Segment.lon1, Segment.lat1, Segment.lon2, Segment.lat2, Segment.dip, Segment.lDep, Segment.bDep)   
      else
         return;
      end      
      
      
   % Change longitude labeling
   case 'Seg.dispMeridian'
      fprintf(GLOBAL.filestream, '%s\n', option);
      value                        = get(findobj(gcf, 'Tag', 'Seg.dispMeridian'), 'Value');
      if value == 1
         set(gca, 'XTickLabel', deblank(strjust(num2str(zero22pi(transpose(get(gca, 'XTick')))), 'center')));
      elseif value == 2
         set(gca, 'XTickLabel', deblank(strjust(num2str(npi2pi(transpose(get(gca, 'XTick')))), 'center')));
      else
         return;
      end

  
   % Load a line file
   case 'Seg.dispPushLine'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'Seg.dispEditLine');
      filename                     = get(ha, 'String');
      if exist(filename, 'file')
         filenameFull              = strcat(pwd, '\', filename);
      else
         [filename, pathname]      = uigetfile({'*'}, 'Load line file');
         if filename == 0
            return;
            set(ha, 'string', '');
         else
            set(ha, 'string', filename);
            filenameFull           = strcat(pathname, filename);
         end
      end
      PlotLine(filenameFull);
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckLine');
      set(hb, 'Value', 1);
   
        
   % Toggle the line file visibility
   case 'Seg.dispCheckLine'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'lineAll');
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckLine');
      if isempty(ha)
         return;
      else
         if get(hb, 'Value') == 0
            set(ha, 'Visible', 'off');
         elseif get(hb, 'Value') == 1
            set(ha, 'Visible', 'on');
         end
      end
     
        
   % Load a xy file
   case 'Seg.dispPushXy'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'Seg.dispEditXy');
      filename                     = get(ha, 'String');
      if exist(filename, 'file')
         filenameFull              = strcat(pwd, '\', filename);
      else
         [filename, pathname]      = uigetfile({'*'}, 'Load xy file');
         if filename == 0
            return;
            set(ha, 'string', '');
         else
            set(ha, 'string', filename);
            filenameFull           = strcat(pathname, filename);
         end
      end
      PlotXy(filenameFull);
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckXy');
      set(hb, 'Value', 1);
     
        
   % Toggle the xy file visibility
   case 'Seg.dispCheckXy'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'xyAll');
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckXy');
      if isempty(ha)
         return;
      else
         if get(hb, 'Value') == 0
            set(ha, 'Visible', 'off');
         elseif get(hb, 'Value') == 1
            set(ha, 'Visible', 'on');
         end
      end
     
      
   % Load a station file
   case 'Seg.dispPushSta'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'Seg.dispEditSta');
      filename                     = get(ha, 'String');
      if exist(filename, 'file')
         filenameFull              = strcat(pwd, '\', filename);
      else
         [filename, pathname]      = uigetfile({'*.sta.data'}, 'Load station file');
         if filename == 0
            return;
            set(ha, 'string', '');
         else
            set(ha, 'string', filename);
            filenameFull           = strcat(pathname, filename);
         end
      end
      Station = PlotSta(filenameFull);

      PlotStaVec(Station)
      setappdata(gcf, 'Station', Station);
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckSta');
      set(hb, 'Value', 1);
      
        
   % Toggle the station file visibility
   case 'Seg.dispCheckSta'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'staAll');
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckSta');
      if isempty(ha)
         return;
      else
         if get(hb, 'Value') == 0
            set(ha, 'Visible', 'off');
         elseif get(hb, 'Value') == 1
            set(ha, 'Visible', 'on');
         end
      end

      
   % Toggle the vector visibility
   case 'Seg.dispCheckStaVec'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'staAllVec');
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckStaVec');
      if isempty(ha)
      	fprintf(GLOBAL.filestream, 'No vectors\n');
         return;
      else
         if get(hb, 'Value') == 0
            set(ha, 'Visible', 'off');
         elseif get(hb, 'Value') == 1
            set(ha, 'Visible', 'on');
         end
      end      

   % Toggle station name visibility
   case 'Seg.dispCheckStaNames'
      fprintf(GLOBAL.filestream, '%s\n', option);
      ha                           = findobj(gcf, 'Tag', 'staNames');
      hb                           = findobj(gcf, 'Tag', 'Seg.dispCheckStaNames');
      if isempty(ha) % plot if it doesn't exit
      	fprintf(GLOBAL.filestream, 'No station names\n');
        % text command goes here
        Station = getappdata(gcf, 'Station');
        text(Station.lon, Station.lat, Station.name, 'Interpreter', ...
             'none', 'FontName', 'Lucida', 'Tag', 'staNames');
        return;
      else
         if get(hb, 'Value') == 0
            set(ha, 'Visible', 'off');
         elseif get(hb, 'Value') == 1
            set(ha, 'Visible', 'on');
         end
      end      
            
      
   case 'Seg.velPushUp'
      ScaleAllVectors(1.1);


   case 'Seg.velPushDown'
      ScaleAllVectors(0.9);


   % Start Navigation Commands
   case 'Seg.navZoomRange'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = GetRangeRbbox(getappdata(gcf, 'Range'));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      Range                        = getappdata(gcf, 'Range');
      Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      
            
   case 'Seg.navZoomIn'
      fprintf(GLOBAL.filestream, '%s\n', option);
      zoomFactor                   = 0.5;
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = (max(Range.lon) - min(Range.lon)) / 2;
      deltaLat                     = (max(Range.lat) - min(Range.lat)) / 2;
      centerLon                    = mean(Range.lon);
      centerLat                    = mean(Range.lat);
      Range.lon                    = [centerLon - zoomFactor * deltaLon, centerLon + zoomFactor * deltaLon];
      Range.lat                    = [centerLat - zoomFactor * deltaLat, centerLat + zoomFactor * deltaLat];
      Range                        = CheckRange(Range);
 		Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);

      
   case 'Seg.navZoomOut'
      fprintf(GLOBAL.filestream, '%s\n', option);
      zoomFactor                   = 2.0;
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = (max(Range.lon) - min(Range.lon)) / 2;
      deltaLat                     = (max(Range.lat) - min(Range.lat)) / 2;
      centerLon                    = mean(Range.lon);
      centerLat                    = mean(Range.lat);
      Range.lon                    = [centerLon - zoomFactor * deltaLon, centerLon + zoomFactor * deltaLon];
      Range.lat                    = [centerLat - zoomFactor * deltaLat, centerLat + zoomFactor * deltaLat];
      Range                        = CheckRange(Range);
      Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul                          = min([size(Range.lonOld, 1)-1 cul+1]);
      st                           = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
            
   case 'Seg.navUpdate'
      fprintf(GLOBAL.filestream, '%s\n', option);
      lonMax                       = str2num(get(findobj(gcf, 'Tag', 'Seg.navEditLonMax'), 'string'));
      lonMin                       = str2num(get(findobj(gcf, 'Tag', 'Seg.navEditLonMin'), 'string'));
      latMax                       = str2num(get(findobj(gcf, 'Tag', 'Seg.navEditLatMax'), 'string'));
      latMin                       = str2num(get(findobj(gcf, 'Tag', 'Seg.navEditLatMin'), 'string'));
      Range                        = getappdata(gcf, 'Range');
      Range.lon                    = Range.lon - translateScale * deltaLon;
      Range.lat                    = Range.lat - translateScale * deltaLat;
      Range                        = CheckRange(Range);
      Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
      
   case 'Seg.navBack'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      RangeLev							  = max([1 cul]);
      Range.lon						  = Range.lonOld(RangeLev, :);
      Range.lat						  = Range.latOld(RangeLev, :);
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      cul 								  = max([1 RangeLev - 1]);

      
   case 'Seg.navSW'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lon                    = Range.lon - translateScale * deltaLon;
      Range.lat                    = Range.lat - translateScale * deltaLat;
      Range                        = CheckRange(Range);
      Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      

   case 'Seg.navS'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lon                    = Range.lon;
      Range.lat                    = Range.lat - translateScale * deltaLat;
      Range                        = CheckRange(Range);
 		Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
      
   case 'Seg.navSE'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lon                    = Range.lon + translateScale * deltaLon;
      Range.lat                    = Range.lat - translateScale * deltaLat;
      Range                        = CheckRange(Range);
 		Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
      
   case 'Seg.navW'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lon                    = Range.lon - translateScale * deltaLon;
      Range.lat                    = Range.lat;
      Range                        = CheckRange(Range);
 		Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
      
   case 'Seg.navC'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lonOld                 = [Range.lonOld(2:end, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(2:end, :) ; Range.lat];
      k                            = waitforbuttonpress;
      point                        = get(gca, 'CurrentPoint');
      point                        = point(1, 1:2);
      Range                        = CheckRange(Range);
 		Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
      
   case 'Seg.navE'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lon                    = Range.lon + translateScale * deltaLon;
      Range.lat                    = Range.lat;
      Range                        = CheckRange(Range);
      Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
      
   case 'Seg.navNW'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lon                    = Range.lon - translateScale * deltaLon;
      Range.lat                    = Range.lat + translateScale * deltaLat;
      Range                        = CheckRange(Range);
      Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
      
   case 'Seg.navN'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lon                    = Range.lon;
      Range.lat                    = Range.lat + translateScale * deltaLat;
      Range                        = CheckRange(Range);
      Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
      
   case 'Seg.navNE'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Range                        = getappdata(gcf, 'Range');
      deltaLon                     = max(Range.lon) - min(Range.lon);
      deltaLat                     = max(Range.lat) - min(Range.lat);
      Range.lon                    = Range.lon + translateScale * deltaLon;
      Range.lat                    = Range.lat + translateScale * deltaLat;
      Range                        = CheckRange(Range);
      Range.lonOld                 = [Range.lonOld(st:cul+1, :) ; Range.lon];
      Range.latOld                 = [Range.latOld(st:cul+1, :) ; Range.lat];
      cul								  = min([size(Range.lonOld, 1)-1 cul+1]);
      st									  = 1 + (cul==(ul-1));
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);
      
        
   % Draw the clean map
   case 'DrawClean'
      fprintf(GLOBAL.filestream, '%s\n', option);
      delete(gca);
      Seg.axHandle                 = axes('parent', gcf, 'units', 'pixels', 'position', [360 80 800 700],     'visible', 'on', 'Tag', 'Seg.axHandle', 'Layer', 'top', 'xlim', [0 360], 'ylim', [-90 90], 'nextplot', 'add');
      WorldHiVectors               = load('WorldHiVectors');
      plot(WorldHiVectors.lon, WorldHiVectors.lat, '-k', 'LineWidth', 0.25, 'visible', 'on', 'tag', 'Seg.coastLow', 'Color', 0.7 * [1 1 1]);
      box on;
      Range.lon                    = [0 360];
      Range.lat                    = [-90 90];
      Range.lonOld                 = repmat(Range.lon, ul, 1);
      Range.latOld                 = repmat(Range.lat, ul, 1);
      setappdata(gcf, 'Range', Range);
      SetAxes(Range);

      
   % Redraw the segments   
   case 'DrawSegments'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Segment                          = getappdata(gcf, 'Segment');
      nSegments                        = length(Segment.lon1);
		segnames									= cellstr([repmat('Segment.', nSegments, 1) strjust(num2str([1:nSegments]'), 'left')]);
		plotsegs 								= line([Segment.lon1'; Segment.lon2'], [Segment.lat1'; Segment.lat2'], 'color', 'k');
		set(plotsegs, {'tag'}, segnames);
    
      
   % Redraw the segments, including update of the name list   
   case 'RedrawSegments'
      fprintf(GLOBAL.filestream, '%s\n', option);
      
      % Delete all old segments
      delete(findobj(gcf, '-regexp', 'tag', '^Segment.\d'));
      
      % Sort new Segment structure and update fault name pulldown
      Segment                          = getappdata(gcf, 'Segment');
      Segment                          = AlphaSortSegment(Segment);
      setappdata(gcf, 'Segment', Segment);
      set(findobj(gcf, 'Tag', 'Seg.modSegList'), 'string', cellstr(strvcat(' ', 'Multiple', Segment.name)));
            
      % Plot segments again
      nSegments                        = length(Segment.lon1);
      segnames									= cellstr([repmat('Segment.', nSegments, 1) strjust(num2str([1:nSegments]'), 'left')]);
      plotsegs 								= line([Segment.lon1'; Segment.lon2'], [Segment.lat1'; Segment.lat2'], 'color', 'k', 'marker', 'none');
      set(plotsegs, {'tag'}, segnames);
	
	% Redraw the blocks, including update of the name list 
	case 'RedrawBlocks'
      fprintf(GLOBAL.filestream, '%s\n', option);
      
      % Delete all old block IPs
      delete(findobj(gcf, '-regexp', 'tag', 'Block.\d'));
      
      % Sort new block structure and update block name pulldown
      Block                          = getappdata(gcf, 'Block');
      Block                          = AlphaSortBlock(Block);
      setappdata(gcf, 'Block', Block);
      set(findobj(gcf, 'Tag', 'Seg.modSegListBlock'), 'string', cellstr(strvcat(' ', 'Multiple', Block.name)));
            
      % Plot blocks again
      nBlocks                          = length(Block.interiorLon);
      blnames									= cellstr([repmat('Block.', nBlocks, 1) strjust(num2str([1:nBlocks]'), 'left')]);
%       plotbls 									= line([Block.interiorLon'; Block.interiorLon'], [Block.interiorLat'; Block.interiorLat'], 'color', 'g', 'marker', '.', 'linestyle', 'none');
      plotbls                      = line([Block.interiorLon'; Block.interiorLon'], [Block.interiorLat'; Block.interiorLat'], 'MarkerFaceColor', 'm', 'MarkerSize', 5, 'marker', 'o', 'linestyle', 'none', 'MarkerEdgeColor', 'k');
      set(plotbls, {'tag'}, blnames);
      
      % Reset the selected block name to blank
      ha                           = findobj(gcf, 'Tag', 'Seg.modSegListBlock');
      set(ha, 'Value', 1);
      
      
   % Print the figure
   case 'Seg.pszPrint'
      fprintf(GLOBAL.filestream, '%s\n', option);
      printdlg(gcf);
   
            
   % Save the figure
   case 'Seg.pszSave'
      fprintf(GLOBAL.filestream, '%s\n', option);
      %%  Get filename
      filename                    = char(inputdlg('Please enter a base filename', 'Base filename', 1));
      if length(filename) > 0
         Zumax(gca);
         SaveCurrentFigure(filename);
         delete(gcf);
      else
         return;
      end
      
      
   % Zumax the figure
   case 'Seg.pszZumax'
      fprintf(GLOBAL.filestream, '%s\n', option);
      Zumax(gca);

end


function Range = GetRangeRbbox(Range)
% GetRangeRbbox
k                            = waitforbuttonpress;
point1                       = get(gca, 'CurrentPoint');
finalRect                    = rbbox;
point2                       = get(gca, 'CurrentPoint');
point1                       = point1(1,1:2);
point2                       = point2(1,1:2);
Range.lon                    = sort([point1(1) point2(1)]);
Range.lat                    = sort([point1(2) point2(2)]);

function SetAxes(Range)
% SetAxes
axis([min(Range.lon) max(Range.lon) min(Range.lat) max(Range.lat)]);
set(findobj(gcf, 'Tag', 'Seg.navEditLonMin'), 'string', sprintf('%7.3f', min(Range.lon)));
set(findobj(gcf, 'Tag', 'Seg.navEditLonMax'), 'string', sprintf('%7.3f', max(Range.lon)));
set(findobj(gcf, 'Tag', 'Seg.navEditLatMin'), 'string', sprintf('%7.3f', min(Range.lat)));
set(findobj(gcf, 'Tag', 'Seg.navEditLatMax'), 'string', sprintf('%7.3f', max(Range.lat)));
yAspect                        = cos(deg2rad(mean(Range.lat)));
daspect([1 yAspect 1]);

if max(Range.lon) == 360
   set(gca, 'XTick', [0 60 120 180 240 300 360]);
   set(gca, 'YTick', [-90 -45 0 45 90]);
else
   set(gca, 'XTickMode', 'auto');
   set(gca, 'YTickMode', 'auto');
end   
SegmentManagerFunctions('Seg.dispMeridian');


function Range = CheckRange(Range)
% CheckRange
Range.lon                    = sort(Range.lon);
Range.lat                    = sort(Range.lat);
Range.lon(Range.lon > 360)   = 360;
Range.lon(Range.lon < 0)     = 0;
Range.lat(Range.lat > 90)    = 90;
Range.lat(Range.lat < -90)   = -90;


function PlotLine(filename)
% *** CHANGED: JPL 7 January 2008 ****
% open the file
fid1 = fopen(filename); frewind(fid1);

in = textscan(fid1, '%s', 'delimiter', '\n', 'whitespace', '');
in = in{1};
in = char(in);
szin = size(in);

fclose(fid1);

% find line separators
in = strjust(in, 'left'); % shift all left
str = in(:, 1)';
pat = '[^-\d]';
blank = regexp(str, pat, 'start');
in(blank, 1:7) = repmat('NaN NaN', length(blank), 1);
in = str2num(in);
plot(zero22pi(in(:, 1)), in(:, 2), '-', 'LineWidth', 0.5, 'Color', 'm', 'Tag', 'lineAll');
% **** 


function PlotXy(filename)
% PlotXy
fileStream      = fopen(filename, 'r');
data            = fgetl(fileStream);
while (isstr(data));
   % Try a conversion to numeric
   vals         = str2num(data);
   if ~isempty(vals)
      plot(zero22pi(vals(1)), vals(2), '.k', 'Tag', 'xyAll');
   end
   
   % Get the next line
   data = fgetl(fileStream);
end
fclose(fileStream);


function Station = PlotSta(filename)
% PlotSta
Station         = GetVelStruct(filename);
on 				 = find(Station.tog);
plot(Station.lon(on), Station.lat(on), 'bo', 'color', [0 0 1], 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'Tag', 'staAll');


function PlotStaVec(Station)
on					 = find(Station.tog);
vecScale = get(findobj(gcf, 'tag', 'Seg.velSlider'), 'value');
quiver(Station.lon(on), Station.lat(on), Station.eastVel(on), Station.northVel(on), 0, 'userdata', vecScale, 'color', [0 0 1], 'tag', 'staAllVec', 'visible', 'off');


% Scale station vectors
function ScaleAllVectors(vecScale)
groups = findobj(gcf, 'type', 'hggroup');
set(groups, 'Udata', vecScale*get(groups, 'UData'));
set(groups, 'Vdata', vecScale*get(groups, 'VData'));


function DeletePropertyLabels
ha              = findobj('Tag', 'propertyLabel');
if isstruct(get(ha))
   delete(ha);
end
title('');


function DeletePropertyLabelsBlock
ha              = findobj('Tag', 'propertyLabelBlock');
if isstruct(get(ha))
   delete(ha);
end
title('');


function ShowPropertyLabels(labels)
DeletePropertyLabels;
Segment          = getappdata(gcf, 'Segment');
lonMid           = (Segment.lon1 + Segment.lon2) / 2;
latMid           = (Segment.lat1 + Segment.lat2) / 2;
text(lonMid, latMid, labels, 'Tag', 'propertyLabel', 'HorizontalAlignment', 'center', 'clipping', 'on', 'FontSize', 8, 'FontName', 'Lucida', 'BackgroundColor', 'none', 'EdgeColor', 'none', 'Interpreter', 'none');


function ShowPropertyLabelsBlock(labels)
DeletePropertyLabelsBlock;
Block            = getappdata(gcf, 'Block');
text(Block.interiorLon, Block.interiorLat, labels, 'Tag', 'propertyLabelBlock', 'HorizontalAlignment', 'center', 'clipping', 'on', 'FontSize', 8, 'FontName', 'Lucida', 'BackgroundColor', 'none', 'EdgeColor', 'none', 'Interpreter', 'none');


function BigTitle(label)
title(label, 'FontSize', 16);
