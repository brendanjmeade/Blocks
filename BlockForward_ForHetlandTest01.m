%Runs BlocksForward for the stations in Mod.Sta...so the results of Mod.Sta
%should be the same as the results of BlocksForward. V_B is the same, V_SD
%is not.

clear all; close all;

%load Hetland stations...we don't actually use these in this test code
[staNames, staLats, staLons] = textread('Izmit_NAF_Station_List_Hetland.csv','%s%f%f', 'delimiter', ',');

folderBlockResults = 'ElasticBlockResults_ForHeltandJune2015/';

inSegFile = sprintf('%s/Mod.Segment',folderBlockResults);
ObsFile = sprintf('%s/Obs.Sta',folderBlockResults);
ModFile = sprintf('%s/Mod.Sta',folderBlockResults);
RotFile = sprintf('%s/Rot.Sta',folderBlockResults);

Obs = ReadStation(ObsFile);
Mod = ReadStation(ModFile);
Rot = ReadStation(RotFile);

Segment = ReadSegmentStruct(inSegFile);

% figfolder = 'HetlandBlocksJune2015Figs';
% mkdir(figfolder);
% 
% [VB, VSD, VT, VS] = BlocksForward_Eileen(staLons, staLats, folderBlockResults, []);
% [V] = BlocksForward_Eileen(staLons, staLats, folderBlockResults, []);

[VB, VSD, VT, VS] = BlocksForward(Mod.lon, Mod.lat, folderBlockResults);
[V] = BlocksForward(Mod.lon, Mod.lat, folderBlockResults);

VBEast = (VB(1:3:end));
VBNorth = (VB(2:3:end));

VSDEast = (VSD(1:3:end));
VSDNorth = (VSD(2:3:end));

VSEast = (VS(1:3:end));
VSNorth = (VS(2:3:end));

VE = (V(1:3:end));
VN = (V(2:3:end));

latlimmin = 37;
latlimmax = 43;
lonlimmin = 25;
lonlimmax = 38;

figure; hold on;
mstruct = defaultm('mercator');
mstruct.origin = [0 90 0];
mstruct = defaultm(mstruct);
% load topo_NAF.xyzObs
lat1 = 37;
lat2 = 43;
lon1 = 25;
lon2 = 38;

% loop over V_B, v_SD, v_Strain (which is zero everywhere here...) and
% plot them all
for ii = 1:4
    if ii == 1
        fieldE = VBEast;
        fieldN = VBNorth;
        titlestr ='$\mathbf{v}_{\mathrm{B}}\,\mathrm{red}=\mathrm{BlocksForward}, \,\mathrm{black}=\mathrm{Blocks}$';
        figName = 'VBHetland.jpg';
        scale(ii) = .04;

    elseif ii == 2
        fieldE = VSDEast;
        fieldN = VSDNorth;
        figName = 'VSDHetland.jpg';
        titlestr ='$\mathbf{v}_{\mathrm{SD}}\,\mathrm{red}=\mathrm{BlocksForward}, \,\mathrm{black}=\mathrm{Blocks}$';
        scale(ii) = .04;
        
    elseif ii == 3
        fieldE = VSEast;
        fieldN = VSNorth; 
        figName = 'Locations.jpg'; % these velocities are just 0, so we just plot station locations in this plot
        titlestr ='Locations, Big Black Circles = Hetland Stations, Small Red Circles = Our Stations';
        scale(ii) = .04;
        
    elseif ii == 4
        fieldE = VE;
        fieldN = VN;
        titlestr='$\mathbf{v}_{\mathrm{I}}, \,$';
        figName = 'VTotalHeltand.jpg';
        scale(ii) = .04;
        
    end
    
    figure('Position', [0 0 1200 1100]); hold on; %eval(mp); eval(mg);
    
    xlim([lon1 lon2]); ylim([lat1 lat2]);
    % lightangle(-45,30)% camlight right
    
    for i = 1:numel(Segment.lon1)
        plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 0.00*[1 1 1], 'linewidth', 6)
    end
    for i = 1:numel(Segment.lon1)
        plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 1.00*[1 1 1], 'linewidth', 3)
    end

    for i = 1:numel(VBEast)
        %                 arrow( llons(i,j),Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
        %   axis(axis)
        if (fieldE(i) ~= 0 && fieldN(i) ~= 0)
            a1 = arrow([Mod.lon(i) Mod.lat(i) 600], [Mod.lon(i)+scale(ii)*fieldE(i), Mod.lat(i)+scale(ii)*fieldN(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'r', 'lineWidth', 3);
        else
            p = plot3(Mod.lon(i), Mod.lat(i), 600, 'ko', 'markersize', 12, 'markerfacecolor', 'k');
            p = plot3(Mod.lon(i), Mod.lat(i), 600, 'ro', 'markersize', 8, 'markerfacecolor', 'r');
            
        end
    end
    
    if ii == 1 % if we're plotting v_B, also plot rotation velocities from Rot.sta
        for i = 1:numel(Mod.lon)
            a1 = arrow([Mod.lon(i) Mod.lat(i) 650], [Mod.lon(i)+scale(ii)*Rot.eastVel(i), Mod.lat(i)+scale(ii)*Rot.northVel(i) 650], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'k', 'edgecolor', 'k', 'lineWidth', 1);
        end
    end
    if ii == 2 % if we're plotting v_SD from BlocksForward, also plot v_SD from Blocks
          fieldSDE = (Rot.eastVel-Mod.eastVel);
          fieldSDN = (Rot.northVel-Mod.northVel);
          for i = 1:numel(Mod.lon)
            a1 = arrow([Mod.lon(i) Mod.lat(i) 650], [Mod.lon(i)+scale(ii)*fieldSDE(i), Mod.lat(i)+scale(ii)*fieldSDN(i) 650], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'k', 'edgecolor', 'k', 'lineWidth', 1);
          end
    end
%     keyboard;
    a = arrow([36 42 600], [36+scale(ii)*20, 42 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 1);
    text(36, 42.2, '20 mm/yr');
    
    set(gca,'xtick', [26 28 30 32 34 36 38], 'ytick', [38 39 40 41 42 43])
    ticks_format('%.0f','%.0f');
    [hx,hy] = format_ticks(gca,'^{\circ}E','^{\circ}N');
    set(gca, 'fontName', 'times');
    title(titlestr, 'interpreter', 'latex');
    
%     export_fig(sprintf('%s/%s', figfolder, figName))
    
end

%plot observed and modeled velocities from blocks
figure('Position', [0 0 1200 1100]); hold on; %eval(mp); eval(mg);
xlim([lon1 lon2]); ylim([lat1 lat2]);

for i = 1:numel(Segment.lon1)
    plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 0.00*[1 1 1], 'linewidth', 6)
end
for i = 1:numel(Segment.lon1)
    plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 1.00*[1 1 1], 'linewidth', 3)
end

for i = 1:numel(Mod.lon)
    %                 arrow( llons(i,j),Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
    %   axis(axis)
    a2 = arrow([Mod.lon(i) Mod.lat(i) 600], [Mod.lon(i)+scale(ii)*Mod.eastVel(i), Mod.lat(i)+scale(ii)*Mod.northVel(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'b', 'edgecolor', 'k', 'lineWidth', 2);
    
    %         set(a, 'facecolor', 'w', 'edgecolor', 'k', 'lineWidth', 1);
    %         m_vec(35, llons(i,j), llats(i,j), VBEast(i,j), VBNorth(i,j), arrowcol, 'shaftwidth', sw, 'headlength', 6.0, 'headangle', 60,'edgecolor', 'k', 'linewidth', 1.5);
end
for i = 1:numel(Mod.lon)
    
    a3 = arrow([Obs.lon(i) Obs.lat(i) 600], [Obs.lon(i)+scale(ii)*Obs.eastVel(i), Obs.lat(i)+scale(ii)*Obs.northVel(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
    
end
a = arrow([36 42 600], [36+scale(ii)*20, 42 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
text(36, 42.2, '20 mm/yr');

set(gca,'xtick', [26 28 30 32 34 36 38], 'ytick', [38 39 40 41 42 43])
ticks_format('%.0f','%.0f');
[hx,hy] = format_ticks(gca,'^{\circ}E','^{\circ}N');
set(gca, 'fontName', 'times');
titlestr = 'Modeled in blue, Observed in red, from Blocks';
title(titlestr, 'interpreter', 'latex');
xlim([lon1 lon2]); ylim([lat1 lat2]);

% export_fig(sprintf('%s/%s', figfolder, 'Obs_Mod_Vels.pdf'))



% plot Blocksforward v_I results with Blocks results on same plot
figure('Position', [0 0 1200 1100]); hold on; %eval(mp); eval(mg);
xlim([lon1 lon2]); ylim([lat1 lat2]);

for i = 1:numel(Segment.lon1)
    plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 0.00*[1 1 1], 'linewidth', 6)
end
for i = 1:numel(Segment.lon1)
    plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 1.00*[1 1 1], 'linewidth', 3)
end

% VE = -VSDEast + VBEast;
% VN = -VSDNorth + VBNorth;
for i = 1:numel(Mod.lon)
    %                 arrow( llons(i,j),Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
    %   axis(axis)
    a2 = arrow([Mod.lon(i) Mod.lat(i) 600], [Mod.lon(i)+scale(ii)*VE(i), Mod.lat(i)+scale(ii)*VN(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'b', 'edgecolor', 'k', 'lineWidth', 2);
    
    %         set(a, 'facecolor', 'w', 'edgecolor', 'k', 'lineWidth', 1);
    %         m_vec(35, llons(i,j), llats(i,j), VBEast(i,j), VBNorth(i,j), arrowcol, 'shaftwidth', sw, 'headlength', 6.0, 'headangle', 60,'edgecolor', 'k', 'linewidth', 1.5);
end
for i = 1:numel(Mod.lon)
    
    a3 = arrow([Obs.lon(i) Obs.lat(i) 600], [Obs.lon(i)+scale(ii)*Mod.eastVel(i), Obs.lat(i)+scale(ii)*Mod.northVel(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
    
end
a = arrow([36 42 600], [36+scale(ii)*20, 42 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
text(36, 42.2, '20 mm/yr');

set(gca,'xtick', [26 28 30 32 34 36 38], 'ytick', [38 39 40 41 42 43])
ticks_format('%.0f','%.0f');
[hx,hy] = format_ticks(gca,'^{\circ}E','^{\circ}N');
set(gca, 'fontName', 'times');
titlestr = 'BlocksForward in blue, Blocks in red';
title(titlestr, 'interpreter', 'latex');
xlim([lon1 lon2]); ylim([lat1 lat2]);

% export_fig(sprintf('%s/%s', figfolder, 'Forward_and_Modeled_Vels.pdf'))


%%


% figure('Position', [0 0 1200 1100]); hold on; %eval(mp); eval(mg);
% xlim([lon1 lon2]); ylim([lat1 lat2]);
% 
% for i = 1:numel(Segment.lon1)
%     plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 0.00*[1 1 1], 'linewidth', 6)
% end
% for i = 1:numel(Segment.lon1)
%     plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 1.00*[1 1 1], 'linewidth', 3)
% end
% 
% 
% for i = 1:numel(Mod.lon)
%     %                 arrow( llons(i,j),Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
%     %   axis(axis)
%     a2 = arrow([Mod.lon(i) Mod.lat(i) 600], [Mod.lon(i)+scale(ii)*VE(i), Mod.lat(i)+scale(ii)*VN(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'b', 'edgecolor', 'k', 'lineWidth', 2);
%     
%     %         set(a, 'facecolor', 'w', 'edgecolor', 'k', 'lineWidth', 1);
%     %         m_vec(35, llons(i,j), llats(i,j), VBEast(i,j), VBNorth(i,j), arrowcol, 'shaftwidth', sw, 'headlength', 6.0, 'headangle', 60,'edgecolor', 'k', 'linewidth', 1.5);
% end
% for i = 1:numel(Mod.lon)
%     
%     a3 = arrow([Obs.lon(i) Obs.lat(i) 600], [Obs.lon(i)+scale(ii)*Mod.eastVel(i), Obs.lat(i)+scale(ii)*Mod.northVel(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
%     
% end
% a = arrow([36 42 600], [36+scale(ii)*20, 42 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
% text(36, 42.2, '20 mm/yr');
% 
% set(gca,'xtick', [26 28 30 32 34 36 38], 'ytick', [38 39 40 41 42 43])
% ticks_format('%.0f','%.0f');
% [hx,hy] = format_ticks(gca,'^{\circ}E','^{\circ}N');
% set(gca, 'fontName', 'times');
% titlestr = 'BlocksForward in blue, Blocks in red';
% title(titlestr, 'interpreter', 'latex');
% xlim([lon1 lon2]); ylim([lat1 lat2]);
% 
% export_fig(sprintf('%s/%s', figfolder, 'Forward_and_Modeled_Vels.pdf'))

% 
% outname = 'BlockModelVels.sta.data';
% S.lon = staLons; S.lat = staLats; %S.name = staNames;
% S.eastVel = VE; S.northVel = VN;
% S.eastSig = zeros(size(VE)); S.northSig = zeros(size(VN));
% S.corr = zeros(size(VN)); S.other1 = zeros(size(VN));
% S.tog = zeros(size(VN));
% % keyboard;
% WriteStation(outname, S.lon, S.lat, S.eastVel, S.northVel, S.eastSig, S.northSig, S.corr, S.other1, S.tog, cell2mat(staNames))

% end



