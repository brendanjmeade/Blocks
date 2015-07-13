%Runs BlocksForward for a grid of stations

clear all; close all;

[staNames, staLats, staLons] = textread('Izmit_NAF_Station_List_Hetland.csv','%s%f%f', 'delimiter', ',');

latlimmin = 37;
latlimmax = 43;
lonlimmin = 25;
lonlimmax = 38;

staLatsi = linspace(latlimmin, latlimmax, 20);
staLonsi = linspace(lonlimmin, lonlimmax, 25);
[staLons, staLats] = meshgrid(staLonsi, staLatsi);
staLons = staLons(:);
staLats = staLats(:);
folderBlockResults = 'ElasticBlockResults_ForHeltandJune2015/';
% folderBlockResults ='0000000001test/';
inSegFile = sprintf('%s/Mod.Segment',folderBlockResults);
ObsFile = sprintf('%s/Obs.Sta',folderBlockResults);
ModFile = sprintf('%s/Mod.Sta',folderBlockResults);

Obs = ReadStation(ObsFile);
Mod = ReadStation(ModFile);
RotFile = sprintf('%s/Rot.Sta',folderBlockResults);
Rot = ReadStation(RotFile);

Segment = ReadSegmentStruct(inSegFile);

% figfolder = 'HetlandBlocksJune2015Figs';
% mkdir(figfolder);
%
% [VB, VSD, VT, VS] = BlocksForward_Eileen(staLons, staLats, folderBlockResults, []);
%
% [V] = BlocksForward_Eileen(staLons, staLats, folderBlockResults, []);

[VB, VSD, VT, VS] = BlocksForward(staLons, staLats, folderBlockResults);

[V] = BlocksForward(staLons, staLats, folderBlockResults);


% [VB, VSD, VT, VS] = BlocksForward(Mod.lon, Mod.lat, folderBlockResults);
% [V] = BlocksForward(Mod.lon, Mod.lat, folderBlockResults);


VBEast = (VB(1:3:end));
VBNorth = (VB(2:3:end));

VSDEast = (VSD(1:3:end));
VSDNorth = (VSD(2:3:end));

VSEast = (VS(1:3:end));
VSNorth = (VS(2:3:end));

VE = (V(1:3:end));
VN = (V(2:3:end));

% VE = VBEast+VSDEast;%(V(1:3:end));
% VN = VBNorth+VSDNorth;%(V(2:3:end));

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
% [htopox, htopoy] = mfwdtran(mstruct, topo_NAF(:, 2), topo_NAF(:, 1));
% topo_NAF(:, 1) = htopox*6378;
% topo_NAF(:, 2) = htopoy*6378;
% N = 1000;
% % Interpolate topography onto a regular grid.
% lonvec = linspace(lon1, lon2, N);
% latvec = linspace(lat1, lat2, N);
% [xvec, yvec] = mfwdtran(mstruct, latvec, lonvec);
% xvec=xvec*6378;
% yvec=yvec*6378;
% [xmat ymat] = meshgrid(xvec, yvec);
%
% tic;
% zmat = griddata(topo_NAF(:,1), topo_NAF(:,2), topo_NAF(:,3), xmat, ymat);

%%
for ii = 1:4
    if ii == 1
        fieldE = VBEast;
        fieldN = VBNorth;
        titlestr ='$\mathbf{v}_{\mathrm{B}}$';
        figName = 'VBHetland.jpg';
        scale(ii) = .04;
        
    elseif ii == 2
        fieldE = VSDEast;
        fieldN = VSDNorth;
        figName = 'SDHetland.jpg';
        titlestr ='$\mathbf{v}_{\mathrm{SD}}$';
        scale(ii) = .04;
        
    elseif ii == 3
        fieldE = VSEast;
        fieldN = VSNorth;
        figName = 'Locations.jpg';
        titlestr ='Locations, Black Circles = Hetland Stations, Red Circles = Our Stations';
        scale(ii) = .04;
        
    elseif ii == 4
        fieldE = VE;
        fieldN = VN;
        titlestr='$\mathbf{v}_{\mathrm{I}}$';
        figName = 'VTotalHeltand.jpg';
        scale(ii) = .04;
        
        %     elseif ii == 5
        %         fieldE = VE_East;
        %         fieldN = VE_North;
        %         titlestr='$\mathbf{v}_{\mathrm{VE}}$';
        %         figName = 'VEWithTopoJet.jpg';
        %     elseif ii == 6
        %         fieldE =  VE_meanEast;
        %         fieldN =  VE_meanNorth;
        %         figName = 'VEMeanWithTopoJet.jpg';
        %         titlestr='$\overline{\mathbf{v}}_{\mathrm{VE}} $';
        %
    end
    
    figure('Position', [0 0 1200 1100]); hold on; %eval(mp); eval(mg);
    
    % cmap2 = (colormap_cpt('19 hue sat light1'));
    % cmap2 = flipud(cmap2(1:end-30,:));
    % cmap2 = cmap2(1:20:end,:);
    % colormap(cmap2);
    % tic;
    % [latMat, lonMat] = meshgrid(latvec, lonvec);
    %
    % zmat = griddata(topo_NAF(:,1), topo_NAF(:,2), topo_NAF(:,3), xmat, ymat);
    % % field
    % magsBnew = zeros(size(lonMat));%griddata(llons, llats, fieldf, lonMat, latMat);
    % %     zmat2 = zmat;
    % %     zmat2 = medfilt2(zmat, [5 5], 'zeros');
    % zmatorig = zmat;
    % nanIdx = find(zmat<0);
    % zmat(nanIdx) = 0;
    %surface(lonMat, latMat, h', magsBnew);
    % h=surface(lonMat, latMat, zmat'/10000, magsBnew);
    % % contour3(lonMat, latMat, zmatorig, 0, 'edgecolor', 'k');
    % c = contour3(lonMat, latMat, zmatorig'/10000, [0 0], 'k');
    % %     c = contour3(lonMat, latMat, magsBnew+100, linspace(lim0+100, lim1+100, 12), 'w');
    % caxis([0 10]);
    % shading interp
    % % material shiny
    % % shading flat;
    % axis equal;
    % %camlight(90, 0)
    % %camlight(45, 0)
    % % camlight left
    % lightangle(-45,30)% camlight right
    xlim([lon1 lon2]); ylim([lat1 lat2]);
    % lightangle(-45,30)% camlight right
    
    for i = 1:numel(Segment.lon1)
        plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 0.00*[1 1 1], 'linewidth', 6)
    end
    for i = 1:numel(Segment.lon1)
        plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 1.00*[1 1 1], 'linewidth', 3)
    end
    %
    
    % for i = 1:numel(Segment.lon1)
    %     p = plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'k');
    %     set(p, 'linewidth', 8);
    % end
    % keyboard;
    % ax1 = axes('position', pos, 'visible', 'off')
    if ii == 3
        for i = 1:numel(Mod.lon)
            p = plot3(Mod.lon(i), Mod.lat(i), 600, 'ro', 'markersize', 12, 'markerfacecolor', 'r');
        end
        for i = 1:numel(VBEast)
            p = plot3(staLons(i), staLats(i), 650, 'ko', 'markersize', 8, 'markerfacecolor', 'k');
            
        end
%         keyboard;
    end
    for i = 1:numel(VBEast)
        %                 arrow( llons(i,j),Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
        %   axis(axis)
        if (fieldE(i) ~= 0 && fieldN(i) ~= 0)
            a1 = arrow([staLons(i) staLats(i) 600], [staLons(i)+scale(ii)*fieldE(i),staLats(i)+scale(ii)*fieldN(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
            %         set(a, 'facecolor', 'w', 'edgecolor', 'k', 'lineWidth', 1);
            %         m_vec(35, llons(i,j), llats(i,j), VBEast(i,j), VBNorth(i,j), arrowcol, 'shaftwidth', sw, 'headlength', 6.0, 'headangle', 60,'edgecolor', 'k', 'linewidth', 1.5);
        end
    end
    if ii == 2
        fieldSDE = (Rot.eastVel-Mod.eastVel);
        fieldSDN = (Rot.northVel-Mod.northVel);
        for i = 1:numel(Mod.lon)
            a1 = arrow([Mod.lon(i) Mod.lat(i) 650], [Mod.lon(i)+scale(ii)*fieldSDE(i), Mod.lat(i)+scale(ii)*fieldSDN(i) 650], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'k', 'edgecolor', 'k', 'lineWidth', 1);
            %         set(a, 'facecolor', 'w', 'edgecolor', 'k', 'lineWidth', 1);
            %         m_vec(35, llons(i,j), llats(i,j), VBEast(i,j), VBNorth(i,j), arrowcol, 'shaftwidth', sw, 'headlength', 6.0, 'headangle', 60,'edgecolor', 'k', 'linewidth', 1.5);
        end
    end
    a = arrow([36 42 600], [36+scale(ii)*20, 42 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
    text(36, 42.2, '20 mm/yr');
    
    set(gca,'xtick', [26 28 30 32 34 36 38], 'ytick', [38 39 40 41 42 43])
    ticks_format('%.0f','%.0f');
    [hx,hy] = format_ticks(gca,'^{\circ}E','^{\circ}N');
    set(gca, 'fontName', 'times');
    title(titlestr, 'interpreter', 'latex');
    
%     export_fig(sprintf('%s/%s', figfolder, figName))
    % keyboard;
    
end




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
titlestr = 'Modeled in blue, Observed in red';
title(titlestr, 'interpreter', 'latex');
xlim([lon1 lon2]); ylim([lat1 lat2]);

% export_fig(sprintf('%s/%s', figfolder, 'Obs_Mod_Vels.pdf'))






figure('Position', [0 0 1200 1100]); hold on; %eval(mp); eval(mg);
xlim([lon1 lon2]); ylim([lat1 lat2]);

for i = 1:numel(Segment.lon1)
    plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 0.00*[1 1 1], 'linewidth', 6)
end
for i = 1:numel(Segment.lon1)
    plot3([Segment.lon1(i) Segment.lon2(i)], [Segment.lat1(i) Segment.lat2(i)], [400 400], 'color', 1.00*[1 1 1], 'linewidth', 3)
end
fieldE = VE;
fieldN = VN;
for i = 1:numel(staLons)
    %                 arrow( llons(i,j),Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
    %   axis(axis)
    a2 = arrow([staLons(i) staLats(i) 600], [staLons(i)+scale(ii)*VE(i), staLats(i)+scale(ii)*VN(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'b', 'edgecolor', 'k', 'lineWidth', 2);
    
    %         set(a, 'facecolor', 'w', 'edgecolor', 'k', 'lineWidth', 1);
    %         m_vec(35, llons(i,j), llats(i,j), VBEast(i,j), VBNorth(i,j), arrowcol, 'shaftwidth', sw, 'headlength', 6.0, 'headangle', 60,'edgecolor', 'k', 'linewidth', 1.5);
end
for i = 1:numel(Mod.lon)
    
    a3 = arrow([Mod.lon(i) Mod.lat(i) 600], [Mod.lon(i)+scale(ii)*Mod.eastVel(i), Mod.lat(i)+scale(ii)*Mod.northVel(i) 600], 'length', 10, 'baseangle', 90, 'tipangle', 40, 'width', 8, 'facecolor', 'r', 'edgecolor', 'k', 'lineWidth', 2);
    
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


% outname = 'BlockModelVels.sta.data';
% S.lon = staLons; S.lat = staLats; %S.name = staNames;
% S.eastVel = VE; S.northVel = VN;
% S.eastSig = zeros(size(VE)); S.northSig = zeros(size(VN));
% S.corr = zeros(size(VN)); S.other1 = zeros(size(VN));
% S.tog = zeros(size(VN));
% % keyboard;
% WriteStation(outname, S.lon, S.lat, S.eastVel, S.northVel, S.eastSig, S.northSig, S.corr, S.other1, S.tog, cell2mat(staNames))

% end


