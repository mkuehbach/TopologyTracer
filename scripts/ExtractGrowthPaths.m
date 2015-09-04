%% IMPORT TOPOTRACER RESULTS
clear;
clc;

%import timemapping <time>,<ngrains>,<gblength>
fname = 'TimeNGrainsGBLength.txt';
rawData1 = importdata(fname);
[~,name] = fileparts(fname);
newData1.(genvarname(name)) = rawData1;
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
clearvars i name newData1 rawData1 vars

% synchronize timesteps
nrows = 24; %how many grains survived?
ncols = 1001; %how many timesteps?

%% reading binarys
suffix = '2015.first.0.offset.10.last.10000.bin';
fileID = fopen(['topotrace.volume.' suffix]);
VOLUME = fread(fileID,[nrows,ncols],'double');
fclose(fileID);

fileID = fopen(['topotrace.nfaces.' suffix]);
NFACES = fread(fileID,[nrows,ncols],'uint32');
fclose(fileID);

fileID = fopen(['topotrace.gbtype.' suffix]);
GBTYPE = fread(fileID,[nrows,ncols],'double');
fclose(fileID);

%% import time scheme
%evol = 'PopulationStatistics.dat';
%newData1 = importdata(evol);
% Create new variables in the base workspace from those fields.
%data <time;fid;nfaces;area
%vars = fieldnames(newData1);
%for i = 1:length(vars)
%    assignin('base', vars{i}, newData1.(vars{i}));
%end

IDS = zeros(nrows,ncols);
for gr=1:nrows
    for c=1:ncols
       IDS(gr,c) = TimeNGrainsGBLength(c,1); %c %data(c,1); %for real time
    end;
end

%% plot population evolution

figure;
plot(TimeNGrainsGBLength(:,1),TimeNGrainsGBLength(:,2),'b','LineWidth',2,'Color',[0 0 0]+0.0 );
grid('on');
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',20,'FontName','Times')
xlabel({'Annealing time (s)'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times')
ylabel({'Total number of grains'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times')
print -depsc Riso2015SCubeBiCrystalPopulationStatistics.eps
print -dpng Riso2015SCubeBiCrystalPopulationStatistics.png


%% calculate face switching matrix
FACESWITCH=zeros(nrows,ncols);
for gr=1:nrows
    for s=2:ncols
        FACESWITCH(gr,s) = NFACES(gr,s) - NFACES(gr,s-1);
    end
end

%% count potential for reduced information content
totalswitches = zeros(1,ncols);
totalneighborrels = zeros(1,ncols);
for c=1:ncols
    totalswitches(1,c) = sum(abs(FACESWITCH(:,c)));
    totalneighborrels(1,c) = sum(NFACES(:,c));
end
plot(totalneighborrels,totalswitches,'.','MarkerSize',14,'Color',[0 0 0]+0);
xlim([min(totalneighborrels(1,:))*0.95 max(totalneighborrels(1,:))*1.05]);
ylim([0 30]);
grid('on');
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',20,'FontName','Times')
xlabel({'Sum of faces in the system (w/o double counting)'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times')
ylabel({'Total number of switches per time step'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times')
print -depsc Riso2015SCubeBiCrystalTopologicalChanges.eps
print -dpng Riso2015SCubeBiCrystalTopologicalChanges.png



%% plotting
finalav=mean(VOLUME(:,ncols));
maxsize=max(VOLUME(:,ncols));
scale=40.0;
i=0;
for gr=1:nrows
    if (VOLUME(gr,ncols) > scale*finalav)
        i=i+1;
    end
end






%% main plot routine
figure
%xlim([0 7705])
maxx = 0.0;
maxid = -1;
view(-45,14);
grid('on');

for gr=1:nrows
    %if (VOLUME(gr,ncols) >= scale*finalav)    
    %if (VOLUME(gr,ncols) <= scale*finalav)
    
    %    if ( max(VOLUME(gr,:)) > maxx )
    %        maxid = gr;
    %    end
    %if (NFACES(gr,1) == 5)
        hold on
       
        %plot only time over volume but mark switches
%         plot(IDS(gr,:),VOLUME(gr,:),'Color',[0 0 1])
%         for s=1:ncols
%             if ( FACESWITCH(gr,s) == -1 )
%                 hold on
%                 plot(IDS(gr,s),VOLUME(gr,s),'.','Color',[1 0 0]) %[Red Green Blue]
%             end
%             if ( FACESWITCH(gr,s) == 1 )
%                 hold on
%                 plot(IDS(gr,s),VOLUME(gr,s),'.','Color',[0 0 0])
%             end
%         end
         
    %plot3(IDS(gr,:),NFACES(gr,:),VOLUME(gr,:),'b','LineWidth',2,'Color',[0 0 0]+(VOLUME(gr,ncols)/maxsize)*0.6 ); %,'.'); %( % 1) fac(:)],'.'),%'.')
    %plot3(IDS(gr,:),NFACES(gr,:),GBTYPE(gr,:),'b','LineWidth',2,'Color',[0 0 0]+(VOLUME(gr,ncols)/maxsize)*0.6 ); %,'.'); %( % 1) fac(:)],'.'),%'.')
    plot3(IDS(gr,:),VOLUME(gr,:),GBTYPE(gr,:),'b','LineWidth',2,'Color',[0 0 0]+(VOLUME(gr,ncols)/maxsize)*0.6 ); %,'.'); %( % 1) fac(:)],'.'),%'.')
     
    %end
end
%set(gca,'xscale','log');
xlim([min(IDS(1,:)) max(IDS(1,:))]);

set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Times')
%NFACES VOLUME
xlabel({'Time (s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
ylabel({'Number of faces'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
zlabel({'Normalized grain area'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')

%NFACES GBTYPE
xlabel({'Time (s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
ylabel({'Number of faces'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
zlabel({'Mobility*Energy ($m^2$/s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')

%VOLUME GBTYPE
xlabel({'Time (s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
ylabel({'Normalized grain area'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
zlabel({'Mobility*Energy ($m^2$/s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')

view(-45,14);
print -dpng SCubeBicrystal.VolTopoTime.png
print -dpng SCubeBicrystal.NFacesMobEnergyTime.png
print -dpng SCubeBicrystal.VolumeMobEnergyTime.png

view(0,0);
print -dpng SCubeBicrystal.VolTime.png
print -dpng SCubeBicrystal.MobEnergyTime.png
print -dpng SCubeBicrystal.MobEnergyTime.png

view(-90,0);
print -dpng SCubeBicrystal.VolTopo.png
print -dpng SCubeBicrystal.MobEnergyNFaces.png
print -dpng SCubeBicrystal.MobEnergyVolume.png

view(0,90);
print -dpng SCubeBicrystal.TopoTime.png
print -dpng SCubeBicrystal.NFacesTime.png


grid('on');
view(-45,14);
print -depsc Riso2015SCubeBicrystalAtBndVolTopoTime.10000.3D.eps
view(0,0);
print -depsc Riso2015SCubeBicrystalAtBndVolTime.10000.3D.eps
view(0,90);
print -depsc vepsRiso2015SCubeBicrystalAtBndTopoTime.10000.3D.eps




figure;
plot3(IDS(maxid,:),NFACES(maxid,:),VOLUME(maxid,:),'b','LineWidth',2,'Color',[0 0 1]); %,'.'); %( % 1) fac(:)],'.'),%'.')
grid('on');


%view(-35,14);
view(0,0); %vol=f(time)
view(-90,0); %vol=f(topo)
view(0,90); %topo=f(time)

figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',12,'fontName','Calibri'); %'fontWeight','bold')
%axes('fontName','Calibri','fontSize',12);
axesHandle = gca;
set(gca,'fontName','Calibri','fontSize',12);

print -dpng AllLargerThan3TimesTheAverage.png
%print -deps fname.eps
%print -depsc fname.eps



%% temporal evolution classes still buggy####
[hc,hw] = hist(NFACES(:,2));
nbins = length(hc);
faceevolution = zeros(nbins,ncols);
figure
for c=1:20:ncols
    [hc,hw] = hist(NFACES(:,c));
    %for b=1:nbins
        faceevolution(b,c) = hc(b);
    %end
    hold on
    plot(hw, hc,'Color',[0 0 c/ncols])
end
    
%% sorting
%HERE IS A BUG AT THE MOMENT ...
AS = sortrows(VOLUME,ncols,'descend');
A = sort(A,ncols,'ascend');


%% isocontour - still debuggy
figure
row = linspace(1,ncols-1,ncols-1);
col = linspace(1,ncols-1,ncols-1);
for r=1:ncols-1
    for c=1:ncols-1
        chi(r,c)=NFACES(r,c);
    end;
end;
pcolor(row,col,NFACES);

%[C, h] = contourf(ee0, eee, chi)
%set(h,'LineColor','none')
pcolor(ee0,eee,chi)
shading flat
%colormap(map)
colormap gray
colorbar
xlabel(['<110>'])
ylabel(['<111>'])
%mark sphere
%hold on
%plot(sq2,sq3,'*','MarkerSize',8,'Color',[1 1 1])
hold on
plot(0.7957,0.6494,'x','MarkerSize',10,'Color',[0 0 0])
%plot(0.7957,0.6494,'o','MarkerSize',8,'Color',[1 1 1])