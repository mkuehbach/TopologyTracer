%% IMPORT FOLLOWING GRAINS RESULTS
% allows following the evolution of an individual grain in response to its
% neighbors
clear;
clc;

nrows = 257; %how many disjoint neighbors ever?
ncols = 1001; %how many timesteps?
offset = 1;

size=nrows*ncols*8;

id = 25843;

%% reading binary
prefix = ['followingtarget.' num2str(id) '.'];
postfix = '.2015.first.0.offset.10.last.10000.bin';

fileID = fopen([prefix 'volume' postfix]);
VOLUME = fread(fileID,[nrows,ncols],'double');
fclose(fileID);

fileID = fopen([prefix 'gbtype' postfix]);
GBTYPE = fread(fileID,[nrows,ncols],'double');
fclose(fileID);


fileID = fopen([prefix 'nfaces' postfix]);
NFACES = fread(fileID,[nrows,ncols],'uint32');
fclose(fileID);

%% import time scheme
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

IDS = zeros(nrows,ncols);
for gr=1:nrows
    for c=1:ncols
        %%TimeNGrains logs for fixed offset of the simulation but
        %%analysis can exclude time steps
        cid = (c-1)*offset + 1;
       IDS(gr,c) = TimeNGrainsGBLength(cid,1);
    end;
end

%% load the target grain itself
id = 25843;
tgrfile = ['followinggrain.' num2str(id) '.Evolution.2015.csv'];
newData1 = importdata(tgrfile);
%data <time;fid;nfaces;area
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
ntlen = length(data(:,1));
IDTARGET = zeros(1,ntlen);
NFACESTARGET = zeros(1,ntlen);
VOLTARGET = zeros(1,ntlen);
GBTARGET = zeros(1,ntlen);
for c=1:ntlen
    IDTARGET(1,c) = IDS(1,c); %%ONLY WORKS IF TIMENGRAINS AND FOLLOW START AT 0!!!
    NFACESTARGET(1,c) = data(c,3);
    VOLTARGET(1,c) = data(c,4);
    GBTARGET(1,c) = data(c,5);
end;

%% main plot routine
figure
n=ncols;
for gr=1:nrows
    hold on
     
    %plot3(IDS(gr,:), NFACES(gr,:),VOLUME(gr,:),'b','LineWidth',2,'Color',[(gr/nrows) (gr/nrows) 1]); %,'.'); %( % 1) fac(:)],'.'),%'.')
    
    %NFACES-VOLUME
    %plot3(IDS(gr,:), NFACES(gr,:),VOLUME(gr,:),'-','LineWidth',2,'Color',[0 0 0]+((gr/nrows)*0.8) ); %,'.'); %( % 1) fac(:)],'.'),%'.')
    %NFACES-GBTYPE
    %DUMMY TOGGLE NOT WORKING
%     st = 1;
%     se = 0;
%     for cc=1:n
%         if ( NFACES(gr,cc) <= 0 ) 
%             st = cc;
%         end
%         if (se <= 0)
%             if ( NFACES(gr,cc) > 0 )
%                 se = 0;
%             else
%                 se = cc;
%             end
%         end
%     end
%         
    plot3(IDS(gr,1:n), NFACES(gr,1:n),GBTYPE(gr,1:n),'-','LineWidth',2,'Color',[0 0 0]+((gr/nrows)*0.8) ); %,'.'); %( % 1) fac(:)],'.'),%'.')
    %VOLUME-GBTYPE
    %plot3(IDS(gr,:), VOLUME(gr,:),GBTYPE(gr,:),'-','LineWidth',2,'Color',[0 0 0]+((gr/nrows)*0.8) ); %,'.'); %( % 1) fac(:)],'.'),%'.')
    
    %plot3(GBTARGET(gr,:),NFACES(gr,:),VOLUME(gr,:),'-','LineWidth',2,'Color',[0 0 0]+((gr/nrows)*0.8) ); %,'.');
end

figure
color = [1 0 0];
hold on
%add the target grain
%plot3(IDTARGET(1,:), NFACESTARGET(1,:), VOLTARGET(1,:),'b','LineWidth',4,'Color',[1 0.5 0]);
plot3(IDTARGET(1,1:n), NFACESTARGET(1,1:n), GBTARGET(1,1:n),'b','LineWidth',4,'Color',[0 1 0]);



plot3(GBTARGET(1,:),NFACESTARGET(1,:), VOLTARGET(1,:),'b','LineWidth',4,'Color',color);
hold on
plot3(GBTARGET(1,1),NFACESTARGET(1,1), VOLTARGET(1,1),'o','MarkerSize',10,'Color',color); %'LineWidth',4,'Color',[1 0.5 0]);
hold on
plot3(GBTARGET(1,n),NFACESTARGET(1,n), VOLTARGET(1,n),'.','MarkerSize',40,'Color',color); %'LineWidth',4,'Color',[1 0.5 0]);



xlim([1 max(IDS(1,:))]);
%set(gca,'xscale','ln');

%NFACES-VOLUME
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Times')
xlabel({'Time (s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
ylabel({'Number of faces'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
zlabel({'Normalized grain area'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')

%NFACE-GBTYPE
xlabel({'Time (s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
ylabel({'Number of faces'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
zlabel({'$m \gamma$ ($m^2$/s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
print -dpng SCubeBicrystal.25843.AgainstNBorsNFacesTime.png

%VOLUME-GBTYPE
xlabel({'Time (s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
ylabel({'Number of faces'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
zlabel({'Mobility * Energy ($m^2$/s)$m \gamma$ ($m^2$/s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')

%GBFACE,NFACES,VOLUME
xlabel({'$m \gamma$ ($m^2$/s)'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
ylabel({'Number of faces'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')
zlabel({'Normalized grain area'},'FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times')


view(-45,14);
grid('on');
print -dpng SCubeBicrystal.FastestSlowest.3D.png

view(0,0); %vol=f(time)

%DUMMY
view(-45,14);
print -dpng SCubeBicrystal.25843.VolTopoTime.png
view(0,0);
print -dpng  SCubeBicrystal.25843.VolTime.png
view(-90,0); %vol=f(topo)
print -dpng  SCubeBicrystal.25843.VolTopo.png
view(0,90);
print -dpng  SCubeBicrystal.25843.TopoTime.png

epsfname = ['Rios2015SCubeBicrystal' num2str(id) 'VolTime.eps'];
pngfname = ['Rios2015SCubeBicrystal' num2str(id) 'VolTime.png'];

print -depsc epsfname
view(0,90);
print -dpng pngfname


%zlim([1e-12 1e-9])
%set(gca,'Units','normalized','YTick',1:1:15,'XTick',0:5:ncols,'Position',[.15 .2 .75 .7],'FontUnits','points','interpreter','latex','FontWeight','normal','FontSize',9,'FontName','Times')
% LATEX compatible output
grid('on');



xlabel('Timestep');
ylabel('Number of faces');
zlabel('Normalized volume');
grid('on');
view(0,0);
view(-45,14);
%view(-35,14);
view(-90,0); %vol=f(topo)
view(0,90); %topo=f(time)

figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',12,'fontName','Calibri'); %'fontWeight','bold')
%axes('fontName','Calibri','fontSize',12);
axesHandle = gca;
set(gca,'fontName','Calibri','fontSize',12);


% DEBUG
for r=1:nrows  
for c=1:ncols 
    if NFACES(r,c) > 20  
        NFACES(r,c) = 0;  
    end
end
end
