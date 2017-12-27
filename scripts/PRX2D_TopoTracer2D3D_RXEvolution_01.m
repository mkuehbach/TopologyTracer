%% TopologyTracer2D3D
% PRX2D_RXEvolution
% Author: m.kuehbach (at) mpie.de, 09/08/2017
% MinRequirement: Matlab v2017a
% Purpose: visualize evolution of area fraction survivors versus all = f(t)
clear;
clc;
format long;
digits(32);

%% user interaction
prefix_track = 'E:\LongRangePaperFINAL\PRX2D\2DataAnalysis\TrackingParallel\TrackingFW_Resolution_1\105200\';
realtime_fn = 'E:\LongRangePaperFINAL\PRX2D\1Coarsening\NrGrains&EnergyStatistics.txt';
first = 10;
offset = 1;
last = 5200;
simid = 105200;
PhysDomainSize = 4.2037e-3; %meter
HAGBEnergy = 1.000; %J/m^2
Dimension = 2;
DomainSize = (PhysDomainSize/1.0e-6)^Dimension;
LargestGrainID = 9900000;

ncols = 5191;
nrows = 1; %because struct of 5 contiguous doubles

%% end of user interaction

%% load RX
prgnm = ['TopoTracer2D3D.SimID.' num2str(simid) '.'];
suffix = ['.F.' num2str(first) '.O.' num2str(offset) '.L.' num2str(last) '.NC.' ...
    num2str(ncols) '.NR.1.bin'];

%% load backtracking grain history in memory
fileID = fopen([prefix_track prgnm 'RXEVO' suffix]);
RX = fread(fileID,[5,ncols],'double');
fclose(fileID);
'Binary raw data successfully read'

%% load realtime data
delimiter = '\t';
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(realtime_fn,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
REALTIME = [dataArray{1:end-1}];
clearvars delimiter formatSpec fileID dataArray ans;


for fid=first:offset:last
    RX(6,(fid-first)/offset+1) = REALTIME(fid,1);
end

X=RX(3,:);
T=RX(6,:);
%% Avrami linearization
% lnT = log(T);
% lnX = log(-1.0 .* log(1.0-X));
% gradient(lnT,lnX)
% max(lnX)
% (lnX(1,1000)-lnX(1,999))/(lnT(1,1000)-lnT(1,999))
% plot(T(1,:),X(1,:),'.')

%% plot "kinetics"
fontsz = 22;
fontnm = 'Calibri Light';
figure('Position',[100 100 1200 1000])
hold('on')
plot(T,X,'-','Color',[0 0 1],'LineWidth',2);
set(gca, 'XScale', 'log')
xlabel({'Annealing time (s)'},'FontSize',fontsz,'FontName',fontnm);
ylabel({'\Sigma_{srv} A_i / \Sigma_{all} A_j \approx X'},'FontSize',fontsz,'FontName',fontnm);
%projections
grid('on');
box('on');
set(gca,'FontSize',fontsz,'FontName',fontnm,'LineWidth',1.5);
set(gcf,'PaperUnits','Inches');
set(gcf,'PaperSize',[30 30]);
pbaspect([1 1 1]);
set(gcf,'color','w');
%colorbar;
%xlim([-0.02 1.02]);
%xt = [0:0.1:1.0];
%xticks(xt);
ylim([-0.02 1.02]);
yt = [0:0.1:1.0];
yticks(yt);
print(gcf,['PRX2D_RXEvolution_5200_01.png'],'-dpng','-r500');
close(gcf);

ORIGINXT(:,1) = reshape(T,[ncols,1]);
ORIGINXT(:,2) = reshape(X,[ncols,1]);
dlmwrite('PRX2D_RXEvolution_5200_01.dat',ORIGINXT,'delimiter','\t');

save('PRX2D_RXEvolution_5200_01.mat');

