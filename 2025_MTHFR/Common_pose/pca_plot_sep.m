close all; clear; clc
format shortG
load results_WT.mat

siz1=size(drt,2);
gs=30;
dat = [drt(1,:)' drt(2,:)' drt(3,:)' drt(4,:)'];
%dat için histogram elde etti

[N, edgesX, edgesY] = histcounts2(dat(:,1), dat(:,2), [gs gs]);     nn1=N/siz1;
[Xgrid, Ygrid] = meshgrid(edgesX(1:end-1), edgesY(1:end-1));
% Convert bin edges to centers for plotting purposes
centersX = edgesX(1:end-1) + diff(edgesX) / 2;
centersY = edgesY(1:end-1) + diff(edgesY) / 2;


figure(1)
h1=imagesc(centersX, centersY, nn1'); hold on
set(gca,'YDir','normal')
set(gca,'FontSize',15,'LineWidth',1.1)
colormap jet
cb=colorbar;
% shading interp
xlabel('PC1','FontName','Times'); ylabel('PC2','FontName','Times');


%% defining of individual runs 

run1=size(t_tot1,1);            %wtrun1
run2=run1+size(t_tot2,1);         %2
run3=run2+size(t_tot3,1);         %3
% run4=run3+size(t_tot4,1);       %mutant1 run1
% run5=run4+size(t_tot5,1);         %2
% run6=run5+size(t_tot6,1);         %3

%plot(dat(1:run1,1),dat(1:run1,2),'square','MarkerSize',2,'Color','y')
%plot(dat(run1:run2,1),dat(run1:run2,2),'square','MarkerSize',2,'Color','m')
%plot(dat(run2:run3,1),dat(run2:run3,2),'square','MarkerSize',2,'Color','g')



%% Find frames present of certain pose
lx=-16; %% Define the x/y bounds
ux=-13;
ly=-6;
uy=-2;

% Selection for PC1 % PC2
[row col]= find(0.015 < nn1); % Define the proportion that is the minima at this pose
[h1, g1]=find(centersX(row)>lx & centersX(row)< ux);
[h2, g2]=find(centersY(col)>ly & centersY(col)< uy);
PC12=intersect(g1,g2);
PC12n=[];
for k=1:size(PC12,2)
    [i1, j1]=find(dat(:,1)>edgesX(row(PC12(k))) & dat(:,1)< edgesX(row(PC12(k))+1));
    [i2, j2]=find(dat(:,2)>edgesY(col(PC12(k))) & dat(:,2)< edgesY(col(PC12(k))+1));
    PC12n=[PC12n; intersect(i1,i2)];
end

PC12n=sort(PC12n);

figure(1);
plot(dat(PC12n,1),dat(PC12n,2),'<','MarkerSize',4,'Color','m')

%selection of frames for individual trajectorys
selPC=PC12n;
szi=[size(t_tot1(1:end,:),1); size(t_tot2(1:end,:),1); size(t_tot3(1:end,:),1);]; %all dcdcs length 
szi1=[szi(1); szi(1)+szi(2); szi(1)+szi(2)+szi(3);]; %all data

rn1a=find(selPC<szi1(1)+1);                 run1=((selPC(rn1a))-1)';
rn2a=find(selPC>szi1(1) & selPC<szi1(2)+1); run2=((selPC(rn2a)-(szi(1)))-1)';
rn3a=find(selPC>szi1(2) & selPC<szi1(3)+1); run3=((selPC(rn3a)-((szi(1)+szi(2))))-1)';



run1=sort(run1);
run2=sort(run2);
run3=sort(run3);

disp('run1');
fprintf('%d ', run1);
fprintf('\n');
disp('run2');
fprintf('%d ', run2);
fprintf('\n');
disp('run3');
fprintf('%d ', run3);
fprintf('\n');

%% ok
[numRows, numCols] = size(N);

for i = 1:numRows
    for j = 1:numCols
        % Get the count for the current cell
        count = N(j, i);

        % Place the text at the center coordinates
        % X_coords and Y_coords are matched to N(i, j)
        text(Xgrid(i, j)+0.5, Ygrid(i, j)+0.5, num2str(count), ...
             'HorizontalAlignment', 'Center', ...
             'VerticalAlignment', 'Middle', ...
             'Color', 'white','FontSize',10); % Optional: makes text white for better contrast
    end
end
