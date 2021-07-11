%% Figure 5 for Case and Kingslake 2021

% plots accumulation rate and elevation across the divide for Korff and
% Skytrain
% Requires ExtractAccumulation.m and output from REMA (in REMAdivides)

%%
locs = {'ST','KF'};
for iloc = 1:length(locs)
    
fs = locs{iloc};
RhoI = 906;
nI = 1.78;
nArrFun = @(rho) 1 + (nI -1)*rho/RhoI;

switch fs
    case 'FP'
        core = load_FP_core;
        RhoS = interp1(core.d,core.rho,0.0,'linear','extrap');
        load NXYwBHe_FP_pRES
        pRES_data = NXYwBHe_FP_pRES;  
        loctitle='Fletcher Promontory';
        
    case 'ST'
        core = load_ST_core;
        RhoS = interp1(core.d,core.rho,0.0,'linear','extrap');
        load NXYwBHe_STpRES_for_Carlos_110719
        load ST_REMAdivide.mat
        pRES_data = NXYwBHe_STpRES;
        loctitle = 'Skytrain Ice Rise';
        cc = '#b84300'; %match QGIS color
        fig = figure('Units','centimeters','Position',[10,10,8.5,10]);
        subplot(211)
        text(.95,.95,'a','Units','Normalized','FontSize',10)
        
    case 'KF'
        load KF_pRES_DxS_09_04_20.mat
        load KF_REMAdivide.mat
        RhoS = [stake_pRES_DxS.RhoS];
        load('NXYwBHe_KorffpRES_for_Carlos_220515.mat');
        pRES_data= NXYwBHe_KorffpRES;        
        loctitle = 'Korff Ice Rise';
        subplot(212)
        text(.95,.95,'b','Units','Normalized','FontSize',10)
        cc = '#b8a000';
end

x = [divide.dist];
clear pRES_data

% exctract bI
bI = ExtractAccumulation(fs);
Lave = mean(bI(x < -500));
Rave = mean(bI(x > 500));

hold on;
h1 = plot(x(x>500),Rave*ones([1,length(x(x>500))]), 'k-','HandleVisibility','off');
h2 = plot(x(x<-500),Lave*ones([1,length(x(x<-500))]), 'k-','HandleVisibility','off');


for i = 1:length(x)-1
    yyaxis right;
    h3 = plot(x(i),divide(i).elev,'o','MarkerFaceColor','#A3BED1','HandleVisibility','off','MarkerSize',4);
    ylim([min([divide.elev]) max([divide.elev])])
    ylabel('Elevation (m)','color','#A3BED1')
    set(h3, 'color', '#A3BED1')
    
    yyaxis left
    h4 = plot(x(i),bI(i),'^','Color',[0,0.0314,0.2314],'HandleVisibility','off','MarkerSize',5)
    ylabel('Accumulation rate (m i.e. a^{-1})')
    xlabel('Distance from divide (m)')
    
    set(gca, 'SortMethod', 'depth')
 
end
yyaxis right
h5 = plot(x(end),divide(end).elev,'o','MarkerFaceColor','#A3BED1','HandleVisibility','off','MarkerSize',4);
ylim([min([divide.elev]) max([divide.elev])])
ylabel('Elevation (m)')
set(h5, 'color', '#A3BED1')

yyaxis left
h6 = plot(x(end),bI(end),'>','Color',cc,'MarkerFaceColor',cc,'LineWidth',1,'HandleVisibility','off','MarkerSize',5);

if strcmp(fs,'KF')
    xlim([-2500 2500])
end
legend('off')
set(gca, 'FontSize', 8,'Box','on')
title(loctitle,'FontSize',10)
end
text(0.6,1.5,'— avg acc. rate of flank','Margin',1.2,'Background','white','Units','Normalized','FontSize',6,'EdgeColor','black')
text(0.6,0.09,'— avg acc. rate of flank','Margin',1.2,'Background','white','Units','Normalized','FontSize',6,'EdgeColor','black')


