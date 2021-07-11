%generates Figure 4 from Case and Kingslake 2021

locations = {'KF','ST'};
figure('Units','centimeters','Position',[10,10,17.9,10]);

for iloc = 1:length(locations)

FS = locations{iloc};
type = 'pres'; %'pres' for pRES-defined DxS or 'fullmin' for 4-parameter inversion

switch FS
    case 'KF'
       
        load KF_pRES_DxS_12_17_20.mat;
        load KF_min_DxS_12_17_20.mat;
        load ST_REMAdivide.mat
        
        loctitle = 'Korff Ice Rise';
        iplot = [1,8,22,30];
        
    case 'ST'
        load ST_min_DxS_04_05_21_allsites.mat;
        load ST_pRES_DxS_04_05_21_allsites.mat;      
        load ST_REMAdivide.mat
        
        loctitle = 'Skytrain Ice Rise';
        iplot = [1,5,19,29];
end

switch type
    case 'pres'
        stake = stake_pRES_DxS;
    case 'fullmin'
        stake = stake_min_Dxs;
end
        
StakeOK=unique([stake.number]);
cmap = colormap(bone(length(StakeOK)+5));
cmap2 = {'#2a9d8f','#a44a3f','#e59500','#264653'};
n = 0;
nCol = 1;

for i = StakeOK
    n = 1 + n;
    

    %define vars
    wFirnPres = stake(i).wFirnInv_InvDensity;
    d = stake(i).d;

    switch FS
        case 'KF'
            subplot(221);
            
        case 'ST'
            subplot(222);
            
    end
    
    hold on;
    scatter(stake(i).wI, stake(i).depthI,12, 'o', 'filled','MarkerFaceColor', cmap(n,:,:))                
    title('Total vertical velocity')
    xlabel('Velocity (m a^{-1})')
    ylabel('Depth {\it\zeta} (m)')
    ylim([0 200])
    axis ij
    alpha 0.5
    set(gca,'FontSize',8,'Box','on')
    switch FS
        case 'KF'
            subplot(223);
            
        case 'ST'
            subplot(224);
            
    end    
    
    hold on;
    scatter(wFirnPres, d, 12,'o', 'filled','MarkerFaceColor', cmap(n,:,:))                
    title('Compaction velocity')
    xlabel('Velocity (m a^{-1})')
    ylabel('Depth {\it\zeta} (m)')
    ylim([0 200])
    axis ij 
    alpha 0.5
    set(gca,'FontSize',8,'Box','on')
    
end


    switch FS
        case 'KF'
            subplot(221);
            text(.95,.95,'a','Units','Normalized','FontSize',10)
            xlim([0 0.33])
            box on
        case 'ST'
            subplot(222);
            xlim([0 0.33])
            text(.95,.95,'b','Units','Normalized','FontSize',10)
            box on
    end
    
    switch FS
        case 'KF'
            subplot(223);
            text(.95,.95,'c','Units','Normalized','FontSize',10)
            xlim([0 0.33])
            box on
        case 'ST'
            subplot(224);
            text(.95,.95,'d','Units','Normalized','FontSize',10)
            xlim([0 0.33])
            box on
    end        
end

hc = colorbar('Location','southoutside','Ticks',[]);
set(hc, 'Position', [.44,.01,.15,.025])
hct = get(hc, 'Title');
set(hct, 'String', {'\it Location along transect'})

h1 = annotation('textbox',[0.22 0.96 0.2 0.06],...
    'string','Korff Ice Rise');
h2 = annotation('textbox',[0.638 0.96 0.2 0.06],...
    'string','Skytrain Ice Rise');
set([h1 h2], 'fitboxtotext','on',...
    'edgecolor','none','FontSize',12)

text(.07,-.26,'$\triangleright$','Units','Normalized','Interpreter','latex','FontSize',12);
