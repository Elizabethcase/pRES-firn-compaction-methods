%% outputs Figure 2 for Case and Kingslake. (2021). Phase Sensitive Radar as a Tool for Measuring Firn Compaction. Journal of Glaciology.
% updated on July 10, 2021

% saves results in 'method 1/method1_output/'
% does not automatically save figure
% note: run Figure2.m in folder 'deliver/method 1/code'
%% reset workspace

clear
%close all 

%% load data

locations = {'FP','ST'};
s = 1; %1 to save

for iloc = 1:length(locations)
fs = locations{iloc};

switch fs
    case 'FP'
        %load('FletcherInv_012820.mat'); %inversion results
        load('../../model/output/FP_model_pres_0416212109.mat'); %model results
        pRESModelOutput = modelOutput;
        
        load('../../model/output/FP_model_density_0416212114.mat'); %model results
        densityModelOutput = modelOutput;
        
        loc = 'Fletcher Promontory';
        locshort = 'FP';
        core = load_FP_core;

        load('NXYwBHe_FP_pRES.mat');
        NXYwBHe_FP_pRES = flipud(NXYwBHe_FP_pRES);       

        pRES(1).stake = 1;
        pRES(1).dI = NXYwBHe_FP_pRES(:,3);
        pRES(1).wI = NXYwBHe_FP_pRES(:,4);
        pRES(1).eI = NXYwBHe_FP_pRES(:,7);
        
        stakeOK = [1];
        
        ub = 160; 
        lb = 309; 
        
        clear NXYwBHe_FP_pRES;    
        
    case 'ST'
        
        load('../../model/output/ST_model_pres_0415211128.mat')
        pRESModelOutput = modelOutput;
        
        load('../../model/output/ST_model_density_0415211145.mat')
        densityModelOutput = modelOutput;
        
        loc = 'Skytrain Ice Rise';
        locshort = 'ST';
        core = load_ST_core;
        
        load('NXYwBHe_STpRES_for_Carlos_110719')
        
        stakeOK = [18];
        
        ub = 122;
        lb = 304;
        
        
        for i = 1:length(stakeOK)
            
            I = NXYwBHe_STpRES(:,1)==stakeOK(i);
            
            pRES(i).stake = stakeOK(i);
            
            dI = NXYwBHe_STpRES(I,3);
            wI = NXYwBHe_STpRES(I,4);
            e = NXYwBHe_STpRES(I,7);

            pRES(i).dI = flipud(dI);
            pRES(i).wI = flipud(wI);
            pRES(i).eI = flipud(e);
        
        end
        clear NXYwBHe_STpRES

end

nI = 1.78;
rhoI = 906;
spy = 365.25*24*3600;
c = 2.983e8;
f = 305e6;


%% core preparation

%find depth of 550 and 830
[value, i550] = min(abs(core.rho-550));
[value, i830] = min(abs(core.rho-830));

td550 = core.d(i550);
td830 = core.d(i830);

    
% interpolate core

switch fs
    case 'FP'
        core.rhoS = 435;
        % interpolate first chunk of firn/snow
        dif = 0.8; % not always constant in core, but most consistent in fletcher core
    case 'ST'
        core.rhoS = 400;
        % interpolate first chunk of firn/snow
        dif = 0.8; % not always constant in core, but most consistent in fletcher core

end
rhoSadd = linspace(core.rhoS, min(core.rho)-dif, round((min(core.d)-dif)/dif));
dadd = linspace(0, min(core.d)-dif, round((min(core.d)-dif)/dif));
core.rhoAdd = [rhoSadd core.rho];
core.dAdd = [dadd core.d];

% translate core into T
n = @(rho) 1 + (nI -1)*rho/rhoI;

diffd = diff(core.dAdd);
midrho = (core.rhoAdd(1:end-1) + core.rhoAdd(2:end))/2;
core.T = cumsum(diff(core.dAdd).*n(midrho)./c);


%% prepare pRES data

for i=1:length(stakeOK)
    
    % translate icewz into T
    pRES(i).T = pRES(i).dI.*nI/c;
    
    %interpolate rho onto T of pRES measurements
    pRES(i).rho = interp1(core.T,midrho,pRES(i).T);
    
    %core goes deeper than pRES
    iGood = ~isnan(pRES(i).rho); 
    %% we want to translate dT/wI into w, and T into z

    % true depth of pRES measurements

    pRES(i).T = pRES(i).T(iGood);
    pRES(i).wI = pRES(i).wI(iGood);
    pRES(i).rho = pRES(i).rho(iGood);
    pRES(i).eI = pRES(i).eI(iGood);
    pRES(i).dI = pRES(i).dI(iGood);
    rhotemp = [[core.rhoS] [pRES(i).rho]'];
    midrho = (rhotemp(1:end-1) + rhotemp(2:end))/2;
    
%    d = cumtrapz([0 pRES(i).T']',c./n([435 pRES(i).rho'])');
    d = cumsum(diff([0 pRES(i).T']).*c./n(midrho));
    pRES(i).d = d';
    % true vertical velocity of pRES
    pRES(i).w = pRES(i).wI.*nI./n(pRES(i).rho);
    pRES(i).e = pRES(i).eI.*nI./n(pRES(i).rho);


%% remove ice dynamic component

    spy = 365.25*24*3600;

    pRES(i).wFirn = firnVelo(pRES(i).d, pRES(i).w, ub, lb);
    
    if s == 1
        save(['../method1_output/' fs '_method1Results_' datestr(now,'mmddyyHH')],'pRES');
    end
%% plot

    colorM1 = '#193764';%[0.0353,0.5686,0.4980];
    colorM1edge = '#C8A219';
    colorCore = '#B1B1B1';
    colorModelpRES = '#6787C8';%#955E42';%[0.1843,0.3294,0.3882];%[0.5686,0.7451,0.7882];
    colorModelcore = '#2A3A17';%[0.6392,0.4471,0];

    %velocities
    switch fs
        case 'FP'
            figure('Units','centimeters','Position',[5,10,17.9,10]);
            subplot(221)
        case 'ST'
            subplot(222)
    end
    hold on;
     
    errorbar(pRES(i).wFirn,pRES(i).d,pRES(i).e,'horizontal','s','DisplayName','Method 1','LineWidth',1,'Color',colorM1,'MarkerSize',8);

    plot(densityModelOutput.wFirn,densityModelOutput.d,':','DisplayName','Model tuned to core{\it \rho}','LineWidth',2,'Color',colorModelcore)    
    plot(pRESModelOutput.wFirn,pRESModelOutput.d,'--','DisplayName','Model tuned to{\it w_c}','LineWidth',2,'Color',colorModelpRES)

    plot(linspace(-.1,0.5,100),td550*ones([1 100]),'--','Color','#C6C6C6','HandleVisibility','off') 
    plot(linspace(-.1,0.6,100),td830*ones([1 100]),'-.','Color','#C6C6C6','HandleVisibility','off') 

    xlabel('Vertical velocity (m a^{-1})')
    ylabel('Depth (m)')
    title(['Compaction velocity'])
    axis ij
    legend('Location','southeast','FontSize',8)
    ylim([0 200])
    xlim([-.05 0.4])
    set(gca,'FontSize',8)
    box on
    
    switch fs
        case 'FP'
            text(.95,.9,'a','Units','Normalized','FontSize',10)
        case 'ST'
            text(.95,.9,'b','Units','Normalized','FontSize',10) 
    end
    %%
    %densities
    switch fs
        case 'FP'
            subplot(223)
            box on
            
        case 'ST'
            subplot(224)
            box on
    end

    
    splot=scatter(core.rho,core.d,'filled','DisplayName',['Core density'],'MarkerFaceColor',colorCore,'MarkerEdgeColor',colorCore);
    alpha(splot,0.3)
    hold on;
    plot(densityModelOutput.rho, densityModelOutput.d,':','Displayname','Model tuned to core{\it \rho}','LineWidth',2,'Color',colorModelcore)
    plot(pRESModelOutput.rho,pRESModelOutput.d,'--','DisplayName','Model tuned to{\it w_c}','LineWidth',2,'Color',colorModelpRES)
    plot(linspace(300,1200,100),td550*ones([1 100]),'--','Color','#C6C6C6','HandleVisibility','off') 
    plot(linspace(300,1200,100),td830*ones([1 100]),'-.','Color','#C6C6C6','HandleVisibility','off') 
    
    ylim([0 200])
    xlim([375,925])
    title('Density')
    legend('Location','southwest','FontSize',8,'Color','none')
    xlabel('Density (kg m^{-3})')
    ylabel('Depth (m)')
    set(gca,'FontSize',8)
    axis ij
    box on
    switch fs
        case 'FP'
            text(0.95,0.9,'c','Units','Normalized','FontSize',10)
        case 'ST'
            text(0.95,0.9,'d','Units','Normalized','FontSize',10) 
            %text(-.34,-.17,sprintf('— — ρ = 550 kg m^{-3} \n—·— ρ = 830 kg m^{-3}'),'Units','Normalized','FontSize',6)
            text(0.6,2.342,'{\it\rho} = 550 kg m^{-3}','Margin',0.1,'Background','white','Units','Normalized','FontSize',6,'Color',[0.6,0.6,0.6])
            text(0.6,2.122,'{\it\rho} = 830 kg m^{-3}','Margin',0.1,'Background','white','Units','Normalized','FontSize',6,'Color',[0.6,0.6,0.6])
    end  

end
 
h1 = annotation('textbox',[0.18 0.96 0.2 0.06],...
    'string','Fletcher Promontory');
h2 = annotation('textbox',[0.64 0.96 0.2 0.06],...
    'string','Skytrain Ice Rise');
set([h1 h2], 'fitboxtotext','on',...
    'edgecolor','none','FontSize',12)

end

function wFirn = firnVelo(z, w, ub, lb)

    I = z>ub & z < lb;

    strain_fit = polyfit(z(I) ,w(I),1);
    strain = polyval(strain_fit,z);
    wFirn = w - strain;
    
end