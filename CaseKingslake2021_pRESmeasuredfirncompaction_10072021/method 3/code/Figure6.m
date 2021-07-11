%% Creates Figure 6 from Case and Kingslake 2021 from Method 3

% uses model output from ModelFigure6.m
% and includes functions for transltaing w, wI, d, dI into TWTT and dTWTT
% transforms w / wI / d into TWTT & dTWTT

% w = array of *ice-equivalent* total or firn compaction velocities
% rho = array of density values
%
% dataType = 'pRES' or 'model' depending on the data type
% 
        
%%function [TWTT, dTWTT] = WIntoTWTT(w,rho,z,loc,dataType)


clear all;
close all;

%% load data


locs = {'FP','ST','KF'};

for ilocs = 1:length(locs)
    loc = locs{ilocs};
    
    max_depth = 400;
    nI = 1.78;
    rhoI = 906;
    secpy = 365.25*24*3600;
    c = 2.983e8;
    f = 305e6;
    n = @(rho) 1 + (nI -1)*rho/rhoI;
    
    switch loc

        case 'FP'
            figure('Units','centimeters','Position',[10,10,17.9,20.4]);
            %subplot(311)
            dt = 353/365;
            
            tt = 'Fletcher Promontory';
            
            load('FP_method3_0416212257.mat')
            load('NXYwBHe_FP_pRES.mat')
            flipud('NXYwBHe_FP_pRES');
            
            pRES(1).stake = 1;
            pRES(1).dI = NXYwBHe_FP_pRES(:,3);
            pRES(1).wI = NXYwBHe_FP_pRES(:,4);
            pRES(1).eI = NXYwBHe_FP_pRES(:,7);
            pRES(1).H = NXYwBHe_FP_pRES(:,6);
            
            MO.x = unique(NXYwBHe_FP_pRES(:,2));
            stakeOK = [1];
            MO.stakeOK = stakeOK;
            
            clear NXYwBHe_FP_pRES

        case 'ST'
            %subplot(312)
            dt = 395/365;
            
            load('ST_method3_0415211425.mat')
                	
            load('NXYwBHe_STpRES_for_Carlos_110719.mat')
            stakeOK = unique(NXYwBHe_STpRES(:,1));
            tt = 'Skytrain Ice Rise';
            for i = 1:length(stakeOK)

                I = NXYwBHe_STpRES(:,1)==stakeOK(i);

                pRES(i).stake = stakeOK(i);

                dI = NXYwBHe_STpRES(I,3);
                wI = NXYwBHe_STpRES(I,4);
                eI = NXYwBHe_STpRES(I,7);
                H = NXYwBHe_STpRES(I,6);
                x = unique(NXYwBHe_STpRES(I,2));
                
                
                pRES(i).dI = flipud(dI(dI < max_depth));
                pRES(i).wI = flipud(wI(dI < max_depth));
                pRES(i).eI = flipud(eI(dI < max_depth));
                pRES(i).H = unique(H);
                MO(i).x = x;
                MO(i).stakeOK = stakeOK(i);
                

            end
            clear NXYwBHe_STpRES

        case 'KF'
            
            %subplot(313)
            dt = 350/365;
            tt = 'Korff Ice Rise';
            load('KF_method3_0415211301.mat') % load model data

            load('NXYwBHe_KorffpRES_for_Carlos_220515.mat') % load pres data
            stakeOK = unique(NXYwBHe_KorffpRES(:,1));
            %stakeOK = unique(NXYwBHe_KorffpRES(abs(NXYwBHe_KorffpRES(:,2)) > 500,1));


            for i = 1:length(stakeOK)

                I = NXYwBHe_KorffpRES(:,1)==stakeOK(i);
                

                pRES(i).stake = stakeOK(i);

                dI = NXYwBHe_KorffpRES(I,3);
                wI = NXYwBHe_KorffpRES(I,4);
                eI = NXYwBHe_KorffpRES(I,7);
                H = NXYwBHe_KorffpRES(I,6);
                x = unique(NXYwBHe_KorffpRES(I,2));

                pRES(i).dI = flipud(dI(dI < max_depth));
                pRES(i).wI = flipud(wI(dI < max_depth));
                pRES(i).eI = flipud(eI(dI < max_depth));
                pRES(i).H = unique(H);
  
                MO(i).x = x;
                MO(i).stakeOK = stakeOK(i);

            end
           

           clear NXYwBHe_KorffpRES
                    
    end

% make figure 3

cmap = {'#2a9d8f','#a44a3f','#e59500','#264653'};
nCol = 1;
nCol1 = 1;
for i = 1:length(MO)

    disp([num2str(i),'***',num2str(MO(i).stakeOK)])
    d = MO(i).d;
    dI = MO(i).dI;
    wFirnI = MO(i).wFirnI;
    rho = MO(i).rho;
    
    
    % add ice flow strain back in
    strain_fit = strainCalc(loc,pRES(i).dI, pRES(i).wI);
    strainI = polyval(strain_fit,dI);
    wIModel = wFirnI + strainI; %adds horizontal strain rate-induced verticla velocity back into model


    % convert to two way travel time
    TModel = dItoT(dI);
    dTModel = wItodT(wIModel,dt);

    Tpres = dItoT(pRES(i).dI);
    dTpres = wItodT(pRES(i).wI,dt);
    eTpres = wItodT(pRES(i).eI,dt);
    
    
    switch loc
        case 'FP'
            subplot(321)
            iplot = 1;
        case 'ST'
            subplot(323)
            iplot = [1,5,19,29];
        case 'KF'
            subplot(325)
            iplot = [1,8,22,30];
    end
    
    if ismember(i,iplot)
        color = cmap{nCol};
        plot(dTModel*1e9,TModel*1e6, 'DisplayName','model','Color',color,'DisplayName','Model','LineWidth',2)
        hold on;

        errorbar(dTpres*1e9,Tpres*1e6,eTpres,'horizontal', 'o','DisplayName','pres','Color',color,'DisplayName','pRES','LineWidth',2)
        hold on;
        axis ij;
        title('Modelled and Measured Travel Time')
        xlabel('Change in two-way travel time,{\it \DeltaT} (ns)')
        ylabel('Two-way travel time,{\it T} (ms)')
        ylim([0 max(TModel*1e6)])
        
        nCol = nCol + 1;
        
        set(gca,'FontSize',8,'Box','on')
    end
    
    
    if i==1
    switch loc
        case 'FP'
            legend('Position',[0.467,0.0267,0.1,0.05],'Units','normalized','FontSize',8)
            text(0.95,0.05,'a','Units','normalized','FontSize',10)
        case 'ST'
            text(0.95,0.05,'c','Units','normalized','FontSize',10)
        case 'KF'
            text(0.95,0.05,'e','Units','normalized','FontSize',10)
    end
    end
    
    switch loc
        case 'FP'
            subplot(322)
        case 'ST'
            subplot(324)
        case 'KF'
            subplot(326)
    end
    
    dTModelInterp = interp1(TModel,dTModel,Tpres);
    diff(i).dpm = dTpres-dTModelInterp;
    diff(i).ppm = (diff(i).dpm)./dTpres;


    plot(diff(i).ppm*100,Tpres*1e6,'s','Color','#a5a58d')
    xlabel('Percent (%)')
    ylabel('Two-way Travel Time, {\it T} (ms)')
    xlim([-20 20])
    axis ij
    hold on;
    title('Model & observed difference in{\it \DeltaT}')
    
    if ismember(i,iplot)
        color = cmap{nCol1};
        plot(diff(i).ppm*100,Tpres*1e6,'s','Color',color,'MarkerFaceColor',color)
        nCol1 = nCol1 + 1;
    end


    

end
    switch loc
        case 'FP'
            text(0.95,0.05,'b','Units','normalized','FontSize',10)
            text(-.47,1.15,'Fletcher Promontory','Units','normalized','FontSize',12)
        case 'ST'
            text(0.95,0.05,'d','Units','normalized','FontSize',10)
            text(-.43,1.15,'Skytrain Ice Rise','Units','normalized','FontSize',12)

        case 'KF'
            text(0.95,0.05,'f','Units','normalized','FontSize',10)
            text(-.38,1.15,'Korff Ice Rise','Units','normalized','FontSize',12)
    
    end

set(gca,'FontSize',8,'Box','on')
end



%% vary ub and lb, find std dev of firn compaction in pres, plot
% only works for FP right now

% n = 50;
% 
% stakeOK = [1,10,19,30];
% 
% for j = stakeOK
%     
%     figure();
%     
%     dIpres = pRES(j).dI;
%     wIpres = pRES(j).wI;
%     dIModel = MO(j).dI;
%     wFirnIModel = MO(j).wFirnI;
%     
%     for i = 1:n
% 
% 
%         ub = 0.33 + .1*(2*rand(1)-1);
%         lb = 0.66 + .1*(2*rand(1)-1);
%         strain = strainCalcRandom(dIpres,wIpres,ub,lb);
%         randBounds(i).ub = ub;
%         randBounds(i).lb = lb;
%         %randBounds(i).wFirnI = pRES.wI - strain;
%         wFirnIRand(:,i) = wIpres - strain;
% 
% %         plot(strain,dIpres);
% %         hold on;
% % 
% %         drawnow;
%     end
% 
%     for k = 1:size(wFirnIRand,1)
%         
%         [mi midx(k)] = min(wFirnIRand(k,:));
%         [ma madx(k)] = max(wFirnIRand(k,:));
% 
%         maxWFirnIRand(k) = ma;
%         minWFirnIRand(k) = mi;
%     end
%     
%         % add ice flow strain back in
% 
%     strainIi = strainCalc(loc,dIpres,wIpres);
%     strainI = interp1(dIpres, strainIi, dI);
%     wIModel = wFirnIModel + strainI;
% 
%     % convert to two way travel time
%     TModel = dItoT(dIModel);
%     dTModel = wItodT(wIModel,dt);
% 
%     Tpres = dItoT(dIpres);
%     dTpres = wItodT(wIpres,dt);
% 
%     plot(maxWFirnIRand,dIpres)
%     hold on;
%     plot(minWFirnIRand,dIpres)
%     x2 = [dIpres', fliplr(dIpres')];
%     inBetween = [maxWFirnIRand, fliplr(minWFirnIRand)];
%     fill(inBetween, x2, [0.3020, 0.7451, 0.9333],'FaceAlpha',0.5,'DisplayName','diff line fits');
%     plot(wFirnIModel,dIModel,'DisplayName','model','Color',[0.8510,0.3255,0.0980],'LineWidth',2)
%     axis ij
%     ylim([0 250])
%     legend('Location','se')
%     title([loc 'at stake ' num2str(j)])
% 
%     
%     clear wFirnIRand maxWFirnIRand minWFirnIRand
% end
% 


%%

function strain = strainCalcRandom(d, w, ub, lb)
    
    I = d > ub & d < lb;

    strain_fit = polyfit(d(I) ,w(I),1);
    strain = polyval(strain_fit,d);

end

function strain_fit = strainCalc(loc, di, wi)
    
    switch loc
        case 'FP'
            ub = 200;
            lb = 300;
        case 'ST'
            ub = 130;
            lb = 280;
        case 'KF'
            ub = 175;
            lb = 300;
    end
    
    I = di > ub & di < lb;
    strain_fit = polyfit(di(I) ,wi(I), 1);

end

function T = dItoT(dI)
    c=3e8;
    nI = 1.78;
    T = 2*nI*dI/c;
end

function dT = wItodT(wI,dt)
    
    nI = 1.78;
    c = 2.98e8;
    dT = wI.*dt.*nI/c; %(m/yr * yr * s/m) = T in units of s
    
end