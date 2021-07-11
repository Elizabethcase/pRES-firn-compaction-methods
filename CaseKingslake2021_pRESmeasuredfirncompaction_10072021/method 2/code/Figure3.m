%generates Figure 3 from Case and Kingslake 2021


locations = {'FP','ST'};

for iloc = 1:length(locations)
    
    fs = locations{iloc};

    switch fs
        case 'FP'
            load FP_min_DxS_04_16_21.mat %loads "stake_min_DxS", from MinimizationMethod2.m
            load FP_pRES_DxS_04_16_21.mat; %loads "stake_pRES_DxS", from MinimizationMethod2.m
            load FP_method1Results_04162121.mat; %loads "pRES", from Figure2.m
            load FP_model_pres_0416212109.mat; %loads "modelOutput", from model output
            StakeOK = 1;

            core=load_FP_core;
            zRho = core.d;
            Rhod = core.rho;

            iplot=1;
            loctitle='Fletcher Promontory';

        case 'ST'

            load ST_min_DxS_04_05_21_site18.mat; %loads "stake_min_DxS"
            load ST_pRES_DxS_04_05_21_site18.mat; %loads "stake_pRES_DxS"
            load ST_method1Results_05182117.mat; %loads "pRES"
            load ST_model_pres_0415211128.mat; %loads "modelOutput"

            StakeOK = unique([stake_min_DxS.number]);
            StakeOK = 18;
            pRES = pRES([pRES.stake]==StakeOK);

            core=load_ST_core;
            zRho = core.d;
            Rhod = core.rho;

            iplot = [1,5,19,29];
            loctitle='Skytrain Ice Rise';

    end


    %% functions and global vars

    max_depth = 400;
    nI = 1.78;
    RhoI = 905;
    secpy = 365.25*24*3600;
    c = 2.983e8;
    f = 305e6;

    colorpRES = '#193764';%[0.2,0.2,0.2];
    color4Min = [.223,.349,.6];
    color3Min = [.023,.102,.259];
    colorCore = '#B1B1B1';

    cmapM3 = {'#2a9d8f','#a44a3f','#e59500','#264653'};

    ShouldBeZero=@(d,dI,L,RhoS) -dI+d+L*(nI-1)/nI*(1-RhoS)*(exp(-d/L)-1);
    TrueDepthfun=@(RhoS,L,dI) bisection(@(d) ShouldBeZero(d,dI,L,RhoS),0,max(dI)+100);

    %normalized by RhoI
    RhoFun=@(RhoS,L,dI) 1 + (RhoS-1).*exp(-TrueDepthfun(RhoS,L,dI)/L);
    nExpFun=@(RhoS,L,dI) 1/nI+(1-1/nI).*(1 + (RhoS-1)*exp(-TrueDepthfun(RhoS,L,dI)/L));
    nArrFun = @(rho) 1 + (nI -1)*rho/RhoI;

    %find depth of 550 and 830
    [value, i550] = min(abs(core.rho-550));
    [value, i830] = min(abs(core.rho-830));
    
    td550 = core.d(i550);
    td830 = core.d(i830);
    %% plot Fletcher or Skytrain single points

    cmap = brewermap(4,'Dark2');

    %cycle through DX options
    DX = {'pres','none'};

    for j=1:length(DX)

        dx = DX{j};

        switch dx
            case 'pres'
                stake = stake_pRES_DxS;

                for i = StakeOK

                %%%load data%%%
                    zetaTop = stake(i).zetaTop;
                    d = stake(i).d;
                    wFirnInvRhoInv = stake(i).wFirnInv_InvDensity;
                    Rho = stake(i).Rho;
                
               %%%plot velocities%%%
                     switch fs
                        case 'FP'
                            figure('Units','centimeters','Position',[10,10,17.9,10]);
                            subplot(221);
                            mn = -0.03;
                            mx = 0.47;
                            axy = [0 200];
                            axx = [-.05 0.5];
   
                        case 'ST'
                            subplot(222);
                            mn = -0.2;
                            mx = 0.8;
                            axy = [0 200];
                            axx = [-.05 0.5];
                     end
                    hold on
                    plot(linspace(mn,mx,100),td550*ones([1 100]),'--','Color','#C6C6C6','HandleVisibility','off') 
                    plot(linspace(mn,mx,100),td830*ones([1 100]),'-.','Color','#C6C6C6','HandleVisibility','off') 
                     
                    plot(wFirnInvRhoInv, d, 'ko',...
                        'Color', color3Min, 'DisplayName', 'Method 2','LineWidth',1.5)

                    axis ij
                    xlabel('Vertical velocity (m a^{-1})')
                    ylabel('Depth {\it\zeta} (m)')
                    ylim([0 200])
                    xlim(axx)
                    legend('Location','se','FontSize',8);
                    title('Compaction velocity')

               %%%plot densities%%%
               
                    switch fs
                        case 'FP'; subplot(223);
                        case 'ST'; subplot(224);
                    end
                    hold on;
                    plot(linspace(300,1220,100),td550*ones([1 100]),'--','Color','#C6C6C6','HandleVisibility','off') 
                    plot(linspace(300,1220,100),td830*ones([1 100]),'-.','Color','#C6C6C6','HandleVisibility','off') 
                    
                    
                    plot(Rhod,zRho,'s',...
                        'Color',colorCore,'DisplayName','Core density')
                    plot(Rho,d,'o',...
                        'Color',color3Min,'DisplayName','Method 2','LineWidth',1.5)

                    xlabel('Density \rho (kg m^{-3})')
                    ylabel('Depth {\it\zeta} (m)')
                    axis ij;
                    ylim([0 200])
                    legend('Location','sw','FontSize',8,'Color','none');
                    title('Density')
                    xlim([375,925])
                end

            case 'none'

                stake = stake_min_DxS;

                for i=StakeOK

                    %load data
                    zetaTop = stake(i).zetaTop;
                    d = stake(i).d;
                    wFirnInvRhoInv = stake(i).wFirnInv_InvDensity;
                    Rho = stake(i).Rho;

                    switch fs
                        case 'FP'
                            subplot(221)
                            text(.95,.95,'a','Units','Normalized','FontSize',10)
                        case 'ST'
                            subplot(222)
                    end
                    set(gca,'FontSize',8,'Box','on')
                    hold on;
                    
                    plot(wFirnInvRhoInv, d, '^',...
                        'Color', color4Min, 'DisplayName', 'Method 2 (4-param)','LineWidth',1.5)                
                    errorbar(pRES.wFirn, pRES.d, pRES.e,'horizontal','s',...
                        'Color', colorpRES, 'DisplayName', 'Method 1')                
                    plot(modelOutput.wFirn,modelOutput.d,'--','DisplayName','Model tuned to{\it w_c}','LineWidth',2,'Color','#6787C8')
                    
                    %densities
                    switch fs
                        case 'FP'
                            subplot(223)
                            text(.95,.95,'c','Units','Normalized','FontSize',10)

                        case 'ST'
                            subplot(224)

                    end
                    
                    
                    hold on;
                    plot(Rho,d,'^',...
                        'Color',color4Min,'DisplayName','Method 2 (4-param)','LineWidth',1.5)
                    
                    plot(modelOutput.rho,modelOutput.d,'--',...
                        'LineWidth',2, 'Color','#6787C8', 'DisplayName','Model tuned to{\it w_c}')
                    
                    set(gca,'FontSize',8,'Box','on')
                    xlim([375,925])
                end
        end

    end

    set(gca,'FontSize',8)
    switch fs
        case 'FP'
            subplot(221);
            
        case 'ST'
            subplot(222);
            text(.95,.95,'b','Units','Normalized','FontSize',10)
    end
    
    switch fs
        case 'FP'
            subplot(223);
            
        case 'ST'
            subplot(224);
            text(.95,.95,'d','Units','Normalized','FontSize',10)
    end  
    
    h1 = annotation('textbox',[0.18 0.96 0.2 0.06],...
        'string','Fletcher Promontory');
    h2 = annotation('textbox',[0.64 0.96 0.2 0.06],...
        'string','Skytrain Ice Rise');
    set([h1 h2], 'fitboxtotext','on',...
        'edgecolor','none','FontSize',12)

end
text(0.6,2.342,'{\it\rho} = 550 kg m^{-3}','Margin',0.1,'Background','white','Units','Normalized','FontSize',6,'Color',[0.6,0.6,0.6])
text(0.6,2.122,'{\it\rho} = 830 kg m^{-3}','Margin',0.1,'Background','white','Units','Normalized','FontSize',6,'Color',[0.6,0.6,0.6])
    
