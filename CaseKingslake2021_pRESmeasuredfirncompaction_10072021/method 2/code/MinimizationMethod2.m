%% Performs minimization at each of the sites


%% House keeping

clear all
%close all
%% Location and choice of DxS
FSopt = {'FP'};%options: 'ST','FP','KF'
DXopt = {'pRES','none'}; %options: {'pRES','none','GPSavg','GPSind'}; 
    % 'pRES sets DxS value from pRES obs (see Section 2.3.2), 'none' includes DxS in
    %       optimization (see discussion)
    % note that 'pRES and 'none' can be run at the same time as {'pRES','none'}; GPS options must be run separately (i.e. {'GPSavg'})
    
for iFS = 1:length(FSopt)
    for iDX = 1:length(DXopt)
        
clearvars -except FSopt DXopt iFS iDX

FS = FSopt{iFS}; %FP, KF, ST;
DX = DXopt{iDX}; %pRES; %none; %GPSavg (only available for 'KF'); %GPSind (only available for 'KF', very noisy beware)

RhoI=906; %kg/m3
c=299792458; %m/s speed of light
nI=1.78; %Refractive index of Ice

%% Define functions

%non normalized
% RhoExpFun=@(RhoS,L,d) RhoI - (RhoI-RhoS).*exp(-d/L);
% nFun=@(RhoS,L,d) 1+(nI-1).*RhoExpFun(RhoS,L,d);


ShouldBeZero=@(d,dI,L,RhoS) -dI+d+L*(nI-1)/nI*(1-RhoS)*(exp(-d/L)-1);
TrueDepthfun=@(RhoS,L,dI) bisection(@(d) ShouldBeZero(d,dI,L,RhoS),0,max(dI)+100);
nArrFun = @(rho) 1 + (nI -1)*rho/RhoI;

%normalized by RhoI
RhoFun=@(RhoS,L,dI) 1 + (RhoS-1).*exp(-TrueDepthfun(RhoS,L,dI)/L);
%normalized by nI
nFun=@(RhoS,L,dI) 1/nI+(1-1/nI).*(1 + (RhoS-1)*exp(-TrueDepthfun(RhoS,L,dI)/L));

%RhoTrue=@(RhoS,L,d) 1 + (RhoS-1).*exp(-d/L);
%% Load pRES & GPS data

switch FS
    case 'KF'
        load('NXYwBHe_KorffpRES_for_Carlos_220515.mat');
        pRES_data= NXYwBHe_KorffpRES;
        clear NXYwBHe_KorffpRES
        label0=pRES_data(:,1);
        StakeOK=[1:6,8:32,34:max(label0)]';
        
        %Load GPS
        load('Xu_KorffpRES_for_Carlos.mat');
        xGPS=Xu_KorffpRES(:,1);
        usGPS=Xu_KorffpRES(:,2);
        uerrGPS = Xu_KorffpRES(:,3);
        dudx = diff(usGPS)./diff(xGPS);
        ix = uerrGPS<.1;
        dudx_avg = mean(diff(usGPS(ix))./diff(xGPS(ix))); %uses only locations where GPS error is low
        clear Xu_KorffpRES ix;
        loctitle = 'Korff Ice Rise';
        
        depthComp=150; %depth where compaction has no effect
        depthDxUniform=300; %depth where we consider Dx to be still uniform

    case 'FP'
        load NXYwBHe_FP_pRES
        pRES_data = NXYwBHe_FP_pRES;
        clear NXYwBHe_FP_pRES
        StakeOK = 1;

        core=load_FP_core;
        zRho = core.d;
        Rhod = core.rho;
        
        loctitle='Fletcher Promontory';
        
        load('FP_method1Results_04162123.mat');
        method1 = pRES;
        clear pRES
        
        depthComp=160; %depth where compaction has no effect
        depthDxUniform=309; %depth where we consider Dx to be still uniform

    case 'ST'
        load NXYwBHe_STpRES_for_Carlos_110719
        pRES_data = NXYwBHe_STpRES;
        StakeOK = unique(NXYwBHe_STpRES(:,1));
        %StakeOK = [17,18,19];
        %StakeOK = [18];
        clear NXYwBHe_STpRES
        
        core=load_ST_core;
        zRho = core.d;
        Rhod = core.rho;
        
        loctitle='Skytrain Ice Rise';
        
        load('ST_method1Results_04152112.mat');
        method1 = pRES;
        clear pRES
        
        depthComp=122; %depth where compaction has no effect
        depthDxUniform=304; %depth where we consider Dx to be still uniform
       
end

% Read data
label0=pRES_data(:,1);
x0=pRES_data(:,2);
depth0=pRES_data(:,3);
w0=pRES_data(:,4);
H0=pRES_data(:,5);
hs0=pRES_data(:,6);
errw0=pRES_data(:,7);
clear pRES_data;

% Go stake by stake, read data into Stake structure
for i=StakeOK'
    
    xi=mean(x0(label0==i));
    depthi=depth0(label0==i);
    hsi=mean(hs0(label0==i));
    Hi=mean(H0(label0==i));
    zi=hsi-depthi;
    wi=w0(label0==i);
    errwi = errw0(label0==i);
    
    stake(i).number=i;
    stake(i).x=xi;    
    stake(i).depthI=depthi(end:-1:1);
    stake(i).zI=zi(end:-1:1);
    stake(i).hs=hsi;
    stake(i).H=Hi;
    stake(i).wI=wi(end:-1:1);
    stake(i).errw=errwi(end:-1:1);
    stake(i).loc = FS;
    
 
    % Get us from GPS
    switch FS
        case 'KF'
            stake(i).us=interp1(xGPS,usGPS,stake(i).x);
            stake(i).uerr = uerrGPS(i);
            switch DX
                case 'GPSind'
                    stake(i).DxS = dudx(i);
                case 'GPSavg'
                    stake(i).DxS = dudx_avg;
            end
            
    end
end

clear x0 depth0 hs0 H0 w0 dudx;

%

%% Using different estimates (or no estimate) for DxS, do inversion
%figure();
switch DX
    case 'pRES'
        
        for i=StakeOK'
            depthi=stake(i).depthI; %depth measured by pres
            wi=stake(i).wI; %i.e. vertical velocity measured by pres
            zetaComp=(depthi>depthComp&depthi<depthDxUniform); % below compaction, above bed effects
            zetaTop=depthi<depthDxUniform; % above bed effects


        % First we estimate Dx from solid ice 
        % sum( (w_pRES - (w_intercept + ex * depth))^2 )
        % v(1) = w_intercept; v(2) = strain rate
            GetStrainRate=@(v) sum((wi(zetaComp)...
                - (v(1)+v(2)*depthi(zetaComp))).^2);

            [v1,~,~]= fminsearch(GetStrainRate,[0.4 -9e-4]);
            
            stake(i).DxS=v1(2);
            DxS=v1(2);
            intDxS = v1(1);
    
        %Now we go for the rest

            MinimizeForRhoSWsL=@(v) trapz(TrueDepthfun(v(3),v(2),depthi(zetaTop)),...
                (wi(zetaTop).*RhoFun(v(3),v(2),depthi(zetaTop))./nFun(v(3),v(2),depthi(zetaTop))...
                - (v(1)*v(3)/nFun(v(3),v(2),0.0)...
                +DxS*(...
                TrueDepthfun(v(3),v(2),depthi(zetaTop))...
                +(1-v(3))*v(2).*(exp(-TrueDepthfun(v(3),v(2),depthi(zetaTop))/v(2))-1)))...
                ).^2);

            [v,fval,exitflag,~]= fminsearch(MinimizeForRhoSWsL,[0.2 90 0.3]);%,options);

            fprintf('wsi=%g Dx=%g LRho=%g RhoS=%g Stake=%g\n',v(1),DxS,v(2),v(3)*RhoI,i)

       %get true and ice equivalent firn velocities
       
            %define vars from inversion
            wsi = v(1);
            L = v(2);
            RhoS_norm = v(3); %normalized by RhoI
            RhoS = v(3)*RhoI;
            Rho_norm = RhoFun(v(3),v(2),depthi); %normalized by RhoI
            Rho=RhoI*RhoFun(v(3),v(2),depthi);
            trued = TrueDepthfun(v(3),v(2),depthi);
            ns = nArrFun(RhoS);
            ws = wsi.*nI./ns;
            
            
            %Equation 15 in paper (w ice equivalent)
             wInv = (wsi*RhoS_norm/(nFun(RhoS_norm,L,0.0))...
                    +DxS.*(...
                    trued...
                    +(1-RhoS_norm)*L.*(exp(-trued/L)-1)))./Rho_norm.*nFun(RhoS_norm,L,depthi);
                
            
            %equation 13 in firn paper (w true)
            wTrue = RhoS.*ws./Rho + (DxS./Rho).*...
                    (RhoI .* trued-...
                    L.*(RhoI - RhoS).*...
                    (1-exp(-trued./L)));
           
            GetStrainRate=@(v) sum((wTrue(zetaComp)...
                - (v(1)+v(2)*trued(zetaComp))).^2);

            [v1,~,~]= fminsearch(GetStrainRate,[0.4 -9e-4]);
            
            %get true firn velocity
            wInvDueToIceFlow = v1(1)+v1(2)*trued;
            wFirnInv_InvDensity = wTrue - wInvDueToIceFlow;
             

        % save data in struct

            stake(i).wFirnInv_InvDensity = wFirnInv_InvDensity;            
            stake(i).wInvTrue = wTrue;
            stake(i).wInv = wInv;
            stake(i).ws=v(1);
            stake(i).RhoI=RhoI;
            stake(i).RhoS=v(3)*RhoI;
            stake(i).L=v(2);
            stake(i).Rho=Rho;
            stake(i).depthComp = depthComp;
            stake(i).depthDxUniform = depthDxUniform;
            stake(i).zetaTop = zetaTop;
            stake(i).d = trued;
            stake(i).TypeOfInversion = 'pRES-derived DxS held constant';
            stake(i).FileOrigin = 'Minimizationmethod2.m';
            stake(i).exitflag = exitflag;
            stake(i).dateCreated = date;
            stake(i).funcVal = fval;
            stake_pRES_DxS = stake;
            
            figure;
            subplot(131)
            plot(wi,depthi,'ko','DisplayName','pRES-observed')
            hold on;
            plot((v1(1)+v1(2)*depthi)...
                 ,depthi,'k-','DisplayName','w from DxS');
    
            plot(wTrue,trued,'*b','DisplayName','True w (pRES DxS)');
            plot(wInv,depthi,'*r','DisplayName','I.e. w (pRES DxS)');
            plot(wInv./nFun(v(3),v(2),depthi),depthi,'*m','DisplayName','I.e. w (pRES DxS)');
            plot(v1(1)+v1(2)*trued, trued,'b-','LineWidth',2); 
            xlabel('Vertical velocity (m yr^{-1}');
            ylabel('Elevation (m)');  
            drawnow;
            legend()
            
            subplot(132)
            hold on;
            plot(wFirnInv_InvDensity,trued,'b-')
            %plot(method1.wFirn,stake(i).hs-method1.d,'mo')
            
            subplot(133)
            if ~strcmp(FS,'KF')
                plot(core.rho,core.d,'o','DisplayName','core','Color',[0.5,0.5,0.5])
            end
            hold on;
            plot(RhoI*RhoFun(v(3),v(2),depthi(zetaTop)),trued(zetaTop),'ob','LineWidth',2,'DisplayName','pRES-derived min');
    
        end

    case 'none'
        
        for i=StakeOK'
            
            depthi=stake(i).depthI;
            wi=stake(i).wI;
            zetaComp=(depthi>depthComp&depthi<depthDxUniform);
            zetaTop=depthi<depthDxUniform;
            
            MinimizeForRhoSWsLDxS=@(v) trapz(TrueDepthfun(v(3),v(2),depthi(zetaTop)),...
            (wi(zetaTop).*RhoFun(v(3),v(2),depthi(zetaTop))./nFun(v(3),v(2),depthi(zetaTop))...
            - (v(1)*v(3)/nFun(v(3),v(2),0.0)...
            +v(4)*(...
            TrueDepthfun(v(3),v(2),depthi(zetaTop))...
            +(1-v(3))*v(2).*(exp(-TrueDepthfun(v(3),v(2),depthi(zetaTop))/v(2))-1)))...
            ).^2);

            [v,fval,exitflag,~]= fminsearch(MinimizeForRhoSWsLDxS,[0.6 40 0.4 -9e-4]);%,options);

            fprintf('ws=%g Dx=%g LRho=%g RhoS=%g Stake=%g\n',v(1),v(4),v(2),v(3)*RhoI,i)

       %get true and ice equivalent firn velocities
       
            %define vars from inversion
            wsi = v(1);
            L = v(2);
            RhoS_norm = v(3); %normalized by RhoI
            RhoS = v(3)*RhoI;
            Rho_norm = RhoFun(v(3),v(2),depthi); %normalized by RhoI
            Rho=RhoI*RhoFun(v(3),v(2),depthi);
            trued = TrueDepthfun(v(3),v(2),depthi);
            ns = nArrFun(RhoS);
            ws = wsi.*nI./ns;
            DxS = v(4);
            
            
            %Equation 15 in paper (w ice equivalent)
             wInv = (wsi*RhoS_norm/(nFun(RhoS_norm,L,0.0))...
                    +DxS.*(...
                    trued...
                    +(1-RhoS_norm)*L.*(exp(-trued/L)-1)))./Rho_norm.*nFun(RhoS_norm,L,depthi);
                
            
            %equation 13 in firn paper (w true)
            wTrue = RhoS.*ws./Rho + (DxS./Rho).*...
                    (RhoI .* trued-...
                    L.*(RhoI - RhoS).*...
                    (1-exp(-trued./L)));
           
            
            GetStrainRate=@(v) sum((wTrue(zetaComp)...
                - (v(1)+v(2)*trued(zetaComp))).^2);

            [v1,~,~]= fminsearch(GetStrainRate,[0.4 -9e-4]);
            
            %get true firn velocity
            wInvDueToIceFlow = v1(1)+v1(2)*trued;
            wFirnInv_InvDensity = wTrue - wInvDueToIceFlow;
            
            %save data in struct
            stake(i).wFirnInv_InvDensity = wFirnInv_InvDensity; 
            stake(i).wInvTrue = wTrue;
            stake(i).wInv = wInv;
            stake(i).ws=v(1);
            stake(i).RhoI=RhoI;
            stake(i).RhoS=v(3)*RhoI;
            stake(i).L=v(2);
            stake(i).Rho=Rho;
            stake(i).DxS = v(4);
            stake(i).zetaTop = zetaTop;
            stake(i).d = trued;
            stake(i).TypeOfInversion = 'minimized for all four variables (wS, rhoS, DxS, L)';
            stake(i).FileOrigin = 'Minimizationmethod2.m';
            stake(i).dateCreated = date;
            stake(i).exitflag = exitflag;
            stake(i).funcVal = fval;
            stake_min_DxS = stake;
            
%             figure 
            subplot(131)
            %hold on;
            %plot(wi,stake(i).hs-depthi,'ko')
            hold on;
            plot(wTrue,trued,'rs','DisplayName','Full min');
            plot(wInv*nI./nArrFun(Rho),trued,'g^');
            plot(wInvDueToIceFlow,trued,'g-');
            drawnow;
            %hold off;
            xlabel('Vertical velocity (m yr^{-1}');
            ylabel('Elevation (m)');
            legend();
            title('Total velocity (m/yr)')
            axis ij
            
            subplot(132)
            
            hold on;
            plot(wFirnInv_InvDensity,trued,'g-','DisplayName','true full-min w')
            title('Firn velocity')
            xlabel('velocity (m/yr or mie/yr)')
            axis ij
            
            subplot(133)
            hold on;
            plot(RhoI*RhoFun(v(3),v(2),depthi),trued,'-g','DisplayName','Full min');
            legend()
            title('Density')
            xlabel('density (kg m^{-3})')
            axis ij            
              
        end

    case 'GPSavg'
        if strcmp(FS,'KF')
            for i=StakeOK'
                depthi=stake(i).depthI;
                wi=stake(i).wI;
                zetaComp=(depthi>depthComp&depthi<depthDxUniform);
                zetaTop=depthi<depthDxUniform;   
                
                DxS = stake(i).DxS;

                erri=@(v) trapz(TrueDepthfun(v(3),v(2),depthi(zetaTop)),...
                (wi(zetaTop).*RhoFun(v(3),v(2),depthi(zetaTop))./nFun(v(3),v(2),depthi(zetaTop))...
                - (v(1)*v(3)/nFun(v(3),v(2),0.0)...
                +DxS*(...
                TrueDepthfun(v(3),v(2),depthi(zetaTop))...
                +(1-v(3))*v(2).*(exp(-TrueDepthfun(v(3),v(2),depthi(zetaTop))/v(2))-1)))...
                ).^2);

                [v,fval,exitflag,~]= fminsearch(erri,[0.6 40 0.4]);

                fprintf('ws=%g Dx=%g LRho=%g RhoS=%g\n ',v(1),DxS,v(2),v(3)*RhoI)

                wInv = (v(1)*v(3)/nFun(v(3),v(2),0.0)...
                        +DxS*(...
                        TrueDepthfun(v(3),v(2),depthi(zetaTop))...
                        +(1-v(3))*v(2).*(exp(-TrueDepthfun(v(3),v(2),depthi(zetaTop))/v(2))-1)))./RhoFun(v(3),v(2),depthi(zetaTop)).*nFun(v(3),v(2),depthi(zetaTop));

%                 subplot(121)
%                 hold on;
%                 plot(wInv,stake(i).hs-depthi(zetaTop),'-m','DisplayName','GPSavg-derived DxS min');
%                 drawnow;
%                 hold off;
%                 xlabel('Vertical velocity (m yr^{-1}');
%                 ylabel('Elevation (m)');
%                 
%                 subplot(122)
%                 plot(RhoI*RhoFun(v(3),v(2),depthi(zetaTop)),stake(i).hs-depthi(zetaTop),'-r','DisplayName','GPSavg-derived DxS min');
%                 xlim([100 1000])
%                 legend()
                
                %get firn velocities
            
                Rho=RhoI*RhoFun(v(3),v(2),depthi);

                GetStrainRate=@(v) sum((wi(zetaComp)...
                    - (v(1)+v(2)*depthi(zetaComp))).^2);

                [v1,~,~]= fminsearch(GetStrainRate,[0.4 -9e-4]);

                wDueToIceFlow = v1(1)+v1(2)*depthi;
                wFirnIPres = wi - wDueToIceFlow;
                wFirnPres_InvDensity = wFirnIPres .* nI./nArrFun(Rho);

                GetStrainRate=@(v) sum((wInv(zetaComp)...
                                - (v(1)+v(2)*depthi(zetaComp))).^2);

                [v1,~,~]= fminsearch(GetStrainRate,[0.4 -9e-4]);

                wDueToIceFlow = v1(1)+v1(2)*depthi(zetaTop);
                wFirnIInv = wInv - wDueToIceFlow;
                wFirnInv_InvDensity = wFirnIInv .*nI./nArrFun(Rho(zetaTop));

                stake(i).wFirnIPres = wFirnIPres;
                stake(i).wFirnIInv = wFirnIInv;
                stake(i).wFirnPres_InvDensity = wFirnPres_InvDensity;
                stake(i).wFirnInv_InvDensity = wFirnInv_InvDensity; 
            
                stake(i).wInv = wInv;
                stake(i).ws=v(1);
                stake(i).RhoI=RhoI;
                stake(i).RhoS=v(3)*RhoI;
                stake(i).L=v(2);
                stake(i).Rho=Rho;
                stake(i).zetaTop = zetaTop;
                stake(i).TypeOfInversion = 'DxS equal to avg of GPS (where err < 0.1)';
                stake(i).FileOrigin = 'Minimizationmethod2.m';
                stake(i).dateCreated = date;
                stake(i).exitflag=exitflag;
                stake(i).funcVal = fval;
                stake_GPSavg_DxS = stake;
                
            end
        else
            disp('GPS data only available for Korff')
        end
        
    case 'GPSind'
        if strcmp(FS,'KF')
            for i=StakeOK'
                
                if stake(i).uerr < 0.1
                    
                    depthi=stake(i).depthI;
                    wi=stake(i).wI;
                    zetaComp=(depthi>depthComp&depthi<depthDxUniform);
                    zetaTop=depthi<depthDxUniform;
                    
                    DxS = stake(i).DxS;

                    erri=@(v) trapz(TrueDepthfun(v(3),v(2),depthi(zetaTop)),...
                    (wi(zetaTop).*RhoFun(v(3),v(2),depthi(zetaTop))./nFun(v(3),v(2),depthi(zetaTop))...
                    - (v(1)*v(3)/nFun(v(3),v(2),0.0)...
                    +DxS*(...
                    TrueDepthfun(v(3),v(2),depthi(zetaTop))...
                    +(1-v(3))*v(2).*(exp(-TrueDepthfun(v(3),v(2),depthi(zetaTop))/v(2))-1)))...
                    ).^2);

                    [v,fval,exitflag,~]= fminsearch(erri,[0.6 40 0.4]);

                    fprintf('ws=%g Dx=%g LRho=%g RhoS=%g\n ',v(1),DxS,v(2),v(3)*RhoI)

                    wInv = (v(1)*v(3)/nFun(v(3),v(2),0.0)...
                            +DxS*(...
                            TrueDepthfun(v(3),v(2),depthi(zetaTop))...
                            +(1-v(3))*v(2).*(exp(-TrueDepthfun(v(3),v(2),depthi(zetaTop))/v(2))-1)))./RhoFun(v(3),v(2),depthi(zetaTop)).*nFun(v(3),v(2),depthi(zetaTop));

%                     subplot(121)
%                     hold on;
%                     plot(wInv,stake(i).hs-depthi(zetaTop),'-g','DisplayName','GPSind-derived DxS min');
%                     drawnow;
%                     hold off;
%                     xlabel('Vertical velocity (m yr^{-1}');
%                     ylabel('Elevation (m)');
%                     
%                     subplot(122)
%                     hold on;
%                     plot(RhoI*RhoFun(v(3),v(2),depthi(zetaTop)),stake(i).hs-depthi(zetaTop),'-r','DisplayName','GPSind-derived DxS min');
%                     legend()
                    
                    %get firn velocities

                    Rho=RhoI*RhoFun(v(3),v(2),depthi);

                    GetStrainRate=@(v) sum((wi(zetaComp)...
                        - (v(1)+v(2)*depthi(zetaComp))).^2);

                    [v1,~,~]= fminsearch(GetStrainRate,[0.4 -9e-4]);

                    wDueToIceFlow = v1(1)+v1(2)*depthi;
                    wFirnIPres = wi - wDueToIceFlow;
                    wFirnPres_InvDensity = wFirnIPres .* nI./nArrFun(Rho);

                    GetStrainRate=@(v) sum((wInv(zetaComp)...
                                    - (v(1)+v(2)*depthi(zetaComp))).^2);

                    [v1,~,~]= fminsearch(GetStrainRate,[0.4 -9e-4]);

                    wDueToIceFlow = v1(1)+v1(2)*depthi(zetaTop);
                    wFirnIInv = wInv - wDueToIceFlow;
                    wFirnInv_InvDensity = wFirnIInv .*nI./nArrFun(Rho(zetaTop));

                    stake(i).wFirnIPres = wFirnIPres;
                    stake(i).wFirnIInv = wFirnIInv;
                    stake(i).wFirnPres_InvDensity = wFirnPres_InvDensity;
                    stake(i).wFirnInv_InvDensity = wFirnInv_InvDensity;  
                    stake(i).wInv = wInv;
                    stake(i).ws=v(1);
                    stake(i).RhoI=RhoI;
                    stake(i).RhoS=v(3)*RhoI;
                    stake(i).L=v(2);
                    stake(i).Rho=Rho;
                    stake(i).zetaTop = zetaTop;
                    stake(i).TypeOfInversion = 'DxS derived from individual measurements from GPS (err < 0.1)';
                    
                else
                    stake(i).wInv = wInv;
                    stake(i).ws=NaN;
                    stake(i).RhoI=RhoI;
                    stake(i).RhoS=NaN;
                    stake(i).L=NaN;
                    stake(i).Rho=NaN;
                    stake(i).zetaTop = zetaTop;
                    stake(i).TypeOfInversion = 'GPS err > 0.1';
                    
                end
                stake(i).FileOrigin = 'Minimizationmethod2.m';
                stake(i).dateCreated = date;
                stake(i).exitflag = exitflag;
                stake(i).funcVal = fval;
                stake_GPSind_DxS = stake;
                
                
                
            end
        else
            disp('GPS data only available for Korff')
        end

    sgtitle(FS)

end


% save output
today = datestr(now,'mm_dd_yy');
% 
switch DX
    case 'pRES'
        save(['../minimization_output/' FS '_pRES_DxS_' today],'stake_pRES_DxS')
    case 'none'
        save(['../minimization_output/' FS '_min_DxS_' today],'stake_min_DxS')
    case 'GPSavg'
        save(['../minimization_output/' FS '_GPSavg_DxS_' today],'stake_GPSavg_DxS')
    case 'GPSind'
        save(['../minimization_output/' FS '_GPSinvert_DxS_' today],'stake_GPSind_DxS')

end

    end
end

