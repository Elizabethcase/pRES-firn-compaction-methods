%% Physics-based firn compaction model based on Arthern and others 2010

%loc = 'ST', 'KF', 'FP' (skytrain, korff, fletcher) -- depending on loc, loads different input variables

%opt_type = 'density' or 'pres' -- optmized to pRES or density profile

%plot_opt = 1 (plot during model run) or 0 (no plot)

%s = 1 (save output) or 0 (don't save)

function modelOutput = Model_041421(loc,opt_type,plot_opt,s,bI,stk)

%define time
    secpy = 3600*24*365.25; % seconds per year
    t_years = 10000.0;
    t_sec = t_years*secpy;
    dt = secpy*.001;
    t = 0:dt:t_sec;
    n = 10000.0;
    
%apres constants
    n_i = 1.78; % index of refraction of ice
    
%model constants
    kg = 1.3e-7 ;%growth constant, m^2 s^-1, default = 1.3e-7
    kc_low = 9.2 * 10^-9 ;%low creep coefficient / compaction (?) kg^-1 s-1 m^3
    kc_high = 3.7 * 10^-9;
    kappa_i = 2.1; %W * m-1 * K-1
    R = 8.31; %gas constant % kg m2 s-2
    g = 9.81; %m s^-2 / m ^ yr 
    rho_i = 906.0; %average density of skytrain core
    c_p = 2009.0; % J kg?1 K?1 %specific heat

    switch loc
        case 'FP'
        
            switch opt_type
                case 'density' %optimized to density
                    Ec = 6.7123e4; %6.6512e4;
                    Eg = 4.8579e4; %4.7974e4;

                case 'pres' %optimized to densification
                    Ec = 6.8721e4; %6.4275e4; %5.8045e4;
                    Eg = 5.0304e4; %4.5769e4; %3.9562e4;
            end
            
            rho_s = 435;
            Tav = 273.15 - 27.1; % fletcher (27.1) from mulvaney 2014
            if ~exist('bI','var')
                b = 280/rho_i/secpy; %ice-equivalent m/s, from Mulvaney personal communications april 2021
            else
                b = bI/secpy;
            end
            
            if ~exist('stk','var')
                stakeOK = 1;
            else
                stakeOK = stk;
            end
            core = load_FP_core;
            da = char(datetime(),'MMddyyHHmm');
            FN = ['../output/FP_model_' [opt_type] '_' da];
            % load radar data
%            load('NXYwBHe_FP_pRES.mat');
%             pRES = flipud(NXYwBHe_FP_pRES);
% 
%             wIE = pRES(:,4);
%             depthIE = pRES(:,3);
%             max_depth = 250;
%             
%             w_cIE_pRES = wFirn_from_wTot('FP', wIE, depthIE);
%             w_cIE_pRES = w_cIE_pRES(depthIE<max_depth);
%             depth_pRES = depthIE(depthIE<max_depth);
%             clear pRES
              load('FP_method1Results_04162121.mat');
            
        case 'ST'
            
            switch opt_type
                case 'density' %optimized to density
                    Ec = 6.6502e4;%5.7572e4;
                    Eg = 4.8506e4;%4.0551e4;

                case 'pres' %optimized to densification
                    Ec = 6.6838e4;%5.8150e4;
                    Eg = 4.9151e4;%4.0293e4;
            end
            
            rho_s = 400;
            Tav = 273.15 - 26.0; %26.0 as reported by Mulvaney et al 2020 (WACSWAIN)
            
            if ~exist('bI','var')
                b = 140/rho_i/secpy; % ice-equivalent m/s
            else
                b = bI/secpy;
            end

            if ~exist('stk','var')
                stakeOK = 18;
            else
                stakeOK = stk;
            end
            core = load_ST_core();
            da = char(datetime(),'MMddyyHHmm');
            FN = ['../output/ST_model_' [opt_type] '_' da];
            % load method 1
            load('ST_method1Results_04162121.mat');
    
        case 'KF'
            if ~exist('bI','var')
                b = 0.21/secpy; % ice-equivalent m/s from RACMO
            else
                b = bI/secpy;
            end
            
            switch opt_type
                case 'density'
                    disp(['No core for Korff'])
                    return
                case 'pres'
                    Ec = 6.3770e4; %6.6215e4;%5 
                    Eg = 4.6580e4; %4.9014e4;%5 
                case {'optr','optw'}
                    Ec = EcIn;
                    Eg = EgIn;
            end
            
            
            rho_s = 400;
            Tav =273.15-24.3;
            
            load('NXYwBHe_KorffpRES_for_Carlos_220515.mat')
            if ~exist('stk','var')
                stakeOK = 5;
            else
                stakeOK = stk;
            end
            pRES = flipud(NXYwBHe_KorffpRES(NXYwBHe_KorffpRES(:,1) == stakeOK,:));
            clear NXYwBHe_KorffpRES

            wIE = pRES(:,4);
            depthIE = pRES(:,3);     
            max_depth = 250;
            
            w_cIE_pRES = wFirn_from_wTot('KF', wIE, depthIE);
            w_cIE_pRES = w_cIE_pRES(depthIE < max_depth);
            depth_pRES = depthIE(depthIE<max_depth);
            
            da = char(datetime(),'MMddyyHHmm');
            FN = ['../output/KF_model_' opt_type '_' da];
             
            
    end
           
    %Setting up grid
    grid_size = 750.0; 
    z = linspace(0,250.0,grid_size); %depth from NXYwBHe_FP_pRES
    dz = mean(diff(z)); %crude
    z_prev = z;
    zs_prev = 0;

    %option 3: start with a line; then evolve
    rho = rho_i + (rho_s-rho_i).*exp(-z./40);
    rho_prev = rho;

    %temperature
    T_s = Tav;
    T_b = Tav + 2; %because of geothermal and viscous heating; 2 increase based on fletcher borehole
    T = linspace(T_s,T_b,length(z));
    T_prev = T;
    
    %grain size
    r2_s = 1e-9;
    r2 = linspace(1e-9,1e-5,length(rho)); % mm ^ 2  [default]
    r2_prev = r2;
    
    %accumulation
    b_avg = b; 

    %initialize vertical velocity
    w = zeros(1,length(z));

    %for plotting
    finished = 0;

%% set up & limitations

%z axes: z is depth (increasing from zero at surface to depth at bottom, depth can change
%w axes: w is positive downwards
%moving down, dw/dz = -, drho/dz = +)

for i=1:length(t)
    
    
    %% updates grid given accumulation, advection 
        dH = (b*rho_i/rho_s - w(1)).* dt;
        %H = H + dH;
        firn_base = z(end) + dH;
        dz_prev = dz;
        z = linspace(0,firn_base,grid_size); %z = 0 at top, depth changes
        dz = mean(diff(z));    
        d_grid = z - (z_prev + dH);  % the change in the location of each grid point in the upcoming interpolation.

    
        %% reinterpolate variables on new grid
 
        % rho
        %rho_prev = interp1(z_prev + zs,rho_prev,z,'linear','extrap');
        rho_prev = interp1Dgridsmalldx(rho_prev,d_grid,dz_prev,dH);
        rho_prev(1) = rho_s;   
        % T     
        %T_prev = interp1(z_prev + zs,T_prev,z,'linear','extrap');         
        T_prev = interp1Dgridsmalldx(T_prev,d_grid,dz_prev,dH);        
        T_prev(1) = T_s;
        % r2
        %r2_prev = interp1(z_prev + zs,r2_prev,z,'linear','extrap');
        r2_prev = interp1Dgridsmalldx(r2_prev,d_grid,dz_prev,dH);        
        r2_prev(1) = r2_s;
        %(w has no t dependence so is already reinitialized every time step)
    
    
    %% update overburden
    sigma = cumtrapz((g.*rho_prev.*dz));
    
    %% update velocity
    w = zeros(1,length(z));

    kc = zeros(size(w));
    kc(rho_prev > 550) = kc_high;
    kc(rho_prev <=550) = kc_low;
    
    int = (1./rho_prev).*kc.*(rho_i-rho_prev).*sigma.*exp(-Ec./(R.*T_prev)).*dz./r2_prev;
    
    Lz = length(z) ;
    w = zeros(1,Lz);
    w(end) = b_avg*rho_i/rho_prev(Lz); %assumes advection out base at average accumulation rate
    w = w(end) + cumsum([int(1:Lz-1) 0],'reverse');
    

    
    %% update temperature
    rho_stag = (rho_prev(1:end-1) + rho_prev(2:end))/2;
    kappa_stag = kappa_i * (rho_stag./rho_i).^2; 
    alpha = 1./(rho_prev.*c_p);
   

    T(2:end-1) = T(2:end-1) +...
           + dt .*...
                (alpha(2:end-1) .* ...
                    (kappa_stag(2:end).*T(3:end) -...
                     T(2:end-1).*(kappa_stag(2:end)+kappa_stag(1:end-1)) +...
                     kappa_stag(1:end-1).*T(1:end-2))./dz^2 -...
           w(2:end-1).*(T(3:end)-T(1:end-2))./(2*dz) ) ;
       

%    BCs
       %T_s2 = 15*sin(2*pi*t(i)/(lambda*secpy))+Tav;
       T(1) = T_s;
       T(end) = T(end-1)+.0034*dz; %T(end-1) + (G/kappa_i * dz);
%     
   
   %% update grain size
       r2(2:end) = r2_prev(2:end)  + dt *((w(1:end-1)+w(2:end))/2.*(r2_prev(1:end-1)-r2_prev(2:end))./dz + kg*exp(-Eg./(R.*T(2:end))) )   ;

    
    %% update density
    
        rho(2:end) = rho_prev(2:end) + dt.*((rho_prev(1:end-1) .* w(1:end-1)) - (rho_prev(2:end) .* w(2:end)))./dz;
        rho(1) = rho_s; %right now, assumes there is ALWAYS new accumulation/the top layer never evolves 


    %% save variables
    
    %save all data
    %dump....
    
    
    %convergence criteria -- sum of squares and ignores first 50 ( takes a
    %minute to deviate from input)
    if (sum((rho-rho_prev).^2)) < 1e-8 && i > 50 || i >= (length(t)-1)
            
        n_rho = 1+(n_i-1).*rho./rho_i;
        n_rho_mid = mean([n_rho(1:end-1);n_rho(2:end)]);

        wi = (w-b_avg).*secpy.*n_rho./n_i;

        zi = cumsum(diff(z).*n_rho_mid./n_i);
        zi = [0 zi];
        
        modelOutput.b = b;
        modelOutput.d = z;
        modelOutput.w = w; 
        modelOutput.wFirn = (w-w(end)).*secpy;
        modelOutput.r2 = r2;
        modelOutput.T = T;
        modelOutput.rho = rho;
        modelOutput.Ec = Ec;
        modelOutput.Eg = Eg;
        modelOutput.wFirnI = wi;
        modelOutput.dI = zi;
        modelOutput.loc = loc;
        modelOutput.stake = stakeOK;
        modelOutput.opt_type = opt_type;
        modelOutput.model_file = '/model/code/Model_041421.m';
        
        if s
            
            save(FN,'modelOutput');
        end
        
        finished = 1;

    end
    
    rho_prev = rho;
    T_prev = T;
    z_prev = z;
    r2_prev = r2;
    

    
    %% plotting
    if rem(i,n) == 0
        disp(['making progress... ' num2str(i)]);
    end
    
    %if plot_opt == 1 && rem(i,n) == 0 || finished == 1
    if plot_opt == 1 && finished ==1          
        model_c = [0/255 73/255 73/255];
        
        figure(50)
        
        set(gcf,'Pos',[ 0.0010    0.0410    1.5360    0.7488]*1e3)
        
        subplot(2,3,1)

        plot(rho,z,'Color',model_c,'LineWidth',2);
        hold on;
        plot(550*ones(1,length(z)), z,'--');
        switch loc
            case {'FP','ST'}
                plot(core.rho,core.d,'o','Color',[130/255,130/255,130/255,.6]);
            case 'KF'
        end        
        hold off;
        set(gca,'YDir','reverse')
        xlim([300,1000])
        ylim([0,max(z)])
        xlabel('\rho [kg/m^3]','FontSize',20)
        ylabel('z [m]','FontSize',20)
        title('Density profile at steady state','FontSize',20)
        set(gca, 'LineWidth',2,'Fontsize',20)
        legend('core','Arthern (2010)', 'Herron & Langway (1988)','Location','sw')
        
        switch loc
            case {'FP','ST'}
                rho_compare = interp1(z,rho,core.d,'linear','extrap');
                rho_compare(1) = rho_s;

                subplot(2,3,2)
                plot(core.rho-rho_compare,core.d,'Color',[143/255,0,0,.7],'LineWidth',2);
                set(gca,'YDir','reverse')
                ylim([0,max(z)])
                xlabel('\rho-fletcher density [kg/m^3]','FontSize',20)
                ylabel('z [m]','FontSize',20)
                title('Density profile mismatch at steady state','FontSize',20)
                set(gca, 'LineWidth',2,'Fontsize',20)
            case 'KF'
        end



    
    
        subplot(2,3,3)
        hold on;
        %plot(w_cIE_pRES, depth_pRES, 'ko','DisplayName','pRES')
        switch loc
            case {'ST','FP'}
                plot((w-b_avg)*secpy,z,'Color',model_c,'LineWidth',2,'DisplayName','model');
                plot(pRES.wFirn,pRES.d,'ko','DisplayName','PRES')
                
                xlim([0 max(w)+.1])
                ylim([0 max(z)])
                xlabel('w_c [m/yr]','FontSize',20)                
            case {'KF'}
                
                Iz = z<max_depth;
                n_rho = 1+(n_i-1).*rho(Iz)./rho_i;
                n_rho_mid = mean([n_rho(1:end-1);n_rho(2:end)]);

                wi = (w(Iz)-b_avg).*secpy.*n_rho./n_i;

                zi = cumsum(diff(z(Iz)).*n_rho_mid./n_i);
                zi = [0 zi];
                    
                plot(wi,zi,'Color',model_c,'LineWidth',2,'DisplayName','model');
                plot(w_cIE_pRES,depth_pRES,'ko','DisplayName','PRES') 
                
                xlim([0 max(wi)+.1])
                ylim([0 max(zi)])
                xlabel('w_c [m i.e./yr]','FontSize',20)                   
        end
        set(gca,'YDir','reverse')

        title('Compaction velocity','FontSize',20)
        set(gca, 'LineWidth',2,'Fontsize',20)
        legend();
        hold off;
        
        subplot(2,3,4)
        plot(r2,z,'Color',model_c,'LineWidth',2);
        set(gca,'YDir','reverse')
        ylim([0 max(z)])
        xlabel('r^2 [m^2]','FontSize',20)
        title('Grain size','FontSize',20);
        set(gca, 'LineWidth',2,'Fontsize',20)
        
        subplot(2,3,5)
        plot(T,z,'Color',model_c,'LineWidth',2);
        hold off;
        set(gca,'YDir','reverse')
        ylim([0 max(z)])
        xlim([246 - 20, 246 + 20])
        xlabel('T [K]','FontSize',20)
        title('Temperature','FontSize',20);
        set(gca, 'LineWidth',2,'Fontsize',20)

        subplot(2,3,6)
        plot(max(z),b*secpy,'o');
        hold on;
        set(gca,'YDir','reverse')
        xlabel('thickness [m]','FontSize',20)
        ylabel('b [m/yr]','FontSize',20)
        title('Accumulation vs thickness','FontSize',20);
        set(gca, 'LineWidth',2,'Fontsize',20)
        

        drawnow;
        
    end
    
   if finished==1
       break;
   end
end

end
