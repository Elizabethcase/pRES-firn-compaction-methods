% extracts strain rates and accumulation rates from FP, KF, ST, based on
% method 3 "ExtractAccumulation"

function [bI, strainRate] = ExtractStrainAccumulation(loc)

    max_depth = 350; 
    
    switch loc
        
        case 'FP'
            
            load('NXYwBHe_FP_pRES.mat');
            NXYwBHe_FP_pRES = flipud(NXYwBHe_FP_pRES);

            pRES(1).stake = 1;
            pRES(1).dI = NXYwBHe_FP_pRES(:,3);
            pRES(1).wI = NXYwBHe_FP_pRES(:,4);
            pRES(1).e = NXYwBHe_FP_pRES(:,6);

            stakeOK = [1];
            
            %pRES = pRES(pRES.dI < max_depth);

            clear NXYwBHe_FP_pRES;
            
            ub = 200;
            lb = 300;
        
        case 'ST'
        
            load('NXYwBHe_STpRES_for_Carlos_110719')
            stakeOK = unique(NXYwBHe_STpRES(:,1));
            
            for i = 1:length(stakeOK)
            
                I = NXYwBHe_STpRES(:,1)==stakeOK(i);

                pRES(i).stake = stakeOK(i);

                dI = NXYwBHe_STpRES(I,3);
                wI = NXYwBHe_STpRES(I,4);
                e = NXYwBHe_STpRES(I,7);
                H = NXYwBHe_STpRES(I,6);

                pRES(i).dI = flipud(dI(dI < max_depth));
                pRES(i).wI = flipud(wI(dI < max_depth));
                pRES(i).e = flipud(e(dI < max_depth));
                pres(i).H = unique(H);

                
            end
                
            ub = 150;
            lb = 280;
                
            clear NXYwBHe_STpRES
        
        case 'KF'
            
            load('NXYwBHe_KorffpRES_for_Carlos_220515.mat')
            stakeOK = unique(NXYwBHe_KorffpRES(:,1));

    
           for i = 1:length(stakeOK)
            
                I = NXYwBHe_KorffpRES(:,1)==stakeOK(i);

                pRES(i).stake = stakeOK(i);

                dI = NXYwBHe_KorffpRES(I,3);
                wI = NXYwBHe_KorffpRES(I,4);
                e = NXYwBHe_KorffpRES(I,7);
                H = NXYwBHe_KorffpRES(I,6);

                pRES(i).dI = flipud(dI(dI < max_depth));
                pRES(i).wI = flipud(wI(dI < max_depth));
                pRES(i).e = flipud(e(dI < max_depth));
                pres(i).H = unique(H);
                                
            
           end
           
           ub = 180;
           lb = 280;
           
           clear NXYwBHe_KorffpRES
            
    end

    nI = 1.78;
    rhoI = 906;
    spy = 365.25*24*3600;
    c = 2.983e8;
    f = 305e6;

%%

    for i = 1:length(pRES)

        dI = pRES(i).dI;
        wI = pRES(i).wI;

        pRES(i).bI = accumulation(dI, wI, ub, lb);

    end
    
    bI = [pRES.bI];

end


function bI = accumulation(z, w, ub, lb)

    I = z>ub & z < lb;

    strain_fit = polyfit(z(I) ,w(I),1);
    strain = polyval(strain_fit,z);
    bI = max(strain);

end