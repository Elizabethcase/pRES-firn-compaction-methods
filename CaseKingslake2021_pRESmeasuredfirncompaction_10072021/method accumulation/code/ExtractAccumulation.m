% extracts ice-equivlent accumulation rates from FP, KF, ST ice-equivalent data, saves to a matfile
% for use in model & Figure5

function bI = ExtractAccumulation(loc)

    max_depth = 400;
    
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

        pRES(i).bI = GetAccumulationRate(loc,dI, wI);

    end
    
    bI = [pRES.bI];

end


