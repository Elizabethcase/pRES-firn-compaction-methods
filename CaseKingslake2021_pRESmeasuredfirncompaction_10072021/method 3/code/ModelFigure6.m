% Creates model output for use in Figure6.m

% runs model using accumulation rates extracted from pRES and saves to a structure

locOpt = {'FP'};%{'FP','KF','ST'};


for iLoc = 1:length(locOpt)
    loc = locOpt{iLoc};%'KF';
    
%     switch loc
%         case {'FP','ST'}
%             opttype = {'density','pres'};
%         case 'KF'
%             opttype = {'pres'};
%     end

opttype = {'pres'};
    
    for iOpt = 1:length(opttype)



bI = ExtractAccumulation(loc);

opt_type = opttype{iOpt};%'pres';

plot_opt = 0;

s = 0;

da = char(datetime(),'MMddyyHHmm');

switch loc
    case 'FP'
        stakeOK = 1;
        FN = ['../model_output/FP_method3_' da '.mat'];

    case 'ST'
        load('NXYwBHe_STpRES_for_Carlos_110719.mat')
        stakeOK = unique(NXYwBHe_STpRES(:,1)); 
        FN = ['../model_output/ST_method3_' da '.mat'];
        
    case 'KF'
        load('NXYwBHe_KorffpRES_for_Carlos_220515.mat') % load pres data
        stakeOK = unique(NXYwBHe_KorffpRES(:,1));
        FN = ['../model_output/KF_method3_' da '.mat'];

end

for i=length(bI):-1:1
    
    %Korff run with activation energies from optimizing to stake 5 (pres)
    %Skytrain runs with activation energies from optimizing to stake 18
    %(pres)
    %Fletcher run with activation energies from optimizing to stake 1
    %(pres)
    
    MO(i) = Model_041421(loc,opt_type,plot_opt,s,bI(i),stakeOK(i));
    
end

save(FN,'MO')
clearvars 'MO'
    end
end