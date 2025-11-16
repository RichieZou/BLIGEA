function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,Fitness,iter,sv,fv)
% The variation operator for BLIGEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [N,D]       = size(ParentDec);
    Parent1Dec  = ParentDec(1:floor(end/2),:);
    Parent2Dec  = ParentDec(floor(end/2)+1:floor(end/2)*2,:);
    Parent1Mask = ParentMask(1:N/2,:);
    Parent2Mask = ParentMask(N/2+1:end,:);
    
    %% Crossover for mask
    OffMask = Parent1Mask;
    numberOfGroups = max(1,floor(D*iter));
    for i = 1 : N/2
        if rand < 0.5
            index = find(Parent1Mask(i,:)&~Parent2Mask(i,:));
            index = index(TS(ceil(mean(sv)*D),-Fitness(index)));
            OffMask(i,index) = 0;
        else
            index = find(~Parent1Mask(i,:)&Parent2Mask(i,:));
            index = index(TS(ceil(mean(sv)*D),Fitness(index)));
            OffMask(i,index) = Parent2Mask(i,index);
        end
    end

      
   %% Group-based mutation for mask (Innovation 2)
    [MaskGroupIndex] = MaskGroup(iter,OffMask,numberOfGroups,fv);

    if rand > 0.5
        for i = 1 : N/2
            half = ceil(max(2,max(MaskGroupIndex))/2);
            chosengroups = [randperm(half,1),randperm(half,1)+floor(max(2,max(MaskGroupIndex))/2)];
            
            Site1 = MaskGroupIndex == chosengroups(:,1);
            Site2 = MaskGroupIndex == chosengroups(:,2);
            sum1 = sum(Site1.*Fitness);
            sum2 = sum(Site2.*Fitness);
            if rand < iter % Flap to 0
                if sum1 < sum2
                    index = Site2;
                else
                    index = Site1;
                end
                ind = index==1;
                OffMask(i,ind) = 0;
            else % Flap to 1
                if sum1 < sum2
                    index = Site1;
                else
                    index = Site2;
                end
                ind = index==1;
                OffMask(i,ind) = 1;
            end
        end
    else
        for i = 1 : N/2
            half = ceil(max(2,max(MaskGroupIndex))/2);
            chosengroups = [randperm(half,1),randperm(half,1)+floor(max(2,max(MaskGroupIndex))/2)];
            
            Site1 = MaskGroupIndex == chosengroups(:,1);
            Site2 = MaskGroupIndex == chosengroups(:,2);
            sum1 = sum(Site1.*(1./(0.1+sv)));
            sum2 = sum(Site2.*(1./(0.1+sv)));
            if rand < iter 
                if sum1 < sum2
                    index = Site2;
                else
                    index = Site1;
                end
                ind = index==1;
                OffMask(i,ind) = 0;
            else 
                if sum1 < sum2
                    index = Site1;
                else
                    index = Site2;
                end
                ind = index==1;
                OffMask(i,ind) = 1;
            end
        end
    end
    
    %% Crossover and mutation for dec
    if any(Problem.encoding~=4)
        OffDec = GLP_OperatorGAhalf(Problem,iter,Parent1Dec,Parent2Dec,4,OffMask,sv); % 4 -- numberofgroups
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(N/2,D);
    end
end

function index = TS(numberOfGroups,Fitness)
%% Tournament selection helper function
    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(numberOfGroups,1,Fitness);
    end
end