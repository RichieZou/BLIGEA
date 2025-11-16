function [Offspring] = GLP_OperatorGAhalf(Problem,iter,Parent1,Parent2, numberOfGroups,OffMask,sv)
% Genetic operators for decision variables based on groups

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [proC,disC,~,disM] = deal(1,20,1,20);
    [N,D]   = size(Parent1);
    
    %% Genetic operators for real encoding (SBX and Polynomial Mutation)
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;
    
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    [outIndexList,~] = CreateGroups(iter,numberOfGroups,Offspring,D,OffMask,sv); 
    chosengroups = randi(numberOfGroups,size(outIndexList,1),1);    
    Site = outIndexList == chosengroups;
    
    mu    = rand(N,1);
    mu    = repmat(mu,1,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));    
end


function [outIndexArray,numberOfGroupsArray] = CreateGroups(iter,numberOfGroups, xPrime, numberOfVariables,OffMask,sv)
%% Create groups by ordered grouping
% Creat groups by ordered grouping

    noOfSolutions = size(xPrime,1);
    outIndexArray = [];
    numberOfGroupsArray = [];
    
    index = sv>0 & sv<1;
    count = sum(index);
    varsPerGroup = floor(count/numberOfGroups);
    
    Ind = -ones(1,numberOfVariables);
    outIndexList = -ones(1,numberOfVariables);
    vars = mean(xPrime);
    [~,I] = sort(vars(index),'descend');
    Ind(index) = I;
    for i = 1:numberOfGroups
        outIndexList(Ind>=(varsPerGroup*(i-1)+1) & Ind<varsPerGroup*i) = i;
    end
    outIndexArray = repmat(outIndexList,noOfSolutions,1);
    numberOfGroupsArray = [numberOfGroupsArray;numberOfGroups];

end