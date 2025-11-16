function [MaskGroupIndex] = MaskGroup(iter,ParentMask,numberOfGroups,fv)
% Grouping strategy for the mask based on fitness and variable sums

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    noOfSolutions = size(ParentMask,1);
    numberOfVariables = size(ParentMask,2);
    varsPerGroup = floor(numberOfVariables/numberOfGroups);
    
    %% Sort variables based on fv and sum of masks
    vars = sum(ParentMask);
    [~,I] = sortrows([fv;vars]','descend');
    
    %% Assign group indices
    outIndexList = -ones(1,numberOfVariables);
    for i = 1:numberOfGroups-1
        outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
    end
    outIndexList(I(((numberOfGroups-1)*varsPerGroup)+1:end)) = numberOfGroups;
    
    MaskGroupIndex = outIndexList;
end