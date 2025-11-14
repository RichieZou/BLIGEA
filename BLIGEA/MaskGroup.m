function [MaskGroupIndex] = MaskGroup(iter,ParentMask,numberOfGroups,fv)


% MaskGroupIndex=[];

noOfSolutions = size(ParentMask,1);
numberOfVariables=size(ParentMask,2);
varsPerGroup = floor(numberOfVariables/numberOfGroups);

%varsPerGroup=max(ceil(numberOfVariables/2-((ceil(iter/0.1)-1)*(numberOfVariables/2)/9)),2);
%numberOfGroups=floor(numberOfVariables/varsPerGroup);

vars=sum(ParentMask);
[~,I] = sortrows([fv;vars]','descend');
outIndexList = -ones(1,numberOfVariables);
for i = 1:numberOfGroups-1
    outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
end
outIndexList(I(((numberOfGroups-1)*varsPerGroup)+1:end)) = numberOfGroups;
%MaskGroupIndex=repmat(outIndexList,noOfSolutions,1);
MaskGroupIndex=outIndexList;

end

