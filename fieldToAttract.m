function attractField = fieldToAttract(aField,isAlate)

global veeR nuVec
%disp(nuVec)
% Transforming field matrix into weighted landing probabilities (i.e. accounting for attraction)
sf=size(aField);
attractField=aField;
attractField(aField==0)=zeros(length(find(aField==0)),1);
attractField(aField==1)=ones(length(find(aField==1)),1);

for i=1:length(nuVec), attractField(aField==((i-1)*10)+2)=nuVec(i)*ones(length(find(aField==((i-1)*10)+2)),1); end

attractField(aField==3)=veeR*ones(length(find(aField==3)),1);
l1=sum(sum(aField==1));
l2=sum(sum(mod(aField,10)==2));
l3=sum(sum(aField==3));
if and(isAlate,((l1+l2+l3)~=sf(1)*sf(2))), disp('ERROR There were plants that were not accounted for'); end
end

