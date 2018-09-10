function out = qClump(field, qType1, qType2)
% For now this expects a vector of plant types for which to generate 
% homogeneous clumping measures, eventually some if clauses could append 
% additional heterogeneous clumping if qType2 differs from qType1

% e.g. qClump(maty, 2, 1)

considerDeads=0;

% pair density is rho_{II}=P_I*P_I|I
% currently returning P_I|I
locType1=find(field==qType1);
if isempty(locType1)~=1
    numType2=zeros(1,length(locType1));
    for i=1:length(locType1)
        %if field(locType1(i))~=qType1, disp('Encountered Error in qClump indexing of qType !!'); break; end
        %disp(field)
        %disp(locType1(i))
        nNbor=makeNearNborField(field,locType1(i));
        % NOTE: accumaray is allegedly fast but doesnt like zeros therefore
        % add 1 but bear in mind when extracting totals
        
        accumForm=accumarray(((nNbor(:)==qType2)+1),1); %nnz(nNbor(:)==qType2);
        locType2temp=accumForm(2);
        if considerDeads==1
            accumForm=accumarray(((nNbor(:)==0)+1),1); %nnz(nNbor(:)==qType2);
            locType3temp=accumForm(2);
        else
            locType3temp=0;
        end
        %locType2temp=nnz(nNbor(:)==qType2);
        %locType3temp=nnz(nNbor(:)==0);
        
        %locType2temp=sum(nNbor(:)==qType2);
        %locType3temp=sum(nNbor(:)==0);
        
        %locType2temp=sum(nNbor(:)==qType2);
        %locType3temp=sum(nNbor(:)==0);
        
        if qType1==qType2, locType2temp=locType2temp-1; else locType3temp=locType3temp-1; end % deleting the centre count
        if considerDeads==1 
           numType2(i)=locType2temp/(4-(locType3temp-4)); % denom is based on 4 cell nhood, problem if different nhood! denom removes 4 corner zeros from the dead tally as they are not dead they are just outside nhoof
        else
           numType2(i)=locType2temp/4; 
        end
    end
    out=mean(numType2);
else
    out=0;
end
