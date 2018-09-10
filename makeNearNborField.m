function subFieldP = makeNearNborField(fieldP,dispIndex)

%   Take in a matrix of plant status' representing a field
%   Take in a single digit index (for matlab this corresponds to repeated vertical counts down and then across columns)
%   Return a sub-matrix that consists of the indexed plant as the focal
%   plant and 8 immediate neighbors
%   When a neigboring plant site is outside the field treat this as a
%   vacant site 
numNbor=4;
sizF=size(fieldP); 

[rowNum,colNum]=ind2sub(sizF,dispIndex); 

% NOTE: Includes toroidal boundary conditions as default
superFieldP=zeros(size(fieldP)+2);
superFieldP(2:end-1,2:end-1)=fieldP; 
superFieldP(1,:)=superFieldP(end-1,:); % make top boundary periodic
superFieldP(end,:)=superFieldP(2,:); % make bottom boundary periodic
superFieldP(:,1)=superFieldP(:,end-1); % make left boundary periodic
superFieldP(:,end)=superFieldP(:,2); % make right boundary perioduc

subFieldP=superFieldP(rowNum:(rowNum+2),colNum:(colNum+2));

if numNbor==4
    subFieldP(1,1)=0;
    subFieldP(1,3)=0;
    subFieldP(3,1)=0;
    subFieldP(3,3)=0;
end

end

