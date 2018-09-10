function fieldPnew = apterousProbingChain(fieldP, dispIndex)

global pfeedS sigA2 Pacq Pinoc epsVec Pmutate returnToSender


% Note: plant attractiveness factors are implemented in fieldToAttract.m

% Example call: out=apterousProbingChain(fieldP=[1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4], dispIndex=1);
%
% fieldP is the matrix of plant status'
% dispIndex is the source plant of the aphid
%
% can include resistant plants but default is susc only
%
% dispIndex=-99 means that an immigrant has just entered the dispersal
% chain... 

strains=length(epsVec);
infected=0;   % i.e., always start out dispersals without virus for NPT (extend here if need to allow immigrating with virus)
dimf=size(fieldP);
mRows=dimf(1);
mCols=dimf(2);

activInd=dispIndex;

feeds=0;  % binary 1 for in state of feeding
fieldPnew=fieldP;
dispR(['new dispersing aptarae is on plant # ' num2str(activInd)]);
twoDind=zeros(2,1);

while(feeds==0)
    subrand1=rand();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Plant selection stage % calculating where the aphid goes next(activInd) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % selecting next plant (activInd)
    [twoDind(1),twoDind(2)]=ind2sub([mRows mCols],activInd);
    subFieldP = makeNearNborField(fieldP,activInd);
    subFieldPattr=fieldToAttract(subFieldP,0);
    if returnToSender==0, subFieldPattr(5)=0; end % i.e. aphid cannot disperse back to source plant
    if sum(sum(subFieldPattr))==0, dispR('looks like aptarae was on a dead plant surrounded by dead plants: Work Around - kill aptarae'); break; end
    fracSubFieldP=subFieldPattr/sum(sum(subFieldPattr));
    scaledCumSum=(1-sigA2)*cumsum(fracSubFieldP(:));
    outTrue=find(subrand1<scaledCumSum);
    if isempty(outTrue)~=1  % Adjusting aptarae position in overall field according to nearest neighbour subfield
        [L,R]=ind2sub([3 3],outTrue(1));
        L=L-2;
        R=R-2;
        twoDind(1)=twoDind(1)+L;
        twoDind(2)=twoDind(2)+R;
    else
        activInd=-98;            % aphid death/emigration
        dispR(['dispersing aptarae mortality' num2str(activInd)]);
        break;
    end
    if twoDind(1)==0, twoDind(1) = mRows; end  % Toroidal boundary conditions
    if twoDind(2)==0, twoDind(2) = mCols; end
    if twoDind(1)==(mRows+1), twoDind(1) = 1; end
    if twoDind(2)==(mCols+1), twoDind(2) = 1; end
    activInd=sub2ind([mRows mCols],twoDind(1),twoDind(2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Probing and feeding stage given activInd %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subrand2=rand();
    if mod(fieldP(activInd),10)==2 % PLANT IS INFECTED
        indAdj=floor(fieldP(activInd)/10); % determine which strain for multi-strain sims/evo
        if subrand2<(pfeedS*epsVec(indAdj+1))
            feeds=1; 
            dispR('feeds infected plant'); 
            infected=0;    % SETTLE
        else % DON'T SETTLE
            subrand3=rand();
            if subrand3<Pacq
            dispR('probes infected, departs with virus');                  % ACQUISITION
            infected=(floor(fieldP(activInd)/10)+1);
            else
            dispR('probes infected, but does not depart with virus');     
            end
        end
    elseif fieldP(activInd)==1 % PLANT IS SUSCEPTIBLE
        if infected>0                                                      % INOCULATION
            subrand4=rand();
            if subrand4<Pinoc
                subrand5=rand() ; % mutation
                dispR(['Old infection was ' num2str(2+10*(infected-1))])
                if subrand5<(Pmutate/2) % mutate left                          % MUTATION
                    infected=infected-1*(infected~=1);
                    dispR(['was there a left mutation?' num2str(infected~=1)])
                elseif subrand5<Pmutate % mutate right
                    infected=infected+1*(infected~=strains);
                    dispR(['was there a right mutation?' num2str(infected~=strains)])
                end
                fieldPnew(activInd)=2+10*(infected-1);
                dispR(['new infection is ' num2str(fieldPnew(activInd)) ' at location ' num2str(activInd)]);
            else
                dispR('visited susceptible with virus but no inoculation ocurred');
            end
        end
        if subrand2<pfeedS
            feeds=1;
            dispR('feeds susceptible');
            infected=0;
        else
            dispR('probes susceptible, departs');
            infected=0;
        end
        %elseif fieldP(activInd)==3 % FOR IF NEED TO INCORP. RESISTANT PLANTS 
        %    if subrand2<pfeedR
        %        feeds=1;
        %        dispR('feeds resistant');
        %    else
        %        dispR('probes resistant, departs');
        %        if probelosevirus, infected=0; end
        %    end
    else
        disp('Error! Aphid seems to be on a vacant site... dumping current plant code and position');
        probnum=fieldP(activInd);
        disp(probnum);
        break;
    end
end

% Appending final position to field array, recover field matrix via out-(out>100)*101 and settling index via find(out>100)
if activInd>0, fieldPnew(activInd)=fieldPnew(activInd)+5; end % adding 5 to record the position of the settled aphid

end
