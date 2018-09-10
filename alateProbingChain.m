function fieldPnew = alateProbingChain(fieldP, dispIndex)

global pfeedS sigA1 Pacq Pinoc epsVec Pmutate

% Note: plant attractiveness factors are implemented in fieldToAttract.m
% Example call: out=alateProbingChain(fieldP=[1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4], dispIndex=1);
%
% fieldP is the matrix of plant status'
% dispIndex is the source plant of the aphid
%
% Can includes resistant-susceptible plant mixtures, default is only
% susceptible
%
% dispIndex=-99 means that an immigrant has just entered the dispersal
% chain... if immigrant is alate it can go anywhere... if immigrant is
% apterous it will first encounter the border rows (see apterousProbingChain.m for case of apterous aphids)

% probelosevirus=1;
strains=length(epsVec);
infected=0;   % i.e. always start out dispersals without virus for NPT !!! Situation for immigrants requires treatment of this, current default is no immigrants
dimf=size(fieldP);
mRows=dimf(1);
mCols=dimf(2);
activInd=dispIndex;
feeds=0;  % binary 1 for in state of feeding
fieldPnew=fieldP;
dispR(['new dispersing alate is on plant # ' num2str(activInd)]); 
dispR(['new dispersing aptarae is on plant # ' num2str(activInd)]); 
twoDind=zeros(2,1);

while(feeds==0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Plant selection stage % calculating where the aphid goes next(activInd) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Alate, selecting next plant
    fieldPattr=fieldToAttract(fieldP,1);
    % transform fieldP to take account of attraction
       
    % converting to a vector with final element = 1
    matchIt=find(rand()*sum(fieldPattr(:))<(1-sigA1)*cumsum(fieldPattr(:)));

    if isempty(matchIt)==1
        if infected>0
            dispR('alate that had virus emigration/death');   
            oldActivInd=activInd;
            activInd=-96;   % NOTE: '-96' means alate emigration/mortality with virus! 
                            % ENTER a coded blimp into the output field in
                            % the last plant the lost aphid was on i.e.
                            % fieldP(oldActivInd)=fieldP(oldActivInd)-1000;
                            % In Sto_spatial look for negative values of
                            % fieldP and add 1000
        else
            dispR('alate emigration/death');  
            activInd=-98;
        end
        break; 
    else
        activInd=matchIt(1);
    end
    dispR(['dispersing alate is on plant # ' num2str(activInd)]);
    dispR(['corresponding plant type is # ' num2str(fieldP(activInd))])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Probing and feeding stage given activInd %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subrand2=rand();
    if mod(fieldP(activInd),10)==2 % PLANT IS INFECTED
        indAdj=floor(fieldP(activInd)/10); % which strain
        if subrand2<(pfeedS*epsVec(indAdj+1))  % SETTLE
            feeds=1; 
            dispR('feeds infected plant'); 
            infected=0;
        else
            subrand3=rand();
            % DON'T SETTLE 
            if subrand3<Pacq                                  % ACQUISITION
                dispR('probes infected, departs with virus');
                infected=(floor(fieldP(activInd)/10)+1);
            else
                dispR('probes infected, but does not acquire virus');
                infected=0;
            end
        end
    elseif fieldP(activInd)==1 % PLANT IS SUSCEPTIBLE
        if infected>0                                            % MUTATION
            subrand4=rand() ;
            if subrand4<Pinoc
                subrand5=rand() ;
                dispR(['Old infection was ' num2str(2+10*(infected-1))])
                if subrand5<(Pmutate/2) % mutate left
                    infected=infected-1*(infected~=1);
                    dispR(['was there a left mutation?' num2str(infected~=1)])
                elseif subrand5<Pmutate % mutate right
                    infected=infected+1*(infected~=strains);
                    dispR(['was there a right mutation?' num2str(infected~=strains)])
                end
                fieldPnew(activInd)=2+10*(infected-1);        % INOCULATION
                dispR(['new infection is ' num2str(fieldPnew(activInd)) ' at location ' num2str(activInd)]);
            else
                dispR('visits susceptible with virus but no inoculation');
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
        %elseif fieldP(activInd)==3 % PLANT IS RESISTANT
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
if activInd==-96, fieldPnew(oldActivInd)=fieldPnew(oldActivInd)-1000; end % subtracting 1000 to mark/record plant on which viruliferous winged aphid mortality/emigration occurred
end
