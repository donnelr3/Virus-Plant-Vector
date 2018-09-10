function returnSet=Sto_spatial(fullTimeIn, inocTimeIn, pfeedSin, inFolder1, inFolder, dimsIn, evoNums, numTpoints, bSS, PacqIn, mIn, aphidTypeIn,printPams,R0runIn)

global vee veeR pfeedS pfeedR sigA1 sigA2 Pacq Pinoc rePlant nuVec epsVec Pmutate returnToSender

tic

% Ruairi Donnelly, copyright, 1/9/2018 University of Cambridge
%
% Equal weighting of alates and aptarae for morph determination (default)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Setting up the field%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Array #1 : the status of the plant in each cell fieldP
% Array #2 : the number of alate aphids in each cell fieldA1
% Array #3 : the number of apterous aphids in each cell fieldA2
% Plant status codes:
% 0 dead or absent
% 1 susceptible
% 2 infected   (12,22 etc. for multiple strains)
% 3 resistant

% By default aphid carrying capacity is constant (i.e. not perdiodic), but code
% facilitates periodicity

% Default is single virus strain with infection coded as 2 (for multi-strains or strain evolution 2+m*(10) denotes (m+1)'th strain)

% Default is single susceptible host variety... but code written to also
% incorporate resistant variety for resistant-susc mixtures

R0run=R0runIn;
outputClump=0;
outputMetrics=1;
outputSnapshot=1;
printCounts=0;
plotBurnIn=1;
flagI1=0;
flagI2=0;

blah=[];
startSeed=evoNums(end);
nuVec=evoNums(1:floor(length(evoNums)/2))
epsVec=evoNums((floor(length(evoNums)/2))+1:(end-1))
strains=length(epsVec);
cstring=cell(1,length(epsVec));
for i=1:length(epsVec)
    cstring{i}=num2str(epsVec(i));
end

Pmutate=mIn;
returnToSender=1;
plantDeathDisperses=0;
plantDeathKillsAphids=0; % By default aphid population continues on instantaneously re-planted plant

aphidType=aphidTypeIn; % Codes for aphid population: 1-Concurrent (alate and apterous);2-alate only ; 3-aptarae only

animate=0;
extinct=0;
sampleField=1;

if fullTimeIn==-1
    fullTime=1000;
    useSSmat=0;
    disp('Not using any stored aphid densities')
else
    fullTime=fullTimeIn;
    useSSmat=1;
    disp('Using stored aphid densities!')
end

% markerSnaps is for setting frequency to output spatial snapshots of populations
markerSnaps=0; 
desiredNum=20; % of snap outputs
snapFileFreq=round(fullTime/desiredNum);

rePlant=1;
k=dimsIn; % field dimensions
inocTime=inocTimeIn; 

% several parameters to govern condition for steady-state vector population dynamics having been
% reached... Use for virus invasion into environment with endemic vector 
tolAphidAvr=0.01;
margAphidAvr=5;
movingAphidAvr=0;
countSS=0;
threshSS=6; 

%%%%%%%%% Ignore, variables used for testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countOne=0;
countTwo=0;
countThree=0;
countFour=0;
countFive=0;
countSeven=0;
countSix=0;
countEight=0;
countNine=0;
countTen=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=0;
seasonL=40;
k1_base=10; 
k2_base=10; 
k1_fluc=0;
k1_prd=2*seasonL; 
k2_fluc=0;
k2_prd=2*seasonL;

vee=1;     % attraction to infected plant
veeR=1;     % attraction to resistant plant
pfeedS=pfeedSin; % prob of accepting probed suc plant
pfeedR=0.75; % prob of accepting probed resis plant (by default there will be no resistants)

% Details for wigned aphid production response curve (to on-plant crowding)
apterWgt=1;
apterSatConst=(2/4)*(k1_base); 
expo=6;

%%%%%%%%%%%%%%%
% plant rates %
%%%%%%%%%%%%%%%
Pacq=PacqIn;
Pinoc=PacqIn;
bS=bSS; 
bI=bS;
bR=bS;
al=bS; % disease death / roguing
gam=bS; % recovery
rate_vecP=[bS bI bR al gam];
tot_rate_vecP=rate_vecP; 
%%%%%%%%%%%%%%%
% aphid rates %
%%%%%%%%%%%%%%%
sigA1=20/100; %probability of loss of alate aphid per flight;
sigA2=20/100; % probability of loss of apterous aphid per apterous journey between plants;
bA1=1/50; % rate of on-plant alate mortality
bA2=1/50; % rate of on-plant apterous mortality
if aphidType==1
    immA1=0; %JOINT:4; %10; % 2; % #immigrants per field i.e. deliberately not density-dependent     # 1 per day?
    immA2=0;
elseif aphidType==2
    immA1=0; %10; % 2; % #immigrants per field i.e. deliberately not density-dependent
    immA2=0;
elseif aphidType==3
    immA2=0; %10; % 1; % #immigrants per field
    immA1=0;
end
phiA1=2; % rate of alate dispersal
phiA2=2; % rate of apterous dispersal
rA1=2; % rate of alate reproduction
rA2=2; % rate of apterous reproduction
rate_vecA1=[immA1 phiA1 rA1 bA1];
rate_vecA2=[immA2 phiA2 rA2 bA2];
tot_rate_vecA1=zeros(size(rate_vecA1)); 
tot_rate_vecA2=zeros(size(rate_vecA2)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating a field of plants
numInocula=1;                               % code to generate 'inocula' x random positions for the inoculum
fieldP=ones(k,k, 'int32');
inocvec=ceil(rand(numInocula,2)*k);         % random x_1 coordinate for inocula
% Generating a field of alate aphid populations (each cell corresponds to above plant)
if (useSSmat~=1)||(fullTimeIn==-1)
    if aphidType==3
        Pdie=1-(1-sigA2)*pfeedS/(1-(1-sigA2)*(1-pfeedS));      % distribute theoretical aphid steady-state over all plants
        AA00=ceil(k2_base*(1-(1/rA2)*(bA2+phiA2*Pdie)))-1;
        fieldA2=AA00*ones(k,k, 'int16');
        fieldA1=zeros(k,k, 'int16');
    elseif aphidType==2
        Pdie=1-(1-sigA1)*pfeedS/(1-(1-sigA1)*(1-pfeedS));      % distribute theoretical aphid steady-state over all plants
        AA00=ceil(k1_base*(1-(1/rA1)*(bA1+phiA1*Pdie)));
        fieldA1=AA00*ones(k,k, 'int16');
        fieldA2=zeros(k,k, 'int16');
    else
        fieldA1=zeros(k,k, 'int16');
        indRand=ceil(rand*(k*k));
        fieldA1(indRand)=1;            
        fieldA2=zeros(k,k, 'int16');
    end
    disp('GENERATED NEW INITIAL APHID DENSITIES !!!!')
else
    disp([inFolder1 'SSaphidDistribution.csv'])
    mat=csvread([inFolder1 'SSaphidDistribution.csv']);
    disp('OLD APHID DENSITIES READ-IN OK !!!!')
    sm=size(mat);
    if aphidType==3
        fieldA2=int16(reshape(mat(randperm(numel(mat))),[sm(1),sm(2)]));
        fieldA1=zeros(k,k, 'int16');
    elseif aphidType==2
        fieldA1=int16(reshape(mat(randperm(numel(mat))),[sm(1),sm(2)]));
        fieldA2=zeros(k,k, 'int16');
    else
        matt1=mat(1:(sm(1)/2),:);
        matt2=mat((sm(1)/2)+1:end,:);
        fieldA1=int16(reshape(matt1(randperm(numel(matt1))),[(sm(1)/2),sm(2)]));
        fieldA2=int16(reshape(matt2(randperm(numel(matt2))),[(sm(1)/2),sm(2)]));
    end
end

disp(['Alate row #1 ' num2str(fieldA1(1,1:10))])
disp(['Aptarae row #1 ' num2str(fieldA2(1,1:10))])
%if animate, MovA(indx)=getframe(gcf); end 

% BLOCK allocating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalAlivePlants=-1*ones(1,1000000, 'int16');
totalInfectedPlants=-1*ones(1,1000000, 'int16');
totalAlates=-1*ones(1,1000000, 'int32');
totalAptarae=-1*ones(1,1000000, 'int32');
timeVec=-1*ones(1,1000000);
inocProd=-1*ones(1,1000000, 'int32');
alatesOnSusc=-1*ones(1,1000000, 'int32');
aptOnSusc=-1*ones(1,1000000, 'int32');
alateStrains=-1*ones(1000000,strains, 'int32');
aptStrains=-1*ones(1000000,strains, 'int32');
strainClump =-1*ones(1000000,strains);
iStrainHosts=zeros(1,strains, 'int16');
sStrainHosts=-1*ones(1,1000000, 'int16');
totalStrains=-1*ones(1000000,strains, 'int16');
K1=-1*ones(1,1000000);
K2=-1*ones(1,1000000);
marker=0;
inoculated=0;
counter=0;
outputFile=[inFolder 'EpidemicMetrics_' inFolder(end-1) '.csv'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newInfR0A1=[];
newInfR0A2=[];
%%% OUTPUT PARAMETER SET
if printPams>0
    nameOutVec={'psai1','psai2','vee', 'pfeedS', 'sigA1', 'sigA2', 'immA1', 'immA2', 'phiA1', 'phiA2', 'rA1', 'rA2', 'bA1', 'bA2', 'Pmutate', 'bS', 'bI', 'al', 'gam', 'Pacq', 'Pinoc', 'rePlant', 'returnToSender', 'startSeed', 'aphidType', 'apterSatConst', 'expo', 'tolAphidAvr', 'margAphidAvr', 'threshSS', 'rePlant', 'numInocula', 'nuVec(:)', 'epsVec(:)'};
    valsOutVec=[k1_base k2_base vee pfeedS sigA1 sigA2 immA1 immA2 phiA1 phiA2 rA1 rA2 bA1 bA2 Pmutate bS bI al gam Pacq Pinoc rePlant returnToSender startSeed aphidType apterSatConst expo tolAphidAvr margAphidAvr threshSS rePlant numInocula  nuVec epsVec];
    fidel=fopen([inFolder 'parameterSet.csv'],'wt');
    disp([inFolder 'parameterSet.csv']);
    fprintf(fidel,'%s,',nameOutVec{1:end-1});
    fprintf(fidel,'%s\n',nameOutVec{end});
    fprintf(fidel,'%d,',valsOutVec(1:end-1));
    fprintf(fidel,'%d\n',valsOutVec(end));
    fclose(fidel);
end

timeFilt=linspace(inocTime,fullTime,numTpoints);
zNew=1;
% START Of GILLESPIE While Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while t<fullTime
    counter=counter+1;
    inocProd(counter)=0;
    % Block Allocating
    if mod(counter,1000000)==0
        totalAlivePlants=[totalAlivePlants -1*ones(1,1000000, 'int16')];
        totalInfectedPlants=[totalInfectedPlants -1*ones(1,1000000, 'int16')];
        totalAlates=[totalAlates -1*ones(1,1000000, 'int32')];
        totalAptarae=[totalAptarae -1*ones(1,1000000, 'int32')];
        timeVec=[timeVec -1*ones(1,100000)];
        sStrainHosts=[sStrainHosts -1*ones(1,1000000, 'int16')];
        inocProd=[inocProd -1*ones(1,1000000, 'int32')];
        alatesOnSusc=[alatesOnSusc -1*ones(1,1000000, 'int32')];
        aptOnSusc=[aptOnSusc -1*ones(1,1000000, 'int32')];
        alateStrains=[alateStrains; -1*ones(1000000,strains, 'int32')];
        aptStrains=[aptStrains; -1*ones(1000000,strains, 'int32')];
        totalStrains=[totalStrains; -1*ones(1000000,strains, 'int16')];
        strainClump=[strainClump; -1*ones(1000000,strains)];
        K1=[K1 -1*ones(1,1000000)];
        K2=[K2 -1*ones(1,1000000)];
    end
    
    % Seeding Infection
    if ((t>inocTime)&&(inoculated==0))&&(fullTimeIn~=(-1))
        for l=1:numInocula, fieldP(inocvec(l,1),inocvec(l,2))=2+10*(startSeed-1); end
        if sum(fieldP==2+10*(startSeed-1))>numInocula, disp('!!! Error! Seeded too many'); break; end
        if sum(fieldP==2+10*(startSeed-1))<numInocula, disp('!!! Error! Seeded too many'); break; end
        disp('Dropping in Infection to Commence Epidemic!!!');
        inoculated=1;
    end
    
    % Calculating Census' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sites_s=find(fieldP==1);
    sites_i=find(mod(fieldP,10)==2); % Using mod counts all the viruses if multi-strain
    sites_r=find(fieldP==3);
    if (inoculated==1)&&(isempty(sites_i)), disp('Epidemic Has Terminated'); extinct=1; break; end
    if (R0run==2)&&(length(sites_i)>=50)   %%% Invasion threshold: use this for studies of invasion success
        disp('Epidemic Invasion Threshold Has Been Reached');
        extinct=0;
        break;
    end
    for i=1:strains, iStrainHosts(i)=length(find(fieldP==2+(i-1)*10)); end
    sStrainHosts(counter)=length(sites_s);
    K1(counter)=k1_base+k1_fluc*sin((2*pi*t/k1_prd)+pi/8);
    K2(counter)=k2_base+k2_fluc*sin((2*pi*t/k2_prd)+pi/8);
    sites_alive=[sites_s; sites_i; sites_r];
    %vhosts=length(sites_v);
    shosts=length(sites_s);
    ihosts=length(sites_i);
    rhosts=length(sites_r);
    ahosts=length(sites_alive);
    A1aphids=sum(sum(fieldA1));
    A2aphids=sum(sum(fieldA2));
    % Alate counts on each plant type
    A1_1=int16((fieldP)==1).*(fieldA1);
    A1_2=int16((mod(fieldP,10)==2)).*(fieldA1);
    alateStrTemp=ones(1,strains)*(-1);
    % Aptarae counts on each plant type
    A2_1=int16((int16(fieldP)==1)).*(fieldA2);
    A2_2=int16((mod(int16(fieldP),10)==2)).*(fieldA2);
    aptStrTemp=ones(1,strains)*(-1);
    for i=1:strains
        isTemp=(fieldP==(2+(i-1)*10));
        alateStrTemp(i)=sum(sum(int16(isTemp).*(fieldA1)));
        aptStrTemp(i)=sum(sum(int16(isTemp).*(fieldA2)));
    end
    
    alatesOnSusc(counter)=sum(sum(A1_1));
    aptOnSusc(counter)=sum(sum(A2_1));
    alateStrains(counter,:)=alateStrTemp;   % total # of aphids on plants bearing a particular infection
    aptStrains(counter,:)=aptStrTemp;
    if sum(alateStrTemp+aptStrTemp)~=sum(sum(A1_2+A2_2)), disp('Error!!! Strain specific aphid count does not match total!'); disp(epsVec); disp(fieldP);  disp('!!!! NEW CHECK !!!!!!'); disp(alateStrTemp); disp(aptStrTemp); disp(sum(sum(A1_2))); disp(sum(sum(A2_2))); break; end
    % myCheck(A2_0, A2_1, A2_2, A2_3, A2aphids) 
    
    % Defining Events
    % We define 5 plant events
    tot_rate_vecP(1)=shosts*bS;              % mortality of susceptible variety hosts
    tot_rate_vecP(2)=ihosts*bI;              % mortality of infected hosts (assuming all virus strains have the same mortality)
    tot_rate_vecP(3)=rhosts*bR;              % mortality of resistant hosts
    tot_rate_vecP(4)=ihosts*al;              % additional infection mortality (assuming all virus strains have the same mortality)
    tot_rate_vecP(5)=ihosts*gam;             % recovery from infection (assuming all virus strains have the same mortality)
    
    if aphidType==1
        transformedFieldA1=ddTransform(rA1,double(fieldA1),double(fieldA1)+double(apterWgt*fieldA2),K1(counter));  % transform the field of alate densities into density dependent birth rates
    else
        transformedFieldA1=ddTransform(rA1,double(fieldA1),double(fieldA1),K1(counter));  % transform the field of alate densities into density dependent birth rates
    end
    transformedFieldA1=transformedFieldA1.*double((1-(transformedFieldA1<0)));
    if sum(sum(transformedFieldA1<0))>0
        disp('Sign A1 Error Break !!!!!');
    end
    if aphidType==1
        transformedFieldA2=ddTransform(rA2,double(fieldA2),double(fieldA1)+double(apterWgt*fieldA2),K2(counter));  % transform the field of apterae densities into density dependent birth rates
    else
        transformedFieldA2=ddTransform(rA2,double(fieldA2),double(fieldA2),K2(counter));  % transform the field of apterae densities into density dependent birth rates
    end
    transformedFieldA2=transformedFieldA2.*double((1-(transformedFieldA2<0)));
    if sum(sum(transformedFieldA2<0))>0
        disp('Sign A2 Error Break !!!!!');
    end
    % We define 4 alate events partita
    tot_rate_vecA1(1)=immA1;
    tot_rate_vecA1(2)=A1aphids*phiA1;              % rate of alate dispersal
    tot_rate_vecA1(3)=sum(sum(transformedFieldA1)); % figure out to the total rate of alate birth over all plants (which are individually density depdendent)
    tot_rate_vecA1(4)=A1aphids*bA1;              % rate of on-plant alate mortality
    % We define 4 aptarae events
    tot_rate_vecA2(1)=immA2;
    tot_rate_vecA2(2)=A2aphids*phiA2;             % rate of apterous dispersal
    tot_rate_vecA2(3)=sum(sum(transformedFieldA2)); % figure out to the total rate of aptarae birth over all plants (which are individually density depdendent)
    tot_rate_vecA2(4)=A2aphids*bA2;              % rate of on-plant apterous mortality
    tot_rate_vec=[tot_rate_vecP tot_rate_vecA1 tot_rate_vecA2];
    tot_rate=sum(tot_rate_vec);
    if sum(tot_rate<0)>0
        disp('Sign Error Break !!!!!');
    end
    %if tot_rate_vecA2(3)==0  % TESTING
        %disp('NEW ERROR !!!!!!!!!!!!!!!!'); disp(fieldA2); disp(rA2); disp(K2(counter)); disp(ddTransform(rA2,double(fieldA2),K2(counter))); 
    %end
    
    % Stochastic Jump and Stochastic Event Assign
    rand1=rand();
    timeto=-(1/tot_rate)*log(rand1);
    % jump to time
    t=t+timeto;
    timeVec(counter)=t;
   
    numDefinedRates=13; % TESTING
    if length(tot_rate_vec)~=numDefinedRates, disp('Error in rate sizings!!!'); break; end
    partita=cumsum(tot_rate_vec)./tot_rate;

    rand2=rand();
    if rand2<partita(1)                                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plant events %%%%%%%%%%%%%%%%%%
        dispO('natural susc death');  
        rand3=ceil(rand()*shosts);
        %if plantDeathDisperses
        %    [fieldP,fieldA1,fieldA2]=plantDeath(sites_s(rand3),fieldP,fieldA1,fieldA2);  % make a plant die and its aphids disperse
        %elseif plantDeathKillsAphids
        %    if rePlant, fieldP(sites_s(rand3))=1; else fieldP(sites_s(rand3))=0; end
        %    fieldA1(sites_s(rand3))=0;
        %    fieldA2(sites_s(rand3))=0;
        %else
        if rePlant, fieldP(sites_s(rand3))=1; else fieldP(sites_s(rand3))=0; end
        %end
    elseif rand2<partita(2)
        dispO('natural inf death');                                             % what happens to the plant's aphid pop when it dies? Disperse them...
        rand4=ceil(rand()*ihosts);
        %if plantDeathDisperses
        %    [fieldP,fieldA1,fieldA2]=plantDeath(sites_i(rand4),fieldP,fieldA1,fieldA2);
        %elseif plantDeathKillsAphids
        %    if rePlant, fieldP(sites_i(rand4))=1; else fieldP(sites_i(rand4))=0; end
        %    fieldA1(sites_i(rand4))=0;
        %    fieldA2(sites_i(rand4))=0;
        %else
        dispInv('Nat Inf death!')
        if rePlant, fieldP(sites_i(rand4))=1; else fieldP(sites_i(rand4))=0; end
        %end
    elseif rand2<partita(3)
        dispO('natural res death');
        rand5=ceil(rand()*rhosts);
        %if plantDeathDisperses
        %    [fieldP,fieldA1,fieldA2]=plantDeath(sites_r(rand5),fieldP,fieldA1,fieldA2);
        %elseif plantDeathKillsAphids
        %    if rePlant, fieldP(sites_r(rand5))=1; else fieldP(sites_r(rand5))=0; end
        %    fieldA1(sites_r(rand5))=0;
        %    fieldA2(sites_r(rand5))=0;
        %else
        if rePlant, fieldP(sites_r(rand5))=1; else fieldP(sites_r(rand5))=0; end
        %end
    elseif rand2<partita(4)
        dispO('disease death');
        rand6=ceil(rand()*ihosts);
        %if plantDeathDisperses
        %    [fieldP,fieldA1,fieldA2]=plantDeath(sites_i(rand6),fieldP,fieldA1,fieldA2);
        %elseif plantDeathKillsAphids
        %    if rePlant, fieldP(sites_i(rand6))=1; else fieldP(sites_i(rand6))=0; end
        %    fieldA1(sites_i(rand6))=0;
        %    fieldA2(sites_i(rand6))=0;
        %else
        dispInv('Disease death!')
        if rePlant, fieldP(sites_i(rand6))=1; else fieldP(sites_i(rand6))=0; end
        %end
    elseif rand2<partita(5)
        dispO('recovery');
        rand7=ceil(rand()*ihosts);
        fieldP(sites_i(rand7))=1;
        dispInv('Recovery!')
    elseif rand2<partita(6)      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Alate events %%%%%%%%%%%%%%%%%%
        dispR('alate immigration');
        out=alateProbingChain(double(fieldP), -99);   % -99 is code for dispersing immigrant (note it is possible that the immigrant might have the virus initially)
        
        seekVirEmm=find(out<0);% 1st up look for emigration/dispersal mortality of aphid, coded as a negative in the output field (-1000 was the marker)
        if isempty(seekVirEmm)~=1, inocProd(counter)=inocProd(counter)+1; out(seekVirEmm)=out(seekVirEmm)+1000; end
        
        outy=out;
        if any(any((mod(out,10)>=5)))==1, outy((mod(outy,10)>=5))=outy(mod(outy,10)>=5)-ones(sum(sum(mod(outy,10)>=5)),1)*5; end
        sorta=sort(find(mod(fieldP,10)==2));
        sortb=sort(find(mod(outy,10)==2));
        newInf=sortb(find(ismember(sortb,sorta)==0));
        if any(any((mod(out,10)>=5)))==1
            fieldA1(mod(out,10)>=5)=fieldA1(mod(out,10)>=5)+1; %updating the aphid count on the settled plant
        end
        if isempty(newInf)==0
            dispR(['New infections alate imm routine ' num2str(newInf(:)')]);
            fieldP(newInf)=outy(newInf);
        end
    elseif rand2<partita(7)
        dispR('alate dispersal');
        findField=find(rand()*(sum(fieldA1(:)))<cumsum(fieldA1(:)));
        indA1=findField(1);
        fieldA1(indA1)=fieldA1(indA1)-1;
        out=alateProbingChain(double(fieldP), indA1);
        
        seekVirEmm=find(out<0);% 1st up look for emigration/dispersal mortality of aphid, coded as a negative in the output field (-1000 was the marker)
        if isempty(seekVirEmm)~=1, inocProd(counter)=inocProd(counter)+1; out(seekVirEmm)=out(seekVirEmm)+1000; end
        
        blah=[blah find(mod(out,10)>=5)];
        outy=out;
        if any(any((mod(out,10)>=5)))==1, outy((mod(outy,10)>=5))=outy(mod(outy,10)>=5)-ones(sum(sum(mod(outy,10)>=5)),1)*5; end
        sorta=sort(find(mod(fieldP,10)==2));
        sortb=sort(find(mod(outy,10)==2));
        newInf=sortb(find(ismember(sortb,sorta)==0));
        if (R0run==1)&&(isempty(newInf)==0)
            dispInv('New Infections from single Alate Dispersal -----------------------------------------------')
            dispInv(newInf)
            dispInv('---------------------------------------------------------------------------------------------')
        end
        if any(any((mod(out,10)>=5)))==1
            fieldA1(mod(out,10)>=5)=fieldA1(mod(out,10)>=5)+1;
        end
        if isempty(newInf)==0
            dispR(['New infections alate disp rountine ' num2str(newInf(:)')]);
            if R0run==1
                newInfR0A1=[newInfR0A1 newInf'];
            else
                fieldP(newInf)=outy(newInf);
            end
            if mod(newInf,10)==1, countOne=countOne+1; end %%%%%%%%%%%%%%%% Testing variables, ignore !! %%%%%%%%%%%%%%%
            if mod(newInf,10)==2, countTwo=countTwo+1; end
            if mod(newInf,10)==3, countThree=countThree+1; end
            if mod(newInf,10)==4, countFour=countFour+1; end
            if mod(newInf,10)==5, countFive=countFive+1; end
            if mod(newInf,10)==6, countSix=countSix+1; end
            if mod(newInf,10)==7, countSeven=countSeven+1; end
            if mod(newInf,10)==8, countEight=countEight+1; end
            if mod(newInf,10)==9, countNine=countNine+1; end
            if mod(newInf,10)==0, countTen=countTen+1; end
        end
    elseif rand2<partita(8) % use transformed dd alate matrix to choose aphid index for event
        dispR('alate reproduces');
        findField=find(rand()*(sum(transformedFieldA1(:)))<cumsum(transformedFieldA1(:)));
        indA1=findField(1);
        if aphidType==1  % code 1 is normal mixed aphid dyamics
            morphIsAlate=morphDetermine2(fieldA1(indA1),fieldA2(indA1),[apterWgt apterSatConst expo]);
        elseif aphidType==2
            morphIsAlate=1;  % code 2 is alate only
        elseif aphidType==3
            morphIsAlate=0;  % code 3 is aptarae only
        else
            disp('Error! incorrect aphidType code');
            break;
        end
        if morphIsAlate, fieldA1(indA1)=fieldA1(indA1)+1; dispR('alate reproduction produces alate'); else fieldA2(indA1)=fieldA2(indA1)+1; dispR('alate reproduction prduces aptarae'); end
    elseif rand2<partita(9)
        dispR('alate mortality');
        findField=find(rand()*(sum(fieldA1(:)))<cumsum(fieldA1(:)));
        indA1=findField(1);
        fieldA1(indA1)=fieldA1(indA1)-1;
    elseif rand2<partita(10)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Aptarae events %%%%%%%%%%%%%%%%%%
        dispR('aptarae immigration');
        out=apterousProbingChain(double(fieldP), ceil(rand()*(k*k))); % instead of -99 just starting dispersal from a random plant, this seems quite reasonable
        outy=out;
        if any(any((mod(out,10)>=5)))==1, outy((mod(outy,10)>=5))=outy(mod(outy,10)>=5)-ones(sum(sum(mod(outy,10)>=5)),1)*5; end
        sorta=sort(find(mod(fieldP,10)==2));
        sortb=sort(find(mod(outy,10)==2));
        newInf=sortb(find(ismember(sortb,sorta)==0));
        if any(any((mod(out,10)>=5)))==1
            fieldA2(mod(out,10)>=5)=fieldA2(mod(out,10)>=5)+1;
        end
        if isempty(newInf)==0
            dispR(['New infections apt imm rountine ' num2str(newInf(:)')]);
            fieldP(newInf)=outy(newInf);
        end
    elseif rand2<partita(11)
        dispR('aptarae dispersal');
        findField=find(rand()*(sum(fieldA2(:)))<cumsum(fieldA2(:)));
        indA2=findField(1);
        fieldA2(indA2)=fieldA2(indA2)-1;
        out=apterousProbingChain(double(fieldP), indA2);
        blah=[blah find(mod(out,10)>=5)];
        outy=out;
        if any(any(mod(out,10)>=5))==1, outy((mod(outy,10)>=5))=outy(mod(outy,10)>=5)-ones(sum(sum(mod(outy,10)>=5)),1)*5; end
        sorta=sort(find(mod(fieldP,10)==2));
        sortb=sort(find(mod(outy,10)==2));
        newInf=sortb(find(ismember(sortb,sorta)==0));
        if (R0run==1)&&(isempty(newInf)==0)
            dispInv('New Infections from single Apterous Dispersal -----------------------------------------------')
            dispInv(newInf)
            dispInv('---------------------------------------------------------------------------------------------')
        end
        if any(any(mod(out,10)>=5))==1
            fieldA2(mod(out,10)>=5)=fieldA2(mod(out,10)>=5)+1;
        end
        if isempty(newInf)==0
            dispInv(['New infections Apterous disp rountine !!!!!!!!!' num2str(newInf(:)')]);
            if R0run==1
                newInfR0A2=[newInfR0A2 newInf'];
            else
                fieldP(newInf)=outy(newInf);
                tempCount=length(find(mod(fieldP,10)==2));
                dispInv(['Just update infected pop, new tot infs= ' num2str(tempCount)])
            end
            if mod(newInf,10)==1, countOne=countOne+1; end %%%%%%%%%%%%%%%% variables used for testing, ignore %%%%%%%%%%%%%%%%%%%%%%
            if mod(newInf,10)==2, countTwo=countTwo+1; end
            if mod(newInf,10)==3, countThree=countThree+1; end
            if mod(newInf,10)==4, countFour=countFour+1; end
            if mod(newInf,10)==5, countFive=countFive+1; end
            if mod(newInf,10)==6, countSix=countSix+1; end
            if mod(newInf,10)==7, countSeven=countSeven+1; end
            if mod(newInf,10)==8, countEight=countEight+1; end
            if mod(newInf,10)==9, countNine=countNine+1; end
            if mod(newInf,10)==0, countTen=countTen+1; end
        end
    elseif rand2<partita(12)
        dispR('aptarae reproduces');
        findField=find(rand()*(sum(transformedFieldA2(:)))<cumsum(transformedFieldA2(:)));
        indA2=findField(1);
        if aphidType==1  % code 1 is normal mixed aphid dyamics
            morphIsAlate=morphDetermine2(fieldA1(indA2),fieldA2(indA2),[apterWgt apterSatConst expo]);  % check the cell density for morph determination
        elseif aphidType==2
            morphIsAlate=1;  % code 2 is alate only
        elseif aphidType==3
            morphIsAlate=0;  % code 3 is aptarae only
        else
            disp('Error! incorrect aphidType code');
            break;
        end
        if morphIsAlate, fieldA1(indA2)=fieldA1(indA2)+1; dispR('aptarae reproduction prduces alate'); else fieldA2(indA2)=fieldA2(indA2)+1; dispR('aptarae reproduction prduces aptarae'); end
    elseif rand2<partita(13)
        dispR('aptarae mortality');
        findField=find(rand()*(sum(fieldA2(:)))<cumsum(fieldA2(:)));
        indA2=findField(1);
        fieldA2(indA2)=fieldA2(indA2)-1;
    else
        disp('Error Break !!!!!');
    end
    
    % Some Census Updates To Arrays Where Needed
    totalStrains(counter,:)=iStrainHosts;
    totalAlivePlants(counter)=ahosts;
    totalInfectedPlants(counter)=ihosts;
    totalAlates(counter)=A1aphids;
    totalAptarae(counter)=A2aphids;
    
    % find post update sites_i
    sites_i_post=find(mod(fieldP,10)==2);
    ihosts_post=length(sites_i_post);
    
    
    %if animate
    %    subplot(2,2,1);
    %    imagesc(fieldP);
    %    AX1=gca;
    %    ylabel('Plant Row','fontsize',10);
    %    title('Plants','fontsize', 12)
    %    subplot(2,2,2)
    %    plot(timeVec(1:sum(timeVec>=0)),totalInfectedPlants(1:sum(totalInfectedPlants>=0))./totalAlivePlants(1:sum(totalAlivePlants>=0)));
    %    hold on;
    %    plot(timeVec(1:sum(timeVec>=0)),totalAlates(1:sum(totalAlates>=0)));
    %    plot(timeVec(1:sum(timeVec>=0)),totalAptarae(1:sum(totalAptarae>=0)));
    %    AX4=gca;
    %    axis([0 seasonL 0 1]);
    %    ylabel('Incidence, alates, aptarae','fontsize',10);
    %    title('Time-Series','fontsize', 12)
    %    subplot(2,2,3)
    %    imagesc(fieldA1);
    %    AX2=gca;
    %    ylabel('Plant Row','fontsize',10);
    %    xlabel('Plant Column','fontsize',10);
    %    title('Alates','fontsize', 12)
    %    subplot(2,2,4)
    %    imagesc(fieldA2);
    %    AX3=gca;
    %    ylabel('Plant Row','fontsize',10);
    %    xlabel('Plant Column','fontsize',10);
    %    title('Aptarae','fontsize', 12)
    %    shading interp;
    %    set(gcf,'color', 'w');
    %    set([AX1 AX2 AX3], ...
    %        'ytick', 0.5:k+1,  ...
    %        'xtick',  0.5:k+1, ...
    %        'yticklabel', [],   ...
    %        'xticklabel', [],   ...
    %        'fontweight'  ,'n'       , ...
    %        'FontSize'    , 9       , ...
    %        'Box'         , 'off'     , ...
    %        'TickDir'     , 'out'     , ...
    %        'TickLength'  , [.05 .05] , ...
    %        'XMinorTick'  , 'off'      , ...
    %        'YMinorTick'  , 'off'      , ...
    %        'YGrid'       , 'off'      , ...
    %        'XColor'      , [.3 .3 .3], ...
    %        'YColor'      , [.3 .3 .3], ...
    %        'LineWidth'   , 0.2         );
    %    MovA(indx) = getframe(gcf);    % put frame
    %    set(gcf, 'Units','centimeters', 'Position',[5 3 18*2 7*3])
    %end
    
    
    %if and(and(sampleField==1,mod(t,snapFileFreq)<0.2),markerSnaps==0)
    %    disp(['Outputting snapshot t= ' num2str(t)]);
    %    markerSnaps=1;
    %    % VECTORISE Qclumping (just here as expensive) %%%
    %    strainClumpTemp=zeros(1,strains);
    
    % IF dryrun check equilibrium has been reached
    if fullTimeIn==(-1)
        if mod(counter,1000)==0
            if or(or(aphidType==1,aphidType==2),aphidType==3) % then mixed... so check alate has reached steady-state
                totAphids=totalAlates+totalAptarae;
                if totAphids(counter)>0
                    if countSS==0
                        disp('...Calculating first avr alates over 1000 steps...')
                        movingAphidAvr=mean(totAphids(counter-999:counter));
                        countSS=1;
                    else
                        newAphidAvr=mean(totAphids(counter-999:counter));
                        if (abs(newAphidAvr-movingAphidAvr)<(tolAphidAvr*movingAphidAvr))||(abs(newAphidAvr-movingAphidAvr)<margAphidAvr)
                            disp(['...Calculating (' num2str(countSS) ') avr alates over 1000 steps...'])
                            countSS=countSS+1;
                            if countSS>=threshSS
                                disp('Alate Run: Steady-State Criterion met!!!!')
                                break;
                            end
                        else
                            countSS=0;
                        end
                    end
                end
            elseif aphidType==2 % then alate
                if countSS==0
                    movingAphidAvr=mean(totAphids(counter-999:counter));
                    countSS=1;
                else
                    newAphidAvr=mean(totAphids(counter-999:counter));
                    if (abs(newAphidAvr-movingAphidAvr)<(tolAphidAvr*movingAphidAvr))||(abs(newAphidAvr-movingAphidAvr)<margAphidAvr)
                        countSS=countSS+1;
                        if countSS>=threshSS
                            disp('Alate Run: Steady-State Criterion met!!!!')
                            break;
                        end
                    else
                        countSS=0;
                    end
                end
            elseif aphidType==3 % then aptarae
                if countSS==0
                    movingAphidAvr=mean(totalAptarae(counter-999:counter));
                    countSS=1;
                else
                    newAphidAvr=mean(totalAptarae(counter-999:counter));
                    if (abs(newAphidAvr-movingAphidAvr)<(tolAphidAvr*movingAphidAvr))||(abs(newAphidAvr-movingAphidAvr)<margAphidAvr)
                        countSS=countSS+1;
                        if countSS>=threshSS
                            disp('Aptarae Run: Steady-State Criterion met!!!!')
                            break;
                        end
                    else
                        countSS=0;
                    end
                end
            end
        end
    end
    
    
    
    if t>timeFilt(zNew)
        disp(t);
        zNew=zNew+1;
        if outputClump        %%%%%%%%%%%%%%%%%%%%% Recording clumping %%%%%%%%%%%%%
            for i=1:strains, strainClumpTemp(i)=qClump(fieldP, (2+(i-1)*10), (2+(i-1)*10)); end % this is the conditional density rather than pair density
            strainClump(counter,:)=strainClumpTemp;
        end
        
    end
    if outputSnapshot
        if or((ihosts_post>=1)*(1-flagI1),(ihosts_post>=k*k/2)*(1-flagI2))
            if (ihosts_post>=1)*(1-flagI1), flagI1=1; end
            if (ihosts_post>=k*k/2)*(1-flagI2), flagI2=1; end
            
            disp(['Outputting snapshot @ t=' num2str(timeFilt(zNew-1)) 'I' num2str(length(ihosts_post))])
            outArray=[fieldP;fieldA1;fieldA2];
            snapFile=[inFolder 'SnapshotI' num2str(length(sites_i)) '.csv'];
            fido = fopen(snapFile,'w');
            arrSz=size(outArray);
            numRows=arrSz(1); numCols=arrSz(2);
            fprintf(fido,',');
            for l=1:numCols
                if l~=numCols
                    fprintf(fido,'%s', ['Col_' num2str(l)]);
                    fprintf(fido,',');
                else
                    fprintf(fido,'%s\n', ['Col_' num2str(l)]);
                end
            end
            for j=1:numRows
                if ceil(j/k)==1
                    fprintf(fido,'%s', ['plant_Type_Row_' num2str(mod(j,(k+1)))]);
                elseif ceil(j/k)==2
                    fprintf(fido,'%s', ['alate_Density_Row_' num2str(j-k)]);
                else
                    fprintf(fido,'%s', ['aptarae_Density_Row_' num2str(j-2*k)]);
                end
                fprintf(fido,',');
                for i=1:numCols-1
                    fprintf(fido,'%s', num2str(outArray(j,i)));
                    fprintf(fido,',');
                end
                fprintf(fido,'%s\n', num2str(outArray(j,end)));
            end
            fclose(fido);
        end
        
    end
end
% END OF GILLESPIE WHILE LOOP
fMis=find(abs(inocProd)>1); % Testing: 'inocProd' should be 1's, 0's or -1's so test for other values
if isempty(fMis)~=1, disp('ERROR !!!!! abs(inocProd) was non-binary!!!!!'); return; end
inocProd=cumsum(inocProd);
% Locating Positive Indices From The Block Allocation
K1pos=find(K1>=0);
K2pos=find(K2>=0);
simInds=max(K1pos(end),K2pos(end));
totalAlivePlants=totalAlivePlants(1:simInds);
timeVec=timeVec(1:simInds);
sStrainHosts=sStrainHosts(1:simInds);
totalStrains=totalStrains(1:simInds,:);
totalAlates=totalAlates(1:simInds);
totalAptarae=totalAptarae(1:simInds);
inocProd=inocProd(1:simInds);
alatesOnSusc=alatesOnSusc(1:simInds);
aptOnSusc=aptOnSusc(1:simInds);
alateStrains=alateStrains(1:simInds,:); % currently in terms of total # aphids... to get per plant average: ./ totalStrains
aptStrains=aptStrains(1:simInds,:);     % ditto...
strainClump=strainClump(1:simInds,:);
K1=K1(1:simInds);
K2=K2(1:simInds);

% Plotting The Simulation Outputs
%if animate==0
%[~,posInfsC]=find(totalStrains>0);
%posInfs=min(posInfsC);
%    seed=1000000;
%    for l=1:strains
%        seed=min(min([(find(totalStrains(:,l)>0)') seed]));
%    end
%    figure;
%    strainMatrix=(totalStrains(seed:end,:)./(repmat(totalAlivePlants(seed:end)',1,strains)))';    % i.e. incidence
%    sVec=sStrainHosts(seed:end)./totalAlivePlants(seed:end);
%    plot(timeVec(seed:end),strainMatrix);
%    if sum(sum(strainMatrix))>0
%        for j=1:strains, text(0.95*timeVec(end),0.01+1.025*strainMatrix(j,end),['Strain #' num2str(j)]); end
%        hold on;
%        plot(timeVec(seed:end),sum(strainMatrix,1),'linewidth',2);
%        hold on;
%        plot(timeVec(seed:end),sVec,'r--')
%        title('strainMatrix');   % plotting time series for each virus strain
%        cstring3=cstring;
%        cstring3{end+1}='S';
%        legend(cstring3);
%    end
%    if aphidType~=3
%        figure;
%        plot(timeVec,totalAlates);
%        hold on;
%        plot(inocTime*ones(1,2),[0 max(totalAlates)*2],'r:')
%        text(inocTime,max(totalAlates)*0.25,'Epi. Seed');
%        title('totalAlates');    % plotting alate time series
%    end
%    if aphidType~=2
%        figure;
%        plot(timeVec,totalAptarae);
%        hold on;
%        plot(inocTime*ones(1,2),[0 max(totalAptarae)*2],'r:')
%        text(inocTime,max(totalAptarae)*0.25,'Epi. Seed');
%        title('totalAptarae');   % plotting aptarae time series
%        plot(timeVec,aptStrains)
%        hold on;
%        size(timeVec)
%        size(totalAptarae')
%        size(aptStrains)
%        plot(timeVec,totalAptarae'-sum(aptStrains,2))
%    end
%    if or(k1_fluc>0,k2_fluc>0)
%        figure;
%        plot(timeVec,K1);        % plotting alate carrying capacity (per plant)
%        title('K1');
%        figure;
%        plot(timeVec,K2);
%        title('K2');             % plotting aptarae carrying capacity (per plant)
%    end
%    figure;
%    if sum(sum(strainMatrix))>0
%        strainClumpOut=strainClump(seed:end,:)';
%        timeVecOut=timeVec(seed:end);
%        %strainClumpOutInds=find(strainClumpOut>0);
%        plot(timeVecOut,strainClumpOut);
%        title('strainClump');    % plotting the pair density for a particular type infection type ! note convert to 'relatedness'
%        legend(cstring);
%        posSusc=find(sStrainHosts>0);
%        if aphidType~=3
%            figure;
%            alatesPerStrainPl=alateStrains./totalStrains;   % i.e. / # infected plants of each strain
%            plot(timeVec(seed:end),alatesPerStrainPl(seed:end,:)')
%            hold on;
%            plot(timeVec(posSusc),alatesOnSusc(posSusc)./sStrainHosts(posSusc),'--b')
%            title('alatesPerStrainPl'); % plotting the average number of alates on plants with a particular infection type
%        end
%        if aphidType~=2
%            figure;
%            aptPerStrainPl=aptStrains./totalStrains;
%            plot(timeVec((seed:end)),aptPerStrainPl(seed:end,:)')
%            %for j=1:strains, text(0.95*timeVec(end),0.01+1.025*aptPerStrainPl(j,end),['Strain #' num2str(j)]); end
%            hold on;
%            plot(timeVec(posSusc),aptOnSusc(posSusc)./sStrainHosts(posSusc),'r--')
%            title('aptaraePerStrainPl'); % plotting the average number of aptarae on plants with a particular infection type
%            cstring2=cstring;
%            cstring2{end+1}='S';
%            legend(cstring2);
%        end
%    end
%end

% Some Final Consistency Checks
disp(['Checking that there are no alates on dead plants, total such aphids:' num2str(sum(fieldA1(fieldP==0)))])
disp(['Checking that there are no aptarae on dead plants, total such aphids:' num2str(sum(fieldA2(fieldP==0)))])


if (fullTimeIn==(-1))*aphidType*plotBurnIn==1
    figure; 
    plot(timeVec,totalAlates); 
    hold on; 
    plot(timeVec,totalAptarae);
    legend('Alates','Aptarae')
elseif (fullTimeIn==(-1))*aphidType*plotBurnIn==2
    figure; plot(timeVec,totalAlates);
elseif (fullTimeIn==(-1))*aphidType*plotBurnIn==3
    disp('time'); disp(t); disp('size(timeVec)'); disp(size(timeVec)); disp('timeVec(end)'); disp(timeVec(end));
    figure; plot(timeVec,totalAptarae);
end
% Animation Outputs
%if animate, close(gcf); movie2avi(MovA, 'C:\Users\donnelr3\OneDrive - University Of Cambridge\Archive_I_Beans\deleteThisVideo.avi'); end

try
    if fullTimeIn==(-1)
        % overwrite filename
        disp(['aphid type is ' num2str(aphidType)])
        outItFile=[inFolder 'SSaphidDistribution.csv'];
        if aphidType==2
            outMat=fieldA1;
        elseif aphidType==3
            outMat=fieldA2;
        else
            outMat=[fieldA1; fieldA2];
        end
        disp(outItFile)
        disp(outMat)
        disp('HERE! Line 917 after print output array but before save')
        disp('matrix type?')
        class(outMat)
        dlmwrite(outItFile,outMat);
        disp('after save')
    elseif outputMetrics>0
        disp('Warning !!!!!! in outputMetrics')
        % Formatting Simulation Output For Data Save
        outMat=zeros(length(timeFilt)-1,1+strains+1+1+1+1+1+1+1+strains,'double');
        lastTime=find(timeFilt<t);
        disp(t)
        disp(timeFilt(lastTime(end)))
        for z=1:(lastTime(end)-1)
            tmpFnd=find(timeVec>timeFilt(z));
            outMat(z,:)=[timeVec(tmpFnd(1)) double(totalStrains(tmpFnd(1),:)) double(totalAlates(tmpFnd(1))) double(totalAptarae(tmpFnd(1))) double(alatesOnSusc(tmpFnd(1))) double(aptOnSusc(tmpFnd(1))) double(inocProd(tmpFnd(1))) K1(tmpFnd(1)) K2(tmpFnd(1)) strainClump(tmpFnd(1),:)];
        end
        disp(outputFile)
        disp(outMat)
        dlmwrite(outputFile,outMat);
        disp('Trying within outputmetrics!!! line 942')
    end
catch
    disp(totalStrains((end-10):end,:))
    disp(size((repmat(totalAlivePlants',1,strains))));
    disp('Catch!!!!!!!!!!!!!!!!')
    disp('print cd')
    currentFolder = pwd
end
totA1=totalAlates(end);
totA2=totalAptarae(end);
disp(['sum(sum(fieldA1))' num2str(totA1)])
disp(['sum(sum(fieldA2))' num2str(totA2)])

disp(['fieldA1(1,:)' num2str(fieldA1(1,1:10))])
disp(['fieldA2(1,:)' num2str(fieldA2(1,1:10))])

if printCounts==1   %%%%%%%%%%%%%%%%%%%%%%  For testing, ignore %%%%%%%%%%%
    disp(['countOne = ' num2str(countOne)]);
    disp(['countTwo = ' num2str(countTwo)]);
    disp(['countThree = ' num2str(countThree)]);
    disp(['countFour = ' num2str(countFour)]);
    disp(['countFive = ' num2str(countFive)]);
    disp(['countSix = ' num2str(countSix)]);
    disp(['countSeven = ' num2str(countSeven)]);
    disp(['countEight = ' num2str(countEight)]);
    disp(['countNine = ' num2str(countNine)]);
    disp(['countTen = ' num2str(countTen)]);
    temp=mod(blah,10);
    temp1=find(temp==1);
    temp2=find(temp==2);
    temp3=find(temp==3);
    temp4=find(temp==4);
    temp5=find(temp==5);
    temp6=find(temp==6);
    temp7=find(temp==7);
    temp8=find(temp==8);
    temp9=find(temp==9);
    temp10=find(temp==0);
    disp(['blah 1 = ' num2str(length(temp1))]);
    disp(['blah 2 = ' num2str(length(temp2))]);
    disp(['blah 3 = ' num2str(length(temp3))]);
    disp(['blah 4 = ' num2str(length(temp4))]);
    disp(['blah 5 = ' num2str(length(temp5))]);
    disp(['blah 6 = ' num2str(length(temp6))]);
    disp(['blah 7 = ' num2str(length(temp7))]);
    disp(['blah 8 = ' num2str(length(temp8))]);
    disp(['blah 9 = ' num2str(length(temp9))]);
    disp(['blah 10 = ' num2str(length(temp10))]);
end

%dlmwrite(outputFile,outMat','-append');
%end
%csvwrite(outputFile,outMat','-append');
%myVideo = VideoWriter('C:\Users\donnelr3\OneDrive - University Of Cambridge\Archive_I_Cambridge\deleteThisVideo.avi');
%myVideo.FrameRate = 15;  % Default 30
%myVideo.Quality = 50;    % Default 75
%open(myVideo);
%movie(Mov,1,1000);
%close(gcf);
%writeVideo(Mov);
%close(myVideo);

finalInfCount=length(find(mod(fieldP,10)==2));

clearvars -except extinct totA1 totA2 aphidType R0run newInfR0A1 newInfR0A2 finalInfCount t
disp(['WithinRoutinePrint' num2str(extinct) 'timeIS' num2str(t) 'checkInf' num2str(finalInfCount)])
disp(['WithinRoutinePrint_new_infectionsA1' num2str(newInfR0A1)])
disp(['WithinRoutinePrint_new_infectionsA2' num2str(newInfR0A2)])
if aphidType==1
    if R0run==1
        returnSet=[length(unique(newInfR0A1)) length(unique(newInfR0A2))];
    else
        returnSet=[extinct totA1 totA2];
    end
elseif aphidType==2
    if R0run==1
        returnSet=length(unique(newInfR0A1));
    else
        returnSet=[extinct totA1];
    end
elseif aphidType==3
    if R0run==1
        returnSet=length(unique(newInfR0A2));
    else
        returnSet=[extinct totA2];
    end
end
toc