clearvars

feedVec=1/5; 

epsVec=[1/4 1/3.5] % 1/3 0.41667 0.45833 0.5 0.54167 2/3 1 3/2 1.8462 2 1/0.45833 1/0.41667 3 3.5 4];
nuFac=[0.5 1 2]; 
seedVec=1:length(epsVec);

bSvec=1/20; 
reps=1;

mutRate=0; % i.e., by default simulations do not feature evolution, but code is extendable for evolution/multi-pathogen strain dynamics

%%%%% 1-Normal (alate and apterous);2-alate only ; 3-aptarae only
fieldDims=20;
firstRun=0;
for aphidType=1
    baseInT=30*2; 
    cd('C:\Users\donnelr3\Desktop\Files4GitHub\')
    newName1=['Sep10_p02_ConcInv\'];
    inT=baseInT;
    
    if firstRun
        mkdir(newName1)
        numTpoints=inT*20;
        Sto_spatial(-1,0,feedVec,[],newName1, fieldDims, [nuFac epsVec 6], numTpoints,bSvec, 0.5,mutRate,aphidType, 0,0)
        firstRun=0;
    end
    
    for i=1:length(nuFac)
        for j=1:length(seedVec) 
            nuVec=nuFac(i)*ones(1,length(epsVec));
            if seedVec(j)==1, inT=inT; end
            clearvars -except nuVec feedVec epsVec seed reps i j bSvec bb newName1 mutRate aphidType seedVec nuFac burnIn inT fieldDims 
            numTpoints=inT*20;
            cd('C:\Users\donnelr3\Desktop\Files4GitHub\')
            newName2=['nu' num2str(nuFac(i)) 'eps' num2str(epsVec(seedVec(j)))];
            for s=(1:reps)
                disp(['Rep # is: ' num2str(s)])
                newName3=['rep' num2str(s)];
                mkdir([newName1 newName2 '\' newName3])
            end
            
            parfor k=(1:reps)
                disp(['Rep # is: ' num2str(k)])
                newName4=['rep' num2str(k)];
                ff=[newName1 newName2 '\' newName4 '\'];
                disp(ff)
                disp(ff(end-6)) % inT, 20
                Sto_spatial(inT+1,1,feedVec,newName1,ff, fieldDims, [nuVec epsVec seedVec(j)], numTpoints, bSvec, 0.5, mutRate,aphidType,1,0)
            end            
        end
    end
end
