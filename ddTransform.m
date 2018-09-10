function transField = ddTransform(R,reproField,totalField,K)
% Function to read in a matrix where each entry represents number of aphids
% on corresponing plant in spatial array
Kadj=(1-(K<0))*K;
Kadj2=Kadj+(Kadj==0)*0.000001;
transField=R*reproField.*(1-(totalField/Kadj2));
transField=double((1-(transField<0))).*transField;
end

