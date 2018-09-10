function boolAlate = morphDetermine2(alatesOnPlant,aptaraeOnPlant,cons)

aptWgt=cons(1);
morphSatCons=cons(2);
expon=cons(3);

a=(alatesOnPlant+aptWgt*aptaraeOnPlant)^expon;
b=(alatesOnPlant+aptWgt*aptaraeOnPlant)^expon;
c=(morphSatCons^expon);

Palate=double(a)/(double(b)+double(c));

if Palate>0.999, disp(['alatesOnPlant aptWgt aptaraeOnPlant expon morphSatCons a b c Palate']); disp([alatesOnPlant aptWgt aptaraeOnPlant expon morphSatCons a b c Palate]); end
randnew=rand();
if randnew<Palate, boolAlate=1; else boolAlate=0; end

end

