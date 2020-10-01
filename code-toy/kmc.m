function [Ytr_p,Yte_p,Ytrp_i,Ytep_i,Rate,Tau]=kmc(K1,K2,K3,K4,Xtr,Ytr,Xte,rk,lambda)
Rate=[];
ER=[];
ERR=[];
OBJ=[];
Ycom=[];
Ycom_s=[];
Err=ones(1000,1);
[Ytr_p,Yte_p,Ytrp_i,Ytep_i,Utr,Vtr,Vte,Err,Obj] = nlmc_bcd(K1,K2,K3,K4,Xtr,Ytr,Xte,rk,lambda);
Ycom=Ytr_p';
Ycom_s=sort(Ycom,'descend');
dY=diff(Ycom_s);
for j=1:size(dY,2)
    [id,fv(j)]=max(abs(dY(:,j)));
    Tau(j)=Ycom_s(fv(j)+1,j);
end
Rate=Err;
end