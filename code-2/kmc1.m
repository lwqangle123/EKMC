function [Ytr_p,Yte_p,Tau]=kmc1(K1,K2,K3,K4,Xtr,Ytr,Xte,rk,lambda)
ER=[];
ERR=[];
OBJ=[];
Ycom=[];
Ycom_s=[];
Err=[];
[Ytr_p,Yte_p,Ytrp_i,Ytep_i,Utr,Vtr,Vte,Err,Obj] = nlmc_bcd(K1,K2,K3,K4,Xtr,Ytr,Xte,1,lambda);
 Ycom=Ytr_p';
 Ycom_s=sort(Ycom,'descend');
 dY=diff(Ycom_s);
 for j=1:size(dY,2)
     [id,fv]=max(abs(dY(:,j)));
     Tau(j)=Ycom_s(fv+1,j);
 end
