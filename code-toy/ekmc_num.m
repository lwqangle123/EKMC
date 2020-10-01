clear
load('matrix_hr.mat');

Ind=1:15;
w=ones(10,1);
r=1;
lambda=1e-6;
max_num=6;
rk=1;
thr=0.01

for k=1:max_num

     XX(1:10,1:9)=diag(w)*XX(1:10,1:9);    
     M_s=ones(15,10);
     M_s(11:15,10)=0;
     Xtr=XX(1:10,1:9)';
     Ytr=XX(1:10,10)';
     Xte=XX(11:15,1:9)';
     K1=rbf1(Xtr,Xtr,Ind(1:10),Ind(1:10),r);
     K2=rbf1(Xtr,Xte,Ind(1:10),Ind(11:15),r);
     K3=rbf1(Xte,Xtr,Ind(11:15),Ind(1:10),r);
     K4=rbf1(Xte,Xte,Ind(11:15),Ind(11:15),r);
     
    [Ytr_p,Yte_p,Ytrp_i,Ytep_i,Rate,Tau]=kmc(K1,K2,K3,K4,Xtr,Ytr,Xte,rk,lambda);

    for kk=1:size(Ytr,2)
        if abs(Ytr(kk)-Ytr_p(kk))<thr
           er(kk)=0;
        else
           er(kk)=1;
        end
        %er(kk)=w(kk)*abs(XX(kk,10)-Ytrn_p(kk));
    end
    erro=sum(er)/sum(w);
    if erro==0
       break
    end    
%     if k>2&abs(sum(Ytr_p-Ytrn_pf(:,k-1)'))<0.011
%         break
%     end
    Bta(k)=log(abs((1-erro)/erro));
    for kk=1:size(Xtr,2)
        w(kk)=w(kk)*exp(Bta(k)*abs(XX(kk,10)-Ytr_p(kk)));
    end
    id=isnan(w);
    w(id)=max(w);
    Ytst_pf(:,k)=Yte_p;
    Ytrn_pf(:,k)=Ytr_p;
    Tau_F(k)=Tau;
    W(:,k)=w;
end
%%

Alp=Bta/sum(Bta);
Tauf=sum(Alp.*Tau_F);
YTEPF=Alp(1)*Ytst_pf(:,1);
for j=1:k
    YTEPF(:,j)=Alp(j)*Ytst_pf(:,j);
    YTRPF(:,j)=Alp(j)*Ytrn_pf(:,j);
end
%YTEPF=alpha.*Ytst_pf;
