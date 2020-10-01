function [Ytr_p,Yte_p,Utr,Vtr,Vte,Err,Obj] = nlmc_bcd(K1,K2,K3,K4,Xtr,Ytr,Xte,r,lambda)

[I1,J1]=size(Xtr);
[I2,J2]=size(Ytr);
[I3,J3]=size(Xte);
Err=[];

Vtr=rand(J2,r)-0.5;
Vte=rand(J3,r)-0.5;
Utr=Ytr*Vtr*pinv(Vtr'*Vtr+2*lambda*eye(r,r));
LE=lambda*eye(r,r);

for iter=1:5000
    Utr_n=Ytr*Vtr*inv(Vtr'*Vtr+2*LE);
    err1=norm(Utr_n-Utr)/norm(Utr);
    Utr=Utr_n;
    C1=inv(Vtr'*Vtr+Vte'*Vte+2*LE); % r times r
    C2=Vtr'*K1*Vtr+Vte'*K3*Vtr; % r times r
    C3=Vtr'*K2*Vte+Vte'*K4*Vte; % r times r
    C4=inv(C1*(C2+C3)*C1+Utr'*Utr+2*LE); %r times r
    Vtr_n=(Ytr'*Utr+(K1*Vtr+K2*Vte)*C1)*C4;
    err2=norm(Vtr_n-Vtr)/norm(Vtr);
    Vtr=Vtr_n;
    Vte_n=(K3*Vtr+K4*Vte)*C1*(C1*(C2+C3)*C1+2*LE);
    err3=norm(Vte_n-Vte)/norm(Vte);
    Vte=Vte_n;
    err=max(err1,err2);
    err=max(err,err3);
    Err(iter)=err;
    Obj(iter)=norm(Ytr-Utr*Vtr')/norm(Ytr);
    if err<1e-4
        break;
    end

end
Ytr_p=Utr*Vtr';
Yte_p=Utr*Vte';
end
