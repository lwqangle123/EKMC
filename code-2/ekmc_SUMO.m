%% construct the independent variable and dependent variable
L=60; % lag L=30, 60, 120
PH=10; % horizon, PH=1,10,60
X1=[];
X2=[];
y1=[];
y2=[];
for ii=1:(86400-L-PH+1)
    for k1=1:L
        X1(1+(k1-1)*32:k1*32,ii)=Det1(:,ii-1+k1);
        X2(1+(k1-1)*32:k1*32,ii)=Det2(:,ii-1+k1);
    end
    y1(:,ii)=Det1(:,ii+L+PH-1);
    y2(:,ii)=Det2(:,ii+L+PH-1);
end
X1=X1';
X2=X2';
y1=y1';
y2=y2';
XY1=[X2,y2];
XY2=[X1,y1];
Ind=1:86361;
%% ensemble KMC
i=1:120;
Ytst_pf=[];
Ytepf=[];
Ytrpf=[];
Ytru=[];
Yth=[];
Alp=[];
error=[];
Tauf=[];
w=ones(1140,1);
r=5; % rbf paremeter
lambda=0.4; 
max_num=50;
rk=15;
thr=0.05;
for ii=1:100
    for j=1:32
        Xx1=X1(1+600*(i(ii)-1):i(ii)*600,:);
        Yy1=y1(1+600*(i(ii)-1):i(ii)*600,j);
        Xx2=X2(1+600*(i(ii)-1):i(ii)*600,:);
        Yy2=y2(1+600*(i(ii)-1):i(ii)*600,j);
        XX1=[Xx1,Yy1];
        XX2=[Xx2,Yy2];
        XX=[XX2;XX1];
        for k=1:max_num         
            XX(1:1140,1:1920)=diag(w)*XX(1:1140,1:1920);    
            M_s=ones(1200,size(XX,2));
            M_s(1141:1200,1921:end)=0;
            Xtr=XX(1:1140,1:1920)';
            Ytr=XX(1:1140,1921:end)';
            Xte=XX(1141:1200,1:1920)';
            K1=rbf1(Xtr,Xtr,[Ind(1:600),Ind(1:540)],[Ind(1:600),Ind(1:540)],r);
            K2=rbf1(Xtr,Xte,[Ind(1:600),Ind(1:540)],Ind(541:600),r);
            K3=rbf1(Xte,Xtr,Ind(541:600),[Ind(1:600),Ind(1:540)],r);
            K4=rbf1(Xte,Xte,Ind(541:600),Ind(541:600),r);
            [Ytr_p,Yte_p,Tau]=kmc1(K1,K2,K3,K4,Xtr,Ytr,Xte,rk,lambda);
            Tauf(k,j,ii)=Tau;
            Ytepf{k}(j,:,ii)=Yte_p;
            Ytrpf{k}(j,:,ii)=Ytr_p;
            Ytru{k}(j,:,ii)=XX(1141:1200,(J+1):end);
            for kk=1:size(Ytru{k}(j,:,ii),2)
                if abs(Ytru{k}(j,kk,ii)-Ytr_p(kk))<thr
                   er(kk)=0;
                else
                   er(kk)=1;
                end
            end
            erro(k,j,ii)=sum(er)/sum(w);
            if erro(k,j,ii)==0
               break
            end    
            if k>2&abs(sum(Ytrpf{k}(j,:,ii)-Ytrpf{k-1}(j,:,ii)))<0.01
               break
            end
            Bta(k,j,ii)=log(abs(1-erro(k,j,ii))/erro(k,j,ii));  
            for kk=1:size(Xtr,2)
                w(kk)=w(kk)*exp(Bta(k,j,ii)*abs(Ytr(kk)-Ytr_p(kk)));
            end
            id=isnan(w);
            w(id)=max(w);
            Alp(k,j,ii)=Bta(k,j,ii);
        end
        for kk=1:k-1
            Alp(kk,j,ii)=Alp(kk,j,ii)/sum(Alp(:,j,ii));
        end
    end
    for j=1:32
        YTEPF(j,:,ii)=Alp(1,j,ii)*Ytepf{1}(j,:,ii);
        YTRPF(j,:,ii)=Alp(1,j,ii)*Ytrpf{1}(j,:,ii);
        TAUF(j,ii)=Alp(1,j,ii)*Tauf(1,j,ii);        
        for kk=2:k-1
            YTEPF(j,:,ii)=YTEPF(j,:,ii)+Alp(kk,j,ii)*Ytepf{kk}(j,:,ii);
            YTRPF(j,:,ii)=YTRPF(j,:,ii)+Alp(kk,j,ii)*Ytrpf{kk}(j,:,ii);
            TAUF(j,ii)=TAUF(j,ii)+Alp(kk,j,ii)*Tauf(kk,j,ii);
        end
        Yth(j,:,ii)=double(YTEPF(j,:,ii)>TAUF(j,ii));
    end    
end
