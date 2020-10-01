
function K = rbf1(coord1,coord2,index1,index2,Sig)
    n1=size(coord1,2);
    n2=size(coord2,2);
    n3=size(coord1,1);
    p=2*pi/90;
    for i=1:n1
        for j=1:n2
           dp(i,j)=min(abs(mod(index1(i)-index2(j),90)),90-abs(mod(index1(i)-index2(j),90)));
           K2(i,j)=exp(-dp(i,j)^2);
           K1(i,j)=exp(-norm(coord1(:,i)-coord2(:,j))^2/((Sig)^2));
           %K2(i,j)=alpha^(n3-norm(coord1(:,i)-coord2(:,j))^2)*(1-alpha)^(norm(coord1(:,i)-coord2(:,j))^2);
           %K3(i,j)=exp((sin(norm(coord1(:,i)-coord2(:,j))/p)-1)/(Sig)^2);
           K(i,j)=K2(i,j)*K1(i,j);
        end
    end   
end