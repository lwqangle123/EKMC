%% load SUMO data
clear
fidin=fopen('outloop1.txt');
data1=textscan(fidin, '%s');
[m,n]=size(data1{1});
for i=1:m/11
Ocp1{i}=data1{1,1}{7+(i-1)*11,1};
end
D1=[];
for i=1:size(Ocp1,2)
if strcmp('occupancy="0.00"',Ocp1{i})
D1(i)=0;
else
D1(i)=1;
end
end
Det1=zeros(32,floor(size(D1,2)/32));
for i=1:size(Det1,2)
Det1(:,i)=D1(1+(i-1)*32:i*32);
end
fidin=fopen('outloop2.txt');
data2=textscan(fidin, '%s');
[m,n]=size(data2{1});
for i=1:m/11
Ocp2{i}=data2{1,1}{7+(i-1)*11,1};
end
D2=[];
for i=1:size(Ocp2,2)
if strcmp('occupancy="0.00"',Ocp2{i})
D2(i)=0;
else
D2(i)=1;
end
end
Det2=zeros(32,floor(size(D2,2)/32));
for i=1:size(Det2,2)
Det2(:,i)=D2(1+(i-1)*32:i*32);
end
Dt1=norm(Det1-Det2); % see the difference between two day data;
m=[];
n=[];
% add random '0' for diversity
for i=1:32
    [m,n]=find(Det1(i,:)==1);
    rd=randperm(length(n),round(length(n)/8));
    Det1(i,n(rd))=0;
end
Dt2=norm(Det1-Det2);
%% Occupancy analysis
for j=1:24
    for i=1:32
    Det1s(i,j)=sum(Det1(i,1+3600*(j-1):j*3600))/3600;
    Det2s(i,j)=sum(Det2(i,1+3600*(j-1):j*3600))/3600;
    end
end
Dets=[Det2s,Det1s];


