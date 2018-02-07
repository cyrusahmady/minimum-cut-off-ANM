# minimum-cut-off-ANM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       This program  makes  rmsd of different ANM  comparablr to ptools rmsca      written by cyrus                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

KB=1.3806488;T=297;KBT=KB*T;A=8*(pi^2)/3;M=12;KBT_M=KBT/M;
ani1=input('type name of protein *.pdb  for the unbound protein : \n','s');
name1=sprintf('%s%s',ani1,'ca.pdb');
fid1=fopen(name1);
A1=textscan(fid1,'%s');
fclose(fid1);

ani2=input('type name of protein *.pdb  for the bound protein : \n','s');
name2=sprintf('%s%s',ani2,'ca.pdb');
%Nmode=input('please input the number of modes to be considered: \n');
tic
fid=fopen(name2);
A2=textscan(fid,'%s');
fclose(fid);
k=0;
D2=A1{1,1};
[a2,b2]=size(D2);
D=A2{1,1};
[a1,b2]=size(D);

for j=1:5:a2
       k=k+1;
      x1(k,1)=str2double(D2(j+1));
      y1(k,1)=str2double(D2(j+2));
      z1(k,1)=str2double(D2(j+3));
     B_exp(k,1)=str2double(D2(j+4)); 
end
%%%%%%%%%%%%%
%%%%%%%%%%%%%%
xyz1=[x1 y1 z1];
[b,boo]=size(x1);
Nca=b;
xyz_u=xyz1;
xyz_u=xyz_u';
for i=1:3*Nca
xyz_u_LINED(i,1)=xyz_u(i);
end

%fid=fopen('as.xyz','w+');
%for i=1:b
%fprintf(fid,'%3f17 %3f17 %3f17\n ',x(i,1),y(i,1),z(i,1));
%end
%fclose(fid);
nodes=xyz1;
%%%%%%%%%%%%%%%%%%%%%%%                                       cutoff10
cutoff=10;
links=[];Nlinks=0;
for i=1:b
    for j=i+1:b
        if sqrt((x1(i)-x1(j))^2 +(y1(i)-y1(j))^2+(z1(i)-z1(j))^2) < cutoff
          Nlinks=Nlinks+1;
          links=[links ; [i,j]];
        end
    end
end
length=zeros(Nlinks,1);xcl=length;
ycl=xcl;zcl=xcl;cx=xcl;cy=xcl;cz=xcl;
stiffnessANM12=zeros(3*b, 3*b);stiffnessANM11=zeros(3*b, 3*b);
stiffnessANM13=zeros(3*b, 3*b);stiffnessANM14=zeros(3*b, 3*b);
stiffnessANM15=zeros(3*b, 3*b);stiffnessANM16=zeros(3*b, 3*b);
stiffnessANM21=zeros(3*b, 3*b);stiffnessANM22=zeros(3*b, 3*b);
stiffnessANM23=zeros(3*b, 3*b);stiffnessANM31=zeros(3*b, 3*b);
stiffnessANM32=zeros(3*b, 3*b);stiffnessANM33=zeros(3*b, 3*b);
for i=1:Nlinks
    %%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%
    j=links(i,1);
    k=links(i,2);
    xcl(i)=x1(k)-x1(j);
    ycl(i)=y1(k)-y1(j);
    zcl(i)=z1(k)-z1(j);
    length(i,1)=sqrt((x1(k)-x1(j))^2 +(y1(k)-y1(j))^2+(z1(k)-z1(j))^2);
    length2(i,1)=((x1(k)-x1(j))^2 +(y1(k)-y1(j))^2+(z1(k)-z1(j))^2);
    %fb(i,1)=exp(- length2(i,1)/16);
   % b0(i,1)=2*fb(i,1)/ length2(i,1);
    cx(i) = xcl(i)/length(i,1); 
    cy(i) = ycl(i)/length(i,1); 
    cz(i) = zcl(i)/length(i,1); 
    %ba cutoff 15 kolan bejaye 2*exp((-length(i,1)^2)/16)
    %benevisim (1/ length(i,1)^2)*
     
    stif1=[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
    
         stif2=(1/ length(i,1))*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
   
         stif3=(1/ length(i,1)^2)*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
  
         stif4=2*exp((-length(i,1)^2)/9)*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
   
         stif5=2*exp((-length(i,1)^2)/16)*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
    
         stif6=2*exp((-length(i,1)^2)/25)*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
    
    %disp(i);
    %disp(stif1);
   % if k-j == 1
      %  stif=10*stif;
    %end
    for r=1:3
        for s=1:3
               stiffnessANM11(3*(j-1)+r,3*(j-1)+s)=stiffnessANM11(3*(j-1)+r,3*(j-1)+s)+stif1(r,s); 
               stiffnessANM11(3*(j-1)+r,3*(k-1)+s)=stiffnessANM11(3*(j-1)+r,3*(k-1)+s)+stif1(r,s+3); 
               stiffnessANM11(3*(k-1)+r,3*(j-1)+s)=stiffnessANM11(3*(k-1)+r,3*(j-1)+s)+stif1(r+3,s); 
               stiffnessANM11(3*(k-1)+r,3*(k-1)+s)=stiffnessANM11(3*(k-1)+r,3*(k-1)+s)+stif1(r+3,s+3); 
        end
    end
    for r=1:3
        for s=1:3
               stiffnessANM12(3*(j-1)+r,3*(j-1)+s)=stiffnessANM12(3*(j-1)+r,3*(j-1)+s)+stif2(r,s); 
               stiffnessANM12(3*(j-1)+r,3*(k-1)+s)=stiffnessANM12(3*(j-1)+r,3*(k-1)+s)+stif2(r,s+3); 
               stiffnessANM12(3*(k-1)+r,3*(j-1)+s)=stiffnessANM12(3*(k-1)+r,3*(j-1)+s)+stif2(r+3,s); 
               stiffnessANM12(3*(k-1)+r,3*(k-1)+s)=stiffnessANM12(3*(k-1)+r,3*(k-1)+s)+stif2(r+3,s+3); 
        end
    end
        for r=1:3
        for s=1:3
               stiffnessANM13(3*(j-1)+r,3*(j-1)+s)=stiffnessANM13(3*(j-1)+r,3*(j-1)+s)+stif3(r,s); 
               stiffnessANM13(3*(j-1)+r,3*(k-1)+s)=stiffnessANM13(3*(j-1)+r,3*(k-1)+s)+stif3(r,s+3); 
               stiffnessANM13(3*(k-1)+r,3*(j-1)+s)=stiffnessANM13(3*(k-1)+r,3*(j-1)+s)+stif3(r+3,s); 
               stiffnessANM13(3*(k-1)+r,3*(k-1)+s)=stiffnessANM13(3*(k-1)+r,3*(k-1)+s)+stif3(r+3,s+3); 
        end
    end
        for r=1:3
        for s=1:3
               stiffnessANM14(3*(j-1)+r,3*(j-1)+s)=stiffnessANM14(3*(j-1)+r,3*(j-1)+s)+stif4(r,s); 
               stiffnessANM14(3*(j-1)+r,3*(k-1)+s)=stiffnessANM14(3*(j-1)+r,3*(k-1)+s)+stif4(r,s+3); 
               stiffnessANM14(3*(k-1)+r,3*(j-1)+s)=stiffnessANM14(3*(k-1)+r,3*(j-1)+s)+stif4(r+3,s); 
               stiffnessANM14(3*(k-1)+r,3*(k-1)+s)=stiffnessANM14(3*(k-1)+r,3*(k-1)+s)+stif4(r+3,s+3); 
        end
    end
        for r=1:3
        for s=1:3
               stiffnessANM15(3*(j-1)+r,3*(j-1)+s)=stiffnessANM15(3*(j-1)+r,3*(j-1)+s)+stif5(r,s); 
               stiffnessANM15(3*(j-1)+r,3*(k-1)+s)=stiffnessANM15(3*(j-1)+r,3*(k-1)+s)+stif5(r,s+3); 
               stiffnessANM15(3*(k-1)+r,3*(j-1)+s)=stiffnessANM15(3*(k-1)+r,3*(j-1)+s)+stif5(r+3,s); 
               stiffnessANM15(3*(k-1)+r,3*(k-1)+s)=stiffnessANM15(3*(k-1)+r,3*(k-1)+s)+stif5(r+3,s+3); 
        end
    end
        for r=1:3
        for s=1:3
               stiffnessANM16(3*(j-1)+r,3*(j-1)+s)=stiffnessANM16(3*(j-1)+r,3*(j-1)+s)+stif6(r,s); 
               stiffnessANM16(3*(j-1)+r,3*(k-1)+s)=stiffnessANM16(3*(j-1)+r,3*(k-1)+s)+stif6(r,s+3); 
               stiffnessANM16(3*(k-1)+r,3*(j-1)+s)=stiffnessANM16(3*(k-1)+r,3*(j-1)+s)+stif6(r+3,s); 
               stiffnessANM16(3*(k-1)+r,3*(k-1)+s)=stiffnessANM16(3*(k-1)+r,3*(k-1)+s)+stif6(r+3,s+3); 
        end
    end




end
%%%%%%%%%%%%%%%%%%%%%%%%

cutoff=15;
links=[];Nlinks=0;
for i=1:b
    for j=i+1:b
        if sqrt((x1(i)-x1(j))^2 +(y1(i)-y1(j))^2+(z1(i)-z1(j))^2) < cutoff
          Nlinks=Nlinks+1;
          links=[links ; [i,j]];
        end
    end
end
length=zeros(Nlinks,1);xcl=length;
ycl=xcl;zcl=xcl;cx=xcl;cy=xcl;cz=xcl;
stiffnessANM21=zeros(3*b, 3*b);stiffnessANM22=zeros(3*b, 3*b);
stiffnessANM23=zeros(3*b, 3*b);
for i=1:Nlinks
    %%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%
    j=links(i,1);
    k=links(i,2);
    xcl(i)=x1(k)-x1(j);
    ycl(i)=y1(k)-y1(j);
    zcl(i)=z1(k)-z1(j);
    length(i,1)=sqrt((x1(k)-x1(j))^2 +(y1(k)-y1(j))^2+(z1(k)-z1(j))^2);
    length2(i,1)=((x1(k)-x1(j))^2 +(y1(k)-y1(j))^2+(z1(k)-z1(j))^2);
    %fb(i,1)=exp(- length2(i,1)/16);
   % b0(i,1)=2*fb(i,1)/ length2(i,1);
    cx(i) = xcl(i)/length(i,1); 
    cy(i) = ycl(i)/length(i,1); 
    cz(i) = zcl(i)/length(i,1); 
    %ba cutoff 15 kolan bejaye 2*exp((-length(i,1)^2)/16)
    %benevisim (1/ length(i,1)^2)*
     
    stif31=[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
    
         stif32=(1/ length(i,1))*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
   
         stif33=(1/ length(i,1)^2)*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
    for r=1:3
        for s=1:3
               stiffnessANM31(3*(j-1)+r,3*(j-1)+s)=stiffnessANM31(3*(j-1)+r,3*(j-1)+s)+stif31(r,s); 
               stiffnessANM31(3*(j-1)+r,3*(k-1)+s)=stiffnessANM31(3*(j-1)+r,3*(k-1)+s)+stif31(r,s+3); 
               stiffnessANM31(3*(k-1)+r,3*(j-1)+s)=stiffnessANM31(3*(k-1)+r,3*(j-1)+s)+stif31(r+3,s); 
               stiffnessANM31(3*(k-1)+r,3*(k-1)+s)=stiffnessANM31(3*(k-1)+r,3*(k-1)+s)+stif31(r+3,s+3); 
        end
    end
    for r=1:3
        for s=1:3
               stiffnessANM32(3*(j-1)+r,3*(j-1)+s)=stiffnessANM32(3*(j-1)+r,3*(j-1)+s)+stif32(r,s); 
               stiffnessANM32(3*(j-1)+r,3*(k-1)+s)=stiffnessANM32(3*(j-1)+r,3*(k-1)+s)+stif32(r,s+3); 
               stiffnessANM32(3*(k-1)+r,3*(j-1)+s)=stiffnessANM32(3*(k-1)+r,3*(j-1)+s)+stif32(r+3,s); 
               stiffnessANM32(3*(k-1)+r,3*(k-1)+s)=stiffnessANM32(3*(k-1)+r,3*(k-1)+s)+stif32(r+3,s+3); 
        end
    end
        for r=1:3
        for s=1:3
               stiffnessANM33(3*(j-1)+r,3*(j-1)+s)=stiffnessANM33(3*(j-1)+r,3*(j-1)+s)+stif33(r,s); 
               stiffnessANM33(3*(j-1)+r,3*(k-1)+s)=stiffnessANM33(3*(j-1)+r,3*(k-1)+s)+stif33(r,s+3); 
               stiffnessANM33(3*(k-1)+r,3*(j-1)+s)=stiffnessANM33(3*(k-1)+r,3*(j-1)+s)+stif33(r+3,s); 
               stiffnessANM33(3*(k-1)+r,3*(k-1)+s)=stiffnessANM33(3*(k-1)+r,3*(k-1)+s)+stif33(r+3,s+3); 
        end
    end
    
    
end
%%%%%
%d.setValue(0.40);
cutoff=10;
links=[];Nlinks=0;
for i=1:b
    for j=i+1:b
        if sqrt((x1(i)-x1(j))^2 +(y1(i)-y1(j))^2+(z1(i)-z1(j))^2) < cutoff
          Nlinks=Nlinks+1;
          links=[links ; [i,j]];
        end
    end
end
length=zeros(Nlinks,1);xcl=length;
ycl=xcl;zcl=xcl;cx=xcl;cy=xcl;cz=xcl;
stiffnessANM21=zeros(3*b, 3*b);stiffnessANM22=zeros(3*b, 3*b);
stiffnessANM23=zeros(3*b, 3*b);
for i=1:Nlinks
    %%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%
    j=links(i,1);
    k=links(i,2);
    xcl(i)=x1(k)-x1(j);
    ycl(i)=y1(k)-y1(j);
    zcl(i)=z1(k)-z1(j);
    length(i,1)=sqrt((x1(k)-x1(j))^2 +(y1(k)-y1(j))^2+(z1(k)-z1(j))^2);
    length2(i,1)=((x1(k)-x1(j))^2 +(y1(k)-y1(j))^2+(z1(k)-z1(j))^2);
    %fb(i,1)=exp(- length2(i,1)/16);
   % b0(i,1)=2*fb(i,1)/ length2(i,1);
    cx(i) = xcl(i)/length(i,1); 
    cy(i) = ycl(i)/length(i,1); 
    cz(i) = zcl(i)/length(i,1); 
    %ba cutoff 15 kolan bejaye 2*exp((-length(i,1)^2)/16)
    %benevisim (1/ length(i,1)^2)*
     
    stif21=[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
    
         stif22=(1/ length(i,1))*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
   
         stif23=(1/ length(i,1)^2)*[cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i),-cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i);
        cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i),-cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i);
        cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i),-cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i);
        -cx(i)*cx(i),-cx(i)*cy(i),-cx(i)*cz(i),cx(i)*cx(i),cx(i)*cy(i),cx(i)*cz(i);
        -cy(i)*cx(i),-cy(i)*cy(i),-cy(i)*cz(i),cy(i)*cx(i),cy(i)*cy(i),cy(i)*cz(i);
        -cz(i)*cx(i),-cz(i)*cy(i),-cz(i)*cz(i),cz(i)*cx(i),cz(i)*cy(i),cz(i)*cz(i)];
    for r=1:3
        for s=1:3
               stiffnessANM21(3*(j-1)+r,3*(j-1)+s)=stiffnessANM21(3*(j-1)+r,3*(j-1)+s)+stif21(r,s); 
               stiffnessANM21(3*(j-1)+r,3*(k-1)+s)=stiffnessANM21(3*(j-1)+r,3*(k-1)+s)+stif21(r,s+3); 
               stiffnessANM21(3*(k-1)+r,3*(j-1)+s)=stiffnessANM21(3*(k-1)+r,3*(j-1)+s)+stif21(r+3,s); 
               stiffnessANM21(3*(k-1)+r,3*(k-1)+s)=stiffnessANM21(3*(k-1)+r,3*(k-1)+s)+stif21(r+3,s+3); 
        end
    end
    for r=1:3
        for s=1:3
               stiffnessANM22(3*(j-1)+r,3*(j-1)+s)=stiffnessANM22(3*(j-1)+r,3*(j-1)+s)+stif22(r,s); 
               stiffnessANM22(3*(j-1)+r,3*(k-1)+s)=stiffnessANM22(3*(j-1)+r,3*(k-1)+s)+stif22(r,s+3); 
               stiffnessANM22(3*(k-1)+r,3*(j-1)+s)=stiffnessANM22(3*(k-1)+r,3*(j-1)+s)+stif22(r+3,s); 
               stiffnessANM22(3*(k-1)+r,3*(k-1)+s)=stiffnessANM22(3*(k-1)+r,3*(k-1)+s)+stif22(r+3,s+3); 
        end
    end
        for r=1:3
        for s=1:3
               stiffnessANM23(3*(j-1)+r,3*(j-1)+s)=stiffnessANM23(3*(j-1)+r,3*(j-1)+s)+stif23(r,s); 
               stiffnessANM23(3*(j-1)+r,3*(k-1)+s)=stiffnessANM23(3*(j-1)+r,3*(k-1)+s)+stif23(r,s+3); 
               stiffnessANM23(3*(k-1)+r,3*(j-1)+s)=stiffnessANM23(3*(k-1)+r,3*(j-1)+s)+stif23(r+3,s); 
               stiffnessANM23(3*(k-1)+r,3*(k-1)+s)=stiffnessANM23(3*(k-1)+r,3*(k-1)+s)+stif23(r+3,s+3); 
        end
    end
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;
%d.setValue(0.50);

for j=1:5:a1
       k=k+1;
      x(k,1)=str2double(D(j+1));
      y(k,1)=str2double(D(j+2));
      z(k,1)=str2double(D(j+3));

end
%%%%%%%%%%%%%
%%%%%%%%%%%%%%
xyz=[x y z];

[b,boo]=size(x);
%fid=fopen('as.xyz','w+');
%for i=1:b
%fprintf(fid,'%3f17 %3f17 %3f17\n ',x(i,1),y(i,1),z(i,1));
%end
%fclose(fid);
xyz_b=xyz;
xyz_b=xyz_b';
for i=1:3*b
xyz_b_LINED(i,1)=xyz_b(i);
end
nodes=xyz;
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%

Dif1 = -(xyz1 - xyz);  % second structure subtract first structure
D2= Dif1';
for i = 1 : 3 *b 
Dif3(i,1) = D2(i);
end
Dif=Dif3/norm(Dif3);
%stifnese unbound
Dif_LINED = (xyz_b_LINED - xyz_u_LINED);  % second structure subtract first structure

%for i = 1:3*Nca
%Dif1(i,1) = D2(i);
%end
norm_Dif_LINED = norm(Dif_LINED);
for i = 1 : 3 *b 
Dif3(i,1) = D2(i);
end
Dif=Dif3/norm(Dif3);
%stifnese unbound
                        %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%
                        %%                                    overlap11  rc 10                                     %
                        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%
 % model 1 radius 1 means 10 angstrom                       
[c,c1]=eig(stiffnessANM11);
[ANM_Values11,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
%for i=1:3*b
   % ANM_Vectors(:,i)=ANM_Vectors(:,i)./norm(ANM_Vectors(:,i));
%end
%%%%%%%%
for i=1:3*b
   ANM_Vectors11(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM11(i,1)=abs((Dif)' * ANM_Vectors11(:,i));
end
maxoverlapvalue11 = max (Dif_ANM11(1:16));
for i=1:16
    if Dif_ANM11(i,1) == maxoverlapvalue11
        maxoverlapindex11= i - 6;
    end 
    
end
   % fprintf('rc = 18 network:const mode: %d max ovlp: %2.3f \n',maxoverlapindex11,maxoverlapvalue11);
    %% %%%%%%%%%%%%%%%soheil
    for Nmode = 1:3*Nca-6
     ANM_Vectors11_Nmode = ANM_Vectors11(:,7:6+Nmode);
   parameters11=ANM_Vectors11_Nmode\Dif_LINED;
  deltrtaghrib11=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib11=deltrtaghrib11+parameters11(i,1)*ANM_Vectors11_Nmode(:,i);
    end
   xyz_p_LINED11 = xyz_u_LINED + deltrtaghrib11;
   Dif_LINED_final11= xyz_p_LINED11-xyz_b_LINED;
   norm_Dif_LINED_final11 = norm(Dif_LINED_final11);
   RMSD_10_const(Nmode,1)= norm(Dif_LINED_final11)/sqrt(Nca);
   %%%
    deltrtaghrib_norm11= norm(deltrtaghrib11);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%
%%                            rc 10          overlap12                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%
%d.setValue(0.60);
[c,c1]=eig(stiffnessANM12);
[ANM_Values12,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
%for i=1:3*b
   % ANM_Vectors(:,i)=ANM_Vectors(:,i)./norm(ANM_Vectors(:,i));
%end
%%%%%%%%
for i=1:3*b
   ANM_Vectors12(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM12(i,1)=abs((Dif)' * ANM_Vectors12(:,i));
end
maxoverlapvalue12 = max (Dif_ANM12(1:16));
for i=1:16
    if Dif_ANM12(i,1) == maxoverlapvalue12
        maxoverlapindex12= i - 6;
    end 
    
end
  %  fprintf('rc = 18 network:1/r mode: %d max ovlp: %2.3f \n',maxoverlapindex12,maxoverlapvalue12);
       ANM_Vectors12_Nmode = ANM_Vectors12(:,7:6+Nmode);
   parameters12=ANM_Vectors12_Nmode\Dif_LINED;
  deltrtaghrib12=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib12=deltrtaghrib12+parameters12(i,1)*ANM_Vectors12_Nmode(:,i);
    end
   xyz_p_LINED12 = xyz_u_LINED + deltrtaghrib12;
   Dif_LINED_final12= xyz_p_LINED12-xyz_b_LINED;
   norm_Dif_LINED_final12 = norm(Dif_LINED_final12);
   RMSD_10_1r(Nmode,1)= norm(Dif_LINED_final12)/sqrt(Nca);
  
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%
%%                               rc 10          overlap13                                  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%
[c,c1]=eig(stiffnessANM13);
[ANM_Values13,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
for i=1:3*b
   ANM_Vectors13(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM13(i,1)=abs((Dif)' * ANM_Vectors13(:,i));
end
maxoverlapvalue13 = max (Dif_ANM13(1:16));
for i=1:16
    if Dif_ANM13(i,1) == maxoverlapvalue13
        maxoverlapindex13= i - 6;
    end 
    
end
    %fprintf('rc = 18 network:1/r^2 mode: %d max ovlp: %2.3f \n',maxoverlapindex13,maxoverlapvalue13);
        ANM_Vectors13_Nmode = ANM_Vectors13(:,7:6+Nmode);
   parameters13=ANM_Vectors13_Nmode\Dif_LINED;
  deltrtaghrib13=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib13=deltrtaghrib13+parameters13(i,1)*ANM_Vectors13_Nmode(:,i);
    end
   xyz_p_LINED13 = xyz_u_LINED + deltrtaghrib13;
   Dif_LINED_final13= xyz_p_LINED13-xyz_b_LINED;
   norm_Dif_LINED_final13 = norm(Dif_LINED_final13);
   RMSD_10_1r2(Nmode,1)= norm(Dif_LINED_final13)/sqrt(Nca);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 14 
    [c,c1]=eig(stiffnessANM14);
[ANM_Values14,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
for i=1:3*b
   ANM_Vectors14(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM14(i,1)=abs((Dif)' * ANM_Vectors14(:,i));
end
maxoverlapvalue14 = max (Dif_ANM14(1:16));
for i=1:16
    if Dif_ANM14(i,1) == maxoverlapvalue14
        maxoverlapindex14= i - 6;
    end 
    
end
   % fprintf('rscale = 3 network:exp mode: %d max ovlp: %2.3f \n',maxoverlapindex14,maxoverlapvalue14);
  %d.setValue(0.70);
  
      ANM_Vectors14_Nmode = ANM_Vectors14(:,7:6+Nmode);
   parameters14=ANM_Vectors14_Nmode\Dif_LINED;
  deltrtaghrib14=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib14=deltrtaghrib14+parameters14(i,1)*ANM_Vectors14_Nmode(:,i);
    end
   xyz_p_LINED14 = xyz_u_LINED + deltrtaghrib14;
   Dif_LINED_final14= xyz_p_LINED14-xyz_b_LINED;
   norm_Dif_LINED_final14 = norm(Dif_LINED_final14);
   RMSD_exp_3(Nmode,1)= norm(Dif_LINED_final14)/sqrt(Nca);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%15
  [c,c1]=eig(stiffnessANM15);
[ANM_Values15,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
for i=1:3*b
   ANM_Vectors15(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM15(i,1)=abs((Dif)' * ANM_Vectors15(:,i));
end
maxoverlapvalue15 = max (Dif_ANM15(1:16));
for i=1:16
    if Dif_ANM15(i,1) == maxoverlapvalue15
        maxoverlapindex15= i - 6;
    end 
    
end
   % fprintf('rscale = 4 network:exp mode: %d max ovlp: %2.3f \n',maxoverlapindex15,maxoverlapvalue15);
    
        ANM_Vectors15_Nmode = ANM_Vectors15(:,7:6+Nmode);
   parameters15=ANM_Vectors15_Nmode\Dif_LINED;
  deltrtaghrib15=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib15=deltrtaghrib15+parameters15(i,1)*ANM_Vectors15_Nmode(:,i);
    end
   xyz_p_LINED15 = xyz_u_LINED + deltrtaghrib15;
   Dif_LINED_final15= xyz_p_LINED15-xyz_b_LINED;
   norm_Dif_LINED_final15 = norm(Dif_LINED_final15);
   RMSD_exp_4(Nmode,1)= norm(Dif_LINED_final15)/sqrt(Nca);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%16
    [c,c1]=eig(stiffnessANM16);
[ANM_Values16,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
for i=1:3*b
   ANM_Vectors16(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM16(i,1)=abs((Dif)' * ANM_Vectors16(:,i));
end
maxoverlapvalue16 = max (Dif_ANM16(1:16));
for i=1:16
    if Dif_ANM16(i,1) == maxoverlapvalue16
        maxoverlapindex16= i - 6;
    end 
    
end
    %fprintf('rscale = 5 network:exp mode: %d max ovlp: %2.3f \n',maxoverlapindex16,maxoverlapvalue16);
        ANM_Vectors16_Nmode = ANM_Vectors16(:,7:6+Nmode);
   parameters16=ANM_Vectors16_Nmode\Dif_LINED;
  deltrtaghrib16=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib16=deltrtaghrib16+parameters16(i,1)*ANM_Vectors16_Nmode(:,i);
    end
   xyz_p_LINED16 = xyz_u_LINED + deltrtaghrib16;
   Dif_LINED_final16= xyz_p_LINED16-xyz_b_LINED;
   norm_Dif_LINED_final16 = norm(Dif_LINED_final16);
   RMSD_exp_5(Nmode,1)= norm(Dif_LINED_final16)/sqrt(Nca);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%21
 [c,c1]=eig(stiffnessANM11);
[ANM_Values11,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
%for i=1:3*b
   % ANM_Vectors(:,i)=ANM_Vectors(:,i)./norm(ANM_Vectors(:,i));
%end
%%%%%%%%

[c,c1]=eig(stiffnessANM21);
[ANM_Values21,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
%for i=1:3*b
   % ANM_Vectors(:,i)=ANM_Vectors(:,i)./norm(ANM_Vectors(:,i));
%end
%%%%%%%%
for i=1:3*b
   ANM_Vectors21(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM21(i,1)=abs((Dif)' * ANM_Vectors21(:,i));
end
maxoverlapvalue21 = max (Dif_ANM21(1:16));
for i=1:16
    if Dif_ANM21(i,1) == maxoverlapvalue21
        maxoverlapindex21= i - 6;
    end 
    
end
   % fprintf('rc = 10 network:const mode: %d max ovlp: %2.3f \n',maxoverlapindex21,maxoverlapvalue21);
    
        ANM_Vectors21_Nmode = ANM_Vectors21(:,7:6+Nmode);
   parameters21=ANM_Vectors21_Nmode\Dif_LINED;
  deltrtaghrib21=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib21=deltrtaghrib21+parameters21(i,1)*ANM_Vectors21_Nmode(:,i);
    end
   xyz_p_LINED21 = xyz_u_LINED + deltrtaghrib21;
   Dif_LINED_final21= xyz_p_LINED21-xyz_b_LINED;
   norm_Dif_LINED_final21 = norm(Dif_LINED_final21);
   RMSD_15_const(Nmode,1)= norm(Dif_LINED_final21)/sqrt(Nca);
   %%%%%%%%%%%%%%%%%%%%%%%%
   %d.setValue(0.80);
   [c,c1]=eig(stiffnessANM22);
[ANM_Values22,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
%for i=1:3*b
   % ANM_Vectors(:,i)=ANM_Vectors(:,i)./norm(ANM_Vectors(:,i));
%end
%%%%%%%%
for i=1:3*b
   ANM_Vectors22(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM22(i,1)=abs((Dif)' * ANM_Vectors22(:,i));
end
maxoverlapvalue22 = max (Dif_ANM22(1:16));
for i=1:16
    if Dif_ANM22(i,1) == maxoverlapvalue22
        maxoverlapindex22= i - 6;
    end 
    
end
   % fprintf('rc = 10 network:1/r mode: %d max ovlp: %2.3f \n',maxoverlapindex22,maxoverlapvalue22);
        ANM_Vectors22_Nmode = ANM_Vectors22(:,7:6+Nmode);
   parameters22=ANM_Vectors22_Nmode\Dif_LINED;
  deltrtaghrib22=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib22=deltrtaghrib22+parameters22(i,1)*ANM_Vectors22_Nmode(:,i);
    end
   xyz_p_LINED22 = xyz_u_LINED + deltrtaghrib22;
   Dif_LINED_final22= xyz_p_LINED22-xyz_b_LINED;
   norm_Dif_LINED_final22 = norm(Dif_LINED_final22);
   RMSD_15_1r(Nmode,1)= norm(Dif_LINED_final22)/sqrt(Nca);
    %%%%%%%%%%%%%
[c,c1]=eig(stiffnessANM23);
[ANM_Values23,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
%for i=1:3*b
   % ANM_Vectors(:,i)=ANM_Vectors(:,i)./norm(ANM_Vectors(:,i));
%end
%%%%%%%%
for i=1:3*b
   ANM_Vectors23(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM23(i,1)=abs((Dif)' * ANM_Vectors23(:,i));
end
maxoverlapvalue23 = max (Dif_ANM23(1:16));
for i=1:16
    if Dif_ANM23(i,1) == maxoverlapvalue23
        maxoverlapindex23= i - 6;
    end 
    
end
    %fprintf('rc = 10 network:1/r^2 mode: %d max ovlp: %2.3f \n',maxoverlapindex23,maxoverlapvalue23);
        ANM_Vectors23_Nmode = ANM_Vectors23(:,7:6+Nmode);
   parameters23=ANM_Vectors23_Nmode\Dif_LINED;
  deltrtaghrib23=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib23=deltrtaghrib23+parameters23(i,1)*ANM_Vectors23_Nmode(:,i);
    end
   xyz_p_LINED23 = xyz_u_LINED + deltrtaghrib23;
   Dif_LINED_final23= xyz_p_LINED23-xyz_b_LINED;
   norm_Dif_LINED_final23 = norm(Dif_LINED_final23);
   RMSD_15_1r2(Nmode,1)= norm(Dif_LINED_final23)/sqrt(Nca);
    %d.setValue(0.850);
   %%%%%%%%%%%%%%%%%%%%%%%%31
    [w,w2]=eig(stiffnessANM31);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM31=t2/(t3*t4);
[c,c1]=eig(stiffnessANM31);
[ANM_Values31,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
for i=1:3*b
   ANM_Vectors31(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM31(i,1)=abs((Dif)' * ANM_Vectors31(:,i));
end
maxoverlapvalue31 = max (Dif_ANM31(1:16));
for i=1:16
    if Dif_ANM31(i,1) == maxoverlapvalue31
        maxoverlapindex31= i - 6;
    end 
    
end
    %fprintf('rc = 15 network:const mode: %d max ovlp: %2.3f \n',maxoverlapindex31,maxoverlapvalue31);
        ANM_Vectors31_Nmode = ANM_Vectors31(:,7:6+Nmode);
   parameters31=ANM_Vectors31_Nmode\Dif_LINED;
  deltrtaghrib31=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib31=deltrtaghrib31+parameters31(i,1)*ANM_Vectors31_Nmode(:,i);
    end
   xyz_p_LINED31 = xyz_u_LINED + deltrtaghrib31;
   Dif_LINED_final31= xyz_p_LINED31-xyz_b_LINED;
   norm_Dif_LINED_final31 = norm(Dif_LINED_final31);
   RMSD_18_const(Nmode,1)= norm(Dif_LINED_final31)/sqrt(Nca);
    
                  %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %%                                    overlap32  rc 18                                %%
                  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%   
[w,w2]=eig(stiffnessANM32);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM32=t2/(t3*t4);
[c,c1]=eig(stiffnessANM32);
[ANM_Values32,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
for i=1:3*b
   ANM_Vectors32(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM32(i,1)=abs((Dif)' * ANM_Vectors32(:,i));
end
maxoverlapvalue32 = max (Dif_ANM32(1:16));
for i=1:16
    if Dif_ANM32(i,1) == maxoverlapvalue32
        maxoverlapindex32= i - 6;
    end 
    
end
   % fprintf('rc = 15 network:1/r mode: %d max ovlp: %2.3f \n',maxoverlapindex32,maxoverlapvalue32); 
    
        ANM_Vectors32_Nmode = ANM_Vectors32(:,7:6+Nmode);
   parameters32=ANM_Vectors32_Nmode\Dif_LINED;
  deltrtaghrib32=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib32=deltrtaghrib32+parameters32(i,1)*ANM_Vectors32_Nmode(:,i);
    end
   xyz_p_LINED32 = xyz_u_LINED + deltrtaghrib32;
   Dif_LINED_final32= xyz_p_LINED32-xyz_b_LINED;
   norm_Dif_LINED_final32 = norm(Dif_LINED_final32);
   RMSD_18_1r(Nmode,1)= norm(Dif_LINED_final32)/sqrt(Nca);
   %%%%%%%%%%%%%%%%%%%%%%%%%33
   [w,w2]=eig(stiffnessANM11);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM11=t2/(t3*t4);
%%%%
   [w,w2]=eig(stiffnessANM12);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM12=t2/(t3*t4);
   %%%%%%%%
      [w,w2]=eig(stiffnessANM13);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM13=t2/(t3*t4);
   %%%%%%%%%
    [w,w2]=eig(stiffnessANM14);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM14=t2/(t3*t4);  
   %%%%%%%%%
     [w,w2]=eig(stiffnessANM15);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM15=t2/(t3*t4); 
   %%%%%%%%%%
      [w,w2]=eig(stiffnessANM16);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM16=t2/(t3*t4);
   %%%%%%
      [w,w2]=eig(stiffnessANM21);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM21=t2/(t3*t4);
   %%%%%%%.
      [w,w2]=eig(stiffnessANM22);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM22=t2/(t3*t4);
   %%%%%%%
      [w,w2]=eig(stiffnessANM23);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM23=t2/(t3*t4);
   %%%%%%%%%%%
[w,w2]=eig(stiffnessANM33);values=diag(w2);
for i=7:3*Nca
    w1(:,i)=w(:,i)./norm(w(:,i));
end
for i=7:3*Nca
    for j=0:Nca-1
        n=3*j+1;
        uij2(j+1,i)=((w1(n,i)^2)+(w1(n+1,i)^2)+(w1(n+2,i)^2));
    end
end
for i=7:3*Nca
    uij2_omeg(:,i)=uij2(:,i)./(values(i,1));
end
for i=1:Nca
    delta_R(i,1)=KBT_M*sum(uij2_omeg(i,7:3*Nca));B_cal(i,1)=A*delta_R(i,1);
end
h=sum(B_exp)/Nca;h1=sum(B_cal)/Nca;g=B_exp(:,1)-h;g1=B_cal(:,1)-h1;
t1=g(:,1).*g1(:,1);t2=sum(t1);t3=norm(g);t4=norm(g1);CCANM33=t2/(t3*t4);
[c,c1]=eig(stiffnessANM33);
[ANM_Values33,Eig_number]=sort(diag(c1),'ascend');
for i=1:3*b
    cd(:,i)=c(:,Eig_number(i));
end
for i=1:3*b
   ANM_Vectors33(:,i)=cd(:,i)./norm(cd(:,i));
end
for i=1:3*b
    Dif_ANM33(i,1)=abs((Dif)' * ANM_Vectors33(:,i));
end
maxoverlapvalue33 = max (Dif_ANM33(1:16));
for i=1:16
    if Dif_ANM33(i,1) == maxoverlapvalue33
        maxoverlapindex33= i - 6;
    end 
    
end
%d.setValue(0.90);
    %fprintf('rc = 15 network:1/r^2 mode: %d max ovlp: %2.3f \n',maxoverlapindex33,maxoverlapvalue33);
        ANM_Vectors33_Nmode = ANM_Vectors33(:,7:6+Nmode);
   parameters33=ANM_Vectors33_Nmode\Dif_LINED;
  deltrtaghrib33=zeros(3*Nca,1);
   for i=1:Nmode
       deltrtaghrib33=deltrtaghrib33+parameters33(i,1)*ANM_Vectors33_Nmode(:,i);
    end
   xyz_p_LINED33 = xyz_u_LINED + deltrtaghrib33;
   Dif_LINED_final33= xyz_p_LINED33-xyz_b_LINED;
   norm_Dif_LINED_final33 = norm(Dif_LINED_final33);
   RMSD_18_1r2(Nmode,1)= norm(Dif_LINED_final33)/sqrt(Nca);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%%%%%%%%%% collectivity experimental  %%%%%%%%%%%%%%%%%%%%%%%%   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for i=1:b
    delta_R(i,1)=(Dif1(i,1)^2+Dif1(i,2)^2+Dif1(i,3)^2)^0.5;
end

alpha=1/sum(delta_R.^2);
for i=1:b
aa(i,1)=alpha*delta_R(i,1)^2*log(alpha.*(delta_R(i,1)^2));
end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  name3=sprintf('%s%s','overlap_',ani1,'.txt');
fid=fopen(name3,'w+'); 
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model :   const Rc= 10',maxoverlapvalue21,'(' ,maxoverlapindex21, ')',' B-factor cc:',CCANM21);
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model :   const Rc= 15',maxoverlapvalue31,'(' ,maxoverlapindex31, ')',' B-factor cc:',CCANM31);
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model :    const Rc =18',maxoverlapvalue11,'(' ,maxoverlapindex11, ')','B-factor cc:',CCANM11); 
 
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model : 1/r    Rc = 10',maxoverlapvalue22,'(' ,maxoverlapindex22, ')',' B-factor cc: ',CCANM22);
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model : 1/r    Rc = 15',maxoverlapvalue32,'(' ,maxoverlapindex32, ')',' B-factor cc: ',CCANM32);   
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model : 1/r    Rc = 18',maxoverlapvalue12,'(' ,maxoverlapindex12, ')',' B-factor cc: ',CCANM12);   

fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model : 1/r^2  Rc = 10',maxoverlapvalue23,'(' ,maxoverlapindex23, ')',' B-factor cc: ',CCANM23); 
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model : 1/r^2  Rc = 15',maxoverlapvalue33,'(' ,maxoverlapindex33, ')',' B-factor cc: ',CCANM33);   
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model : 1/r^2  Rc = 18',maxoverlapvalue13,'(' ,maxoverlapindex13, ')',' B-factor cc: ',CCANM13); 


fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f \n ','model : exp    Rs = 3 ',maxoverlapvalue14,'(' ,maxoverlapindex14, ')',' B-factor cc: ',CCANM14); 
fprintf(fid,'%s  %2.3f %s%d%s %s %2.3f\n ','model : exp    Rs = 4 ',maxoverlapvalue15,'(' ,maxoverlapindex15, ')',' B-factor cc: ',CCANM15);
fprintf(fid,'%s  %2.3f  %s%d%s %s %2.3f\n ','model : exp    Rs = 5 ',maxoverlapvalue16,'(',maxoverlapindex16, ')','B-factor cc: ',CCANM16);  
%d.setValue(0.99);
fclose(fid);
    disp('Done');
toc
%d.setVisible(false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gereftane residue displacement

for i=1:b
    delta_R(i,1)=(Dif1(i,1)^2+Dif1(i,2)^2+Dif1(i,3)^2)^0.5;
    
end
for i=1:b
  delta_R_norm(i,1) =delta_R(i,1)./norm(delta_R(:,1));
    
end

    %%%%%%%%%%%%%% total RMSD
    RMSD1 = 0;
    for i = 1:Nca
        DeltaR(i,1) = (Dif1(i,1)^2+Dif1(i,2)^2+Dif1(i,3)^2);
    RMSD1 =RMSD1 +  (DeltaR(i,1) / Nca);
    end
    RMSD_TOTAL =sqrt( RMSD1);
  %%%%%%%%%%      10_const
 tot_rate_10_const=(RMSD_TOTAL-RMSD_10_const)/RMSD_TOTAL * 100;
 RMSD_rate_10_const=(RMSD_TOTAL-RMSD_10_const)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_10_const(i,1)= tot_rate_10_const(i,1)-tot_rate_10_const(i-1,1);
 end
    %%%%%%%%%%      10_1r
 tot_rate_10_1r=(RMSD_TOTAL-RMSD_10_1r)/RMSD_TOTAL * 100;
 RMSD_rate_10_1r=(RMSD_TOTAL-RMSD_10_1r)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_10_1r(i,1)= tot_rate_10_1r(i,1)-tot_rate_10_1r(i-1,1);
 end
    %%%%%%%%%%      10-1/r2
 tot_rate_10_1r2=(RMSD_TOTAL-RMSD_10_1r2)/RMSD_TOTAL * 100;
 RMSD_rate_10_1r2=(RMSD_TOTAL-RMSD_10_1r2)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_10_1r2(i,1)= tot_rate_10_1r2(i,1)-tot_rate_10_1r2(i-1,1);
 end
    %%%%%%%%%%      exp3
 tot_rate_exp_3=(RMSD_TOTAL-RMSD_exp_3)/RMSD_TOTAL * 100;
 RMSD_rate_exp_3=(RMSD_TOTAL-RMSD_exp_3)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_exp_3(i,1)= tot_rate_exp_3(i,1)-tot_rate_exp_3(i-1,1);
 end
    %%%%%%%%%%      exp4
 tot_rate_exp_4=(RMSD_TOTAL-RMSD_exp_4)/RMSD_TOTAL * 100;
 RMSD_rate_exp_4=(RMSD_TOTAL-RMSD_exp_4)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_exp_4(i,1)= tot_rate_exp_4(i,1)-tot_rate_exp_4(i-1,1);
 end
 %%%%%%%%%%      exp5
 tot_rate_exp_5=(RMSD_TOTAL-RMSD_exp_5)/RMSD_TOTAL * 100;
 RMSD_rate_exp_5=(RMSD_TOTAL-RMSD_exp_5)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_exp_5(i,1)= tot_rate_exp_5(i,1)-tot_rate_exp_5(i-1,1);
 end
 %%%%%%%%%%      15_const
 tot_rate_15_const=(RMSD_TOTAL-RMSD_15_const)/RMSD_TOTAL * 100;
 RMSD_rate_15_const=(RMSD_TOTAL-RMSD_15_const)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_15_const(i,1)= tot_rate_15_const(i,1)-tot_rate_15_const(i-1,1);
 end 
    %%%%%%%%%%      15_1/r
 tot_rate_15_1r=(RMSD_TOTAL-RMSD_15_1r)/RMSD_TOTAL * 100;
 RMSD_rate_15_1r=(RMSD_TOTAL-RMSD_15_1r)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_15_1r(i,1)= tot_rate_15_1r(i,1)-tot_rate_15_1r(i-1,1);
 end
    %%%%%%%%%%      15_1/r2
 tot_rate_15_1r2=(RMSD_TOTAL-RMSD_15_1r2)/RMSD_TOTAL * 100;
 RMSD_rate_15_1r2=(RMSD_TOTAL-RMSD_15_1r2)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_15_1r2(i,1)= tot_rate_15_1r2(i,1)-tot_rate_15_1r2(i-1,1);
 end
    %%%%%%%%%%      18_const
 tot_rate_18_const=(RMSD_TOTAL-RMSD_18_const)/RMSD_TOTAL * 100;
 RMSD_rate_18_const=(RMSD_TOTAL-RMSD_18_const)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_18_const(i,1)= tot_rate_18_const(i,1)-tot_rate_18_const(i-1,1);
 end
 %%%%%%%%%%      8_1r
 tot_rate_18_1r=(RMSD_TOTAL-RMSD_18_1r)/RMSD_TOTAL * 100;
 RMSD_rate_18_1r=(RMSD_TOTAL-RMSD_18_1r)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_18_1r(i,1)= tot_rate_18_1r(i,1)-tot_rate_18_1r(i-1,1);
 end
    %%%%%%%%%%      8_1r2
 tot_rate_18_1r2=(RMSD_TOTAL-RMSD_18_1r2)/RMSD_TOTAL * 100;
 RMSD_rate_18_1r2=(RMSD_TOTAL-RMSD_18_1r2)/RMSD_TOTAL * 100;
 for i=2:3*Nca-6
    RMSD_rate_18_1r2(i,1)= tot_rate_18_1r2(i,1)-tot_rate_18_1r2(i-1,1);
 end
