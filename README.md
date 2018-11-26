# GAMS
genetic algorithm in molecular simulation
clc
clear all
close all
format long

MFI1;

%number of adsorbent atoms
nox=length(OX(:,1));
nsi=length(SI(:,1));

%%%number of initial population
minpop=2;
pop=minpop;

%%%critcal point Of Ethane
TC=305.5;   %K
PC=48.72; %bar
w=.099;
R=82.056e-6;  %atm.m3/gmol.K


%%%Temprature & Pressure &...
P=.001;%bar
T=300;%K
K=8.314;%J/gmole.K
KB=1.3806503e-23; %J/K
beta=1/(KB*T);
Rg=8.314;   %J/gmole.K

%parameter of GA
pc=.5;%Probability of cross over
pm=.2;%Probability of mutation

%parameters of LJ potential For H-SIOX
eps=[82.9848 44.2697 9.17;44.2697 23.6164 48.1;90.17 48.1 98];
sigma=[3.3 3.72 3.51;3.72 4.2 3.96;3.51 3.96 3.75];%K

%Structure of Adsorbent
mind=zeros(1,3);
for j=1:3
    maxd(j)=max(OX(:,j));
end

%Accuracy of structure of lattice
accrx=0.001;
accry=0.001;
accrz=0.001;
dx=(maxd(1,1)-mind(1,1))/accrx;
dy=(maxd(1,2)-mind(1,2))/accry;
dz=(maxd(1,3)-mind(1,3))/accrz;

%%%number of bits of X Y Z
nbit1=1;
n1=dx;
while n1>2
   n1=dx/2;
   nbit1=nbit1+1;
   dx=n1;
end

nbit2=1;
n2=dy;
while n2>2
   n2=dy/2;
   nbit2=nbit2+1;
   dy=n2;
end
nbit3=1;
n3=dz;
while n3>2
   n3=dz/2;
   nbit3=nbit3+1;
   dz=n3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Lattice
avag=6.023e23;% avagadro number
Vlattice=1e-30;
for i=1:3
    Vlattice=maxd(i)*Vlattice;%m^3
end
dens=.5e3;%Kg/m^3
%%%%%%%%%%%%%%%%%%%%%%%%%%

nbit=nbit1+nbit2+nbit3;
maxdec1=2^nbit1-1;
maxdec2=2^nbit2-1;
maxdec3=2^nbit3-1;
ii=1;
while P<.011
    deltaY=-1;
    while deltaY<0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Generating initial population%%%%%%%%%
           i=0;
           Iflag=0;
          while i<pop
                  if Iflag==0
                    i=i+1;
                  end
          Iflag=0;
          TEHS=0;
          TEHO=0;
          TEHH=0;
     %generating  random CHROMOSOME 
             for j=1:nbit
                rd=rand();
                   if rd>.5
                   cr(i,j)=1;
                   else 
                   cr(i,j)=0;
                   end
          end
         %Transfer CHROMOSOME to X Y Z
          sx=0;
          for k=1:(nbit1) 
             sx=cr(i,k)*2^((nbit1)-k)+sx;
          end
          x(i)=mind(1,1)+(sx/maxdec1)*(maxd(1,1)-mind(1,1));
          sy=0;
          for k=(nbit1)+1:nbit2+nbit1
             sy=cr(i,k)*2^(nbit2+nbit1-k)+sy;
          end
           y(i)=mind(1,2)+(sy/maxdec2)*(maxd(1,2)-mind(1,2));
          sz=0;
          for k=nbit2+nbit1+1:nbit
             sz=cr(i,k)*2^(nbit-k)+sz;
          end
          z(i)=mind(1,3)+(sz/maxdec3)*(maxd(1,3)-mind(1,3));
    %Coordinate of Hydrogen atoms
     Ethane(i,1)=x(i);
     Ethane(i,2)=y(i);
     Ethane(i,3)=z(i);
        %calculating distanse between ONE atome of adsorbate and other
        %adsorbate and adsorbent atoms &&& ENERGY of ONE atom of adsorbate
        for k=1:nox
         rEthaneox(i,k)=norm(OX(k,:)-Ethane(i,:));
        end
        for j=1:nox
         EEthaneox(i,j)=4*eps(1,3)*((sigma(1,3)/rEthaneox(i,j))^12-(sigma(1,3)/rEthaneox(i,j))^6);
         TEHO=TEHO+EEthaneox(i,j);
        end
        for k=1:nsi
          rEthanesi(i,k)=norm(SI(k,:)-Ethane(i,:));
        end
        for j=1:nsi
          EEthanesi(i,j)=4*eps(2,3)*((sigma(2,3)/rEthanesi(i,j))^12-(sigma(2,3)/rEthanesi(i,j))^6);
          TEHS=TEHS+EEthanesi(i,j);
        end
        for k=1:i
          rEthaneEthane(i,k)=norm(Ethane(k,:)-Ethane(i,:));
          rEthaneEthane(k,i)=rEthaneEthane(i,k);
        end
        for j=1:i
            if i==j
               EEthaneEthane(i,j)=0;
            end
            if i~=j
            EEthaneEthane(i,j)=4*eps(3,3)*((sigma(3,3)/rEthaneEthane(i,j))^12-(sigma(3,3)/rEthaneEthane(i,j))^6);
            end
            TEHH=TEHH+EEthaneEthane(i,j);
        end
        TE(i)=TEHO+TEHS+TEHH;
           if TE(i)>0 
              Iflag=1;
           end
           if Iflag==1 & i==pop
              i=pop-1;
           end
           fitness(i)=TE(i);
           cr(i,nbit+1)=-1*fitness(i);
   
  end
  sumf=sum(cr(:,nbit+1));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  STE1(1)=sumf;
  STEO=sumf;
  maxitr=5;
  kk=1;
  Kflag=0;
  while Kflag==0
   for itr=1:maxitr
     if kk==1
          if itr>1
            STE1(itr)=STE2(itr-1);%%%%%%%TOTAL fitness befor GENETIC operation
          end
     else
          if itr==1
            STE1(itr)=STE2(maxitr);
          else
            STE1(itr)=STE2(itr-1);%%%%%%%TOTAL fitness befor GENETIC operation
          end
     end

     %%%computing individual & cumulative probabilities
      for i=1:pop
       ip(i)=cr(i,nbit+1)/sumf;
      end
      cp(1)=ip(1);
      for k=2:pop
         cp(k)=ip(k)+cp(k-1);
      end
   
      CR=cr;%%%Saving CHROMOSOME from GENETIC operation

      %%%Selecting new generation
      crs=cr;
      for k=1:pop
         rd=rand();
         if rd<cp(1)
           cr(k,:)=crs(1,:);
         else
           i=1;
           while i<pop
              if rd>cp(i) & rd<cp(i+1)
                  cr(k,:)=crs(i+1,:);
                  i=pop;
              else
                  i=i+1;
              end
           end
         end
      end
       
     %cross over
     dd=1;
     for k=1:pop
         rd=rand();
         if rd<pc 
            candid(dd)=k;
            dd=dd+1;
         end
     end
     dd=dd-1;
     if dd>1
        r=mod(dd,2);
        if r==0
           for j=1:2:dd
               rd=1+fix(rand()*(nbit-2));
               l=candid(j);
               u=candid(j+1);
               for i=rd:nbit
                   t(i)=cr(l,i);
               end
               for i=rd:nbit
                   cr(l,i)=cr(u,i);
                   cr(u,i)=t(i);
               end
           end
        else
            for j=1:2:dd-1
               rd=1+fix((rand()*(nbit-2)));
               l=candid(j);
               u=candid(j+1);
               for i=rd:nbit
                   t(i)=cr(l,i);
               end
               for i=rd:nbit
                cr(l,i)=cr(u,i);
                cr(u,i)=t(i);
               end
            end
            t=[];
        end
     end


     %%%mutation
     nm=pop*nbit+1;
     for i=1:nm-1
         rd=rand();
         if rd<pm
            mc=fix(i/(nbit+1))+1;
            mb=i-((mc-1)*(nbit+1))+1;
         if  mb~=nbit+1 
             if cr(mc,mb)==0
                cr(mc,mb)=1;
             else
                cr(mc,mb)=0;
             end
         end
         end
     end
 
    %computing fitness of each CHROMOSOME and TOTAL fitness after GENETIC
    %operation
    Jflag=0; 
    i=0;
    while i<pop
         TEHH=0;
         TEHO=0;
         TEHS=0;
         if Jflag==0
            i=i+1;
         end
          %%%Transfer to X Y Z
          sx=0;
          for k=1:(nbit1) 
             sx=cr(i,k)*2^((nbit1)-k)+sx;
          end
          x(i)=mind(1,1)+(sx/maxdec1)*(maxd(1,1)-mind(1,1));
          sy=0;
          for k=(nbit1)+1:nbit2+nbit1
             sy=cr(i,k)*2^(nbit2+nbit1-k)+sy;
          end
           y(i)=mind(1,2)+(sy/maxdec2)*(maxd(1,2)-mind(1,2));
          sz=0;
          for k=nbit2+nbit1+1:nbit
             sz=cr(i,k)*2^(nbit-k)+sz;
          end
          z(i)=mind(1,3)+(sz/maxdec3)*(maxd(1,3)-mind(1,3));
    
     Ethane(i,1)=x(i);
     Ethane(i,2)=y(i);
     Ethane(i,3)=z(i);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   %%%%%calculating Distanse between atoms & ENERGY
   for j=1:pop
         rEthaneEthane(i,j)=norm(Ethane(i,:)-Ethane(j,:));      
   end
    for j=1:nsi
        rEthanesi(i,j)=norm(SI(j,:)-Ethane(i,:)); 
    end
    for j=1:nox
        rEthaneox(i,j)=norm(OX(j,:)-Ethane(i,:));
    end

      Jflag=0;
      for j=1:pop
        if i==j
             EEthaneEthane(i,j)=0;
        end
        if i~=j
         EEthaneEthane(i,j)=4*eps(3,3)*((sigma(3,3)/rEthaneEthane(i,j))^12-(sigma(3,3)/rEthaneEthane(i,j))^6);
        end
      TEHH=TEHH+EEthaneEthane(i,j);
      end
      EHH(i)=TEHH;
  
      for j=1:nox
       EEthaneox(i,j)=4*eps(1,3)*((sigma(1,3)/rEthaneox(i,j))^12-(sigma(1,3)/rEthaneox(i,j))^6);
       TEHO=TEHO+EEthaneox(i,j);
      end
      EHO(i)=TEHO;
    
      for j=1:nsi
       EEthanesi(i,j)=4*eps(2,3)*((sigma(2,3)/rEthanesi(i,j))^12-(sigma(2,3)/rEthanesi(i,j))^6);
       TEHS=TEHS+EEthanesi(i,j);
      end
      EHS(i)=TEHS;
    
      TE(i)=EHH(i)+EHO(i)+EHS(i);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      if TE(i)>0
         Jflag=1; 
         %%% extra mutation for fixing CHROMOSOME  
         for k=1:nbit
            rd=rand();
               if rd<pm
                  if cr(i,k)==0
                     cr(i,k)=1;
                  else
                  cr(i,k)=0;
                  end
               end
         end
      end
    
      if Jflag==1 & i==pop
         i=pop-1;
      end
      if Jflag==0
         fitness(i)=TE(i);
         cr(i,nbit+1)=-1*fitness(i);
      end
    end 
      sumf=sum(cr(:,nbit+1));
      STE2(itr)=sumf;%%%%%%%TOTAL fitness after GENETIC operation
          %%%Difrence between TOTAL fitness>>>> 
          if STE2(itr)<STE1(itr)
              cr=CR;
              STE2(itr)=STE1(itr);
          end
  end
  STEN=STE1(maxitr);
  if STEN<=STEO
       Kflag=1;
  else
      STEO=STEN;
      STES(kk,:)=STE1(:);
      kk=kk+1;
  end
            
  end
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%peng-robinson equation%%%%%%%%%%
  PP=P/(1.01325);%bar to atm
  apr=(0.457235*(R*TC)^2)/PC;
  bpr=0.077796*R*TC/PC;
  kpr=0.37464+1.54226*w-0.26992*w^2;
  alfapr=(1+kpr*(1-sqrt(T/TC)))^2;

  Apr=apr*alfapr*PP/((R*T)^2);
  Bpr=bpr*(PP/(R*T));

  a1=-(1-Bpr);
  a2=Apr-2*Bpr-3*Bpr^2;
  a3a=-Apr*Bpr+Bpr^2+Bpr^3;

  zzz=solve('zz^3+a1*zz^2+a2*zz+a3a=0');
  zzzz=subs(zzz);

  zg=max(zzzz);

  lnfc=(zg-1)-log(zg-Bpr)-((Apr)/(2*sqrt(2)*Bpr))*log((zg+(1+sqrt(2))*Bpr)/(zg+(1-sqrt(2))*Bpr));
  miux=Rg*T*lnfc;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  miu=miux+KB*T*log(pop);
  Y(pop)=-pop*miu+KB*T*log(factorial(pop))-STE1(maxitr);
  if pop>minpop
     deltaY=Y(pop)-Y(pop-1);
  end
  pop=pop+1
   end  
  MEthane=(pop-1)*2/avag;
  Mlattice=MEthane+dens*Vlattice*1000;
  teta(ii)=(MEthane/Mlattice)*100
  Pressure(ii)=P
  P=P+.001;
  ii=ii+1;
  minpop=pop-1;
  pop=minpop;
end
plot(Pressure,teta);
xlabel('pressure(bar)')
ylabel('teta(%wt)')
