# GAMS
genetic algorithm in molecular simulation

#Lattice Must import to Code

#number of adsorbent atoms

nox=length(OX(:,1));

nsi=length(SI(:,1));

#number of initial population

minpop=2;

pop=minpop;

#critcal point Of Ethane

TC=305.5;   %K

PC=48.72; %bar

w=.099;

R=82.056e-6;%atm.m3/gmol.K


#Temprature & Pressure &...

P=.001;%bar

T=300;%K

K=8.314;%J/gmole.K

KB=1.3806503e-23; %J/K

beta=1/(KB*T);

Rg=8.314;
%J/gmole.K

#parameter of GA

pc=.5;%Probability of cross over

pm=.2;%Probability of mutation


#parameters of LJ potential For H-SIOX

eps=[82.9848 44.2697 9.17;44.2697 23.6164 48.1;90.17 48.1 98];

sigma=[3.3 3.72 3.51;3.72 4.2 3.96;3.51 3.96 3.75];%K

#Structure of Adsorbent

mind=zeros(1,3);

for j=1:3

    maxd(j)=max(OX(:,j));

end


#Accuracy of structure of lattice

accrx=0.001;

accry=0.001;

accrz=0.001;

dx=(maxd(1,1)-mind(1,1))/accrx;

dy=(maxd(1,2)-mind(1,2))/accry;

dz=(maxd(1,3)-mind(1,3))/accrz;

#number of bits of X Y Z

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

#Lattice

avag=6.023e23;% avagadro number

Vlattice=1e-30;

for i=1:3

Vlattice=maxd(i)*Vlattice;%m^3

end

dens=.5e3;%Kg/m^3



nbit=nbit1+nbit2+nbit3;

maxdec1=2^nbit1-1;

maxdec2=2^nbit2-1;

maxdec3=2^nbit3-1;

ii=1;

while P<.011

deltaY=-1;

while deltaY<0

#Generating initial population
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
    
