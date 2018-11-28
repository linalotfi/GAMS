# GAMS
genetic algorithm in molecular simulation

#Lattice Must import to Code

# Determine number of adsorbent atoms 

nox=length(OX(:,1));

nsi=length(SI(:,1));

# specify number of initial population

minpop=2;

pop=minpop;

# an exapmle : critcal point Of Ethane

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

#parameter of Genetic Algoritm can specify by the following codes

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

while P<Pmax

deltaY=-1;

while deltaY<0
 
 Step 1
 
Generating initial population
               
Generating  random CHROMOSOME 
             
Transfer CHROMOSOME to X Y Z

Step 2

Calculate Toltal Energy

Step 3

Calculate chemical potential from an Equation of States

Step 4
  
  Calculate Fitness Function
  
Step 5 
 
 Add new Chromosomes
 
 end
 
  increase Pressure
 
 end
