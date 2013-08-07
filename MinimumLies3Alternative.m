%This script makes use of the plotregion function availbale on the matlab
%exchange
close all;
clear all;


M=5;

Beta=1; %Logit Parameter
R=[1;2;3];% True reward for the 3 alternatives, can be modified
lb=[-4/Beta,-4/Beta,-4/Beta]; %Lower bound for the LP because sufficient conditions
%derived by Gersgorin theorem imply l_a>-4/Beta for all a


%This represent all the inequalities needed to represent the union of the
%24 convex sets.
%We define M=max_(p!=q) |l_p+l_q| so we need to consider independantly all
%the possible cases
%l_{1}+l_{2}><0 ... l_{1}+l_{3}><0 ... l_{2}+l_{3}><0
%This makes 8 different cases, each of which can produce 3 different cases
%because M=l_1+l_2 or l_2+l_3 or l_1+l_3
%We hence have 24 cases to deal with each of which can be written as a
%linear inequality

b1=[0;0;-8/Beta;-8/Beta;-8/Beta;0;0;0]; %Inequalities which must be satisfied for the LP,
% The 2 first element represent the inequality M=l_p+l_q ie we choose 1
% possibility among the3
%The 3 middle element represent -M/2+l_a>-4/Beta
%The last 3 elements deal with the signs of l_1+l_2,l_2+l_3,l_1+l_3

B1=[1,1,0;1,0,1;0,1,1];%l_1+l_2>0, l_1+l_3>0, l_2+l_3>0
B2=[1,1,0;1,0,1;0,-1,-1];%l_1+l_2>0, l_1+l_3>0, -(l_1+l_3)>0
B3=[1,1,0;-1,0,-1;0,1,1];%l_1+l_2>0, -(l_1+l_3)>0, l_2+l_3>0
B4=[1,1,0;-1,0,-1;0,-1,-1];%l_1+l_2>0, -(l_1+l_3)>0, -(l_2+l_3)>0

%Those matrices represent the inequalities who is M and -M/2+l_a>-4/Beta

%Cases A1* deal with l_1+l_2>0, l_1+l_3>0, l_2+l_3>0, each* case represent
%a different maximum M, A11 represents l_1+l_2>l_1+l_3 and l_1+l_2>l_2+l_3,
%A12 represents l_1+l_3>l_1+l_2 and l_1+l_3>l_2+l_3,
%A13 represents l_2+l_3>l_1+l_2 and l_2+l_3>l_1+l_3
A11=[0,1,-1;1,0,-1;1,-1,0;-1,1,0;-1,-1,2];
A12=[0,-1,1;1,-1,0;1,0,-1;-1,2,-1;-1,0,1];
A13=[0,1,-1;1,2,1;1,-1,0;-1,1,0;-1,-1,2];

%Cases A2* deal with l_1+l_2>0, l_1+l_3>0, -(l_2+l_3)>0,
A21=[0,1,-1;1,2,1;1,-1,0;-1,1,0;-1,-1,2];
A22=[0,-1,1;1,1,2;1,0,-1;-1,2,-1;-1,0,1];
A23=[-1,-2,-1;-1,-1,-2;2,1,1;0,3,1;0,1,3];

%Cases A3* deal with l_1+l_2>0, -(l_1+l_3)>0, l_2+l_3>0
A31=[2,1,1;1,0,-1;1,-1,0;-1,1,0;-1,-1,2];
A32=[-2,-1,-1;-1,-1,-2;3,0,1;1,2,1;1,0,3];
A33=[-1,0,1;1,1,2;2,-1,-1;0,1,-1;0,-1,1];

%Cases A4* deal with l_1+l_2>0, -(l_1+l_3)>0, -(l_2+l_3)>0
A41=[2,1,1;1,2,1;1,-1,0;-1,1,0;-1,1,2];
A42=[-2,-1,-1;-1,1,0;3,0,1;1,2,1;1,0,3];
A43=[-1,-2,-1;1,-1,0;2,1,1;0,3,1;0,1,3];

%Cases A5* deal with -(l_1+l_2)>0, l_1+l_3>0, l_2+l_3>0
A51=[-2,-1,-1;-1,-2,-1;3,1,0;1,3,0;1,1,2];
A52=[2,1,1;1,-1,0;1,0,-1;-1,2,-1;-1,0,1];
A53=[1,2,1;-1,1,0;2,-1,-1;0,1,-1;0,-1,1];

%Cases A6* deal with -(l_1+l_2)>0, l_1+l_3>0, -(l_2+l_3)>0
A61=[-2,-1,-1;-1,0,1;3,1,0;1,3,0;1,1,2];
A62=[2,1,1;1,1,2;1,0,-1;-1,2,-1;-1,0,1];
A63=[1,0,-1;-1,-1,-2;2,1,1;0,3,1;0,1,3];

%Cases A7* deal with -(l_1+l_2)>0, -(l_1+l_3)>0, l_2+l_3>0
A71=[0,-1,1;-1,-2,-1;3,1,0;1,3,0;1,1,2];
A72=[0,1,-1;-1,-1,-2;3,0,1;1,2,1;1,0,3];
A73=[1,2,1;1,1,2;2,-1,-1;0,1,-1;0,-1,1];

%Cases A8* deal with -(l_1+l_2)>0, -(l_1+l_3)>0, -(l_2+l_3)>0
A81=[0,-1,1;-1,0,1;3,1,0;1,3,0;1,1,2];
A82=[0,1,-1;-1,1,0;3,0,1;1,2,1;1,0,3];
A83=[1,0,-1;1,-1,0;2,1,1;0,3,1;0,1,3];

%Concatenate the matrix A** with their B* matrix corresponding
A11=[A11;B1];
A12=[A12;B1];
A13=[A13;B1];

A21=[A21;B2];
A22=[A22;B2];
A23=[A23;B2];

A31=[A31;B3];
A32=[A32;B3];
A33=[A33;B3];

A41=[A41;B4];
A42=[A42;B4];
A43=[A43;B4];

A51=[A51;-B4];
A52=[A52;-B4];
A53=[A53;-B4];

A61=[A61;-B3];
A62=[A62;-B3];
A63=[A63;-B3];

A71=[A71;-B2];
A72=[A72;-B2];
A73=[A73;-B2];

A81=[A81;-B1];
A82=[A82;-B1];
A83=[A83;-B1];

%Preallocate Aglob for speed
Aglob=zeros(8,3,24);

%Allocate the corresponding matrices
Aglob(:,:,1)=A11;
Aglob(:,:,2)=A12;
Aglob(:,:,3)=A13;
Aglob(:,:,4)=A21;
Aglob(:,:,5)=A22;
Aglob(:,:,6)=A23;
Aglob(:,:,7)=A31;
Aglob(:,:,8)=A32;
Aglob(:,:,9)=A33;
Aglob(:,:,10)=A41;
Aglob(:,:,11)=A42;
Aglob(:,:,12)=A43;
Aglob(:,:,13)=A51;
Aglob(:,:,14)=A52;
Aglob(:,:,15)=A53;
Aglob(:,:,16)=A61;
Aglob(:,:,17)=A62;
Aglob(:,:,18)=A63;
Aglob(:,:,19)=A71;
Aglob(:,:,20)=A72;
Aglob(:,:,21)=A73;
Aglob(:,:,22)=A81;
Aglob(:,:,23)=A82;
Aglob(:,:,24)=A83;



%Create a meshgrid for Pi1 and Pi2
PiMax=0.95; %Maximum value of Pi1
Pimin=0.05; %Maximum value of Pi2
NumPoint=5; %Number of point in the Mesh
%Building of the Mesh
Pivec=linspace(Pimin,PiMax,NumPoint);
[Pi1grid,Pi2grid]= meshgrid(Pivec,Pivec);



%Bootstrap to stay in the simplex,
%it checks that Pi stays in the simplex i.e. Pi1+Pi2<1
%Note that Pi1<0.95 and so is Pi2, hence when Pi2grid=0.01 we get
%Pi1grid+Pi2grid=0.96
for k=1:NumPoint
    for j=1:NumPoint
        if ((Pi1grid(k,j)+Pi2grid(k,j))>=0.95)
            vec=[Pi1grid(k,j);Pi2grid(k,j)];
            [valmax,indmax]=max(vec);
            if(indmax==1)
                Pi2grid(k,j)=0.01;
            else
                Pi1grid(k,j)=0.01;
            end
        end
    end
end

%Build Pi3 as 1-Pi1-Pi2
Pi3grid=(ones(NumPoint,NumPoint)-Pi2grid-Pi1grid);
%Compute the difference x2-x1 and x3-x1 in vectors
X21grid=(1/Beta)*log(Pi2grid./Pi1grid);
X31grid=(1/Beta)*log(Pi3grid./Pi1grid);

%Preallocation            
valvec=zeros(NumPoint,NumPoint); %Store the value of the objective function
lvecglob=-100*ones(NumPoint,NumPoint,3); %Store the value of the lies obtained through the optimization process
CredValueGlob=10*ones(NumPoint,NumPoint); %Store the value of the credibility, the weighted average of the lies
%by the probability they will be encountered
Aset=zeros(4,24,NumPoint,NumPoint);
AllPossibleLies=zeros(4,24);
%Set the tolerance, after obtaining a lie in the stability region, the second part of the algorithm
%is to check if the lie actually enforce an equilibrium or if it is too
%far.
Tol=1e-3;

for k=1:NumPoint
    for j=1:NumPoint
        %This loop computes the values of the lies that minimize the
        %objective functtion. The objective function measures how far from
        %rest point is the differential equation. p goes from 1 to 24
        %because the admissible region is the union of 24 convex sets.
  for p=1:24  
options = optimset('Algorithm','interior-point');
x0=[X21grid(k,j),X31grid(k,j)];
Pivec=[Pi1grid(k,j);Pi2grid(k,j);1-(Pi1grid(k,j)+Pi2grid(k,j))];
%objfunDeltaX outputs how far from the equilibrium point one is
[lvec,val] = fmincon(@(L) objfunDeltaX( x0,R,L,Beta ),[0;0;0],-Aglob(:,:,p),-b1,[],[],lb,[],options);
valvec(k,j)=val;
AllPossibleLies(1,p)=val;
AllPossibleLies(2:end,p)=lvec;
[ CredValue ] = QuadNudgeExpected( lvec,Pivec );
Aset(:,p,k,j)=[val;lvec];

  end
  %We know check how far stands the differential equation from the rest
  %point, and discards the solution if it superior to Tol previously defined.
  [value,indexvalue]=min(AllPossibleLies,[],2);
  if (value(1)<Tol)
      CredValueGlob(k,j)=sum(Pivec.*AllPossibleLies(2:4,indexvalue(1)).^2);
      lvecglob(k,j,:)=AllPossibleLies(2:4,indexvalue(1));
  else
      %As a convention the credibility value will be -10 when no feasible
      %lies to enforce the equilibrium can be found
      CredValueGlob(k,j)=-10;
  end
        end
end

%We plot the desired data,

%Plot the value of the credibility obtained,
%l1*\pi_{1}^{*}^2+l2*\pi_{2}^{*}^2+l3*\pi_{3}^{*}^2
figure;
surf(Pi1grid,Pi2grid,CredValueGlob)
xlabel('\pi_{1}^{*}');
ylabel('\pi_{2}^{*}');
zlabel('Credibility Value');
title('Value of the Credibility at Equilibrium');

%Plot the values of l1 for obtaining \pi_{1}^{*} and \pi_{2}^{*}
figure;
surf(Pi1grid,Pi2grid,lvecglob(:,:,1))
xlabel('\pi_{1}^{*}');
ylabel('\pi_{2}^{*}');
zlabel('Credibility Value');
title('Value of l_{1} at Equilibrium');

%Plot the value of l2 for obtaining \pi_{1}^{*} and \pi_{2}^{*}
figure;
surf(Pi1grid,Pi2grid,lvecglob(:,:,2))
xlabel('\pi_{1}^{*}');
ylabel('\pi_{2}^{*}');
zlabel('Credibility Value');
title('Value of l_{2} at Equilibrium');

%Plot the value of l3 for obtaining \pi_{1}^{*} and \pi_{2}^{*}
figure;
surf(Pi1grid,Pi2grid,lvecglob(:,:,3))
xlabel('\pi_{1}^{*}');
ylabel('\pi_{2}^{*}');
zlabel('Credibility Value');
title('Value of l_{3} at Equilibrium');
