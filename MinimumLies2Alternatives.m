close all;
clear all;

Beta=1;%Logit Parameter Value
x0=[5;5];%Starting Point for the optimization
EpsTol=1e-2;%Tolerance
%Input here where you want to begin and end and the precision in the
%probability governed by Pistep
Pibegin=0.05;
Piend=0.95;
Pistep=10;%Precision for the Pivector
%Create a meshgrid for P21vecstar
P21VecStar=linspace(Pibegin,Piend,Pistep);
P21vecNormal=linspace(Pibegin,Piend,Pistep);
%Preallocation
Data=zeros(Pistep,Pistep); %Value of the objective function at optimum
lTildeHistoryMat=zeros(Pistep,Pistep); %Value of the transformed lie at equilibrium
R2TildeHistoryMat=zeros(Pistep,Pistep); % Difference between the payoffs value
exitFlagMat=zeros(Pistep,Pistep); % Store the exitflag of the optimization
l1HistoryMat=zeros(Pistep,Pistep); % Value of l1 at equilibrium
l2HistoryMat=zeros(Pistep,Pistep); % Value of l2 at equilibrium

R21ml1 = zeros(Pistep,Pistep); %R21-l1 value
Ftil= zeros(Pistep,Pistep);% Value of the modified objective function
ConEqSatis = zeros(Pistep,Pistep); %Check if the equality constraints are satisfied
ConIneqSatis=zeros(Pistep,Pistep); %Check if the inequality constraints are satisfied

for k=1:Pistep %The index of the normal probability
    for i=1:Pistep %The index of the influenced prob
        options = optimset('Algorithm','interior-point');
        Pistar1=P21VecStar(i);
        %Compute the payoff difference based on the original probability
        R21=(1/Beta)*log(1/P21vecNormal(k)-1);
        X21 = (1/Beta)*log(1/P21VecStar(i)-1);
        %Inequality and equality constraints as defined in the paper for the QP
        A=-[1,1];
        b=4/Beta-0.1;
        beq=R21-X21;
        Aeq=[1-Pistar1,-Pistar1];
        x0=[1;1];
        
        [l,fval,exitflag]=fmincon(@(l)objToMinimize(l,Pistar1),x0,A,b,Aeq,beq,[],[],[],options);
%x0is defined at the beginning of the paper
        Data(k,i)=fval;
       lTildeHistoryMat(k,i)=sum(l);
       R2TildeHistoryMat(k,i)=R21-l(1);
       l1HistoryMat(k,i)=l(1);
       l2HistoryMat(k,i)=l(2);
       exitFlagMat(k,i)=exitflag;
       R21ml1(k,i)=R21-X21;
       Ftil(k,i)=(1-Pistar1)*l(1)-Pistar1*l(2);
      ConIneqSatis(k,i)=l(1)+l(2);
    end
    
end

%Plot the various interesting quantities, the minimum credibility value,
%the convergence check,....

 ConEqSatis=Ftil-R21ml1;
[X,Y]=meshgrid(P21VecStar,P21vecNormal);
% X is the prob without intervention, Y is the Prob with intervention
figure;
surf(X,Y,Data)
xlabel('\pi_{1}^{*}');
ylabel('\pi_{1}');
zlabel('\pi_{1}^{*}*l_{1}^{2}+(1-\pi_{1}^{*})*l_{2}^{2}');
title(' Minimum Credibility for \beta=1')


figure;
surf(X,Y,l1HistoryMat);
xlabel('\pi_{1}^{*}');
ylabel('\pi_{1}');
zlabel('l_{1}');
title('Minimum l_{1} for \beta=1')


figure;
surf(X,Y,l2HistoryMat);
xlabel('\pi_{1}^{*}');
ylabel('\pi_{1}');
zlabel('l_{2}');
title('Minimum l_{2} for \beta=1')


figure;
[C,h] = contour(X,Y,Data);
xlabel('\pi_{1}^{*}');
ylabel('\pi_{1}');
title('Value of the Lying Metric \beta=1')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
colormap cool


lTildeHistoryVec=zeros(Pistep*Pistep,1);
R2TildeHistoryVec=zeros(Pistep*Pistep,1);
for i=1:Pistep
    lTildeHistoryVec((i-1)*Pistep+1:i*Pistep)=lTildeHistoryMat(:,i);
    R2TildeHistoryVec((i-1)*Pistep+1:i*Pistep)= R2TildeHistoryMat(:,i);
end
