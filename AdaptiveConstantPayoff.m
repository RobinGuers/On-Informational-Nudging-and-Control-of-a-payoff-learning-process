clear all;
close all;
Beta_=1; %Logit parameter
ChoiceNumber=3; %Number of choice available to player
PlayerNumber=5; %Number of different players
StepNumber=100; %Number of step for simulating the discrete averaging process,
%the bigger it is the better the convergence is
alpha=1/10; % Parameter of the fast averaging factor gamma_f=1/n^(1/2+alpha)
GlobalPayoffAv=transpose(linspace(0,0.1,ChoiceNumber)); %The average GlobalPayoff,
%it is avergae because it will be perturbed by a martingale
xstar=[-1;0;1];% The payoff estimate desired by the end

Pistar=(1/sum(exp(Beta_*xstar)))*exp(Beta_*xstar);% Compute with the logit rule
PiNormal=(1/sum(exp(Beta_*GlobalPayoffAv)))*exp(Beta_*GlobalPayoffAv);


DvecAllPlayer=zeros(ChoiceNumber,PlayerNumber);
%Initializing the choice attributes, the code works for interrelated payoff
%hence we need to provide a Ca and Tstar which doesn't matter in the
%example we are investigating.Hence do not touch the next 3 lines
Ca = 10*linspace(1,ChoiceNumber,ChoiceNumber);%Do not touch
Tstar = (1/10)*linspace(100,100+10*ChoiceNumber,ChoiceNumber);%Do not touch
Network1 = Network(Ca,Tstar);%Do not touch

%Initializing an array of object Player that have different
%charachteristics
PlayerArray=Player.empty(PlayerNumber,0);
PayoffHistInit = -ones(ChoiceNumber,1);


%Preallocation for speed
Piest=zeros(ChoiceNumber,PlayerNumber);%Contain the estimate of the choice probability
DvecAllPlayer=zeros(ChoiceNumber,PlayerNumber);%Contains the last choice of each player
PiEstimatorPlayer=(1/ChoiceNumber)*ones(ChoiceNumber,PlayerNumber);


DiffEstTrue  = zeros(ChoiceNumber,PlayerNumber,StepNumber); %Error in the estimator
ConvToxstar  = zeros(ChoiceNumber,PlayerNumber,StepNumber); %Error in the payoff convergence
DivFromr = zeros(ChoiceNumber,PlayerNumber,StepNumber); %Divergence from the true reward
ConvToPistar = zeros(ChoiceNumber,PlayerNumber,StepNumber); %Convergence to Pistar
ConvToPiNormal = zeros(ChoiceNumber,PlayerNumber,StepNumber);


%Initialization of players' charachteristics
for p=1:PlayerNumber
    PlayerArray(p).Beta=Beta_;
    PayoffHistInit=zeros(ChoiceNumber,1);
    
    %Generate a random vector that sum to 1
    v=rand(ChoiceNumber,1);
    ProbVecInit = v/(sum(v));
    
    %As a convention for initializing the payoff values, the payoff value 2
    %is assigned to the choice that has highest probability.
    [value,Index]=max(ProbVecInit);
    PayoffHistInit(Index)=2;
    
    %Compute the initial payoff that gets on well wit the initial
    %estimation
    for k=1:ChoiceNumber
        if(k==Index)
            continue;
        else
            PiRatio=ProbVecInit(k)/ProbVecInit(Index);
            BetaDeltax=log(PiRatio);
            PayoffHistInit(k)=(1/Beta_)*BetaDeltax+PayoffHistInit(Index); %As a convention we add
            %the deltax value to the initial payoff
        end
    end
            
%     (1/sum(exp(Beta_*PayoffHistInit)))*exp(Beta_*PayoffHistInit)-ProbVecInit
%Preceeding line ensures that the passage between pi and payoff is OK,
%uncomment if you want to check

%Here we just initialize the values of the object
    PlayerArray(p).PayoffHist=PayoffHistInit;
    PlayerArray(p).ProbVec=ProbVecInit;
    PlayerArray(p).ProbVecDesired=(1/ChoiceNumber)*ones(ChoiceNumber,1);
end

%Preallocation for storing the averaging factors
slowFactorHistory = zeros(1,StepNumber);
fastFactorHistory = zeros(1,StepNumber);

for k=1:StepNumber
    
    %Compute and store averaging factors
    slowWeightFactor = (1/(k+1));
    fastWeightFactor = (1/(k+1))^(1/2+alpha);
    slowFactorHistory(k)=slowWeightFactor;
    fastFactorHistory(k)=fastWeightFactor;
    
    % Perturbation of the payoff by a 0 mean random term
    Errorterm=normrnd(0,1,ChoiceNumber,1);
    GlobalPayoff=GlobalPayoffAv+Errorterm;


    %Here Decision is taken for every user according to their choice
    %probability
    for p=1:PlayerNumber
        DvecAllPlayer(:,p)=PlayerArray(p).pickArc();
        %Note that pick arc also update the last chosen arc
    end
    DvecTot = sum(DvecAllPlayer,2);
    %Now look at the outcome of the game remember the more positive the payoff is the more chance it has to be chosen!!!!!
    

    
    
 
    
   
     %Update the payoff estimate
    for p=1:PlayerNumber
        [NewPayoff] = PlayerArray(p).updatePayoff(PiEstimatorPlayer(:,p),...
            GlobalPayoff,xstar,slowWeightFactor);
        [newProb] = PlayerArray(p).payoffToProb(NewPayoff);
        PlayerArray(p).ProbVec=newProb;
    end
    
    
        %Now Planner estimate the probability choice of each user based on
        %their past choice
    for p=1:PlayerNumber
    PiEstimatorPlayer(:,p)=Network1.estimatingPie(PiEstimatorPlayer(:,p),...
        DvecAllPlayer(:,p),fastWeightFactor);
    %Compute the difference between the true probability and the estimator 
        DiffEstTrue(:,p,k)=abs( PiEstimatorPlayer(:,p)-PlayerArray(p).ProbVec);
        
    %Compute the distance between interesting quantities
        ConvToxstar(:,p,k)=abs( PlayerArray(p).PayoffHist-xstar);
        DivFromr(:,p,k)=abs( PlayerArray(p).PayoffHist-GlobalPayoffAv);
        ConvToPistar(:,p,k) = abs(PlayerArray(p).ProbVec-Pistar);
        ConvToPiNormal(:,p,k) = abs(PlayerArray(p).ProbVec-PiNormal);
    end
    
   
end

%Prealllocation
PayoffFinal = zeros(ChoiceNumber,PlayerNumber);
TruePiFinal= zeros(ChoiceNumber,PlayerNumber);

%Gather the final probability and Payoff estimate ov every user in a single
%matrix

for p=1:PlayerNumber
    PayoffFinal(:,p)=PlayerArray(p).PayoffHist;
    PiFinal(:,p)=PlayerArray(p).ProbVec;
end

Absciss=linspace(1,StepNumber,StepNumber);

%Plot the averagin factor versus time , this is just to check that the two
%time scales assumption is respected
plot(Absciss,slowFactorHistory,Absciss,fastFactorHistory);
legend('Slow Factor','Fast Factor')
xlabel('Step')
ylabel('yvalue')


%Preallocation, we want to look at a global metric over all the choice,
%hence we reduce the difference by player
DiffEstTot=zeros(PlayerNumber,StepNumber);
ConvToxstarTot  = zeros(PlayerNumber,StepNumber);
DivFromrTot=zeros(PlayerNumber,StepNumber);
ConvToPistarTot=zeros(PlayerNumber,StepNumber);
ConvToPiNormalTot =zeros(PlayerNumber,StepNumber);
DiffEstMoy=zeros(1,StepNumber);

    
    for k=1:StepNumber
    DiffEstTot(:,k)=sum(DiffEstTrue(:,:,k),1);%Sum the difference of all choice
    %to create a global metric
    DivFromrTot(:,k)=sum(DivFromr(:,:,k),1);
    ConvToxstarTot(:,k)=sum(ConvToxstar(:,:,k),1);
    ConvToPistarTot(:,k)=sum(ConvToPistar(:,:,k),1);
    ConvToPiNormalTot(:,k)=sum(ConvToPiNormal(:,:,k),1);
    DiffEstMoy(1,k) = (1/PlayerNumber)*sum(DiffEstTot(:,k),1);
    end
    
    %Plot the convergence error of the probability estimate versus time for
    %every player
   figure;
   XConv=linspace(1,StepNumber,StepNumber);
   for p=1:PlayerNumber
       plot(XConv,DiffEstTot(p,:))
       hold on;
       title('Convergence for \alpha=1/10, \beta=1, PiNormal=[0.31,0.33,0.36] and Pistar=[0.1,0.24,0.66]');
       xlabel('Step Number');
       ylabel('Error on the probability estimate');
   end
   
   %Plot the average convergence error for the probability estimate, the average is
   %taken over the players
   figure;
      plot(XConv,DiffEstMoy);
       title('Average Convergence for \alpha=1/10, \beta=1, PiNormal=[0.31,0.33,0.36] and Pistar=[0.1,0.24,0.66]');
       xlabel('Step Number');
       ylabel('Average error on the probability estimate');
       
       
    %Plot the divergence between the true payoff and the actual payoff
    %versus time, the divergence is due to the adaptive nudging strategy 
       figure;
          for p=1:PlayerNumber
       plot(XConv,DivFromrTot(p,:))
       hold on;
       title('Divergence From r for \alpha=1/10, \beta=1, PiNormal=[0.31,0.33,0.36] and Pistar=[0.1,0.24,0.66]');
       xlabel('Step Number');
       ylabel('l1 norm to r');
          end
   
      %Plot the convergence to the desired payoff achieved through the adaptive nudging strategy    
          figure;
      for p=1:PlayerNumber
       plot(XConv,ConvToxstarTot(p,:))
       hold on;
       title('Convergence to xstar for \alpha=1/10, \beta=1, PiNormal=[0.31,0.33,0.36] and Pistar=[0.1,0.24,0.66]');
       xlabel('Step Number');
       ylabel('l1 norm to xstar');
      end
   
     %Plot the convergence to the desired probability achieved through the adaptive nudging strategy 
          figure;
       for p=1:PlayerNumber
       plot(XConv,ConvToPistarTot(p,:),XConv,ConvToPiNormalTot(p,:))
       hold on;
       title('Convergence to PiValues for \alpha=1/10, \beta=1, PiNormal=[0.31,0.33,0.36] and Pistar=[0.1,0.24,0.66]');
       xlabel('Step Number');
       ylabel('l1 norm to \pi');
       legend('Pistar','PiNormal')
          end
