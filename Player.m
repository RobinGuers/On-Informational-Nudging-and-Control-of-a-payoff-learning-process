classdef Player<handle
    %PLAYER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Beta;%Logit Parameter for choosing
        PayoffHist; %what our Currant player earned
        ProbVec;%What is the current Vec Probability
        ProbVecDesired;%The Final Desired Probability Vector
        LastChosenArc;%A zero 1 vector
        
    end
    
    methods
        %Constructor
        function Pla = Player(Beta_,PayoffHist_,ProbVec_,ProbVecDesired_)
            if nargin>0
            Pla.Beta= Beta_;
            Pla.PayoffHist = PayoffHist_;
            Pla.ProbVec = ProbVec_;
            Pla.ProbVecDesired = ProbVecDesired_;
            Pla.LastChosenArc=1;%This Choice is irrelevant since it will be modif
            end
        end
       
        %this function returns (0,0,0,..1,0,0) and chooses according to
        %probability
        function [dvec,arcChosen] = pickArc(Pla)
            val = random('Uniform',0,1,1,1);
            length = max(size(Pla.PayoffHist));
            sumProb = 0;
            count = 0;
            while(sumProb<val && count<length)
                count = count+1;
                sumProb = sumProb+Pla.ProbVec(count);
            end
            %count is the number of the chosen arc
            dvec = zeros(length,1);
            dvec(count) = 1;
            arcChosen = count;
            Pla.LastChosenArc=arcChosen;
        end
        
        %To compute the expected travel time on the arc whenplayer does not
        %take it one has to add 1 car to the true network state.
        
        function [CorrectedDemand] = correctedDemand(Pla,Totaldemand)
            %First Create a vector where there is one where player did not
            %take arc
            Dtemp=ones(size(Pla.ProbVec));
            Dtemp(Pla.LastChosenArc)=0;
            CorrectedDemand = Totaldemand+Dtemp;
            %We artificially augmented by one all the arc that have not
            %been taken to compute Cis.
        end
                
        %This function compute the new payoff as a weighted average of the
        %former and the new one
        function [NewPayoff] = updatePayoff(Pla,Piestimate,r,xstar,weightFactor)
            %[~,arcChosen] = Pla.pickArc()
            %First Compute the lies
            Y=0;
            l=AdaptiveLieDelta( Y,r,Pla.Beta,xstar,Piestimate );
            %l=zeros(size(Pla.ProbVec));
            
            %Do the lie for everyone
            NewPayoff = (1-weightFactor)*Pla.PayoffHist + ...
                weightFactor*(r+l);
            %Compute the one that is not updated through the lie
            NewPayoff(Pla.LastChosenArc)=(1-weightFactor)*Pla.PayoffHist(Pla.LastChosenArc) + ...
                weightFactor*r(Pla.LastChosenArc);
            %Update the Payoff
            Pla.PayoffHist=NewPayoff
            
            
        end
            
            
        
        function [newProb] = payoffToProb(Pla,newPayoff)
            %First Create a vector of one that has the dimension of
            %PayoffHist 
            temp1 = ones(max(size(Pla.PayoffHist)),1);
            temp2 = (Pla.Beta)*newPayoff;
            tempexp = exp(temp2);%create [exp(Bx1);exp(Bx2);...exp(Bxn)]
            weight = sum(tempexp);
            temp = (1/weight)*temp1;%Create [1/(exp(Bx1)+exp(Bx2)+...exp(Bxn)
            newProb = tempexp.*temp;
            %Update the probability
     
        end
        
        
        %For a given Probability Vector Gives the payoff
        function [DesiredPayoff] = getDesiredPayoff(Pla,WardropEq)
            dim = size(Pla.PayoffHist);
            x = sym('x',dim);
            assume(x<0)
            xBeta  = Pla.Beta*x;
            xExp = exp(xBeta);
            xWeight = sum(xExp);
            xProb = (1/xWeight)*xExp;
            %Now define the equation that is at stake
            EqVecExpr = xProb - Pla.ProbVecDesired;
%             EqVecExpr = sum(EqVecTemp.*EqVecTemp+(x+WardropEq).*(x+WardropEq));
            %Transform them into functions
            EqVecFun = matlabFunction(EqVecExpr,'vars',{x});
            %Now define the parameter for fmin
%             A =[];
%             b=[];
%             Aeq=[];
%             beq=[];
%             lb = -2*WardropEq;
%             ub = zeros(size(WardropEq));
            JacVecExpr = jacobian(EqVecExpr, x);
            JacVecFun = matlabFunction(JacVecExpr,'vars',{x});
            %define the objective function that will be fed into the solver
            function [fun1,jacfun1] = OBJ(xp)
            fun1 = EqVecFun(xp);
            jacfun1 = JacVecFun(xp);
            end
            
            %One needs to first guess the solution
            PayOffGuess = -WardropEq;
            
%             options =optimset('Jacobian','on');
              options =optimset('Jacobian','on');
              DesiredPayoff=fsolve(@OBJ,PayOffGuess,options);
        end
        
    end
    
end
