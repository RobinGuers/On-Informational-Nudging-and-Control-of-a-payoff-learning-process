function [ Adlie,ddta,MaxAdlie ] = AdaptiveLieDelta( Y,r,Beta,xstar,Prob_ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%Delta=abs(xMstar-xmstar);
[ xmstar,xMstar,Delta,Pimax,MaxXmr ] = InterestingQuantities( r,xstar,Beta );
%This is to find the min Value of (xa1-xa2)
% Delta=100;
% for a1=1:max(size(xstar))
%     for a2=1:max(size(xstar))
%         if(abs(xstar(a1)-xstar(a2))<=Delta)
%             Delta=abs(xstar(a1)-xstar(a2));
%         end
%     end
% end

if(Y~=0)
[Prob] = payoffToProb(Y,Beta);
else
    Prob=Prob_;
end
 temp = ones(max(size(xstar)),1);
% temp1 =(xmstar-Delta/2)*temp;
% temp1(1)=xMstar;
% [Probtemp1] = payoffToProb(temp1,Beta)
ddta=min([1-Pimax,0.5]);

CurrentProb= Prob;
%Look at the maximum value of the current probability
[MaxCurr,indMaxCurr]=max(CurrentProb);


if((1-ddta)<=MaxCurr)
    CurrentProb(indMaxCurr)=1-ddta;
end

Adlie=(xstar-r)./(temp-CurrentProb);
MaxAdlie=MaxXmr/(1-Pimax)
    


end
