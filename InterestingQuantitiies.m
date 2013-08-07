function [ xmstar,xMstar,Delta,Pimax,MaxXmr ] = InterestingQuantities( r,xstar,Beta )
%INTERESTINGQUANITITIES Summary of this function goes here
%   Detailed explanation goes here
 [xmstar,indxmstar]=min(xstar);
 [xMstar,indxMstar]=max(xstar);
 [rM,indrM]=max(r);
 
 M=max(size(xstar));
 Delta=1e9;
 MaxXmr=-5;
 IndMaxXmr=[0;0];
 for k=1:M
     for j=1:M
         if (abs(xstar(k)-xstar(j))<Delta && (k-j)~=0 )
             Delta=abs(xstar(k)-xstar(j));
         end
         if (abs(xstar(k)-r(j))>MaxXmr)
             IndMaxXmr=[k;j];
             MaxXmr=abs(xstar(k)-r(j));
         end
     end
 end
         

 temp =ones(size(xstar));
 Deltavec2= (Delta/2)*temp;
 %First every entry is updated
 xstarvec=xstar;
 rVec=xstarvec;
 %Then the special one corresponding to M is updated.
 xstarvec(indxMstar)=xMstar;
 rVec(indrM)=rM;
 
 PiM1=exp(Beta*(xMstar))/(sum(exp(Beta*xstarvec)));
 PiM2=exp(Beta*(rM))/(sum(exp(Beta*rVec)));
 Pimax=max([PiM1,PiM2]);
 
end
