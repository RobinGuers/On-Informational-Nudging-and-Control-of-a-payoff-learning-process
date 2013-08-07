function [ value ] = objfunDeltaX( DeltaX,R,L,Beta )
%UNTITLED2 Summary of this function goes here
%The value output by this function is the amount it moves.
%   Detailed explanation goes here
Prob = PayoffToProbDeltaX( DeltaX,Beta );
DeltaR = [R(2)-R(1);R(3)-R(1)];
valtemp=ones(2,1);
%Valtemp is the differential equation driving deltax
valtemp(1)=-DeltaX(1)+DeltaR(1)+L(2)*(1-Prob(2))-L(1)*(1-Prob(1));
valtemp(2)=-DeltaX(2)+DeltaR(2)+L(3)*(1-Prob(3))-L(1)*(1-Prob(1));
value=sum(valtemp.^2);

end
