function [ Prob ] = PayoffToProbDeltaX( DeltaX,Beta )
%UNTITLED Summary of this function goes here
%This function only works for a 3 choice alternative
%DeltaX is given by DeltaX(1)=x2-x1, DeltaX(2)=x3-x1
%   Detailed explanation goes here
            Prob = ones(max(size(DeltaX))+1,1);
            %Use of the logit Rule
            Prob(1)=1/(1+exp(Beta*DeltaX(1))+exp(Beta*DeltaX(2)));
            Prob(2)=1/(1+exp(-Beta*DeltaX(1))+exp(Beta*(DeltaX(2)-DeltaX(1))));
            Prob(3)=1-(Prob(1)+Prob(2));

     

end
