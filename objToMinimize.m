function [ val,Grad ] = objToMinimize( l,pistar1 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

val=pistar1*l(1)^2+(1-pistar1)*l(2)^2;
Grad=[2*pistar1*l(1);2*(1-pistar1)*l(2)];

end
