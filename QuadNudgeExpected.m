function [ CredValue ] = QuadNudgeExpected( lvec,Pivec )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
CredValue=sum((lvec.^2).*Pivec);

end
