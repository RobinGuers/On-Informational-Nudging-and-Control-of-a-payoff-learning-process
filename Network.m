classdef Network<handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
      C %Capacities
      TStar %Free Flow Travel Time
      ArcNumber
      
        
    end
    
    methods
        %Constructor
        function Nt = Network(C_,TStar_)
            if nargin > 0
            Nt.C = C_;
            Nt.TStar = TStar_;
            Stemp = size(C_);
            Nt.ArcNumber = Stemp(1);
            end
        end
        
        %TravelTime function that needs a demand distribution Nu_vec
        function [TrTim] = travelTime(Nt,Nu_Vec) 
        %It computes a BPR function travel time with (1+^2) Note that the
        %output is a vector
            TrTim = Nt.TStar.*(1+(Nu_Vec).*(Nu_Vec)./(Nt.C.*Nt.C));
        end;
        
        %Inverse Travel time function output the vector Nu_vec from travel
        %time
        function [Nu_Vec] = invTravelTime(Nt,T_vec)
            %Compute the inverse of BPR function
        Nu_Vec = Nt.C.*sqrt(T_vec.*Nt.TStar.^(-1)-1);
        
        end
        
%         function [CisEst] = estimatingCis(Nt,CisFormerStep,ArcToken,...
%                 TotalDemand,AveragingFactor)
%             %Arc Token is a vector where there is a one where the arc was
%             %token and 0 elsewhere
%             TrTim = Nt.travelTime(Nt,TotalDemand);
%             %The only Components that are going to be updated
%         end
    end
        methods(Static)
        function [Piest] = estimatingPie(PiestFormerStep,Outcome,AveragingFactor)
            %Outcome represent a vector with zero and one that describes
            %the choice
            %The averaging factor that shall be used is the fast averaging
            %factor
            Outcome;
            Piest=(1-AveragingFactor)*PiestFormerStep+AveragingFactor*Outcome;
            %+AveragingFactor*(PiestFormerStep).*not(Outcome) ;
        end
        
        function [CiPiEst] = estimatingPiCis(PiCisestFormerStep,Outcome,TravelTime,AveragingFactor)
            %Outcome represent a vector with zero and one that describes
            %the choice
            %The averaging factor that shall be used is the fast averaging
            %factor
            %TRAVEL TIME IS POSITIVE HENCE THE MINUS SIGN
            CiPiEst=(1-AveragingFactor)*PiCisestFormerStep-AveragingFactor*Outcome.*TravelTime; ...
                %+AveragingFactor*not(Outcome).*PiCisestFormerStep;
                
        end
        end

    
end
