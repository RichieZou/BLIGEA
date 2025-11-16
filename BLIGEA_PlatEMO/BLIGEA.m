classdef BLIGEA < ALGORITHM
% <multi> <real/binary> <large/none> <constrained/none> <sparse>
% A sparse large-scale multi-objective evolutionary optimization based on
% bi-level interactive grouping

%------------------------------- Reference --------------------------------
% Y. Zou, J. Zou, S. Wang, Y. Liu, and S. Yang. A sparse large-scale
% multi-objective evolutionary optimization based on bi-level interactive
% grouping. Swarm and Evolutionary Computation, 2025, 99: 102209.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting and fitness initialization
            TDec    = [];
            TMask   = [];
            TempPop = [];
            Fitness = zeros(1,Problem.D);
            for i = 1 : 1+4*any(Problem.encoding~=4)
                Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                Dec(:,Problem.encoding==4) = 1;
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                Fitness    = [Fitness ; NDSort([Population.objs,Population.cons],inf)];
            end
            Fitness = median(Fitness);
            
            %% Generate initial population
            Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            Dec(:,Problem.encoding==4) = 1;
            Mask = zeros(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
            end
            Population = Problem.Evaluation(Dec.*Mask);
            sv = zeros(1,Problem.D);
            Last_temp_num = 0;
            t = mean(Fitness);
            [Population,Dec,Mask,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                iter = Problem.FE/Problem.maxFE;

                %% Bi-level fitness and variable update (Innovation 1)
                First_Mask = Mask(FrontNo==1,:);
                [temp_num,~] = size(First_Mask);
                temp_vote = sum(First_Mask,1);
                sv(1,:) = (Last_temp_num/(Last_temp_num+temp_num))*sv(1,:) + (temp_num/(Last_temp_num+temp_num))*(temp_vote/temp_num);
                Last_temp_num = temp_num;
                fv = std(Population(FrontNo==1).decs,0,1);
                Fitness = max(1,Fitness-sqrt(t)*mpower(iter,1)*sv);
                Fitness = min(max(Fitness),Fitness+sqrt(t)*mpower(iter,1)*(1-sv));
                
                %% Mating pool selection and variation
                MatingPool       = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness,iter,sv,fv);
                Offspring        = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
            end
        end
    end
end