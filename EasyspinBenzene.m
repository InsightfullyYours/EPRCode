function EasyspinBenzene(ATa,AMTa,APMTa)
%Sys.s = [1;1/2];
clear FitOpt UserCommand h
clear N15_Ta N15_TaVary N15_TaExp
clear N15_MTa N15_MTaVary N15_MTaExp
clear NH_MTa NH_MTaVary NH_MTaExp
clear N7_MTa N7_MTaVary N7_MTaExp
clear N14_2_MTa N14_2_MTaVary N14_2_MTaExp
clear N15_Ta N15_TaVary N15_TaExp
clear N15_PMTa N15_PMTaVary N15_PMTaExp
clear NH_PMTa NH_PMTaVary NH_PMTaExp
clear N7_PMTa N7_PMTaVary N7_PMTaExp
clear N14_2_PMTa N14_2_PMTaVary N14_2_PMTaExp



options=optimset('Display','iter','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-8,'TolFun',1E-8);


%--------------------------------------------------------------------------

%TaMTa is the pure spin trap in toluene and benzophenone-pure N peak
%(accounts for shifts)
N15_MTa.g = [2.00813];
N15_MTa.Nucs = 'N';
N15_MTa.n=[1];
N15_MTa.A = [mt2mhz(1.5373,N15_MTa.g)];
N15_MTa.lwpp = [.075];

N15_MTaVary.g = [.001];
N15_MTaVary.A = [0];
N15_MTaVary.lwpp = .000;

N15_MTaExp.mwFreq = 9.7560;
N15_MTaExp.CenterSweep = [347.4 10.0];

[BN15_MTa,SN15_MTa]=garlic(N15_MTa,N15_MTaExp);
%20


%Pure spin trap and benzophenone irradiated-N-H interaction peak
NH_MTa.g = [2.0087];
NH_MTa.Nucs = 'N,H';
NH_MTa.n=[1 1];
NH_MTa.A = [mt2mhz(1.426,NH_MTa.g);mt2mhz(.445,NH_MTa.g)];
NH_MTa.lwpp = .03;

NH_MTaVary.g = [.004];
NH_MTaVary.A = [mt2mhz(.1,NH_MTa.g);mt2mhz(.1,NH_MTa.g)];
NH_MTaVary.lwpp = 0.02;

NH_MTaExp.mwFreq = 9.7560;
NH_MTaExp.CenterSweep = [347.4 10.0];

[BNH_MTa,SNH_MTa]=garlic(NH_MTa,NH_MTaExp);
%3.5

%Pure spin trap and benzophenone irradiated-7.5 gauss ketone peak
N7_MTa.g = [2.00889];
N7_MTa.Nucs = 'N';
N7_MTa.n=[1];
N7_MTa.A = [mt2mhz(.8038,N7_MTa.g)];
N7_MTa.lwpp = .07;

N7_MTaVary.g = [0];
N7_MTaVary.A = [mt2mhz(.01,N7_MTa.g)];
N7_MTaVary.lwpp = 0.02;

N7_MTaExp.mwFreq = 9.7560;
N7_MTaExp.CenterSweep = [347.4 10.0];

[BN7_MTa,SN7_MTa]=garlic(N7_MTa,N7_MTaExp);
%12

%Pure spin trap and benzophenone irradiated-mystery peak
N14_2_MTa.g = [2.00829];
N14_2_MTa.Nucs = 'N';
N14_2_MTa.n=[1];
N14_2_MTa.A = [mt2mhz(1.3691,N14_2_MTa.g)];
N14_2_MTa.lwpp = .35;

N14_2_MTaVary.g = [.001];
N14_2_MTaVary.A = [mt2mhz(.1,N14_2_MTa.g);mt2mhz(.2,N14_2_MTa.g)];
N14_2_MTaVary.lwpp = .2;

N14_2_MTaExp.mwFreq = 9.7560;
N14_2_MTaExp.CenterSweep = [347.4 10.0];

[BN14_2_MTa,SN14_2_MTa]=garlic(N14_2_MTa,N14_2_MTaExp);
%3.5

start_point = rand(1, 4);
estimates2(1,:) = patternsearch(@expfun2, start_point,[],[],[],[],[0.02 0 0 0],[1 1 1 1],[],options);
%fitted2(:,1)=estimates2(1,1).*SN15_MTa+estimates2(1,2).*SNH_MTa+estimates2(1,3).*SN7_MTa+estimates2(1,4).*SN14_2_MTa;
fitted2(:,1)=0.*SN15_MTa+estimates2(1,2).*SNH_MTa;%+estimates2(1,3).*SN7_MTa+0.*SN14_2_MTa;


figure
subplot(2,1,1)
plot(BN14_2_MTa,fitted2.*.7,'r')
hold
plot(BN14_2_MTa,AMTa(:,10)./10000)
ylabel('Intensity')
subplot(2,1,2)
plot(BN14_2_MTa,AMTa(:,10)./10000-fitted2.*.7)
xlabel('B (mT)')
title('Subtraction of Fit from Data')
%-------------------------------------------------------------------------


FitOpt.Method = 'genetic int';
[bestsys,bestspc]=esfit('garlic',AMTa(:,10),NH_MTa,NH_MTaVary,NH_MTaExp,[],FitOpt);
bestspc=bestspc';
mhz2mt(bestsys.A,bestsys.g)

    function sse = expfun1(params)
        A = params(1);
        FittedCurve=A.*SN15_Ta';
        ErrorVector = FittedCurve - ATa(:,8)./10000;
        sse = sum(ErrorVector .^ 2);
    end

    function sse = expfun2(params)
        A = 0;
        B = params(2);
        C = params(3);
        D = 0;
        FittedCurve=A.*SN15_MTa'+B.*SNH_MTa'+C.*SN7_MTa'+D.*SN14_2_MTa';
        ErrorVector = FittedCurve - AMTa(:,10)./10000;
        sse = sum(ErrorVector .^ 2);
    end

    function sse = expfun3(params)
        A = params(1);
        B = params(2);
        C = params(3);
        D = params(4);
        FittedCurve=A.*SN15_Ta'+B.*SNH_PMTa'+C.*SN7_PMTa'+D.*SN14_2_PMTa';
        ErrorVector = FittedCurve - APMTa(:,10)./10000;
        sse = sum(ErrorVector .^ 2);
    end

disp('Done')

clear FitOpt UserCommand h
clear N15_Ta N15_TaVary N15_TaExp
clear N15_MTa N15_MTaVary N15_MTaExp
clear NH_MTa NH_MTaVary NH_MTaExp
clear N7_MTa N7_MTaVary N7_MTaExp
clear N14_2_MTa N14_2_MTaVary N14_2_MTaExp
clear N15_Ta N15_TaVary N15_TaExp
clear N15_PMTa N15_PMTaVary N15_PMTaExp
clear NH_PMTa NH_PMTaVary NH_PMTaExp
clear N7_PMTa N7_PMTaVary N7_PMTaExp
clear N14_2_PMTa N14_2_PMTaVary N14_2_PMTaExp
end