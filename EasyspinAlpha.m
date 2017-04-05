function EasyspinAlpha(ATa,AMTa,APMTa)
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

%Ta is the pure spin trap in toluene-Pure N peak
N15_Ta.g = [2.00991];
N15_Ta.Nucs = 'N';
N15_Ta.n=[1];
N15_Ta.A = [mt2mhz(1.5373,N15_Ta.g)];
N15_Ta.lwpp = [.075];

N15_TaVary.g = [1];
N15_TaVary.A = [1];
N15_TaVary.lwpp = 1;

N15_TaExp.mwFreq = 9.7360;
N15_TaExp.CenterSweep = [346.0 10.0];

[BN15_Ta,SN15_Ta]=garlic(N15_Ta,N15_TaExp);

start_point = rand(1, 1);
estimates1(1,:) = patternsearch(@expfun1, start_point,[],[],[],[],0,1,[],options);
fitted1(:,1)=estimates1(1,1).*SN15_Ta;

figure
subplot(2,1,1)
plot(BN15_Ta,fitted1,'r')
hold
plot(BN15_Ta,ATa(:,8)./10000)
ylabel('Intensity')
subplot(2,1,2)
plot(BN15_Ta,ATa(:,8)./10000-fitted1)
xlabel('B (mT)')
title('Subtraction of Fit from Data')


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
NH_MTa.g = [2.00829];
NH_MTa.Nucs = 'N,H';
NH_MTa.n=[1 1];
NH_MTa.A = [mt2mhz(1.3672,NH_MTa.g);mt2mhz(.19,NH_MTa.g)];
NH_MTa.lwpp = .1;

NH_MTaVary.g = [0];
NH_MTaVary.A = [mt2mhz(.1,NH_MTa.g);mt2mhz(.1,NH_MTa.g)];
NH_MTaVary.lwpp = 0.05;

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
fitted2(:,1)=0.*SN15_MTa+estimates2(1,2).*SNH_MTa+estimates2(1,3).*SN7_MTa+estimates2(1,4).*SN14_2_MTa;


figure
subplot(2,1,1)
plot(BN14_2_MTa,fitted2.*.7,'r')
hold
plot(BN14_2_MTa,AMTa(:,8)./10000)
ylabel('Intensity')
subplot(2,1,2)
plot(BN14_2_MTa,AMTa(:,8)./10000-fitted2.*.7)
xlabel('B (mT)')
title('Subtraction of Fit from Data')
%-------------------------------------------------------------------------

% pure spin trap in toluene and benzophenone-pure N peak
%(accounts for shifts)
N15_PMTa.g = [2.00905];
N15_PMTa.Nucs = 'N';
N15_PMTa.n=[1];
N15_PMTa.A = [mt2mhz(1.5375,N15_PMTa.g)];
N15_PMTa.lwpp = [.06185];

N15_PMTaVary.g = [0];
N15_PMTaVary.A = [0];
N15_PMTaVary.lwpp = 0;

N15_PMTaExp.mwFreq = 9.7070;
N15_PMTaExp.CenterSweep = [345.2 10.0];

[BN15_PMTa,SN15_PMTa]=garlic(N15_PMTa,N15_PMTaExp);
%0.9

%Pure spin trap and benzophenone irradiated-N-H interaction peak
NH_PMTa.g = [2.00852];
NH_PMTa.Nucs = 'N,H';
NH_PMTa.n=[1 1];
NH_PMTa.A = [mt2mhz(1.3542,NH_PMTa.g);mt2mhz(.1860,NH_PMTa.g)];
NH_PMTa.lwpp = .075;

NH_PMTaVary.g = [0];
NH_PMTaVary.A = [mt2mhz(.01,NH_PMTa.g);mt2mhz(.002,NH_PMTa.g)];
NH_PMTaVary.lwpp = 0.01;

NH_PMTaExp.mwFreq = 9.7070;
NH_PMTaExp.CenterSweep = [345.2 10.0];

[BNH_PMTa,SNH_PMTa]=garlic(NH_PMTa,NH_PMTaExp);
%8

%Pure spin trap and benzophenone irradiated-7.5 gauss ketone peak
N7_PMTa.g = [2.00908];
N7_PMTa.Nucs = 'N';
N7_PMTa.n=[1];
N7_PMTa.A = [mt2mhz(.7969,N7_PMTa.g)];
N7_PMTa.lwpp = .06007;

N7_PMTaVary.g = [0];
N7_PMTaVary.A = [mt2mhz(.05,N7_PMTa.g)];
N7_PMTaVary.lwpp = 0.02;

N7_PMTaExp.mwFreq = 9.7070;
N7_PMTaExp.CenterSweep = [345.2 10.0];

[BN7_PMTa,SN7_PMTa]=garlic(N7_PMTa,N7_PMTaExp);
%2.6


%Pure spin trap and benzophenone irradiated-mystery peak
N14_2_PMTa.g = [2.00829];
N14_2_PMTa.Nucs = 'N';
N14_2_PMTa.n=[1];
N14_2_PMTa.A = [mt2mhz(1.3691,N14_2_PMTa.g)];
N14_2_PMTa.lwpp = .35;

N14_2_PMTaVary.g = [.001];
N14_2_PMTaVary.A = [mt2mhz(.1,N14_2_PMTa.g)];
N14_2_PMTaVary.lwpp = .2;

N14_2_PMTaExp.mwFreq = 9.7070;
N14_2_PMTaExp.CenterSweep = [345.2 10.0];

[BN14_2_PMTa,SN14_2_PMTa]=garlic(N14_2_PMTa,N14_2_PMTaExp);
%3.5

start_point = rand(1, 4);
estimates3(1,:) = patternsearch(@expfun3, start_point,[],[],[],[],[0 0 0 0],[1 1 1 1],[],options);
fitted3(:,1)=estimates3(1,1).*SN15_PMTa+estimates3(1,2).*SNH_PMTa+estimates3(1,3).*SN7_PMTa+estimates3(1,4).*SN14_2_PMTa;

figure
subplot(2,1,1)
plot(BN14_2_PMTa,fitted3,'r')
hold
plot(BN14_2_PMTa,APMTa(:,10)./10000)
ylabel('Intensity')
subplot(2,1,2)
plot(BN14_2_PMTa,APMTa(:,10)./10000-fitted3)
xlabel('B (mT)')
title('Subtraction of Fit from Data')

% close all
% plot(GPMTa(:,10)./10,APMTa(:,10)./10000)
% hold
% plot(BN14_2_PMTa,SN14_2_PMTa./3.5,'r')
% hold



% FitOpt.Method = 'genetic int';
% [bestsys,bestspc]=esfit('garlic',ATa(:,8),N15_Ta,N15_TaVary,N15_TaExp,[],FitOpt);
% bestspc=bestspc';
% mhz2mt(bestsys.A,bestsys.g)

    function sse = expfun1(params)
        A = params(1);
        FittedCurve=A.*SN15_Ta';
        ErrorVector = FittedCurve - ATa(:,8)./10000;
        sse = sum(ErrorVector .^ 2);
    end

    function sse = expfun2(params)
        A = params(1);
        B = params(2);
        C = params(3);
        D = params(4);
        FittedCurve=A.*SN15_MTa'+B.*SNH_MTa'+C.*SN7_MTa'+D.*SN14_2_MTa';
        ErrorVector = FittedCurve - AMTa(:,8)./10000;
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