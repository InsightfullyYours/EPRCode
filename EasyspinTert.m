%function EasyspinTert(ATa,AMTa,APMTa)
%Sys.s = [1;1/2];
clear FitOpt UserCommand h
clear N15_Tt N15_TtVary N15_TtExp
clear N15_MTt N15_MTtVary N15_MTtExp
clear NH_MTt NH_MTtVary NH_MTtExp
clear N7_MTt N7_MTtVary N7_MTtExp
clear N14_2_MTt N14_2_MTtVary N14_2_MTtExp
clear N15_Tt N15_TtVary N15_TtExp
clear N15_PMTt N15_PMTtVary N15_PMTtExp
clear NH_PMTt NH_PMTtVary NH_PMTtExp
clear N7_PMTt N7_PMTtVary N7_PMTtExp
clear N14_2_PMTt N14_2_PMTtVary N14_2_PMTtExp
clear HH_Tt HH_TtVary HH_TtExp

close all

options=optimset('Display','iter','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-8,'TolFun',1E-8);

%This is the pure spin trap in toluene-Pure N peak
N15_Tt.g = [2.01142];
N15_Tt.Nucs = 'N';
N15_Tt.n=[1];
N15_Tt.A = [mt2mhz(1.38,N15_Tt.g)];
N15_Tt.lwpp = [.25];

N15_TtVary.g = [0.001];
N15_TtVary.A = [mt2mhz(.05,N15_Tt.g)];
N15_TtVary.lwpp = .05;

N15_TtExp.mwFreq = 9.7710;
N15_TtExp.CenterSweep = [347.0 10.0];

[BN15_Tt,SN15_Tt]=garlic(N15_Tt,N15_TtExp);
%0.3

% start_point = rand(1, 1);
% estimates1(1,:) = patternsearch(@expfun1, start_point,[],[],[],[],0,1,[],options);
% fitted1(:,1)=estimates1(1,1).*SN15_Tt;

%Pure spin trap and Toluene irradiated-N-H interaction peak
NH_Tt.g = [2.01135];
NH_Tt.Nucs = 'N,H';
NH_Tt.n=[1 1];
NH_Tt.A = [mt2mhz(1.38,NH_Tt.g);mt2mhz(3.634,NH_Tt.g)];
NH_Tt.lwpp = .12;

NH_TtVary.g = [.001];
NH_TtVary.A = [mt2mhz(.5,NH_Tt.g);mt2mhz(.5,NH_Tt.g)];
NH_TtVary.lwpp = 0.01;

NH_TtExp.mwFreq = 9.7710;
NH_TtExp.CenterSweep = [347.0 10.0];

[BNH_Tt,SNH_Tt]=garlic(NH_Tt,NH_TtExp);
%1


%Pure spin trap and Toluene irradiated-H-H interaction peak
HH_Tt.g = [2.009326];
HH_Tt.Nucs = 'N,H,H';
HH_Tt.n=[1 1 1];
HH_Tt.A = [mt2mhz(1.0,HH_Tt.g);mt2mhz(0.19,HH_Tt.g);mt2mhz(0.19,HH_Tt.g)];
HH_Tt.lwpp = .1;

HH_TtVary.g = [.001];
HH_TtVary.A = [mt2mhz(.05,HH_Tt.g);mt2mhz(.05,HH_Tt.g);mt2mhz(.05,HH_Tt.g)];
HH_TtVary.lwpp = 0.04;

HH_TtExp.mwFreq = 9.7710;
HH_TtExp.CenterSweep = [347.0 10.0];

[BHH_Tt,SHH_Tt]=garlic(HH_Tt,HH_TtExp);
%6

% figure
% subplot(2,1,1)
% plot(BHH_Tt,SN15_Tt'./.3+SNH_Tt'./1+SHH_Tt'./6,'r')
% hold
% plot(BHH_Tt,ATt(:,10)./10000)
% axis([343 351 -.7 .7])
% ylabel('Intensity')
% subplot(2,1,2)
% plot(BHH_Tt,ATt(:,10)./10000-(SN15_Tt'./.3+SNH_Tt'./1+SHH_Tt'./6))
% axis([343 351 -.7 .7])
% ylabel('Intensity')
% xlabel('B (mT)')
% title('Subtraction of Fit from Data')

% %--------------------------------------------------------------------------
% 
%This is the pure spin trap in toluene-Pure N peak
N15_PMTt.g = [2.00882];
N15_PMTt.Nucs = 'N';
N15_PMTt.n=[1];
N15_PMTt.A = [mt2mhz(1.38,N15_Tt.g)];
N15_PMTt.lwpp = [.23];

N15_PMTtVary.g = [0.001];
N15_PMTtVary.A = [mt2mhz(.05,N15_Tt.g)];
N15_PMTtVary.lwpp = .05;

N15_PMTtExp.mwFreq = 9.7560;
N15_PMTtExp.CenterSweep = [347.0 10.0];

[BN15_PMTt,SN15_PMTt]=garlic(N15_PMTt,N15_PMTtExp);
%0.3

% start_point = rand(1, 1);
% estimates1(1,:) = patternsearch(@expfun1, start_point,[],[],[],[],0,1,[],options);
% fitted1(:,1)=estimates1(1,1).*SN15_PMTt;

%Pure spin trap and Toluene irradiated-N-H interaction peak
NH_PMTt.g = [2.00872];
NH_PMTt.Nucs = 'N,H';
NH_PMTt.n=[1 1];
NH_PMTt.A = [mt2mhz(1.38,NH_PMTt.g);mt2mhz(3.634,NH_PMTt.g)];
NH_PMTt.lwpp = .12;

NH_PMTtVary.g = [.001];
NH_PMTtVary.A = [mt2mhz(.5,NH_PMTt.g);mt2mhz(.5,NH_PMTt.g)];
NH_PMTtVary.lwpp = 0.01;

NH_PMTtExp.mwFreq = 9.7560;
NH_PMTtExp.CenterSweep = [347.0 10.0];

[BNH_PMTt,SNH_PMTt]=garlic(NH_PMTt,NH_PMTtExp);
%1


%Pure spin trap and Toluene irradiated-H-H interaction peak
HH_PMTt.g = [2.00677];
HH_PMTt.Nucs = 'N,H,H';
HH_PMTt.n=[1 1 1];
HH_PMTt.A = [mt2mhz(1.0,HH_PMTt.g);mt2mhz(0.19,HH_PMTt.g);mt2mhz(0.19,HH_PMTt.g)];
HH_PMTt.lwpp = .1;

HH_PMTtVary.g = [.001];
HH_PMTtVary.A = [mt2mhz(.05,HH_PMTt.g);mt2mhz(.05,HH_PMTt.g);mt2mhz(.05,HH_PMTt.g)];
HH_PMTtVary.lwpp = 0.04;

HH_PMTtExp.mwFreq = 9.7560;
HH_PMTtExp.CenterSweep = [347.0 10.0];

[BHH_PMTt,SHH_PMTt]=garlic(HH_PMTt,HH_PMTtExp);
%6

%Pure spin trap and Toluene irradiated-N-H interaction peak
NH_2_PMTt.g = [2.00872];
NH_2_PMTt.Nucs = 'N,H';
NH_2_PMTt.n=[1 1];
NH_2_PMTt.A = [mt2mhz(1.35,NH_2_PMTt.g);mt2mhz(2.974,NH_2_PMTt.g)];
NH_2_PMTt.lwpp = .12;

NH_2_PMTtVary.g = [.001];
NH_2_PMTtVary.A = [mt2mhz(.5,NH_2_PMTt.g);mt2mhz(.5,NH_2_PMTt.g)];
NH_2_PMTtVary.lwpp = 0.01;

NH_2_PMTtExp.mwFreq = 9.7560;
NH_2_PMTtExp.CenterSweep = [347.0 10.0];

[BNH_2_PMTt,SNH_2_PMTt]=garlic(NH_2_PMTt,NH_2_PMTtExp);
%.7

figure
subplot(2,1,1)
plot(BHH_PMTt,SNH_2_PMTt'./1.5+SHH_PMTt'./2.2+SNH_PMTt'./1+SN15_PMTt'./.15,'r')
hold
plot(BHH_PMTt,APMTt(:,11)./10000)
axis([343 351 -2 2])
ylabel('Intensity')
subplot(2,1,2)
plot(BHH_Tt,APMTt(:,11)./10000-(SNH_2_PMTt'./1.5+SHH_PMTt'./2.2+SNH_PMTt'./1+SN15_PMTt'./.15),'b')
axis([343 351 -2 2])
ylabel('Intensity')
xlabel('B (mT)')
title('Subtraction of Fit from Data')
% %-------------------------------------------------------------------------
% 
% % pure spin trap in toluene and benzophenone and Polystyrene-pure N peak
% %(accounts for shifts)
% N15_PMTt.g = [2.00905];
% N15_PMTt.Nucs = 'N';
% N15_PMTt.n=[1];
% N15_PMTt.A = [mt2mhz(1.5375,N15_PMTt.g)];
% N15_PMTt.lwpp = [.06185];
% 
% N15_PMTtVary.g = [0];
% N15_PMTtVary.A = [0];
% N15_PMTtVary.lwpp = 0;
% 
% N15_PMTtExp.mwFreq = 9.7070;
% N15_PMTtExp.CenterSweep = [345.2 10.0];
% 
% [BN15_PMTt,SN15_PMTt]=garlic(N15_PMTt,N15_PMTtExp);
% %0.9
% 
% %Pure spin trap and benzophenone irradiated-N-H interaction peak
% NH_PMTt.g = [2.00852];
% NH_PMTt.Nucs = 'N,H';
% NH_PMTt.n=[1 1];
% NH_PMTt.A = [mt2mhz(1.3542,NH_PMTt.g);mt2mhz(.1860,NH_PMTt.g)];
% NH_PMTt.lwpp = .075;
% 
% NH_PMTtVary.g = [0];
% NH_PMTtVary.A = [mt2mhz(.01,NH_PMTt.g);mt2mhz(.002,NH_PMTt.g)];
% NH_PMTtVary.lwpp = 0.01;
% 
% NH_PMTtExp.mwFreq = 9.7070;
% NH_PMTtExp.CenterSweep = [345.2 10.0];
% 
% [BNH_PMTt,SNH_PMTt]=garlic(NH_PMTt,NH_PMTtExp);
% %8
% 
% %Pure spin trap and benzophenone irradiated-7.5 gauss ketone peak
% N7_PMTt.g = [2.00908];
% N7_PMTt.Nucs = 'N';
% N7_PMTt.n=[1];
% N7_PMTt.A = [mt2mhz(.7969,N7_PMTt.g)];
% N7_PMTt.lwpp = .06007;
% 
% N7_PMTtVary.g = [0];
% N7_PMTtVary.A = [mt2mhz(.05,N7_PMTt.g)];
% N7_PMTtVary.lwpp = 0.02;
% 
% N7_PMTtExp.mwFreq = 9.7070;
% N7_PMTtExp.CenterSweep = [345.2 10.0];
% 
% [BN7_PMTt,SN7_PMTt]=garlic(N7_PMTt,N7_PMTtExp);
% %2.6
% 
% 
% %Pure spin trap and benzophenone irradiated-mystery peak
% N14_2_PMTt.g = [2.00829];
% N14_2_PMTt.Nucs = 'N';
% N14_2_PMTt.n=[1];
% N14_2_PMTt.A = [mt2mhz(1.3691,N14_2_PMTt.g)];
% N14_2_PMTt.lwpp = .35;
% 
% N14_2_PMTtVary.g = [.001];
% N14_2_PMTtVary.A = [mt2mhz(.1,N14_2_PMTt.g)];
% N14_2_PMTtVary.lwpp = .2;
% 
% N14_2_PMTtExp.mwFreq = 9.7070;
% N14_2_PMTtExp.CenterSweep = [345.2 10.0];
% 
% [BN14_2_PMTt,SN14_2_PMTt]=garlic(N14_2_PMTt,N14_2_PMTtExp);
% %3.5
% 
% start_point = rand(1, 4);
% estimates3(1,:) = patternsearch(@expfun3, start_point,[],[],[],[],[0.02 0 0 0],[1 1 1 1],[],options);
% fitted3(:,1)=estimates3(1,1).*SN15_PMTt+estimates3(1,2).*SNH_PMTt+estimates3(1,3).*SN7_PMTt+estimates3(1,4).*SN14_2_PMTt;
% 
% figure
% subplot(2,1,1)
% plot(BN14_2_PMTt,fitted3,'r')
% hold
% plot(BN14_2_PMTt,APMTa(:,10)./10000)
% ylabel('Intensity')
% subplot(2,1,2)
% plot(BN14_2_PMTt,APMTa(:,10)./10000-fitted3)
% xlabel('B (mT)')
% title('Subtraction of Fit from Data')
% 
% % close all
% % plot(GPMTa(:,10)./10,APMTa(:,10)./10000)
% % hold
% % plot(BN14_2_PMTt,SN14_2_PMTt./3.5,'r')
% % hold



% FitOpt.Method = 'genetic int';
% [bestsys,bestspc]=esfit('garlic',ATt(:,10)./10000-SNH_Tt'./1-SN15_Tt'./0.3,HH_Tt,HH_TtVary,HH_TtExp,[],FitOpt);
% bestspc=bestspc';
% mhz2mt(bestsys.A,bestsys.g)

%     function sse = expfun1(params)
%         A = params(1);
%         FittedCurve=A.*SN15_Tt';
%         ErrorVector = FittedCurve - ATa(:,8)./10000;
%         sse = sum(ErrorVector .^ 2);
%     end
% 
%     function sse = expfun2(params)
%         A = params(1);
%         B = params(2);
%         C = params(3);
%         D = params(4);
%         FittedCurve=A.*SN15_MTt'+B.*SNH_MTt'+C.*SN7_MTt'+D.*SN14_2_MTt';
%         ErrorVector = FittedCurve - AMTa(:,8)./10000;
%         sse = sum(ErrorVector .^ 2);
%     end
% 
%     function sse = expfun3(params)
%         A = params(1);
%         B = params(2);
%         C = params(3);
%         D = params(4);
%         FittedCurve=A.*SN15_Tt'+B.*SNH_PMTt'+C.*SN7_PMTt'+D.*SN14_2_PMTt';
%         ErrorVector = FittedCurve - APMTa(:,10)./10000;
%         sse = sum(ErrorVector .^ 2);
%     end



clear FitOpt UserCommand h
clear N15_Tt N15_TtVary N15_TtExp
clear N15_MTt N15_MTtVary N15_MTtExp
clear NH_MTt NH_MTtVary NH_MTtExp
clear N7_MTt N7_MTtVary N7_MTtExp
clear N14_2_MTt N14_2_MTtVary N14_2_MTtExp
clear N15_Tt N15_TtVary N15_TtExp
clear N15_PMTt N15_PMTtVary N15_PMTtExp
clear NH_PMTt NH_PMTtVary NH_PMTtExp
clear N7_PMTt N7_PMTtVary N7_PMTtExp
clear N14_2_PMTt N14_2_PMTtVary N14_2_PMTtExp
clear HH_Tt HH_TtVary HH_TtExp
%end