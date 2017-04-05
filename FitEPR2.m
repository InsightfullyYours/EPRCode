function [fitted estimates]=FitEPR(xData,yData,xSim1,Sim1,xSim2,Sim2,xSim3,Sim3)

if size(xData,2)==1
    xData=repmat(xData,1,size(yData,2));
end

for i=1:size(yData,2)
    y=yData(:,i);
    x=xData(:,i);
    xS1=xSim1(:,i);
    S1=Sim1(:,i);
    
    options=optimset('Display','iter','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-8,'TolFun',1E-8);
    
    if nargin==4
        start_point = max(x).*rand(1, 1);
        estimates(i,:) = fminsearch(@expfun, start_point,options);
        fitted(:,i)=estimates(i,1).*S1(:,i);
    elseif nargin==6
        xS2=xSim2(:,i);
        S2=Sim2(:,i);
        start_point = max(x).*rand(1, 2);
        estimates(i,:) = fminsearch(@expfun2, start_point,options);
        fitted(:,i)=estimates(i,1).*S1(:,i)+estimates(i,2).*S2(:,i);
    elseif nargin==8
        xS2=xSim2(:,i);
        S2=Sim2(:,i);
        xS3=xSim3(:,i);
        S3=Sim3(:,i);
        start_point = max(x).*rand(1, 3);
        estimates(i,:) = fminsearch(@expfun3, start_point,options);
        fitted(:,i)=estimates(i,1).*S1(:,i)+estimates(i,2).*S2(:,i)+estimates(i,3).*S3(:,i);
    end
    
    subplot(2,1,1)
    plot(x,y(:,i),'r')
    V=axis
    ylim
    hold
    plot(xS1,fitted(:,i),'b')
    ylabel('Intensity (AU)')
    legend(['Data';'Fit ']);
    axis(V)
    ylim
    
    subplot(2,1,2)
    plot(x,y(:,i)-fitted(:,i),'b')
    ylim(V(3:4));
    axis(V)
    xlabel('Hyperfine Splitting (Gauss)')
    ylabel('Intensity (AU)')
    title('Subtracted')
    subplot(2,1,1)
end


    function sse = expfun(params)
        m = params(1);
        b = params(2);
        FittedCurve=m.*S1;
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end

    function sse = expfun2(params)
        m = params(1);
        b = params(2);
        FittedCurve=m.*S1+b.*S2;
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end

    function sse = expfun3(params)
        m = params(1);
        b = params(2);
        k = params(3);
        FittedCurve=m.*S1+b.*S2+k.*S3;
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end

end