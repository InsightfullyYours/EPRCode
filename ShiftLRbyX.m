function SpectrumOut=ShiftLRbyX(Data,X)

%this just takes data and shifts it by X.  It's intended as a subsunction
%of EPR data analysis programs. Negative is left, positive right

X=round(X); %because it always has to be an integer

if X==0
    SpectrumOut=Data;
elseif X<0
    SpectrumOut=[Data(abs(X)+1:end,:); zeros(abs(X),size(Data,2))];
elseif X>0
    SpectrumOut=[zeros(X,size(Data,2)); Data(1:end-X,:)];
end
if size(Data,1)~=size(SpectrumOut,1)
    keyboard
end
end
    