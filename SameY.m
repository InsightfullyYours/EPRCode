function [y1 y2]=SameY(x1,x2)

test=x1(1);
x2;

if (x1(1)~=0)
    y1=x1-x1(1);
end
if x2(1)~=0
    y2=x2-x2(1);
end
end