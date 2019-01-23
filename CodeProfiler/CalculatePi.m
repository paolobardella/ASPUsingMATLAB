function Pi=CalculatePi
NumSamples=10e6;
R=0.5; Inside=0;
%figure;hold on;
x=rand(NumSamples,1)-R;
y=rand(NumSamples,1)-R;
for k=1:NumSamples
    if sqrt(x(k)*x(k)+y(k)*y(k))<=R
        Inside=Inside+1;
        %plot(x,y,'ko');
    else
        %plot(x,y,'b.');
    end
end
Pi=4*Inside/NumSamples;