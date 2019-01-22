%Goal: maximize the number of considered points per second
tic
NumSamples=10000;
R=0.5; 
Inside=0;
figure; 
hold on;
for k=1:NumSamples
    x=rand()-R;
    y=rand()-R;
    if sqrt(x*x+y*y)<=R
        Inside=Inside+1;
        plot(x,y,'ko');
    else
        plot(x,y,'b.');
    end
    
end
Pi=4*Inside/NumSamples;
fprintf('You analyzed %g points per second\n',NumSamples/toc) 

