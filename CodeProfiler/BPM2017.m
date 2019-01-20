function [FieldMatrix, AverageIndexMatrix,X,Z]=BPM2017Test
%2D Beam Propagation Method witn transparent Boundary Conditions
%Optimized version
%% Input parameters
lambda=1.55;                            %operating wavelength [um]
W=10;                                   %box width [um]
L=150;                                  %box length [um]
DX=0.01;                                %discretization step in X direction	[um]
DZ=0.05;                                %discretization step in Z direction [um]
nref=3.2;                               %reference refractive index []
SubSample=int16(100);                   %Subsample data in the Z direction for final plots [#]

%% Nodes in X dir
Nx=round(W/DX);                         %how many nodes with the supplied DX? [#]
if rem(Nx,2)==0                         %is N even? [#]
    Nx=Nx+1;                            %impose N odd, to have a node in x=0 [#]
end
X=linspace(-W/2,W/2,Nx).';              %column vector of X nodes [um]
DX=X(2)-X(1);                           %X dir. step [um]

%% Nodes in Z dir
Nz=round(L/DZ);                         %how many nodes with the supplied DZ? [#]
Z=linspace(0,L,Nz);                     %row vector of Z nodes [um]
DZ=Z(2)-Z(1);                           %Z dir. step [um]
%% To ensure reproducibility, reinitialize random number generator
rng(0);
%% Initial field in z=0
kprop=2*pi/lambda;                      %free space propagation constant [um^-1]
field=InitialFieldProfile(X);           %initial distribution of the field [A.U.]
FieldMatrix(:,1)=field;                 %save the input field
%% Main loop in z direction
for m=1:Nz-1
    %% Refractive index and losses in fractionary node m+1/2
    AverageIndex=CalculateAverageIndex(X,Z(m)+DZ/2, W, L, nref);      
    %% Matrices construction   
    o=ones(Nx,1);
    oo=ones(Nx-1,1);
    psi=4j*nref*kprop*DX^2/DZ*o;
    Aconst=diag(psi+2,0)-diag(oo,+1)-diag(oo,-1);
    Bconst=diag(psi-2,0)+diag(oo,+1)+diag(oo,-1);

    chi=-kprop^2*DX^2*(AverageIndex.^2-nref^2);
    ABvar=diag(chi);
    A=Aconst+ABvar;
    B=Bconst-ABvar;
    %% Transparent Boundary Condition in x=0
    Exp=1;
    if field(2)~=0
        Exp=field(1)/field(2);
        if imag(log(Exp))>0             
            Exp=conj(Exp);
        end
    end
    %Enforce the updated relation between components 1 and 2
    field(1)=field(2)*Exp;
    %Correct elements (1,1) in A and B
    A(1,1)=A(1,1)-Exp;
    B(1,1)=B(1,1)+Exp;
    
    %% Transparent Boundary Condition in x=L
    Exp=1;
    if field(end-1)~=0
        Exp=field(end)/field(end-1);
        if imag(log(Exp))>0
            Exp=conj(Exp);              
        end
    end
    %Enforce the updated relation between elements end and end-1
    field(end)=field(end-1)*Exp;
    %Correct elements (end,end) in A and B
    A(end,end)=A(end,end)-Exp;
    B(end,end)=B(end,end)+Exp;
    
    %% Calculate the new field profile solving the linear system
    field=A\(B*field);
    %% Save data for final plots
    FieldMatrix(:,m+1)=field+rand(Nx,1)*0.02*max(abs(field));
    AverageIndexMatrix(:,m)=AverageIndex.';
    

        Power=abs(field).^2;
        Status=sprintf('Iteration: %g (%g%%), position: %g um, norm. power: %g',m,m/Nz*100,m*DZ,sum(Power));
        disp(Status);
        %Current field distribution, in ||^2
        subplot(2,1,1);
        plot(X,Power);
        title(Status);
        xlabel('x [\mum]');
        ylabel('|E|^2 [A.U.]');
        %Refractive index
        subplot(2,1,2);
        plot(X,AverageIndex);
        xlabel('x [\mum]');
        ylabel('n [ ]');
        %Force figure redraw
        drawnow;
end
%% final plots
figure;
subplot(2,1,1);
Z=Z(1:SubSample:end);                 
surf(Z,X,abs(FieldMatrix(:,1:SubSample:end)).^2);
xlabel('Z [um]');
ylabel('X [um]');
zlabel('|E|^2');
shading interp;
subplot(2,1,2);
surf(Z,X,AverageIndexMatrix(:,1:SubSample:end));
xlabel('Z [um]');
ylabel('X [um]');
zlabel('Index [\mum^{-1}]');
shading interp;

%% save data to file
hFile=fopen('results.txt','w');
fwrite(hFile,abs(FieldMatrix).^2);
fclose(hFile);

if nargout==0
    clearvars('FieldMatrix', 'AverageIndexMatrix','X','Z');
end
end

function AverageIndex = CalculateAverageIndex(X,z, W, L, nref)    
%Vector of refractive index profile in slice z
%Input: 
%   X:      coordinates of the nodes in the x direction [um]
%   z:      z direction position [um]
%   W:      box width (x dir) [um]
%   L:      box length (z dir) [um]
%Output:
%  AverageIndex: refractive index in each x point at the specified z position []
Wg=1;                                   %Waveguide width [um]
L1=50;                                  %Initial linear part length [um]
L2=100;                                 %start point for the final linear part length [um]
x2=2;                                   %Wavegudie center in the final linear part [um]
dn=0.05;                                %Refractive index in the core is nref+dn

if z<L1                                 %first part: linear, center in x=0   
    AverageIndex=nref+dn*(abs(X)<Wg/2);    
elseif z>L2                             %last part: linear, center in x=x2
    AverageIndex=nref+dn*(abs(X-x2)<Wg/2);
else
    xc=(z-L1)/(L2-L1)*x2;               %center position [um]
    AverageIndex=nref+dn*(abs(X-xc)<Wg/2);    
end  
end

function E=InitialFieldProfile(X)
%Initial (z=0) distribution of the field
%Input: 
%   X:      coordinates of the nodes in the x direction [um]
%Output:
%   E:      vector with electric field distribution [A.U.]

E=exp(-X.^2); %for simplicity: Gaussian profile [A.U.]

end