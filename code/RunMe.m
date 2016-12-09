%Submitted by
%Suhas M,            10033          suhas.msh@gmail.com          
%Mukund Seethamraju, 09969          mukund.seethamraju@gmail.com 


clc
%Set number of iterations.
c=100; 

% Set circuit Parameters
RL=50;
rmu= 3.3504;
cmu=0.6271e-12;
rsigma= 0.0861;
csigma=0.2266e-12;


%%%%%%%%%%%%%%%%%% Chaos Method %%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The program shall first simulate the circuit by Wiener Chaos.') 
fprintf('\n')
disp('Chaos Function is being called to compute the Wiener coefficients.')
tic;
[t,y]=Chaos(rmu,cmu,rsigma,csigma,RL);
toc;
disp('Chaos Function has returned the Wiener Coefficients.')
fprintf('The program shall now compute %d trajectories and then calculate their mean. \n',c)
tic;


for k=1:c 

X= normrnd(0,1);

%Hermite Polynomials
h0=1;
h1=X;
h2=X^2-1;
h3=X^3-3*X;
    

y1 = y(:,1)*h0 + y(:,2)*h1+ y(:,3)*h2 +y(:,4)*h3;
y2 = y(:,5)*h0 + y(:,6)*h1 + y(:,7)*h2 + y(:,8)*h3;
y3 = y(:,9)*h0+ y(:,10)*h1+ y(:,11)*h2 + y(:,12)*h3;
y4 = y(:,13)*h0+ y(:,14)*h1+ y(:,15)*h2 +y(:,16)*h3;

Y=[y1 y2 y3 y4];

if k==1
    A=zeros(10001,4);
end

A=A+Y;

end
A=A/k;
toc


fprintf('Wiener Chaos Simulation has ended. Mean Trajectory has been computed. \n')
fprintf('\n')
%Plot the Graph
figure();
plot(t,A)
title('Wiener Chaos Simulation')
xlabel('Time in s');
legend('I','V_1','V_2','V_3')

%%%%%%%%%%%%%%%%%%%% End of Chaos Simulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Start Monte Carlo Simulation%%%%%%%%%%%%%%%%%%%%%
fprintf('The program shall now simulate by %d iterations of Monte Carlo method. \n',c)

%Make a storage Matrix
Y = zeros(10001,4,c);
tic;
for i=1:c

fprintf('%d \t',i)
if mod(i,10) ==0 
    fprintf('\n')
end 

Crandom = normrnd(cmu,csigma);
while Crandom>(1e-280)==0
  Crandom = normrnd(cmu,csigma);  
end

[T,Y(:,:,i)]=MonteCarlo(normrnd(rmu,rsigma),Crandom,RL);
    
end
fprintf('\n')
toc;

fprintf('\nMonte Carlo Simulation has ended.\n')

B=sum(Y,3)/c;

%Plot the Graph
figure
plot(T,B)
title('Monte Carlo Simulation')
xlabel('Time in s');
legend('I','V_1','V_2','V_3')

%Error
DiffAB = A-B;
DiffABSQR = DiffAB.^2;
MeanSQRDev= (sum(DiffABSQR))/10001;
fprintf('\nThe RMS Deviation between Monte Carlo and Chaos Expansion for [I V1 V2 V3] is: \n')
format shortEng;
((MeanSQRDev).^(1/2))
