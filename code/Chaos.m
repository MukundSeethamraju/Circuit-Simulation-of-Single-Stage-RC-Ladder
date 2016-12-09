%Submitted by
%Suhas M,            10033          suhas.msh@gmail.com          
%Mukund Seethamraju, 09969          mukund.seethamraju@gmail.com 


function [t, y] = Chaos(rmu,cmu,rsigma,csigma,RL)

	tspan = [0:1e-12:1e-8]; %We want to evaluate a period of 10ns, in steps of 1ps.
y0 = [1e-9;1e-9;1e-9;1e-9;1e-9;1e-9;1e-9;1e-9;1e-9;1e-12;1e-9;1e-9;1e-9;1e-9;1e-9;1e-9]; %column vector of initial conditions. y= [I V1 V2 V3 V4]
 
    
    
    

	

	M = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  %y01 y(1)
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %y11 y(2) 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %y21 y(3) 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %y31 y(4) 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %y02 y(5) 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %y12 y(6) 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %y22 y(7) 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %y32 y(8) 
%12 3 4 5 6 7 8 9 0 1 2 3 4 5 6
0 0 0 0 0 0 0 0 cmu*rmu+csigma*rsigma cmu*rsigma+csigma*rmu 2*csigma*rsigma 0 0 0 0 0; %y03 y(9) 
0 0 0 0 0 0 0 0 cmu*rsigma+csigma*rmu cmu*rmu+3*csigma*rsigma 2*cmu*rsigma+2*csigma*rmu 6*csigma*rsigma 0 0 0 0; %y13 y(10) 
0 0 0 0 0 0 0 0 csigma*rsigma cmu*rsigma+csigma*rmu cmu*rmu+5*csigma*rsigma 3*cmu*rsigma+3*csigma*rmu 0 0 0 0; %y23 y(11) 
0 0 0 0 0 0 0 0 0 csigma*rsigma cmu*rsigma+csigma*rmu cmu*rmu+7*csigma*rsigma 0 0 0 0; %y33 y(12) 
0 0 0 0 0 0 0 0 0 0 0 0 cmu*rmu+csigma*rsigma cmu*rsigma+csigma*rmu 2*csigma*rsigma 0; %y04 y(13) 
0 0 0 0 0 0 0 0 0 0 0 0 cmu*rsigma+csigma*rmu cmu*rmu+3*csigma*rsigma 2*cmu*rsigma+2*csigma*rmu 6*csigma*rsigma; %y14 y(14) 
0 0 0 0 0 0 0 0 0 0 0 0 csigma*rsigma cmu*rsigma+csigma*rmu cmu*rmu+5*csigma*rsigma 3*cmu*rsigma+3*csigma*rmu; %y24 y(15) 
0 0 0 0 0 0 0 0 0 0 0 0 0 csigma*rsigma csigma*rmu+cmu*rsigma cmu*rmu+7*csigma*rsigma];%y34 y(16) 


	options = odeset('Mass',M, 'MassSingular','yes','RelTol',1e-8,'AbsTol',1e-8 );



 
	[t, y] = ode15s(@rcDAE, tspan, y0, options); 
	

	function out = rcDAE(t,y)
		out = [ -y(9)+y(5)-(y(1)*rmu)-(y(2)*rsigma);
-y(10)+y(6)-(y(2)*rmu)-(y(1)*rsigma)-(2*rsigma*y(3));
-y(11)+y(7)-(y(3)*rmu)-(y(2)*rsigma)-(3*rsigma*y(4));
-y(12)+y(8)-(y(4)*rmu)-(y(3)*rsigma);
-sin((1e9)*t)+y(5);
-y(6);
-y(7);
-y(8);
-(2*y(9))+y(5)+y(13);
-(2*y(10))+y(6)+y(14);
-(2*y(11))+y(7)+y(15);
-(2*y(12))+y(8)+y(16);
-(y(13))+(y(9))-(y(13)*rmu/RL)-(y(14)*rsigma/RL);
-(y(14))+(y(10))-(y(14)*rmu/RL)-(y(13)*rsigma/RL)-(y(15)*2*rsigma/RL);
-(y(15))+(y(11))-(y(15)*rmu/RL)-(y(14)*rsigma/RL)-(y(16)*3*rsigma/RL);
-(y(16))+(y(12))-(y(16)*rmu/RL)-(y(15)*rsigma/RL)];
    end

end

	

