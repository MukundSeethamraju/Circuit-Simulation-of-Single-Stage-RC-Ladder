%Submitted by
%Suhas M,            10033          suhas.msh@gmail.com          
%Mukund Seethamraju, 09969          mukund.seethamraju@gmail.com 


function [t, y] = MonteCarlo(R,C,RL)

	tspan = [0:1e-12:1e-8];
  
	y0 = [1e-9;0;0;0];
	
	M = [0 0 0 0; 0 0 0 0; 0 0 C 0; 0 0 0 C;];
	M = -M; 
	options = odeset('Mass',M, 'MassSingular','yes','RelTol',1e-8,'AbsTol',1e-8 );
	
	[t, y] = ode15s(@rcDAE, tspan, y0, options);

    
	function out = rcDAE(t,y) 
		out = [ -y(1) + (y(2) - y(3))/R;
		 	y(2) - sin(1e9*t);
			 -(y(2) - y(3))/R + (y(3) - y(4))/R;
			 -(y(3) - y(4))/R + y(4)/RL]; 
    end


end

	

