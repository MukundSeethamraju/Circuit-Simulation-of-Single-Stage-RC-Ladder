%Submitted by
%Suhas M,            10033          suhas.msh@gmail.com          
%Mukund Seethamraju, 09969          mukund.seethamraju@gmail.com 


clc
syms y01 y11 y21 y31 y02 y12 y22 y32 y03 y13 y23 y33 y04 y14 y24 y34 y41 y42 y43 y44
syms Zy01 Zy11 Zy21 Zy31 Zy02 Zy12 Zy22 Zy32 Zy03 Zy13 Zy23 Zy33 Zy04 Zy14 Zy24 Zy34 Zy41 Zy42 Zy43 Zy44

syms h0 h1 h2 h3 h4 rmu rsigma C cmu csigma RL

syms X

%hermite polynomials
h0=1;
h1=X;
h2=X^2-1;
h3=X^3-3*X;
w=exp(-(X^2)/2);

%Chaos expansion
y1 = y01*h0+ y11*h1+ y21*h2+ y31*h3;
y2 = y02*h0+ y12*h1+ y22*h2+ y32*h3;
y3 = y03*h0+ y13*h1+ y23*h2+ y33*h3;
y4 = y04*h0+ y14*h1+ y24*h2+ y34*h3;

%Derivatives
Zy1 = Zy01*h0+ Zy11*h1+ Zy21*h2+ Zy31*h3;
Zy2 = Zy02*h0+ Zy12*h1+ Zy22*h2+ Zy32*h3;
Zy3 = Zy03*h0+ Zy13*h1+ Zy23*h2+ Zy33*h3;
Zy4 = Zy04*h0+ Zy14*h1+ Zy24*h2+ Zy34*h3;

e4 = -((rmu*h0+ rsigma*h1)*y1)+y2-y3;
e5 = y2-C;
e6 = -(y2-y3)+(rmu*h0+ rsigma*h1)*(cmu*h0+ csigma*h1)*Zy3+y3-y4;
e7 = -(y3-y4)+(rmu*h0+ rsigma*h1)*(cmu*h0+ csigma*h1)*Zy4+(y4*(rmu*h0+ rsigma*h1)/RL);

disp('From Equation 4')
int(e4*h0*w,X,-inf,inf) 
int(e4*h1*w,X,-inf,inf)
int(e4*h2*w,X,-inf,inf) 
int(e4*h3*w,X,-inf,inf) 

disp('From Equation 5')
int(e5*h0*w,X,-inf,inf) 
int(e5*h1*w,X,-inf,inf)
int(e5*h2*w,X,-inf,inf) 
int(e5*h3*w,X,-inf,inf) 


disp('From Equation 6')
int(e6*h0*w,X,-inf,inf) 
int(e6*h1*w,X,-inf,inf)
int(e6*h2*w,X,-inf,inf) 
int(e6*h3*w,X,-inf,inf) 

disp('From Equation 7')
int(e7*h0*w,X,-inf,inf) 
int(e7*h1*w,X,-inf,inf)
int(e7*h2*w,X,-inf,inf) 
int(e7*h3*w,X,-inf,inf) 
