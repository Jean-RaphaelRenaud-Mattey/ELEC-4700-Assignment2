%%ELEC 4700 Assignment 2: Finite Difference Method
%% Student Name: Jean-Raphael Renaud-Mattey
%% Student Number: 100791195
%% Date Submitted: February 25 2018
%% Question 1 of Assignment 2
% The purpose of this code is to use the Finite Difference Method to solve
% for the electrostatic potential in the rectangular region L*W shown in
% Figure 1 of the Assignment, where the ratio of L/W is 3/2. First we solve
% the case of V = Vo at x = 0 and V = 0 at x = L. For this case,I treat it
% as a 1-D case and solve it via iterations process. The second part of
% question 1 requires me to now solve it for V=Vo at x=0, x=L and V=0 at y=0,
% y=W, which will then be compared to a analytical solution. Some conclusions
% will be drawn from it and explained at the end of the code for the first
% question. I used iteration process to solve the first question:


%clear all
clearvars
clearvars -GLOBAL
close all
format shortE
%% Question 1 parameters:
% Here we will define all the parameters for question 1,as well the region
% size for both questions. Any change to see different result can be made
% here.
size = 20; %size to establish the Length and Width of the rectangular region
Vo=1; %The voltage set the boundary conditions as well to find the analytical solution
L = 3*size; %Length of the rectangular region
W = 2*size;%Width of the rectangular region
G1 = sparse(1,L); % Vector to solve question 1 a via iteration process
G2 = sparse(W,L); % Matrix to solve question 1 b via iteration process
G3 = sparse(W,L); % Matrix to solve question 1 b analytical sol
BC1 = 0; % Boundary condition for question 1
BC2 = Vo;% Boundary condition for question 1
iterations =100; %iteration limit for the iteration process, can be changed for the different result

%% Question 1.a):
% This is the iteration process for solving the 1-D case. It plots the
% change and update over the number of iterations. This is answer for
% question 1.a)
f1 = figure; 
G1(L)=BC1; % Set the boundary condition for the Vector to be V=0 at x=L 
G1(1)=BC2; % Set the boundary condition for the Vector to be V=Vo at x=0
G1(2:L-1)=deal(0.5*(BC1+BC2)); %populate the rest of the vector with the average of the Boundary Conditions
figure(f1)
for steps =1:iterations
    plot(G1);
    title('1-D finite iteration process')
    pause(0.01);
    for i=2:L-1
        G1(i)= 0.5*(G1(i-1)+G1(i+1));
    end
end

%% Question 1.b):
% This is the code to answer question 1.b of the assignment. First we
% establish the analytical solution for the equation :
%
% $$V(x,y)=\frac{4Vo}{\pi}\sum_{n=1,3,5...}^\infty\frac{1}{n}\frac{\cosh\left(\frac{n\pi x}{a}\right)}{\cosh\left(\frac{n\pi b}{a}\right)}\sin\left(\frac{n\pi y}{a}\right)$$
%
% And then we will compare the analytical solution to the meshing and
% explain some of the advantages and disadvantages.
Nterms = iterations;
f2 = figure;
figure(f2)
for i = 1:L       
    for j = 1:W
        G3(j,i) = analyticalsol2(i,j,Vo,L,W,Nterms); 
    end
end
subplot(3,1,1);
surf(G3);
title('analytical solution');
pause(0.01);

% This is to populate our 2-D matrix where; for the boundary conditions
% are: V=Vo at x=0,x=L and V=0 at y=0, y=W
% For the rest, we populate it with the average of the boundary conditions
[G2(1,1),G2(1,L),G2(W,1),G2(W,L),G2(2:W-1,2:L-1)]=deal( 0.5*(BC1+BC2) );
[G2(1,2:L-1),G2(W,2:L-1)]=deal(BC1);
[G2(2:W-1,1),G2(2:W-1,L)]=deal(BC2);
%This code section is similar to the first question, except instead of
%being the iteration process for 1-D, it will be for a 2-D iteration
%process, so four terms must be used: (j+1,i), (j-1,i), (j,i+1), (j,i-1)
%and must be divided by four instead of 2.
for step = 1:iterations
    subplot(3,1,2);
    surf(G2);
    title('numerical solution');
    subplot(3,1,3);
    surf(log(abs(G2-G3)))
    set(gca,'ZScale','log')
    title('difference between analytical and numerical')
    pause(0.01);
    for i=2:L-1
        for j=2:W-1
            G2(j,i) = 0.25*(G2(j+1,i)+G2(j-1,i)+G2(j,i+1)+G2(j,i-1));
        end
    end
end

%% Conclusion to question 1.b):
% Through iteration process, we can see that the meshing does approach the
% analytical solution around 100 iterations, although they both have
% reading error in the range of -10e-2, which seems to be the approrpiate
% range when to stop the analytical series. Going any higher than 100 and
% the graph barely has any change. This shows that the advantages of
% numerical solution is that for partial differential equations, it is
% easier to use numerical method since it will approximate to the actual
% value for a complicated analytical problem, however since it is simple
% meshing, the process will take a long time. A disadvantage to numerical
% method is that is difficult to use it for non rectangular geometric and
% for non-uniform meshing.



%% Question 2 of Assignment 2
%  For this section, I must use the Finite Difference Method to solve the
%  current flow in the rectangular region but with different resistivity,
%  higher risistivity outside the box compared to inside. From we will plot
%  the sigma(x,y), V(x,y) Ex Ey, and J(x,y). And then will draw some
%  conclusion from varying the mesh density, space of the bottle neck, and
%  changing the conductivity of the box (sigma).

%% Question 2 parameters:
% These are the parameters to determine the bottle neck size. The length of
% the bottle neck is about 1/3 of the total length of the rectangular
% region, while the width is about 2/8 of the total width of the
% rectangular region. Our conductivity inside the boxes is sig = 10e-2.
% These values will change when comparing the initial results with
% different parameters to see the results.

cMap2=ones(L,W); % create the conductivity mapping
LGap = L/3; %Establish the length of the bottle neck
WGap =3*W/8; % Width of the boxes enclosing the bottle neck
sig = 10e-2;
[cMap2(floor((L-LGap)/2):floor((L+LGap)/2),1:floor((W-WGap)/2)),...
    cMap2(floor((L-LGap)/2):floor((L+LGap)/2),floor((W+WGap)/2):W)]= deal(sig);
% The previous two lines finalize the bottle neck for our conductivity map.
% The following code is to create the G matrix and B vector
G4 = sparse(L*W,L*W); 
B3 = zeros(1,L*W);
%equations will be 
%4T_{i,j} -T_{i-1,j} - T_{i+1,j}-T_{i,j-1} - T_{i,j+1}= 0
%T_{1,j} = V0=BC2
%T_{L,j} = V0=BC2
%T_{i,1} = 0=BC1
%T_{i,W} = 0=BC1
%T_{1,1}=T_{L,1}=T_{1,W}=T_{L,W}=V0/2=(BC1+BC2)/2
[B3(2:W-1),B3((L-1)*W+1:L*W-1)]=deal(BC2);% Populate the vector with 1 volt for the values that will interact with G4
[B3(1),B3(W),B3((L-1)*W),B3(L*W)]=deal(0.5*(BC1+BC2)); %Set the boundary conditions in the vector that will interact with G4
for i=1:L
    for j=1:W
        n= j + (i-1) * W;   
        if i==1 || i == L || j == 1 || j == W
            G4(n,n) = 1;
        else
            nxp = j+i*W;
            nxm = j+(i-2)*W;
            nyp = n+1;
            nym = n-1;
            
            rxp = 0.5* (cMap2(i,j)+cMap2(i+1,j));
            rxm = 0.5* (cMap2(i,j)+cMap2(i-1,j));
            ryp = 0.5* (cMap2(i,j)+cMap2(i,j+1));
            rym = 0.5* (cMap2(i,j)+cMap2(i,j-1));
            
            G4(n,n) = -(rxp+rxm+ryp+rym);
            G4(n,nxp)=rxp;
            G4(n,nxm)=rxm;
            G4(n,nyp)=ryp;
            G4(n,nym)=rym;
        end    
    end
end

V2Vec=G4\(B3'); % The matrix neccessary to obtained our V(x,y) from the conductivity map

V2Map= zeros(L,W);
for i=1:L
    for j=1:W
        V2Map(i,j)=V2Vec(j+(i-1)*W);
    end
end
%% Process to obtain the Electric Field values for the vector format of the Electric Field 
% The code lines 187-206 is to obtained our absolute value of the electric
% field, where we create two matrix , E(x), E(y), and depending on the
% location of the matrix, set the appropriate value from gradient equation.
eX= zeros(L,W);
eY= zeros(L,W);
for i=1:L
    for j=1:W
        if i == 1
            eX(i,j) = V2Map(i+1,j)-V2Map(i,j);
        elseif i == L
            eX(i,j) = V2Map(i,j)-V2Map(i-1,j);
        else
            eX(i,j) = (V2Map(i+1,j)-V2Map(i-1,j))*0.5;
        end
        if j == 1
            eY(i,j) = V2Map(i,j+1)-V2Map(i,j);
        elseif j == W
            eY(i,j) = V2Map(i,j)-V2Map(i,j-1);
        else
            eY(i,j) = (V2Map(i,j+1)-V2Map(i,j-1))*0.5;
        end
    end
end

eX=-eX;
eY=-eY;
%% Code responsible for creating the Current Density Vectors
% The following code lines are responsible for making the Current Density
% Vector, by first assigning the value of the current Density to their
% matching position with the Electric Field (since J = sig*E). 
Jx=cMap2.*eX;% This line of code create the X value of the Current Density vector
Jy=cMap2.*eY;% This line of code create the Y value of the Current Density vector

pointSampling = 3;
[x,y]= meshgrid(1:floor(W/pointSampling),1:floor(L/pointSampling));
x=pointSampling*x-pointSampling+1;
y=pointSampling*y-pointSampling+1;

eXSample = zeros(floor(L/pointSampling),floor(W/pointSampling));
eYSample = zeros(floor(L/pointSampling),floor(W/pointSampling));
JxSample = zeros(floor(L/pointSampling),floor(W/pointSampling));
JySample = zeros(floor(L/pointSampling),floor(W/pointSampling));

for i=1:floor(L/pointSampling)
    for j = 1:floor(W/pointSampling)
        eXSample(i,j)=eX(pointSampling*i,pointSampling*j);
        eYSample(i,j)=eY(pointSampling*i,pointSampling*j);
        JxSample(i,j)=Jx(pointSampling*i,pointSampling*j);
        JySample(i,j)=Jy(pointSampling*i,pointSampling*j);
    end
end
%% Plot of the question 2.a)
% This is plots for the conductivity sig(x,y), V(x,y), E(x,y), and J(x,y).
% This will be standard plot that will compare to the following figures,
% where the setting are modified to investigate the changes. To avoid huge
% bulk of code, I have taken the code for answering question 2.a and put it
% in a function file to reduce the size of the main code.
f3=figure;
figure(f3)
subplot(2,2,1)
surface(cMap2');
title('Conductivity');
subplot(2,2,2)
surface(V2Map');
title('Voltage Vx,Vy');
subplot(2,2,3)
quiver(y,x,eXSample,eYSample,10);
title('Electric Field, Ex,Ey');
subplot(2,2,4)
quiver(y,x,JxSample,JySample,10);
title('Current Density, Jx,Jy');
%% Answer of Question 2.b :
% For here I created a function code that depends on inputs such as size,
% conductivity of the box creating the bottle neck, and the point sampling
% ratio, to change the size of the arrow for the quiver function.
%% This is with meshing size being double of the standard set:

Plotconductivity(40,10e-2,3);
%% This is with meshing size being triple of the standard set:

Plotconductivity(80,10e-2,3);

%% This is with meshing size being half of the standard set:

Plotconductivity(10,10e-2,3);

%% Answer of Question 2.c:
% Here are the different version of the a set with different bottle
% neck size, the conductivity will be the same, the meshing size will be
% the same, just the bottle neck dimension will change:
%
%% This is the set the comparison will be made from:

Plotconductivity(100,10e-2,3);
%% This is one with longer bottle neck:
Plotconductivity2(100,10e-2,3,0.5,0.375);
%% This is one with narrower bottle neck:
Plotconductivity2(100,10e-2,3,0.3333,0.5);

%% Answer of Question 2.d:
% In this question, we wil compare the results with different resistivity
% inside the box:
%% This one is with higher resistivity than the original
Plotconductivity(20,10,3);

%% This one is with a lower resistivity than the original
Plotconductivity(20,10e-5,3);