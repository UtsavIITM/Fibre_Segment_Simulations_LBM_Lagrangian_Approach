% Code Developed for elementary Air Filtration Simulations 
% PhD Project (UQIDAR 00218), The University of Queensland (Australia) and IIT Delhi (India), 2021-2024
% Utsav Bhardwaj, Dr Bahni Ray, Prof Dipayan Das, Dr Travis Mitchell, Prof Apurba Das, and Dr Christopher Leonardi
% Lattice Boltzmann Method (LBM) + Lagrangian Approach  
% Contact email address: utsav_bhardwaj@alumni.iitm.ac.in

pii = 3.14;
mua_Dim = 1.825e-5;
nu_Dim = 1.516e-5;
tau = 0.75;

b_Dim = 4e-6; %% fibre diameter
AR = 4;
L_fib_Dim = (AR)*(b_Dim);

dp_Dim = 2e-6; %% Diameter of particle
dx_Dim = (dp_Dim)/(13);
dy_Dim = dx_Dim;
Cl = dx_Dim;

dx = (dx_Dim)/(Cl); 
dy = (dy_Dim)/(Cl);
dt = 1;
CsSqr = (1/3)*(dx/dt)*(dx/dt);
nu=(CsSqr*dt*((tau/dt)-0.5));
dt_Dim = (CsSqr)*(tau - 0.5)*((dx_Dim*dx_Dim)/nu_Dim);
Ct = dt_Dim;
Cu = (Cl/Ct);

uain_Dim = 0.9;
uain = (uain_Dim)/(Cu);

b = round((b_Dim)/(Cl));
L_fibre = (AR)*b;

rho_Dim = 1.204;
C_rho = rho_Dim;
rho_NDim = (rho_Dim)/(C_rho);
rhoOut = rho_NDim;

conv_limit = 5e-16;

Li = round(5*((L_fib_Dim/2)/Cl)); 
Lf = round(5*((L_fib_Dim/2)/Cl)); 
Wi = round(2*((L_fib_Dim/2)/Cl)); 
Wf = round(2*((L_fib_Dim/2)/Cl)); 
nl = Li+1+Lf;
nw = Wi+1+Wf; 

Fibre_Vol_Fraction = (b*(L_fibre))/((Li+Lf)*(Wi+Wf));

Total_time_one_pass = ((Li+Lf)*Cl)/(uain_Dim);
Number_of_passes = 5;
Total_time_all_passes = (Number_of_passes)*(Total_time_one_pass);
maxnt = (Total_time_all_passes)/(dt_Dim); 

%% Initializations
x = zeros(nw,nl);
y = zeros(nw,nl);
% x = 0:dx:((nl-1)*dx);
% y = 0:dy:((nw-1)*dy);

Fluid = zeros(nw,nl);
uax = zeros(nw,nl);
uaxPrev = zeros(nw,nl);
uay = zeros(nw,nl);
uayPrev = zeros(nw,nl);
rho = zeros(nw,nl);
rhoPrev = zeros(nw,nl);
% rhoOut = zeros(nw,nl);

f = zeros(nw,nl,9);
feq = zeros(nw,nl,9);
fcoll = zeros(nw,nl,9);
fneq = zeros(nw,nl,9);

w = [(1/9) (1/9) (1/9) (1/9) (1/36) (1/36) (1/36) (1/36) (4/9)];
Cx = [1 0 -1 0 1 -1 -1 1 0];
Cy = [0 1 0 -1 1 1 -1 -1 0];

k = zeros(9,1);
kb1 = zeros(9,1); 
kb2 = zeros(9,1);

%% Mesh Generation 
for j = 1:1:nl
    y(1,j) = 0;
    for i = 1:1:(nw-1)
        y(i+1,j)=y(i,j)+dy;
    end
end

for i = 1:1:nw
    x(i,1) = 0;
    for j = 1:1:(nl-1)
        x(i,j+1)=x(i,j)+dx;
    end
end
    

% %% Defining Solid and fluid Zones %
% dff = b;
% Lff = (AR)*b;
% 
% theta_degree = 30;
% theta = (pii/180)*(theta_degree);
% dydx = tan(theta);
% 
% a1 = (dff*((cos(theta))^2))/(sin(theta));
% a2 = dff*(sin(theta));
% 
% z1 = dff*(sin(theta))*(cos(theta));
% Hf1 = ((Lff/2)-z1)*sin(theta);
% z2 = (dff*((cos(theta))^3))/(sin(theta));
% Hf2 = ((Lff/2)-z2)*sin(theta);
% 
% Hf3 = (dff)*(cos(theta)); 
% Hf4 = (dff)*(cos(theta));
% 
% a3 = Hf2/(tan(theta));
% a4 = Hf1/(tan(theta));
% 
% for j=1:1:nl
%     for i=1:1:nw
%         if y(i,j)>=(Wi-Hf1) && y(i,j)<=(Wi+Hf2) && x(i,j)>=(((Li-a1)+((1/tan(theta))*(Hf1)))-((1/tan(theta))*(y(i,j)-(Wi-Hf1)))) && x(i,j)<=(((Li+a2)+((1/tan(theta))*(Hf1)))-((1/tan(theta))*(y(i,j)-(Wi-Hf1))))
%             Fluid(i,j)= 0;
%         elseif y(i,j)>=(Wi-(Hf1+Hf3)) && y(i,j)<=(Wi-Hf1) && x(i,j)>=(((Li-a1)+((1/tan(theta))*(Hf1+Hf3)))-((1/tan(theta))*(y(i,j)-(Wi-(Hf1+Hf3))))) && x(i,j)<=(((Li-a1)+((1/tan(theta))*(Hf1+Hf3)))+((tan(theta))*(y(i,j)-(Wi-(Hf1+Hf3))))) 
%             Fluid(i,j)= 0;
%         elseif y(i,j)>=(Wi+Hf2) && y(i,j)<=(Wi+(Hf2+Hf4)) && x(i,j)>=(((Li-a1)-((1/tan(theta))*(Hf2)))+((tan(theta))*(y(i,j)-(Wi+Hf2)))) && x(i,j)<=(((Li+a2)-((1/tan(theta))*(Hf2)))-((1/tan(theta))*(y(i,j)-(Wi+Hf2)))) 
%             Fluid(i,j)= 0;
%         else 
%             Fluid(i,j)= 1;
%         end
%     end
% end


%% Defining Solid and fluid Zones for 0 and 90 Degree %
dff = b;
Lff = (AR)*b;

for j=1:1:nl
    for i=1:1:nw
        if x(i,j)>=(Li-(dff/2)) && x(i,j)<=(Li+(dff/2)) && y(i,j)>=(Wi-(Lff/2)) && y(i,j)<=(Wi+(Lff/2)) 
            Fluid(i,j) = 0;
        else
            Fluid(i,j) = 1;
        end
    end
end


%% Force Calculation %%
 XForce = 0.0;
 YForce = 0.0;
 Fdrag = 0.0;
 Flift = 0.0;


%% Pairs required for free-slip boundary condition (specular reflection) at top and bottom walls of flow channel // 

for k = 1:1:9
    for k1=k:1:9
        if Cx(k1)==Cx(k) && Cy(k1)+Cy(k)==0
            kb1(k)=k1;
            kb1(k1)=k;
        end
    end
end

%% Pairs required for no-slip boundary condition (bounce back) at the fibre-air interface //

for k = 1:1:9
    for k1 = k:1:9
        if Cx(k1)+Cx(k)==0 && Cy(k1)+Cy(k)==0
            kb2(k)=k1;
            kb2(k1)=k;
        end
    end
end

%% Initialization of the PDF (using Eqbm. PDF) %%
nt=0;
for j = 1:1:nl
    for i = 1:1:nw
        for k = 1:1:9
            
            uax(i,j) = 0;
            uay(i,j) = 0;
            rho(i,j) = rho_NDim;
            c_dot_ua = ((Cx(k))*(uax(i,j))) + ((Cy(k))*(uay(i,j)));
            ua_dot_ua = ((uax(i,j))*(uax(i,j))) + ((uay(i,j))*(uay(i,j)));
            feq(i,j,k) = (w(k))*(rho(i,j))*(1 + ((c_dot_ua)/(CsSqr)) + ((c_dot_ua*c_dot_ua)/(2*CsSqr*CsSqr)) - ((ua_dot_ua)/(2*CsSqr)));
            f(i,j,k) = feq(i,j,k);
        end
    end
end

%% Entering the Time Loop //

while nt<=(maxnt-1)
    
    RhoSum = 0;
    RhoAvg = 0;
    
    % % Calculation of pre-collision density and velocities //
    
    for j = 1:1:nl
        for i = 1:1:nw
            
            uaxPrev(i,j) = uax(i,j);
			uayPrev(i,j) = uay(i,j);
			rhoPrev(i,j) = rho(i,j);
			
	% % local lattice_point_level initialization (zeroing) //
			rhosum = 0;
			uaxsum = 0;
			uaysum = 0;
            
            for k = 1:1:9
                
                rhosum = rhosum + f(i,j,k);
				uaxsum = uaxsum + (Cx(k))*(f(i,j,k));
				uaysum = uaysum + (Cy(k))*(f(i,j,k));
            end
            
            rho(i,j) = rhosum;
            uax(i,j) = (uaxsum)/(rho(i,j));
			uay(i,j) = (uaysum)/(rho(i,j));
            
            RhoSum = RhoSum + rho(i,j);
            
        end
    end
    
    RhoAvg = (RhoSum)/((nl)*(nw));
    
    % % Collision Operation at Time = nt; Finding post-collision PDFs (f_star) //
    
    for j  = 1:1:nl
        for i = 1:1:nw
            for k = 1:1:9
                
                c_dot_ua = ((Cx(k))*(uax(i,j))) + ((Cy(k))*(uay(i,j)));
                ua_dot_ua = ((uax(i,j))*(uax(i,j))) + ((uay(i,j))*(uay(i,j)));
                feq(i,j,k) = (w(k))*(rho(i,j))*(1 + ((c_dot_ua)/(CsSqr)) + ((c_dot_ua*c_dot_ua)/(2*CsSqr*CsSqr)) - ((ua_dot_ua)/(2*CsSqr)));
                Omega_alpha = (-1.0)*((f(i,j,k) - feq(i,j,k))/(tau));
                fcoll(i,j,k) = f(i,j,k) + ((Omega_alpha)*(dt));
            end
        end
    end
    
    % % Streaming Operation: Moving from TIME [nt] to [nt+1] //
	nt = nt + 1;
    for j = 1:1:nl
        for i = 1:1:nw
            for k = 1:1:9
                
                % % Streaming and free slip for Fluid domain
				 	if Fluid(i,j)==1
                        ja = j + (Cx(k))*(dt);
                        ia = i + (Cy(k))*(dt);
                        
                        if ja>=1 && ja<=nl && ia>=1 && ia<=nw
                            
                            if Fluid(ia,ja)==0
                                f(i,j,(kb2(k))) = fcoll(i,j,k);  %% Bounce back inside fluid domain (at the external surface of fibre, back into the fluid)
                            else 
                                f(ia,ja,k) = fcoll(i,j,k); %% Streaming inside fluid domian (within the fluid domain)
                            end
                        end
                        
                        if ja>=1 && ja<=nl && ia>nw  
                            f(i,ja,kb1(k)) = fcoll(i,j,k); %% Free slip (Half-way Specular Reflection) at the bottom boundary of the flow channel ("hypothetical")
                        end
                        
                        if ja>=1 && ja<=nl && ia<1  
                            f(i,ja,kb1(k)) = fcoll(i,j,k); %% Free slip (Half-way Specular Reflection) at the top boundary of the flow channel ("hypothetical")
                        end
                    end
                    
     % % Streaming inside Solid
                    
     if Fluid(i,j)==0
         ja = j + (Cx(k))*(dt);
         ia = i + (Cy(k))*(dt);
         
         if Fluid(ia,ja)==1
             f(i,j,kb2(k)) = fcoll(i,j,k); %% bounce back inside solid domain (at the internal surface of fibre, back into the solid (fibre))
         else
             f(ia,ja,k) = fcoll(i,j,k); %% Streaming inside Solid domain (within the fibre)
             
         end
     end
        
            end
        end
    end
    
    
 %% Velocity inlet boundary condition \\
    j=1;
	for i = 1:1:nw
           
		uax(i,j) = uain;
%       uay(i,j) = 0;
        rho(i,j)	= (f(i,j,9)+f(i,j,2)+f(i,j,4)+(2*(f(i,j,3)+f(i,j,6)+f(i,j,7))))/(1-uax(i,j));
		f(i,j,1)	=  f(i,j,3)+((2/3)*rho(i,j)*uax(i,j));
		f(i,j,5)	=  ((1/6)*rho(i,j)*uax(i,j)) + f(i,j,7) - (0.5*(f(i,j,2)-f(i,j,4)));
		f(i,j,8)	=  ((1/6)*rho(i,j)*uax(i,j)) + f(i,j,6) + (0.5*(f(i,j,2)-f(i,j,4)));
    end
    
    
%% Pressure Oulet boundary condition
	j = nl;
	for i = 1:1:nw
        
		 	rho(i,j)    = rhoOut;
		 	uax(i,j)	= ((f(i,j,9)+f(i,j,2)+f(i,j,4)+ (2*(f(i,j,1)+f(i,j,5)+f(i,j,8))))/(rho(i,j)))-1;
		 	f(i,j,3)	=  f(i,j,1)-((2/3)*rho(i,j)*uax(i,j));
		 	f(i,j,6)	=  ((-1/6)*rho(i,j)*uax(i,j))  + f(i,j,8) +(0.5*(f(i,j,4)-f(i,j,2))); 
		 	f(i,j,7)	=  ((-1/6)*rho(i,j)*uax(i,j))  + f(i,j,5) -(0.5*(f(i,j,4)-f(i,j,2)));
    end
            

    % Force 

    
		XForceOld = XForce;
		YForceOld = YForce;
		XForce = 0.0;
		YForce = 0.0;
		
			for j=1:1:nl
                for i=1:1:nw
                    for k=1:1:9

				 	if Fluid(i,j)==1

				 	 ja = j + (Cx(k))*(dt);
                     ia = i + (Cy(k))*(dt);
                     
						if ja>=1 && ja<=nl && ia>=1 && ia<=nw

							if Fluid(ia,ja)==0
					
								Fx(i,j,k) = 2.0*(fcoll(i,j,k) - f(ia,ja,kb2(k)))*Cx(k); %% X component of force at a given location and a given direction
								Fy(i,j,k) = 2.0*(fcoll(i,j,k) - f(ia,ja,kb2(k)))*Cy(k); %% Y component of force at a given location and a given direction
								XForce = XForce + Fx(i,j,k);
								YForce = YForce + Fy(i,j,k);
                            end
                        end
                    end
    				end
    			end
        	end
		
		if (nt>1)

				% Average the force with previous time step value
				XForce = 0.5*(XForce + XForceOld);
		 		YForce = 0.5*(YForce + YForceOld);
		 		% Drag and Lift
		 		Fdrag = XForce;
				Flift = YForce;
				CD = Fdrag/(0.5*uain*uain*(dff));
				CL = Flift/(0.5*uain*uain*(dff));
        end

%% calculation of NEW (Final) density and velocity components for the entire domain	

for j = 1:1:nl
        for i = 1:1:nw
            
	% % local lattice_point_level initialization (zeroing) //
			rhosum = 0;
			uaxsum = 0;
			uaysum = 0;
            
            for k = 1:1:9
                
                rhosum = rhosum + f(i,j,k);
				uaxsum = uaxsum + (Cx(k))*(f(i,j,k));
				uaysum = uaysum + (Cy(k))*(f(i,j,k));
            end
            
            rho(i,j) = rhosum;
            uax(i,j) = (uaxsum)/(rho(i,j));
			uay(i,j) = (uaysum)/(rho(i,j));
            
            RhoSum = RhoSum + rho(i,j);
            
        end
    end
    
    RhoAvg = (RhoSum)/((nl)*(nw));
    
%     if rem(nt,5000)==0
%         for j = 1:1:nw
%             for i = 1:1:nl
%                     z = [x(i);y(j);uax(i,j);uay(i,j)];
%                     fileID = fopen('FlowField.plt','w');
%                     fmt = '%5d %5d %5d %5d\n';
%                     fprintf(fileID,fmt,z);
%                     fclose(fileID);
%             end
%         end
%     end
    
%% Checking the Convergence //
	Res = 0;
    for j = 1:1:nl
        for i = 1:1:nw
            uaxRes = abs(uax(i,j) - uaxPrev(i,j));
            if uaxRes > Res
                Res = uaxRes;
            end
        end
    end
    
    if Res < conv_limit
        break
    end
            
end


%% Evaluation of Dimensionless Drag Force %%
Cf = ((Cl^4)*(C_rho))/(Ct^2);
FD_Dim = (Cf)*(Fdrag);
F_NDim = (FD_Dim)/((((rho_Dim)*((uain_Dim)^2))/2)*(L_fib_Dim)); 

%% Dimensionalization of Coordinate System (x,y) and Air Velocity Field (uax,uay) [Scaling]

for j = 1:1:nl
    for i = 1:1:nw
        x(i,j) = (Cl)*(x(i,j));
        y(i,j) = (Cl)*(y(i,j));
        uax(i,j)=(Cu)*(uax(i,j));
        uay(i,j)=(Cu)*(uay(i,j));
    end
end

%% Saving Data for plotting Streamlines

% for i=1:1:nw
%     for j=1:1:nl
%         xa(j+((i-1)*(nl))) = x(i,j);
%         ya(j+((i-1)*(nl))) = y(i,j);
%         ua(j+((i-1)*(nl))) = uax(i,j);
%         va(j+((i-1)*(nl))) = uay(i,j);
%         Fluid(j+((i-1)*(nl))) = Fluid(i,j);
%     end
% end
% 
% x_SL = xa';
% y_SL = ya';
% u_SL = ua';
% v_SL = va';
% Fluid_SL = Fluid';
% 
% save x_SL.csv x_SL -ascii
% save y_SL.csv y_SL -ascii
% save u_SL.csv u_SL -ascii
% save v_SL.csv v_SL -ascii
% save Fluid_SL.csv Fluid_SL -ascii

     

%% Motion of Particles

dp = (dp_Dim)/Cl; %% Dimensionless Diameter of particle
rp_Dim = (dp_Dim)/2;
rp = (rp_Dim)/Cl; 
rhoP_Dim = 912; %% Density of particle
pii = 3.14;
e = 2.718;

Captured = 0;
Bypassed = 0;

HH = 0; %% Height (altitude) in meter at which air filtration is taking place %%
T = 20 + 273.15; %% Temperature at which air filtration is taking place %%
P = 101325; %% Pressure at which air filtration is taking place %% 
lambda = (8.314*T)/(1.414*3.14*((3.57e-10)^2)*6.022e23*P*(exp(-(0.02897*9.8*HH)/(8.314*288)))); %% Mean free path of air %%
Cc = 1 + ((lambda/dp_Dim)*(2.514 + 0.8*exp(-0.55*(dp_Dim/lambda)))); %% Cunningham slip correction for an aerosol particle %% 
m_p = ((pii*((dp_Dim)^3))/6)*(rhoP_Dim); %% Mass of particle 

K = (3*pii*(mua_Dim)*(dp_Dim))/((m_p)*Cc); %% Constant

% tInitial = 0;
% h = 3.75e-5; %% time step for particle motion
% tFinal = ((120)*(h));


tInitial = 0.0;
TTotal_01 = ((Li-(dff/2))*(Cl))/(uain_Dim);
TTotal_02 = ((dff)*Cl)/(uain_Dim);
NT2 = 3*(dff);
h = ((TTotal_02)/NT2); % time step for particle motion
TTotal_03 = ((Lf-(dff/2))*(Cl))/(uain_Dim);
TTotal = TTotal_01 + TTotal_02 + TTotal_03;
N = 3*(round(TTotal/h)); % Total No. of time steps
tFinal = (N*h);


xpInitial = 0;
upxInitial = uain_Dim;
uaxInitial = uain_Dim;

upyInitial = 0;
uayInitial = 0;

%% for xp and yp
% qq = ((tFinal/h) + 1); 

% t = zeros(qq,1);
% xp = zeros(qq,1);
% upx = zeros(qq,1);
% yp = zeros(qq,1);
% upy = zeros(qq,1);


%% Number and Random Initial Positioning of particles %%
N_particles = (AR)*(500);

rng(0,'twister');

Min_ypInitial = 0 + (dp_Dim/500);
Max_ypInitial = (Wi+Wf)*Cl;
RNs_ypInitial = ((Max_ypInitial) - (Min_ypInitial)).*rand(N_particles,1) + (Min_ypInitial);



for particle = 1:1:(N_particles)
    
    ypInitial = RNs_ypInitial(particle);     

t(1) = 0;

xp(1) = xpInitial;
upx(1) = upxInitial;
uaxc = uaxInitial;

yp(1) = ypInitial;
upy(1) = upyInitial;
uayc = uayInitial;

Particles_Position_x(1,particle) = xpInitial;
Particles_Position_y(1,particle) = ypInitial;

Particles_Velocity_x(1,particle) = upxInitial;
Particles_Velocity_y(1,particle) = upyInitial;

Particles_Time(1,particle) = 0;


for z = 1:1:N
    
    tc = t(z);
    
    xpc = xp(z);
    upxc = upx(z);
    ypc = yp(z);
    upyc = upy(z);
    Fx = @(tc,xpc,upxc)(K*(uaxc-upxc));
    Fy = @(tc,ypc,upyc)(K*(uayc-upyc));
       
    
    dxpc1 = h*(upxc);
    dupxc1 = h*(Fx(tc,xpc,upxc));
    dypc1 = h*(upyc);
    dupyc1 = h*(Fy(tc,ypc,upyc));
        
    dxpc2 = h*(upxc + ((dupxc1)/2));
    dupxc2 = h*(Fx((tc+(h/2)),(xpc+((dxpc1)/2)),(upxc+((dupxc1)/2))));
    dypc2 = h*(upyc + ((dupyc1)/2));
    dupyc2 = h*(Fy((tc+(h/2)),(ypc+((dypc1)/2)),(upyc+((dupyc1)/2))));
    
    dxpc3 = h*(upxc + ((dupxc2)/2));
    dupxc3 = h*(Fx((tc+(h/2)),(xpc+((dxpc2)/2)),(upxc+((dupxc2)/2))));
    dypc3 = h*(upyc + ((dupyc2)/2));
    dupyc3 = h*(Fy((tc+(h/2)),(ypc+((dypc2)/2)),(upyc+((dupyc2)/2))));
    
    dxpc4 = h*(upxc + (dupxc3));
    dupxc4 = h*(Fx((tc+h),(xpc+(dxpc3)),(upxc+(dupxc3))));
    dypc4 = h*(upyc + (dupyc3));
    dupyc4 = h*(Fy((tc+h),(ypc+(dypc3)),(upyc+(dupyc3))));
    
    dxpc = ((dxpc1 + (2*dxpc2) + (2*dxpc3) + dxpc4)/(6));
    dupxc = ((dupxc1 + (2*dupxc2) + (2*dupxc3) + dupxc4)/(6));
    dypc = ((dypc1 + (2*dypc2) + (2*dypc3) + dypc4)/(6));
    dupyc = ((dupyc1 + (2*dupyc2) + (2*dupyc3) + dupyc4)/(6));


    t(z+1) = t(z)+h;
    xp(z+1) = xp(z) + (dxpc);
    upx(z+1) = upx(z) + (dupxc);
    yp(z+1) = yp(z) + (dypc);
    upy(z+1) = upy(z) + (dupyc);


    Particles_Position_x(z+1,particle) = xp(z+1);
    Particles_Position_y(z+1,particle) = yp(z+1);

    Particles_Velocity_x(z+1,particle) = upx(z+1);
    Particles_Velocity_y(z+1,particle) = upy(z+1);
    
    Particles_Time(z+1,particle) = t(z+1);



       
for j = 1:1:nl-1
    for i = 1:1:nw-1
        
        if xp(z+1)==x(i,j) && yp(z+1)==y(i,j)
            uaxc = uax(i,j);
            uayc = uay(i,j);
   
        elseif x(i,j)<xp(z+1) && xp(z+1)<x(i,j+1) && y(i,j)<yp(z+1) && yp(z+1)<y(i+1,j) 
            uaxc = ((uax(i,j) + uax(i,j+1) + uax(i+1,j) + uax(i+1,j+1))/4);
            uayc = ((uay(i,j) + uay(i,j+1) + uay(i+1,j) + uay(i+1,j+1))/4);
            
        elseif x(i,j)<xp(z+1) && xp(z+1)<x(i,j+1) && yp(z+1)==y(i,j) && yp(z+1)==y(i,j+1)  
            uaxc = ((uax(i,j) + uax(i,j+1))/2);
            uayc = ((uay(i,j) + uay(i,j+1))/2);
            
        elseif xp(z+1)==x(i,j) && xp(z+1)==x(i+1,j) && y(i,j)<yp(z+1) && yp(z+1)<y(i+1,j)  
            uaxc = ((uax(i,j) + uax(i+1,j))/2);
            uayc = ((uay(i,j) + uay(i+1,j))/2);
            
        end
    end
end

if xp(z+1)==x(1,nl) && xp(z+1)==x(nw,nl)
    for i=1:1:nw-1      
        if yp(z+1)==y(i,nl)
            uaxc = uax(i,nl);
            uayc = uay(i,nl);
        elseif y(i,nl)<yp(z+1) && yp(z+1)<y(i+1,nl)
            uaxc = ((uax(i,nl) + uax(i+1,nl))/2);
            uayc = ((uay(i,nl) + uay(i+1,nl))/2);
        end
    end
end

if yp(z+1)==y(nw,1) && yp(z+1)==y(nw,nl) 
    for j=1:1:nl-1        
        if xp(z+1)==x(nw,j)
            uaxc = uax(nw,j);
            uayc = uay(nw,j);
        elseif x(nw,j)<xp(z+1) && xp(z+1)<x(nw,j+1)
            uaxc = ((uax(nw,j) + uax(nw,j+1))/2);
            uayc = ((uay(nw,j) + uay(nw,j+1))/2);
        end
    end
end

if xp(z+1)==x(nw,nl) && yp(z+1)==y(nw,nl)
    uaxc = uax(nw,nl);
    uayc = uay(nw,nl);
end


xp_NDim = abs(round((xp(z+1))/Cl));

yp_NDim = abs(round((yp(z+1))/Cl));



if  xp_NDim>0 && xp_NDim<=(Li+Lf) && yp_NDim>0 && yp_NDim<=(Wi+Wf) && Fluid(yp_NDim,xp_NDim) == 0
    Captured = Captured + 1;
    break

elseif xp(z+1)>=(Li+Lf)*(Cl)
    Bypassed = Bypassed + 1;
    break
    
elseif yp(z+1)>=(Wi+Wf)*(Cl)
    Bypassed = Bypassed + 1;
    break

elseif yp(z+1)<=0
    Bypassed = Bypassed + 1;
    break

end

end

end

Filtration_Efficiency = (Captured)/(N_particles);
Press_Drop = (CD)*((0.5)*rho_Dim*(uain_Dim*uain_Dim));

for particle = 1:1:(N_particles)
    if particle ==1
        contourf(x,y,Fluid);
        hold on
    end
    plot(Particles_Position_x(:,particle),Particles_Position_y(:,particle))
    hold on
end


