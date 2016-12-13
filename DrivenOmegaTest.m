clear all;
%
A_o = 0.2; %input('Enter Intial Amplitude: ');
vacsize = input('Enter Vacancy Size: ');
vacnum = 0; %input('Enter Total Number of Vacancies: ');
al = 2.5; %input('Enter Attenuation Coefficient Value (1-10): ');
%approximating/modifying attenuation coefficient to fit the model
alpha = 0.00002*al;
omega = 0.5; %input('Enter Angular Frequency: ');
for vn = 1:10
    vacnum = vacnum+1;
    omega = 0.5;
for om = 1:10 %record increasing angular frequency values for each vacancy number
    omega = omega+1.5;
%Set up the coordinates
res=[102;102];
x = linspace(0,1,res(1));
y = linspace(0,1,res(2));
%Calculate the cell spacing based on the resolution
dx = x(2)-x(1);
dy = y(2)-y(1);
C0 = 2; %Wave Speed must be an arbitrary initial value since simulation is generalized
dt = dx/(2*C0); %Compute timestep resolution, dt
tmax = 1/dt;
%Maximum timestep iteration (based on dt)
%Initialize the Pressure matrix with 3D size dimensions
Pres  = zeros(length(x),length(y),3);   
%set dimensions of vacancies based on user input
xi = 20+vacsize; yi = 3+vacsize; xf = 100-vacsize; yf = xf;
%create matrices to record relevant data during loop(s)
vacindex = zeros(vacnum,3);
LeadAvg = zeros(10,1);
Avgx = zeros(101,1);
Avgy = zeros(101,1);
xL = zeros(vacnum,1); xH = zeros(vacnum,1);
yL = zeros(vacnum,1); yH = zeros(vacnum,1);
%Establish iteration limiters' initial values
first = 1;
disp = 1;
G = 1;
%
Pres(:,:,1) = Pres(:,:,2);%Ensure that tstep-1 timestep is also initialized
for tstep = 1:tmax   %loop hard!
    if G == 1
       for kk = 1:vacnum %Set and record randomly generated locations for vacancies
                vacindex(kk,1) = kk;
                vacindex(kk,2) = randi([xi, xf],1);
                vacindex(kk,3) = randi([yi, yf],1);
                %Set upper and lower boundary on the vacancies based on their size
                xL(kk,1) = vacindex(kk,2)-vacsize; yL(kk,1) = vacindex(kk,3)-vacsize;
                xH(kk,1) = vacindex(kk,2)+vacsize; yH(kk,1) = vacindex(kk,3)+vacsize;
       end
       G = 0; %Set the condition to 0 so it doesn't make new vacancies each loop
    end
    %Step through the x and y directions
    for xstep = 2:length(x)-1
        for ystep = 2:length(y)-1 
            Pres(2,:,2) = A_o*cos(omega*tstep);
            %Solve for the pressure amplitude at each point   
             Pres(xstep,ystep,3) = dt^2*C0^2*((Pres(xstep-1,ystep,2)-2*Pres(xstep,ystep,2)...
             +Pres(xstep+1,ystep,2))/dx^2+(Pres(xstep,ystep-1,2)-2*Pres(xstep,ystep,2)...
             +Pres(xstep,ystep+1,2))/dy^2)+2*Pres(xstep,ystep,2)-Pres(xstep,ystep,1);
             %Multiply by dissipative term w/ attentuation coefficient alpha
             Pres(xstep,ystep,3) = exp(-alpha*disp)*Pres(xstep,ystep,3);
             %
                for LL = 1:vacnum
                    %Set Pressure wave to zero inside the vacancies for Neumann Boundary reflection
                    if (yL(LL,1)<=ystep)&&(ystep<=yH(LL,1))&&(xL(LL,1)<=xstep)&&(xstep<=xH(LL,1))
                        Pres(xstep,ystep,3) = 0; %Vacancies should resemble perfect vacuum; perfect reflection
                    end
                end 
        end 
    end
    %
    tTest = any(any(abs(Pres(101,:,2))>= 0.0005));
    %record the number of time units it took for the wave to reach the far
    %side of the sheet with a significant amplitude (0.0005)
    %significant amplitude may be much smaller with smaller initial amplitude
    if tTest ~= 0 && first <= 10
        if first == 5
            DispTime = tstep*dt;
            first = first+1;
        else
            tscale = tmax-tstep;
            LeadAvg(first) = (1/100)*sum(abs(Pres(101,:,2))); %Record Average on far edge at the wave front
            first = first+1;
        end
    end
    %Calculate the new Pres(xstep,ystep,tstep+1) at each internal grid location
    %surf(x,y,Pres(:,:,3)); %Plot the latest tstep+1 3D surface
    %xlabel('X Direction'); ylabel('Y Direction');
    %zlabel('Pressure Intensity/Amplitude');
    %caxis([-0.095 0.095]);      %Fix the color range.
    %axis([0 1 0 1 -.25 .25]);  %Fix figure axes
    %text(0.1,1,0.22,sprintf('time=%.2f',(tstep-1)*dt)) %Display time
    %Save the current figure in an array of frames
    %frames(tstep) = getframe; 
    %
    Pres(:,:,1) = Pres(:,:,2);
    Pres(:,:,2) = Pres(:,:,3);
    disp = disp+1;
    %timestep tstep becomes the new tstep-1
    %timestep tstep+1 becomes the new tstep
    %(tstep+1 will be overwritten in the next iteration)
    
end
%Find average value at wave front after loops finish and record it
LeadAverage = (1/first)*sum(LeadAvg(:,1));
AvgVal = 100*abs((1/(10000-vacnum*(2*vacsize)^2))*sum(Pres(:,:,2)));
AvgValue = (1/100)*sum(AvgVal); %find the average pressure value across the whole plane at the 
%end of the loops and record it
format shorte;
sprintf('Time to reach (x=102) is: %.2f',DispTime) %Print time to reach far edge on screen
fid = fopen('DrivenOmegaData.txt','a'); %Open the file for recording and stockpiling data
%
if LeadAverage ~= 0  
    %Record to the file if the simulation hasn't bugged out and decided the wave never reached the far side
    %Which happens about 1/200 simulations... probably because my laptop is decidedly mediocre
    fprintf(fid,'%e %e %e %e %e %e %e \r\n',[DispTime A_o vacsize vacnum LeadAverage AvgValue omega]);
    fclose(fid);
end
end
end
%Yay!
%