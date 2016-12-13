clear all;
%
format shorte;
data = importdata('DispTimeData2.txt');
data = data.data;
%
dtime = data(:,1);   Vnum = data(:,4);
alpha = data(:,2);   IntEdge = data(:,5);
Vsize = data(:,3);   IntFin = data(:,6);
%
Vtotal = Vnum.*((2*Vsize).^2);
%
Vval = zeros(256,1);
realval = zeros(38,2);
intEdgebin = zeros(38,1);  intFinbin = zeros(38,1);
save = zeros(5,1);   save2 = zeros(5,1);
yy = 1;
zz = 1;
for tt = 2:2:512
    Vval(tt) = sum(Vtotal==tt);
    if any(Vval(tt)) == 1
        realval(yy,2) = 2*zz;
        realval(yy,1) = Vval(tt);
        yy = yy+1;
    end
    
    zz=zz+1;
end
ind = 1;
for rstep = 1:length(realval(:,1))
   
    for vstep = 1:length(Vtotal)
        if Vtotal(vstep) == realval(rstep,2)
            save(ind) = IntEdge(vstep);
            save2(ind) = IntFin(vstep);
            ind=ind+1;
        end
        
    end
    intEdgebin(rstep) = (1/ind)*sum(save);
    intFinbin(rstep) = (1/ind)*sum(save2);
end

%
IntL = length(IntEdge);
%
figure(1);
plot(realval(:,2),intEdgebin,'k+'); ylabel('Pressure Magnitude at Far Plane Edge at Initial Reflection');...
    xlabel('Total Number of Vacancy Points'); title('Pressure Magnitude on Far Edge at Wave Front vs Vacancies');
figure(2);
plot(realval(:,2),intFinbin,'b+'); ylabel('Average Final Pressure Magnitude Over Plane');...
    xlabel('Total Number of Vacancy Points'); title('Average Final Pressure Magnitude vs Vacancies');
%figure(3);
%plot(IntEdge,Vtotal,'r+');
%figure(4);
%plot(IntFin,Vtotal,'m+');
  %  xlabel('Time for Wave Front to Reach Far Edge');...
   % ylabel('Total Number of Vacancy Points'); title('Time to Reach Edge vs Vacancies');
[xData, yData] = prepareCurveData( realval(:,2), intEdgebin );

% Set up fittype and options.
ft = fittype( 'rat21' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf];
opts.StartPoint = [0.655477890177557 0.171186687811562 0.706046088019609 0.0318328463774207];
opts.Upper = [Inf Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'Vacancies vs Magnitude', 'exp fit', 'Location', 'NorthEast' );
% Label axes
ylabel( 'Pressure Magnitude at Far Plane Edge at Initial Reflection' );
xlabel( 'Total Number of Vacancy Points' ); title('Pressure Magnitude on Far Edge at Wave Front vs Vacancies');
grid on
%
[xData2, yData2] = prepareCurveData( realval(:,2), intFinbin );
opts.StartPoint = [1.49221668551578e+26 5.63433324918348];
%
ft = fittype( 'exp2' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf];
opts.StartPoint = [5.57612715250609e-05 0.000139661052611739 -2.61255607466845e-05 -0.0318621717954263];
opts.Upper = [Inf Inf Inf Inf];
% Fit model to data.
[fitresult, gof] = fit( xData2, yData2, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData2, yData2 );
legend( h, 'Vacancies vs Final Magnitude', 'exp fit', 'Location', 'SouthEast' );
% Label axes
ylabel( 'Average Final Pressure Magnitude Over Plane' );
xlabel( 'Total Number of Vacancy Points' ); title('Average Final Pressure Magnitude vs Vacancies');
grid on

