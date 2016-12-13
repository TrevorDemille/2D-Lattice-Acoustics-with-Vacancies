clear all;
%
format shorte;
data = importdata('DrivenData.txt');
data = data.data;
%
dtime = data(:,1);   Vnum = data(:,4);   omega = data(:,7);
alpha = data(:,2);   IntEdge = data(:,5);
Vsize = data(:,3);   IntFin = data(:,6);
%
Vtotal = Vnum.*((2*Vsize).^2);
%
IntL = length(IntEdge);
E = IntL/10;
R = IntL/2;
avgdtime = zeros(100,1);
avgVtotal = zeros(E,1);
avgIntEdge = zeros(E,1);
avgIntFin = zeros(E,1);
IntEdgebin = zeros(100,1);
IntFinbin = zeros(100,1);
for cc = 1:100
    if cc <= 60
        avgdtime(cc) = (1/3)*(dtime(cc)+dtime(cc+100)+dtime(cc+200));
        IntEdgebin(cc) = (1/3)*(IntEdge(cc)+IntEdge(cc+100)+IntEdge(cc+200));
        IntFinbin(cc) = (1/3)*(IntFin(cc)+IntFin(cc+100)+IntFin(cc+200));
    else
        avgdtime(cc) = (1/2)*(dtime(cc)+dtime(cc+100));
        IntEdgebin(cc) = (1/2)*(IntEdge(cc)+IntEdge(cc+100));
        IntFinbin(cc) = (1/2)*(IntFin(cc)+IntFin(cc+100));
    end
end
for qq = 1:E
    avgVtotal(qq) = Vtotal(qq*10);
    avgIntEdge(qq) = (1/10)*sum((IntEdge(((qq*10)-9):(qq*10))));
    avgIntFin(qq) = (1/10)*sum((IntFin(((qq*10)-9):(qq*10))));
end
%
figure(1);
plot(avgVtotal,avgIntEdge,'k+'); ylabel('Pressure Magnitude at Far Plane Edge at Initial Reflection');...
    xlabel('Total Number of Vacancy Points'); title('Pressure Magnitude on Far Edge at Wave Front vs Vacancies');
figure(2);
plot(avgVtotal,avgIntFin,'b+'); xlabel('Average Final Pressure Magnitude Over Plane');...
    ylabel('Total Number of Vacancy Points'); title('Average Final Pressure Magnitude vs Vacancies');
figure(3);
plot(omega(1:100),IntEdgebin,'r'); xlabel('Angular Frequency');...
    ylabel('Pressure Magnitude at Far Plane Edge at Initial Reflection');...
    title('Angular Frequency vs Pressure Magnitude on Far Edge');
figure(4);
plot(omega(1:100),avgdtime); xlabel('Angular Frequency'); ylabel('Time for Wave Front To Reach Far Edge');...
    title('Angular Frequency vs Time to Reach Far Edge');
%
[xData, yData] = prepareCurveData( avgVtotal, avgIntEdge );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'Magnitude vs Vacancies ', 'Polynomial fit', 'Location', 'NorthEast' );
% Label axes
ylabel('Pressure Magnitude at Far Plane Edge at Initial Reflection');...
    xlabel('Total Number of Vacancy Points'); title('Pressure Magnitude on Far Edge at Wave Front vs Vacancies');
grid on

