

clear all
close all

load BSfludat2.txt % loading data

format long

tt = BSfludat2(:,1); % define array with t-coordinates of the data

yy = BSfludat2(:,2); % define array with y-coordinates of the data

% cumulative case
yc=cumsum(yy);
%tt1=tt/tfr;

figure(1);
plot(tt, yy, tt, yc )
title('Flu borarding school incidence and cum. cases')
legend('raw incidence', 'cum. cases');
xlabel('days')
ylabel('cases')

% guessing the initial values
nm=max(yy)
ind=find(yy==nm, 1, 'first');
tm=tt(ind)
hguess=1.0;
sircoeffs=[1.0;hguess;tm;nm/hguess]

% fit by i
[c2, R2] = isapfit(tt, yy, sircoeffs)

% for plotting
tt1=linspace(0,20,100);

figure(2)
% double-exponentail model
zap2=c2(4)*sirapprox3(tt1-c2(3), c2(1), c2(2));
plot(tt, yy, tt1, zap2)
legend('raw incidence', 'asymp');
xlabel('days')
ylabel('cases')

% fit by r
[c3, R3] = rsapfit(tt, yc, sircoeffs)

figure(3)
% double-exponentail model
zap3=c3(4)*rapprox(tt1'-c3(3), c3(1), c3(2))';
plot(tt, yy, tt, yc, tt1, zap2, tt1, zap3)
legend('raw incidence', 'cum. cases', 'asymp i', 'asymp r');
xlabel('days')
ylabel('cases')

figure(4)
zap21=c2(4)*sirapprox3(tt1-c3(3), c3(1), c3(2));
plot(tt, yy, tt1, zap2,  tt1, zap21)
xlabel('days')
legend('raw incidence', 'asymp i', 'asymp r');
ylabel('cases')