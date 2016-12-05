Fwlc = @(x,L,p) 4.11/p*(1/4*(1-x/L).^-2-1/4+x/L);


%kbt = 4.11 %pN*nm;

PxScale = 0.110906330000000;
L = .34/1000*6000; %um
p = 50;%nm
%F = kbt*L/varX

close all;

figure();

x = linspace(0.1,0.99*L);
F = Fwlc(x,L,p);

plot(x,F);
set(gca,'yscale','log');
xlabel(sprintf('Extension (Lo=%0.2f\\mum)',L));
ylabel('Force [pN]')
title('Expected F v L');

varX = 4.11*x*1000./F;
figure;
plot(x,sqrt(varX)/(PxScale*1000));
xlabel('Extension (µm)')
ylabel('sqrt(var(x)) [px^2]');

