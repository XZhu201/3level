%% 
D=importdata('Et.dat');

tE = D(:,1);
Ex = D(:,2);
Ey = D(:,3);

figure; plot3(tE,Ex,Ey)


%%
D=importdata('odeResults.dat');

t = D(:,1);
a1 = D(:,2)+i*D(:,3);
a2 = D(:,4)+i*D(:,5);
a3 = D(:,6)+i*D(:,7);

figure; plot(t,abs(a1),t,abs(a2),t,abs(a3))

%% analyze the phase
% alignment of the surperposed p state

phi1 = angle(a1);
phi2 = angle(a2);

% phi1 = angle(a1.* exp(-i*Energy1*t)) ;
% phi2 = angle(a2.* exp(-i*Energy2*t)) ;

Dphi = phi2-phi1;

figure; 
subplot(121)
plot(t, (phi1)/pi, t, (phi2)/pi, t, phi1/pi-phi2/pi);
xlabel('time /T')
ylabel('phase /\pi')
% title('phases of coef.')
legend('phase1','phase2','phase1-phase2')

subplot(122)
plot(t, unwrap(phi1)/pi, t, unwrap(phi2)/pi, t, unwrap(phi1)/pi-unwrap(phi2)/pi );
xlabel('time /T')
ylabel('phase /\pi')
% title('phases of coef.')
legend('phase1','phase2','phase1-phase2')


%% analyze the angle
alignment = Dphi/2;
alignment = alignment - pi.*(alignment>pi/2) + pi.*(alignment<-pi/2);     % to [-pi/2, pi/2]

figure;
plot(t,alignment/pi,'.r'); 

% angle of electric field
% thetaE = atan(Ey./Ex) + pi.*( Ex<0 & Ey>0 ) - pi.*( Ex<0 & Ey<0 );
thetaE = atan(Ey./Ex);                  % to [-pi/2, pi/2]
absE = sqrt(Ex.^2+Ey.^2);

hold on;
plot(tE,thetaE/pi,'.');
xlabel('time /T')
ylabel('angle /\pi')
ylim([-0.5 0.5])

yyaxis right;
plot(tE,absE);
ylabel('Electric field /a.u.')

saveas(gcf,'angles','jpg')
saveas(gcf,'angles','fig')

%% population
figure; 
plot(t,abs(a1).^2,t,abs(a2).^2,t,abs(a3).^2);
legend('up','um','us')

saveas(gcf,'population','jpg')
saveas(gcf,'population','fig')



disp('Job done!')
