
function []=testgrad(dat)

graddat=dat(2:11,:)-dat(12:21,:);
jacobianCentre=(graddat(:,3)*ones(1,10))./graddat(:,19:28);
jacobianFwd=((dat(2:11,3)-dat(1,3))*ones(1,10))./graddat(:,19:28);
jacobianBkw=((dat(1,3)-dat(12:21,3))*ones(1,10))./graddat(:,19:28);

figure
subplot(1,3,1)
surf(log10(abs(jacobianCentre)))
zlabel('absolute of jacobian - central dif')

subplot(2,3,2)
surf(log10(abs(jacobianFwd)))
zlabel('absolute of jacobian - Fwd dif')
subplot(2,3,5)
surf(log10(abs(jacobianBkw)))
zlabel('absolute of jacobian - Bkw dif')


subplot(2,3,3)
surf(((graddat(:,19:28))))
zlabel(' desVarChange')

subplot(2,3,6)

l=plot(diag(jacobianCentre));
l(end+1)=plot(diag(jacobianFwd));
l(end+1)=plot(diag(jacobianBkw));
hold on,
l(end+1)=plot(graddat(:,3));
l(end+1)=plot(diag(graddat(:,19:28)));