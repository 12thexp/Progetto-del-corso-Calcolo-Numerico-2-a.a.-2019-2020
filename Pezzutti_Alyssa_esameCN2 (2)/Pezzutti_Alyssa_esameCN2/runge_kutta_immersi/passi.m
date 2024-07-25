function passi(T)
% visualizzazione ampiezze passi
% T vettore dei tempi, contiene anche t iniziale
disp(['numero passi = ',num2str(length(T)-1)])
disp(['passo max = ',num2str(max(diff(T)))])
disp(['passo min = ',num2str(min(diff(T)))])
figure()
plot(T(2:end),diff(T),'o',T(2:end),diff(T),':')
title('ampiezza dei passi')
