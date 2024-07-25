clear all;

%stima dell'errore

n = 400;    %n iniziale
metodo = 'rk4_';
p = 4;
X = load(strcat(metodo,string(n)));
C = 2;  %costante moltiplicativa dei passi

i = 1; %indice di riga della matrice err
while (n <= 51200)
    Y = load(strcat(metodo,string(C*n)));
    err(i,1) = n;
    
    for (k = 1:n+1)
        diff(1:4) = abs((X(k,2:5)-Y(C*(k-1)+1,2:5)));   %C*(k-1)+1 e' l'indice di riga nel file '2*n' corrispondente all'indice k nel file 'n'
        N(k) = norm(diff);
    end
    
    err(i,2) = (C^p/(C^p-1))*max(diff(1));
    err(i,3) = (C^p/(C^p-1))*max(diff(2));
    err(i,4) = (C^p/(C^p-1))*max(norm(diff));
    i = i+1;
    n = C*n;
    X = Y;
end

fprintf('\npassi \t errore1 \t errore2 \t norma err \n\n');   %stampo la matrice err
[m,n] = size(err);
for (i = 1:m)
    fprintf('%d\t %1.4e\t %1.4e\t %1.4e\t \n',err(i,1),err(i,2), err(i,3), err(i,4));
end


% %grafici

s = '51200';
X = load(strcat(metodo,s));

x = X(:,2);
y = X(:,3);

N = length(X(:,1));
z = linspace(0,250,N);

figure();     %grafico componenti
plot(z,x,'b','LineWidth',1.35);
legend({'φ1(t)'},'Location','southeast');
title(strcat(s,' passi: Componente φ1(t)'));

figure();
plot(z,y,'r','LineWidth',1.35);
legend({'φ2(t)'},'Location','southeast');
title(strcat(s,' passi: Componente φ2(t)'));


figure();   %grafico soluzione
plot(x,y,'LineWidth',1.35);
hold on
ylim([-1.5 1.5]);
xlim([-1.5 1.5]);
xL = xlim;
yL = ylim;
line([0 0], yL, 'Color','black');  %y-axis
line(xL, [0 0], 'Color','black');  %x-axis
title(strcat(s,' passi: grafico soluzione'));


