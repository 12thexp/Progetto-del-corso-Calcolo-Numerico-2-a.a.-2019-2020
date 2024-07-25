clear all;

format short e
% N = 100.*2.^[1:6];
n = 200;
X = load(string(n));
i = 1;
C = 2;
p = 4;

while (n <= 51200)
    Y = load(string(C*n));
    err(i,1) = n;
    
    diff1(1:4) = abs(X(2,2:5)-Y(C+1,2:5));  %inizializzo diff1, k=2
    N1 = norm(diff1);
    
    k = 3;
    while (k <= n+1)
        diff2(1:4) = abs((X(k,2:5)-Y(C*(k-1)+1,2:5)));
        N2 = norm(diff2);
        N1 = max([N1 N2]);
        k = k+1;
    end
    
    err(i,2) = (C^p/(C^p-1))*N1;
    
    n = C*n;
    X = Y;
    i = i+1;
    
end

err
