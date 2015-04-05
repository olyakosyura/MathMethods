function  U = calcU()

    function W = W(t)
        t0 = 0.5;
        Q = 10;
        if (0<=t && t < 0.5*t0)
            W = 2*Q*t;
        else
            if (t >=t0)
                W = 0;
            else
                W=2*Q*(t0-t);
            end
        end
    end

    function K = K(u)
        a = 0.5;
        b = 2;
        sigma = 1.25;
    
        K = a + b * u ^ sigma;
    end

    l = 1;
    den = 1;
    c = 2;
    u0 = 0.2;
    T = 1.0;
    
    M = 50;
    N = 50;
    I = 10; 
    
    h = l / N;
    tau = T / M;
    
    
    x = 0 : h : l;
    t = 0 : tau : T;

    U = zeros(M + 1, N + 1);

    U(1, :) = u0;
    
    coef = tau / (den*c*h^2);
    for i = 2 : M + 1
        U(i, :) = U(i - 1, :);
    
        for k = 1 : I
            A = zeros(N + 1, N + 1);
            F = zeros(N + 1, 1);
    
            A(1, 1) = -1;
            A(1, 2) = 1;
            Kp = ( K(U(i, 1)) + K(U(i, 2)) ) / 2;
            F(1) = 0;
            for j = 2 : N
                Km = ( K(U(i, j - 1)) + K(U(i, j)) ) / 2;
                Kp = ( K(U(i, j)) + K(U(i, j + 1)) ) / 2;
                
                A(j, j - 1) = -Km * coef;
                A(j, j)     = ((Km + Kp) * coef + 1);
                A(j, j + 1) = -Kp * coef;
            
                F(j) = U(i - 1, j);
            end

            A(N+1, N) = -1;
            A(N+1, N+1) = 1;
            F(N+1) = W(t(i)) * h / Kp;           

            U(i, :) = A \ F;
            U(N+1)=u0;
        end
    end
    
    [max_t, max_j] = max(U(:,ceil(N/2)));    
    fprintf( 'Max temperature: %6.6f\n', max_t);
    fprintf( '           Time: %6.6f\n', (max_j-1)*tau);
    surf(x, t, U);
    ylabel('x', 'FontSize', 14);
    xlabel('t', 'FontSize', 14);  
    
end
