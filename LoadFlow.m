% Program for Newton-Raphson Load Flow Analysis..
% nbs=input('No. of buses : ');
% nmc=input('No. of Machines : ');
% bus_dat=input('Bus Data Source : ');
% line_dat=input('Line Data Source : ');

% Extracting data from bus_dat
bus_dat = sortrows(bus_dat,2);  % Input bus data is sorted so that Slack bus appears at the last row- COnvention
bus = bus_dat(:,1);            % Bus Number
type = bus_dat(:,2);           % 103-Slack, 102-PV, 101-PQ
V = bus_dat(:,3);              % Voltage magnitude
del = bus_dat(:,4);            % Voltage Angle
Pg = bus_dat(:,5);              
Qg = bus_dat(:,6);             
Pl = bus_dat(:,7);            
Ql = bus_dat(:,8);          

% Extracting data from line_dat

fb = line_dat(:,1);             % From bus number
tb = line_dat(:,2);             % To bus number
r = line_dat(:,3);              
x = line_dat(:,4);              
b = line_dat(:,5);              % Line Charging
a = line_dat(:,6);              % Tap setting value

b = b./2;                       % Half Line Charging
z = r + i*x;                    % Impedance 
y = 1./z;                       % Admittance
b = i*b;                        % Making B imaginary
nl = length(fb);                % No. of branches
Y = zeros(nbs,nbs);             % Initialise YBus

P = Pg - Pl;                
Q = Qg - Ql;                
Psp = P;                    
Qsp = Q;                    
pq = find(type == 101);             % PQ Buses
npq = length(pq);                   % No. of PQ buses
pv = find(type == 102);             % PQ Buses
npv = length(pv);                   % No. of PQ buses

% Formation of Ybus....
% Off- Diagonal Elements
for k = 1:nl
     Y(fb(k),tb(k)) = Y(fb(k),tb(k)) - y(k)/a(k);
     Y(tb(k),fb(k)) = Y(fb(k),tb(k));
 end
 
 % Diagonal Elements....
 for m = 1:nbs
     for n = 1:nl
         if fb(n) == m
             Y(m,m) = Y(m,m) + y(n)/(a(n)^2) + b(n);
         elseif tb(n) == m
             Y(m,m) = Y(m,m) + y(n) + b(n);
         end
     end
 end

%Splitting ybus
G = real(Y);                % Conductance matrix..
B = imag(Y);                % Susceptance matrix..

% Iterations starting
Tolerance = 1;  
Iterations = 0;
while (Tolerance>0.000001)   
    Iterations =Iterations+1;
    P = zeros(nbs,1);
    Q = zeros(nbs,1);
    % Calculate P and Q
    for i = 1:nbs
        for k = 1:nbs
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end
    
    % Calculate change from specified value
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbs
        if type(i) == 101
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(1:nbs-1);
    M = [dP; dQ];       % Mismatch Vector
    
    % Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbs-1,nbs-1);
    for i = 1:(nbs-1)
        m = i;
        for k = 1:(nbs-1)
            n = k;
            if n == m
                for n = 1:nbs
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
        
    % J2 - Derivative of Real Power Injections with V..
    J2 = zeros(nbs-1,npq);
    for i = 1:(nbs-1)
        m = i;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbs
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(npq,nbs-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbs-1)
            n = k;
            if n == m
                for n = 1:nbs
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbs
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    J = [J1 J2; J3 J4];    % Jacobian Matrix..
    
    %     Finding X by J^-1 * M (inverse using Gausian Elimination)
%     
%     [m1,n1] = size(J);
%     A = [J, eye([m1, n1])];
%     for k1 = 1:m1
%       A(k1,k1:end) = A(k1,k1:end)/A(k1,k1);
%       A([1:k1-1,k1+1:end],k1:end) = A([1:k1-1,k1+1:end],k1:end) - A([1:k1-1,k1+1:end],k1)*A(k1,k1:end);
%     end
%     invers = A(:,n1+1:end);
%     
%     X = invers*M ;   % Correction Vector
       
    % Finding X using Gaussian elimination
  N = 2*npq + npv;
  for ld=1:N-1
    for row=ld+1:N 
       multiplier = J(row,ld)/J(ld,ld);
       for col=ld:N 
	         J(row,col) =J(row,col)- J(ld,col)*multiplier;
       end
       M(row,1) =M(row,1)- M(ld,1) * multiplier;
    end
  end
  X=zeros(N,1);% Correction Vector
  
  for row=N:-1:1 
    X(row,1) = M(row,1);
    for col=N:-1:row+1 
      X(row,1) = X(row,1)-J(row,col) * X(col,1);
    end
   
    X(row,1) =X(row,1)/ J(row,row);
  end
 % Correction Vector Completed
    
    dTh = X(1:nbs-1);      % Change in Voltage Angle..
    dV = X(nbs:end);       % Change in Voltage Magnitude..
    
    % Updating State Vectors..
    del(1:nbs-1) = dTh + del(1:nbs-1);    % Voltage Angle..
    k = 1;
    for i = 1:nbs-1
        if type(i) == 101
            V(i) = dV(k) + V(i);        % Voltage Magnitude..
            k = k+1;
        end
    end
   
    Tolerance = max(abs(M));                  % Tolerance..
    
end

Vm = V.*cos(del)+j*V.*sin(del); % Converting polar to rectangular..
Del = del*180/pi;               % Bus Voltage Angles in Degree...
Pl = bus_dat(:,7);                 % PLi..
Ql = bus_dat(:,8);                 % QLi..

Iij = zeros(nbs,nbs);
Sij = zeros(nbs,nbs);
Si = zeros(nbs,1);

% Bus Current Injections..
 I = Y*Vm;
 Im = abs(I);
 Ia = angle(I);
 
%Line Current Flows..
for m = 1:nl
    p = fb(m); q = tb(m);
    Iij(p,q) = -(Vm(p) - Vm(q))*Y(p,q); % Y(m,n) = -y(m,n)..
    Iij(q,p) = -Iij(p,q);
end
Iijm = abs(Iij);
Iija = angle(Iij);

% Line Power Flows..
for m = 1:nbs
    for n = 1:nbs
        if m ~= n
            Sij(m,n) = Vm(m)*conj(Iij(m,n));
        end
    end
end
Pij = real(Sij);
Qij = imag(Sij);
 
% Bus Power Injections..
for i = 1:nbs
    for k = 1:nbs
        Si(i) = Si(i) + conj(Vm(i))* Vm(k)*Y(i,k);
    end
end
Pi = real(Si);
Qi = -imag(Si);
Pg = Pi+Pl;
Qg = Qi+Ql;
 
disp('                                    Load Flow Analysis Using Newton Raphson Method ');
fprintf('\n');
fprintf('                                           Iterations for convergence = %d',Iterations);
fprintf('\n');
fprintf('\n');

disp(' 1. Bus Voltages and Power Flow through the buses');
fprintf('\n');

disp('|  Bus  |   |V|  |  Angle  | Power Injection (p.u.) |');
disp('|  No   |   pu   |  Degree |  Active   |  Reactive  |');
for m = 1:nbs
    fprintf('  %3g', m); fprintf('  %8.4f ', V(m)); fprintf('   %8.4f', Del(m));
    fprintf(' %8.3f       ', Pi(m)); fprintf('%6.3f', Qi(m));  
     fprintf('\n');
end
fprintf('\n');
fprintf('\n');

disp(' 2. Reactive power generation at the P-V buses  ');
fprintf('\n');
disp('Bus No.        Reactive Power Generation (p.u.)');
for m = 1:nbs
    if type(m)==102
        k=Qg(m);
        fprintf('%d                 %f\n ',bus(m), k );
    end
end
  
  fprintf('\n');fprintf('\n');  

disp(' 3. Power Flow Through the Lines  ');
fprintf('\n');
disp('|  From |  To   |  Power Flow   |  From  |  To   |   Power Flow |');
disp('|  Bus  |  Bus  |      pu       |  Bus   |  Bus  |      pu      |');
for m = 1:nl
    p = fb(m); q = tb(m);
    fprintf('  %4g  ', p); fprintf(' %4g  ', q);fprintf('%7.3f', Pij(p,q)); fprintf(' +j'); fprintf('%7.3f', Qij(p,q)); 
    fprintf('  %4g', q); fprintf('  %4g  ', p); fprintf(' %7.3f', Pij(q,p));fprintf(' +j'); fprintf('%7.3f', Qij(q,p));
    
    fprintf('\n');
end

