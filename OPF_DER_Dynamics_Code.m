%% ===== TWO-BUS SCOPF + DER TRANSIENT (hard V2 band, secure final PF) ====
clc; clear; close all; yalmip('clear')

%% -------------------- USER PARAMETERS -----------------------------------
DER_amp_limit = 350;                 % DER current limit [A]
Qg_cap        = 2;                   % |Qg| limit [pu] for the generator
Ts = 0.01; Tf = 0.20;                % sampling & horizon
t_fault = 0.10; t_ctrl = 0.11;       % fault start / control start [s]
V2_low  = 0.93;   V2_high = 0.955;    % recovery band (HARD)
V2_fault_floor = 0;               % voltage floor during fault
RELAX_DER_DYNAMICS = true;          % set true to use ">=" instead of "=="

% Simple operating cost
wPg = 50;  wQg = 10;  wI = 10;

%% -------------------- Network & Base quantities -------------------------
baseMVA = 100; Vbase_kV = 230;
r = 0.01;  x = 0.10;  b = 0.076;
y = 1/(r+1j*x); g = real(y); bb = imag(y); bsh = b/2;

Ibase   = (baseMVA*1e6)/(sqrt(3)*Vbase_kV*1e3); % [A]
Ipu_max = DER_amp_limit / Ibase;                % pu

%% -------------------- Pre-fault loads (nominal) -------------------------
Pd1 = 100/baseMVA; Qd1 = 0;
Pd2 = 200/baseMVA; Qd2 = 50/baseMVA;

%% -------------------- 1) Steady-state OPF (no DER) ---------------------
V1 = sdpvar(1,1); V2 = sdpvar(1,1); d2 = sdpvar(1,1);   % δ1 = 0 (slack)
Pg0 = sdpvar(1,1); Qg0 = sdpvar(1,1);

c12 = cos(-d2); s12 = sin(-d2);
P1 =  V1^2*g        + V1*V2*(-g*c12 - bb*s12);
Q1 = -V1^2*(bb+bsh) + V1*V2*(-g*s12 + bb*c12);
P2 =  V2^2*g        + V1*V2*(-g*c12 + bb*s12);
Q2 = -V2^2*(bb+bsh) + V1*V2*( g*s12 + bb*c12);

conOPF = [ Pg0 - Pd1 == P1,  Qg0 - Qd1 == Q1, ...
           -Pd2 == P2,       -Qd2 == Q2, ...
           0.95 <= V1 <= 1.05, 0.95 <= V2 <= 1.05, ...
           0 <= Pg0 <= 5,    -2 <= Qg0 <= 2, ...
           -pi <= d2 <= pi ];

ip = sdpsettings('solver','ipopt','verbose',0,'usex0',1);
assign([V1 V2 d2 Pg0 Qg0],[1 1 0 3 0]);
optimize(conOPF, 50*Pg0 + 10*Qg0, ip);

V1_ss = value(V1); V2_ss = value(V2);
Pg_ss = value(Pg0); Qg_ss = value(Qg0); d2_ss = value(d2);


%% -------------------- Time grid & fault profile ------------------------
N = max(2, round(Tf/Ts)); time = (1:N)*Ts;
k_fault = min(N, max(1, round(t_fault/Ts)));
k_ctrl  = min(N, max(1, round(t_ctrl /Ts)));

Qd2_fault = Qd2*ones(1,N);
Qd2_fault(k_fault:N) = 1.5*Qd2;      % set =Qd2 for "no fault" test

%% -------------------- DER dynamics (hard equality) ----------------------
% I[k] = a_d * I[k-1] + b_d * V2[k] + c_d
Bdyn = 1e-6;  Acoeff = 1.533759;  Ccoeff = -0.130297;
denom = 1 - Ts/Bdyn;
a_d = 1/denom;
b_d = -(Ts/Bdyn)*Acoeff/denom;
c_d = -(Ts/Bdyn)*Ccoeff/denom;

%% -------------------- Horizon decision variables ------------------------
Vh   = sdpvar(2,N);        % [V1; V2]
angh = sdpvar(2,N);        % [δ1; δ2]
Pgk  = sdpvar(1,N); Qgk = sdpvar(1,N);
Ipu  = sdpvar(1,N);        % DER reactive current (pu)

% global bounds
con = [ 0 <= Pgk <= 5,  -Qg_cap <= Qgk <= Qg_cap, ...
        0.95 <= Vh(1,:) <= 1.05, ...
        0 <= Ipu <= Ipu_max, ...
        -pi <= angh(:) <= pi ];

% initial
con = [ con, ...
    Vh(:,1)==[V1_ss;V2_ss],  angh(:,1)==[0; d2_ss], ...
    Pgk(1)==Pg_ss,           Qgk(1)==Qg_ss, ...
    Ipu(1)==0 ];

% power balance helper
pb = @(k) deal( ...
    Vh(1,k)^2*g + Vh(1,k)*Vh(2,k)*(-g*cos(angh(1,k)-angh(2,k)) - bb*sin(angh(1,k)-angh(2,k))), ...
   -Vh(1,k)^2*(bb+bsh) + Vh(1,k)*Vh(2,k)*(-g*sin(angh(1,k)-angh(2,k)) + bb*cos(angh(1,k)-angh(2,k))), ...
    Vh(2,k)^2*g + Vh(1,k)*Vh(2,k)*(-g*cos(angh(1,k)-angh(2,k)) + bb*sin(angh(1,k)-angh(2,k))), ...
   -Vh(2,k)^2*(bb+bsh) + Vh(1,k)*Vh(2,k)*( g*sin(angh(1,k)-angh(2,k)) + bb*cos(angh(1,k)-angh(2,k))) );

%% -------------------- 7a) Pre-fault segment -----------------------------
if k_fault > 1
  for k = 1:(k_fault-1)
    con = [ con, ...
        Vh(:,k)==[V1_ss;V2_ss],  angh(:,k)==[0; d2_ss], ...
        Pgk(k)==Pg_ss,           Qgk(k)==Qg_ss,  Ipu(k)==0 ];
  end
end

%% -------------------- 7b) Fault segment (DER off, gen fixed) -----------
k_fault_to = max(k_fault, k_ctrl-1);
for k = k_fault:k_fault_to
    [P1k,Q1k,P2k,Q2k] = pb(k);
    con = [ con, ...
        Vh(1,k)==V1_ss, angh(1,k)==0, Pgk(k)==Pg_ss, Qgk(k)==Qg_ss, Ipu(k)==0, ...
        V2_fault_floor <= Vh(2,k) <= 1.05, ...
        Pgk(k)-Pd1==P1k,  Qgk(k)-Qd1==Q1k,  -Pd2==P2k,  -Qd2_fault(k)==Q2k ];
end

%% -------------------- 7c) Recovery segment (hard band) ------------------
obj = 0;
for k = k_ctrl:N-1
    [P1k,Q1k,P2k,Q2k] = pb(k);
    Qder = Vh(2,k)*Ipu(k);

    % balances + reference angle + hard V band
    con = [ con, ...
        Pgk(k)-Pd1==P1k,  Qgk(k)-Qd1==Q1k,  -Pd2==P2k, ...
        -Qd2_fault(k) + Qder == Q2k, ...
        %angh(1,k)==0, ...
        V2_low <= Vh(2,k) <= V2_high ];

    % DER dynamics
    if RELAX_DER_DYNAMICS
        con = [ con, Ipu(k) >= a_d*Ipu(k-1) + b_d*Vh(2,k) + c_d ];
    else
        con = [ con, Ipu(k) == a_d*Ipu(k-1) + b_d*Vh(2,k) + c_d ];
    end

    % accumulate operating cost
    obj = obj + wPg*(Pgk(k)) + wQg*(Qgk(k)) + wI*(Ipu(k));
end

% Final setpoints must not drop below steady-state (security)
%con = [ con, Pgk(N) >= Pg_ss, Qgk(N) >= Qg_ss ];
% No end drop
con = [ con, Vh(2,N) == Vh(2,N-1), Ipu(N) == Ipu(N-1) ];

%% ---- Base-case feasibility block (no DER, nominal loads) ---------------
V1_bc = sdpvar; V2_bc = sdpvar; d2_bc = sdpvar;
c12b = cos(-d2_bc); s12b = sin(-d2_bc);
P1b =  V1_bc^2*g        + V1_bc*V2_bc*(-g*c12b - bb*s12b);
Q1b = -V1_bc^2*(bb+bsh) + V1_bc*V2_bc*(-g*s12b + bb*c12b);
P2b =  V2_bc^2*g        + V1_bc*V2_bc*(-g*c12b + bb*s12b);
Q2b = -V2_bc^2*(bb+bsh) + V1_bc*V2_bc*( g*s12b + bb*c12b);

con = [ con, ...
    0.95 <= V1_bc <= 1.05,  0.95 <= V2_bc <= 1.05,  -pi <= d2_bc <= pi, ...
    Pgk(N) - Pd1 == P1b, ...
    Qgk(N) - Qd1 == Q1b, ...
    -Pd2 == P2b, ...
    -Qd2 == Q2b ];     % all four equalities enforced -> true base-case feasibility

%% -------------------- Solve SCOPF --------------------------------------
assign(Vh,   repmat([V1_ss;V2_ss],1,N));
assign(angh, repmat([0; d2_ss],1,N));
assign(Pgk,  Pg_ss*ones(1,N));
assign(Qgk,  Qg_ss*ones(1,N));
assign(Ipu,  zeros(1,N));
assign([V1_bc V2_bc d2_bc],[V1_ss V2_ss d2_ss]); % good initial guess

ip = sdpsettings(ip,'verbose',1,'usex0',1);
ip.ipopt.max_iter = 5000; ip.ipopt.tol = 1e-6; ip.ipopt.constr_viol_tol = 1e-5;
sol = optimize(con, obj, ip);
if sol.problem, warning('Solver status: %s', sol.info); end

Vh  = value(Vh);  angh = value(angh);
Pgk = value(Pgk); Qgk = value(Qgk); Ipu = value(Ipu);
V1_bc = value(V1_bc); V2_bc = value(V2_bc); d2_bc = value(d2_bc);
IderA = Ipu*Ibase;

%% -------------------- PRINTED RESULTS -----------------------------------
% fault sag (pre-control)
if k_ctrl > k_fault
  [Vm_fault, idx_fault] = min(Vh(2,k_fault:k_ctrl-1));
  t_fault_min = time(k_fault + idx_fault - 1);
else
  Vm_fault = Vh(2,k_fault); t_fault_min = time(k_fault);
end
[Vm_overall, idx_overall] = min(Vh(2,1:N));

fprintf('\n--- Initial Steady OPF (No Disturbance) -----------------------------------\n')
fprintf('Pg %.4f  Qg %.4f   V1 %.4f   V2 %.4f   δ2 %.2f°\n', ...
        Pg_ss, Qg_ss, V1_ss, V2_ss, d2_ss*180/pi)

fprintf('\n---Fault Voltage Analysis ----------------------------------\n')
fprintf('Voltage Drop During Fault: Vmin %.4f pu at t = %.3f s (before control)\n', Vm_fault, t_fault_min)
fprintf('Overall Minimum Voltage : Vmin %.4f pu at t = %.3f s (entire period)\n', Vm_overall, time(idx_overall))
fprintf('Final (recovery) : V2 %.4f  I_DER %.1f A\n', Vh(2,N), IderA(N))

fprintf('\n--- Final Generator Setpoints (for base-case) --------------\n');
fprintf('Pg_final = %.4f pu,  Qg_final = %.4f pu\n', Pgk(N), Qgk(N));

%% -------- Base-case PF (no DER, nominal loads) USING base-case vars ----
c12_pf = cos(-d2_bc); s12_pf = sin(-d2_bc);
P1_pf =  V1_bc^2*g        + V1_bc*V2_bc*(-g*c12_pf - bb*s12_pf);
Q1_pf = -V1_bc^2*(bb+bsh) + V1_bc*V2_bc*(-g*s12_pf + bb*c12_pf);
P2_pf =  V2_bc^2*g        + V1_bc*V2_bc*(-g*c12_pf + bb*s12_pf);
Q2_pf = -V2_bc^2*(bb+bsh) + V1_bc*V2_bc*( g*s12_pf + bb*c12_pf);

misP1 = (Pgk(N) - Pd1) - P1_pf;
misQ1 = (Qgk(N) - Qd1) - Q1_pf;
misP2 = (-Pd2)        - P2_pf;
misQ2 = (-Qd2)        - Q2_pf;

fprintf('\n--- Power Flow Verification (No DER, No Disturbance) -------\n');
fprintf('Pg = %.4f pu,  Qg = %.4f pu\n', Pgk(N), Qgk(N));
fprintf('V1 = %.4f pu,  V2 = %.4f pu,  d2 = %.2f°\n', V1_bc, V2_bc, d2_bc*180/pi);
fprintf('P1 = %.4f pu,  Q1 = %.4f pu,  P2 = %.4f pu,  Q2 = %.4f pu\n', P1_pf, Q1_pf, P2_pf, Q2_pf);
fprintf('PF mismatches: dP1=%.2e, dQ1=%.2e, dP2=%.2e, dQ2=%.2e\n', misP1, misQ1, misP2, misQ2);

%% -------------------- Plots --------------------------------------------
figure('Color','w')
subplot(2,1,1)
plot(time, IderA, 'LineWidth',2); grid on
ylabel('I_{DER} (A)'); title(sprintf('DER Reactive Current  (limit = %.0f A)',DER_amp_limit))
xline(time(k_fault),'r:','Fault Start','LineWidth',1.2)
xlim([time(1) time(end)])

subplot(2,1,2)
plot(time, Vh(2,:), 'r', 'LineWidth',2); grid on; hold on
yline([0.95 1.05],'k--','LineWidth',1.0);
ylim([0.9 1.02]);
%yline([V2_low V2_high],'k-.','LineWidth',1.0);
ylabel('V_{bus2} (pu)',FontSize=12); xlabel('Time (s)', FontSize=12)
title('Bus-2 Voltage: Fault sag and recovery ')
xline(time(k_fault),'r:','Fault Start','LineWidth',1.2)
%legend('V_2','Std limits','Recovery band','Location','best')
xlim([time(1) time(end)])
