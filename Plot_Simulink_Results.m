% 1) Load the results struct
S = load('Ouput_Data_Volt_Var_New.mat');
R = S.results;  % adjust if your nesting is different: e.g. S.results

% 2) Pull out each variable
V = R.Voltage_at_PCC;            % voltages at PCC (3‑phase)
Q_rms = R.I_RMS;    % RMS reactive injection (3‑phase)
Q_sine = R.I_sine;  % Sine reactive injection (3‑phase)

% 3) Time base 
if isa(V,'timeseries')
  tV    = V.Time;
  dataV = V.Data;
  tQ    = Q_rms.Time;
  dataQ = Q_rms.Data;
  tS    = Q_sine.Time;
  dataS = Q_sine.Data;
else
  % otherwise assume simple numeric arrays of size [N×3]
  dataV = V;  
  dataQ = Q_rms;
  dataS = Q_sine;
  tV    = (1:size(dataV,1)).';   
  tQ    = (1:size(dataQ,1)).';
  tS    = (1:size(dataS,1)).';
end

% 4) Define deeper, clearer colors
colV = [0.00 0.10 0.40;   % deep blue
        0.00 0.40 0.00;   % dark green
        0.60 0.00 0.00];  % dark red
colQ = [0.00 0.20 0.60;
        0.00 0.60 0.00;
        0.80 0.00 0.00];
colS = [0.00 0.20 0.80;
        0.00 0.40 0.00;
        0.80 0.00 0.00];

% 5) Plot
figure('Color','w'); hold on; grid on

% Voltage at PCC
subplot(3,1,1)
hV = plot(tV, dataV, 'LineWidth', 2);
hold on
% add horizontal limit line at y = 0.95
yline(0.95, 'k--', 'LineWidth', 1);
hold off
for i = 1:3
    hV(i).Color = colV(i,:);
end
title('Voltage at PCC', FontSize=16)
ylabel('Voltage (pu)', FontSize=18)
legend('Phase A','Phase B','Phase C','Location','northeast', Fontsize=13)
ax=gca;
ax.FontSize=18;
grid on

% Reactive Injection (RMS)
subplot(3,1,2)
hQ = plot(tQ, dataQ, 'LineWidth', 2);
for i = 1:3
    hQ(i).Color = colQ(i,:);
end
title('Reactive Injection (RMS)', FontSize=16)
ylabel('Current (A)', FontSize=18)
legend('Phase A','Phase B','Phase C','Location','northeast', FontSize=13)
ax=gca;
ax.FontSize=18;
grid on

% AC Reactive Injection (Sine)
subplot(3,1,3)
hS = plot(tS, dataS, 'LineWidth', 2);
for i = 1:3
    hS(i).Color = colS(i,:);
end
title('AC Reactive Injection (Sine)', FontSize=16)
xlabel('Time (s)', FontSize=18)
ylabel('Current (A)', FontSize=18)
% move legend to top-right
legend('Phase A','Phase B','Phase C','Location','northeast', FontSize=13)
ax=gca;
ax.FontSize=18;
grid on
box on