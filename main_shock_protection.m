
clear; clc; close all;
if exist('results.txt', 'file'), delete('results.txt'); end
diary('results.txt');

makeFigures = true;     

%% ------------------------------------------------------------------------
%  Baseline parameters (Section 5)
%  ------------------------------------------------------------------------
lambda = 12;                     % surges per year
alpha  = 1/15;                   % protection failure rate (mean life 15 years)
qP = [0.85, 0.13, 0.02];         % (q10, q11, q12): protected regime
qU = [0.50, 0.35, 0.15];         % (q20, q21, q22): unprotected regime
b  = 4;                          % Model 1
b1 = 4; b2 = 2;                  % Model 2

[Q1m1, u1m1, am1] = model1_matrices(qP, b);
[Q2m1, u2m1]      = model1_matrices(qU, b);
[Q1m2, u1m2, am2] = model2_matrices(qP, b1, b2);
[Q2m2, u2m2]      = model2_matrices(qU, b1, b2);

allOK = true;
fprintf('Mixed shock models with a protection block - numerical results\n');
fprintf('Baseline: lambda = %g, alpha = 1/15, q1 = (%.2f, %.2f, %.2f), q2 = (%.2f, %.2f, %.2f)\n', ...
        lambda, qP, qU);
fprintf('          Model 1: b = %d;  Model 2: (b1, b2) = (%d, %d)\n\n', b, b1, b2);

%% ------------------------------------------------------------------------
%  Tables 1 and 2: survival function for varying alpha (Section 5.1)
%  ------------------------------------------------------------------------
tGrid  = [0.5 1 2 3 5 10 15]';
alphas = [0.020 0.067 0.200 0.500];

T1 = zeros(numel(tGrid), numel(alphas));
T2 = zeros(numel(tGrid), numel(alphas));
for k = 1:numel(alphas)
    T1(:, k) = survival_protection(Q1m1, Q2m1, am1, lambda, alphas(k), tGrid);
    T2(:, k) = survival_protection(Q1m2, Q2m2, am2, lambda, alphas(k), tGrid);
end

fprintf('Table 1. Model 1: P(T > t), b = 4\n');
fprintf('   t      alpha=0.020  alpha=0.067  alpha=0.200  alpha=0.500\n');
for i = 1:numel(tGrid)
    fprintf('%5.1f    %10.4f   %10.4f   %10.4f   %10.4f\n', tGrid(i), T1(i, :));
end
fprintf('\nTable 2. Model 2: P(T > t), (b1, b2) = (4, 2)\n');
fprintf('   t      alpha=0.020  alpha=0.067  alpha=0.200  alpha=0.500\n');
for i = 1:numel(tGrid)
    fprintf('%5.1f    %10.4f   %10.4f   %10.4f   %10.4f\n', tGrid(i), T2(i, :));
end

%% ------------------------------------------------------------------------
%  Tables 3 and 4: MTTF for varying run thresholds (Section 5.2)
%  ------------------------------------------------------------------------
bGrid = 2:6;
T3 = zeros(numel(bGrid), 1);
for i = 1:numel(bGrid)
    [Qa, ~, aa] = model1_matrices(qP, bGrid(i));
    Qb          = model1_matrices(qU, bGrid(i));
    T3(i) = mttf_protection(Qa, Qb, aa, lambda, alpha);
end

b1Grid = 3:5; b2Grid = 1:3;
T4 = NaN(numel(b1Grid), numel(b2Grid));
for i = 1:numel(b1Grid)
    for j = 1:numel(b2Grid)
        if b2Grid(j) < b1Grid(i)                 % assumption b2 < b1
            [Qa, ~, aa] = model2_matrices(qP, b1Grid(i), b2Grid(j));
            Qb          = model2_matrices(qU, b1Grid(i), b2Grid(j));
            T4(i, j) = mttf_protection(Qa, Qb, aa, lambda, alpha);
        end
    end
end

fprintf('\nTable 3. Model 1: MTTF (years) for varying b\n');
fprintf('   b     MTTF\n');
for i = 1:numel(bGrid)
    fprintf('%4d   %7.3f\n', bGrid(i), T3(i));
end
fprintf('\nTable 4. Model 2: MTTF (years) for varying (b1, b2)\n');
fprintf('  b1 \\ b2      1          2          3\n');
for i = 1:numel(b1Grid)
    fprintf('%4d   ', b1Grid(i));
    for j = 1:numel(b2Grid)
        if isnan(T4(i, j)), fprintf('%10s ', '--'); else, fprintf('%10.3f ', T4(i, j)); end
    end
    fprintf('\n');
end

%% ------------------------------------------------------------------------
%  Table 5: value of the protection block (Section 5.3)
%  ------------------------------------------------------------------------
scen  = [0, alpha, Inf];          % perfect protection, with protection, none
tSel  = [1 2 5 10];
T5m1 = zeros(1 + numel(tSel), 3);
T5m2 = zeros(1 + numel(tSel), 3);
for k = 1:3
    T5m1(1, k)     = mttf_protection(Q1m1, Q2m1, am1, lambda, scen(k));
    T5m1(2:end, k) = survival_protection(Q1m1, Q2m1, am1, lambda, scen(k), tSel)';
    T5m2(1, k)     = mttf_protection(Q1m2, Q2m2, am2, lambda, scen(k));
    T5m2(2:end, k) = survival_protection(Q1m2, Q2m2, am2, lambda, scen(k), tSel)';
end

rowNames = {'E[T]', 'P(T > 1)', 'P(T > 2)', 'P(T > 5)', 'P(T > 10)'};
fprintf('\nTable 5. MTTF and survival under three protection scenarios\n');
fprintf('                Perfect      With        No\n');
fprintf('Model 1\n');
for i = 1:5
    fprintf('  %-10s %10.4f %10.4f %10.4f\n', rowNames{i}, T5m1(i, :));
end
fprintf('Model 2\n');
for i = 1:5
    fprintf('  %-10s %10.4f %10.4f %10.4f\n', rowNames{i}, T5m2(i, :));
end

%% ------------------------------------------------------------------------
%  Table 6 and Figure: distribution of N (Section 5.4)
%  ------------------------------------------------------------------------
nmax = 20000;                     % long enough for the identities below
[pm1, sm1] = pmf_protection(Q1m1, u1m1, Q2m1, u2m1, am1, lambda, alpha, nmax);
[pm2, sm2] = pmf_protection(Q1m2, u1m2, Q2m2, u2m2, am2, lambda, alpha, nmax);
[~, ENm1] = mttf_protection(Q1m1, Q2m1, am1, lambda, alpha);
[~, ENm2] = mttf_protection(Q1m2, Q2m2, am2, lambda, alpha);

nSel = [1 2 5 10 20 50 100]';
T6 = [pm1(nSel), pm2(nSel)];
fprintf('\nTable 6. Selected pmf values of N and E[N]\n');
fprintf('    n    P(N=n) Model 1   P(N=n) Model 2\n');
for i = 1:numel(nSel)
    fprintf('%5d   %12.4f     %12.4f\n', nSel(i), T6(i, :));
end
fprintf('E[N]   %12.2f     %12.2f\n', ENm1, ENm2);

%% ------------------------------------------------------------------------
%  Table 7: no-protection limit versus prior work (Section 5.5)
%  ------------------------------------------------------------------------
T7present = [T5m1(1, 3); T5m2(1, 3)];
T7prior   = [prior_work_EN('ET2019',    qU(2), qU(3), b)      / lambda;
             prior_work_EN('Ozkut2023', qU(2), qU(3), b1, b2) / lambda];
fprintf('\nTable 7. MTTF at the no-protection limit (alpha -> Inf)\n');
fprintf('            Present formula   Prior-work formula\n');
fprintf('  Model 1      %10.4f        %10.4f   (Eryilmaz and Tekin 2019, Cor. 1)\n', T7present(1), T7prior(1));
fprintf('  Model 2      %10.4f        %10.4f   (Ozkut 2023, Cor. 2)\n', T7present(2), T7prior(2));

%% ------------------------------------------------------------------------
%  Table 8: Polya arrivals (Section 5.6)
%  ------------------------------------------------------------------------
sGrid = [2 3 5 10 50];
T8 = zeros(numel(sGrid) + 1, 4);   % columns: M1 E[T], M1 P(T>5), M2 E[T], M2 P(T>5)
for i = 1:numel(sGrid)
    [T8(i, 1), T8(i, 2)] = polya_reliability(Q1m1, Q2m1, am1, alpha, lambda, sGrid(i), 5);
    [T8(i, 3), T8(i, 4)] = polya_reliability(Q1m2, Q2m2, am2, alpha, lambda, sGrid(i), 5);
end
T8(end, :) = [T5m1(1, 2), T5m1(4, 2), T5m2(1, 2), T5m2(4, 2)];   % Poisson case
fprintf('\nTable 8. Polya arrivals, gamma shape s, E[Lambda] = 12\n');
fprintf('     s       M1 E[T]   M1 P(T>5)    M2 E[T]   M2 P(T>5)\n');
for i = 1:numel(sGrid)
    fprintf('%6d    %9.3f   %9.3f  %9.3f   %9.3f\n', sGrid(i), T8(i, :));
end
fprintf('Poisson   %9.3f   %9.3f  %9.3f   %9.3f\n', T8(end, :));

%% ------------------------------------------------------------------------
%  Numbers quoted in the text of Section 5
%  ------------------------------------------------------------------------
fprintf('\nNumbers quoted in the text\n');
fprintf('  Sec 5.1: P(T > 5), alpha = 0.067:  Model 1 = %.4f,  Model 2 = %.4f\n', T1(5, 2), T2(5, 2));
fprintf('           P(T > 5), alpha = 1/15 :  Model 1 = %.4f,  Model 2 = %.4f\n', T5m1(4, 2), T5m2(4, 2));
fprintf('  Sec 5.2: MTTF ratio b2 = 2 vs b2 = 1:  b1 = 3: %.2f,  b1 = 4: %.2f,  b1 = 5: %.2f\n', ...
        T4(:, 2) ./ T4(:, 1));
fprintf('  Sec 5.3: gain from protection (With / No):  Model 1 = %.2f,  Model 2 = %.2f\n', ...
        T5m1(1, 2) / T5m1(1, 3), T5m2(1, 2) / T5m2(1, 3));
fprintf('           Model 2 / Model 1 MTTF without protection = %.2f\n', T5m2(1, 3) / T5m1(1, 3));
fprintf('  Sec 5.6: Polya s = 2 vs Poisson:  Model 1 +%.0f%%,  Model 2 +%.0f%%\n', ...
        100 * (T8(1, 1) / T8(end, 1) - 1), 100 * (T8(1, 3) / T8(end, 3) - 1));
n999m1 = find(cumsum(pm1) >= 0.999, 1);
n999m2 = find(cumsum(pm2) >= 0.999, 1);
[~, modem1] = max(pm1);
[~, modem2] = max(pm2);
fprintf('  Sec 5.4: smallest n with P(N <= n) >= 0.999:  Model 1 = %d,  Model 2 = %d\n', n999m1, n999m2);
fprintf('           mode of P(N = n):  Model 1 at n = %d,  Model 2 at n = %d\n', modem1, modem2);

%% ------------------------------------------------------------------------
%  Consistency checks
%  ------------------------------------------------------------------------
fprintf('\nConsistency checks\n');

% Remark 2.3: sum_n P(N=n) = 1 and sum_n n P(N=n) = E[N]
fprintf('  Remark 2.3, Model 1: |1 - sum P(N=n)| = %.1e,  |E[N] - sum n P(N=n)| = %.1e\n', ...
        abs(1 - sum(pm1)), abs(ENm1 - (1:nmax) * pm1));
fprintf('  Remark 2.3, Model 2: |1 - sum P(N=n)| = %.1e,  |E[N] - sum n P(N=n)| = %.1e\n', ...
        abs(1 - sum(pm2)), abs(ENm2 - (1:nmax) * pm2));

% Theorem 2.1(ii) against 1 - cumulative sum of (i)
fprintf('  Theorem 2.1(i) vs (ii): max |P(N>n) - (1 - sum_{k<=n} P(N=k))| = %.1e (M1), %.1e (M2)\n', ...
        max(abs(sm1 - (1 - cumsum(pm1)))), max(abs(sm2 - (1 - cumsum(pm2)))));

% E[T] = E[N]/lambda and E[T] = int_0^Inf P(T > t) dt
% (for Model 1, P(T > t) < 1e-40 beyond t = 300, so the range [0, 300] suffices)
fS1 = @(x) survival_protection(Q1m1, Q2m1, am1, lambda, alpha, x, 'exact');
ETint = integral(fS1, 0, 20, 'ArrayValued', true, 'RelTol', 1e-10, 'AbsTol', 1e-12) + ...
        integral(fS1, 20, 300, 'ArrayValued', true, 'RelTol', 1e-10, 'AbsTol', 1e-12);
fprintf('  Model 1: E[N]/lambda = %.6f,  int_0^Inf P(T > t) dt = %.6f\n', ENm1 / lambda, ETint);

% Proposition 2.1: quadrature (Algorithm 1) versus exact block exponential
dS1 = max(max(abs(T1 - cell2mat(arrayfun(@(al) survival_protection(Q1m1, Q2m1, am1, lambda, al, tGrid, 'exact'), ...
               alphas, 'UniformOutput', false)))));
dS2 = max(max(abs(T2 - cell2mat(arrayfun(@(al) survival_protection(Q1m2, Q2m2, am2, lambda, al, tGrid, 'exact'), ...
               alphas, 'UniformOutput', false)))));
fprintf('  Proposition 2.1, quadrature vs exact: max difference %.1e (M1), %.1e (M2)\n', dS1, dS2);

% Remark 4.1: Model 2 with b2 = 1 equals Model 1 with b = b1
dM = 0;
for bb = 2:6
    [QA, uA] = model1_matrices(qU, bb);
    [QB, uB] = model2_matrices(qU, bb, 1);
    dM = max([dM, max(abs(QA(:) - QB(:))), max(abs(uA - uB))]);
end
fprintf('  Remark 4.1, Model 2 (b2 = 1) vs Model 1: max matrix difference %.1e\n', dM);

% Example of Section 4: b1 = 3, b2 = 2
[QE, uE] = model2_matrices(qU, 3, 2);
q0 = qU(1); q1 = qU(2); q2 = qU(3);
QEx = [q0 q1 0 q2; q0 0 q1 q2; q0 0 0 q2; q0 q1 0 0];
uEx = [0; 0; q1; q2];
fprintf('  Section 4 example (b1 = 3, b2 = 2): max difference %.1e\n', ...
        max([abs(QE(:) - QEx(:)); abs(uE - uEx)]));

% rows of [Q, u] sum to one
rs = [sum([Q1m1, u1m1], 2); sum([Q2m1, u2m1], 2); sum([Q1m2, u1m2], 2); sum([Q2m2, u2m2], 2)];
fprintf('  Row sums of [Q_r, u_r]: max |1 - row sum| = %.1e\n', max(abs(1 - rs)));
fprintf('  Spectral radii: Model 1 (%.4f, %.4f), Model 2 (%.4f, %.4f)\n', ...
        max(abs(eig(Q1m1))), max(abs(eig(Q2m1))), max(abs(eig(Q1m2))), max(abs(eig(Q2m2))));



%% ------------------------------------------------------------------------
%  Figures (Section 5)
%  ------------------------------------------------------------------------
if makeFigures
    tPlot = linspace(0, 15, 301)';
    styles = {'-', '--', '-.', ':'};

    figure('Name', 'Survival Model 1');
    hold on;
    for k = 1:numel(alphas)
        plot(tPlot, survival_protection(Q1m1, Q2m1, am1, lambda, alphas(k), tPlot, 'exact'), ...
             styles{k}, 'LineWidth', 1.5);
    end
    hold off; box on; grid on;
    xlim([0 15]); ylim([0 1]);
    xlabel('t (years)'); ylabel('P(T > t)');
    legend(arrayfun(@(x) sprintf('\\alpha = %.3f', x), alphas, 'UniformOutput', false), ...
           'Location', 'northeast');
    title('Model 1: survival function, b = 4');
    print('-dpng', '-r300', 'fig_survival_model1.png');

    figure('Name', 'Survival Model 2');
    hold on;
    for k = 1:numel(alphas)
        plot(tPlot, survival_protection(Q1m2, Q2m2, am2, lambda, alphas(k), tPlot, 'exact'), ...
             styles{k}, 'LineWidth', 1.5);
    end
    hold off; box on; grid on;
    xlim([0 15]); ylim([0 1]);
    xlabel('t (years)'); ylabel('P(T > t)');
    legend(arrayfun(@(x) sprintf('\\alpha = %.3f', x), alphas, 'UniformOutput', false), ...
           'Location', 'southwest');
    title('Model 2: survival function, (b_1, b_2) = (4, 2)');
    print('-dpng', '-r300', 'fig_survival_model2.png');

    figure('Name', 'pmf of N', 'Position', [100 100 1000 380]);
    subplot(1, 2, 1);
    bar(1:n999m1, pm1(1:n999m1), 1, 'EdgeColor', 'none');
    xlim([0 n999m1]); box on;
    xlabel('n'); ylabel('P(N = n)'); title('Model 1, b = 4');
    subplot(1, 2, 2);
    bar(1:n999m2, pm2(1:n999m2), 1, 'EdgeColor', 'none');
    xlim([0 n999m2]); box on;
    xlabel('n'); ylabel('P(N = n)'); title('Model 2, (b_1, b_2) = (4, 2)');
    print('-dpng', '-r300', 'fig_pmf_N.png');

    fprintf('\nFigures saved: fig_survival_model1.png, fig_survival_model2.png, fig_pmf_N.png\n');
end

diary off;
