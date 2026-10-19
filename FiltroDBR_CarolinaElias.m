%{
                 SIMULACAO 1 - FILTRO OPTICO DBR
- Disciplina: SEL0366 - Comunicacoes Opticas
- Aluna: Carolina Elias de Almeida Americo

Este script modela a resposta espectral de um filtro optico DBR
(Distributed Bragg Reflector) centrado em 1550 nm, pelo metodo das
matrizes de transferencia. Analisa-se a transmissividade para os modos
TE e TM sob incidencia obliqua (30 graus) e a evolucao do fator de
qualidade (Q) com o numero de periodos.

O filtro e um espelho de Bragg com um DEFEITO central: N periodos de
(TiO2/SiO2) de quarto de onda, uma camada defeito de meia onda (que abre
a banda de passagem estreita - a ressonancia em 1550 nm), e mais N
periodos. A camada de meia onda funciona como uma microcavidade
Fabry-Perot entre dois espelhos de Bragg.

Consideracoes do modelo: indice de refracao real e constante (sem
dispersao cromatica nem absorcao - meios ideais).

Como a fase acumulada e linear em frequencia, a resposta em FREQUENCIA 
(ultimo grafico da parte 1) sai simetrica. Como consequência, a resposta em 
COMPRIMENTO DE ONDA, que tem relacao inversa (nao linear) com a frequencia, 
nao e perfeitamente simetrica.
%}

clear; clc; close all;

%% 
% PARTE 1 - Resposta do filtro (indice de refracao fixo)

% --- Parametros de entrada ---
lambda0 = 1550e-9;      % comprimento de onda de projeto [m]
theta_inc_deg = 30;     % angulo de incidencia no ar [graus]
n_air = 1.0;            % indice do ar
n_H   = 2.45;           % TiO2 (alto indice)
n_L   = 1.44;           % SiO2 (baixo indice)
eta0  = 377;            % impedancia do vacuo [Ohm]
num_periodos = 5;       % periodos de cada lado do defeito
% (com N=5, Q ~ 1.5e3 e FWHM ~ 1 nm: o pico fica resolvido pela grade de
%  0.08 nm. Com N=20, Q ~ 5e10 e FWHM ~ 3e-8 nm - nenhuma grade resolve o
%  pico, que apareceria so no ponto de 1550 nm inserido no vetor.)

% --- Vetor de comprimentos de onda ---
lambda_base = linspace(1200e-9, 2000e-9, 10000);
lambda_vec  = unique(sort([lambda_base, lambda0]));  % garante 1550 nm no vetor

% --- Angulos internos (lei de Snell) e espessuras ---
theta_inc = deg2rad(theta_inc_deg);
theta_H = asin(n_air*sin(theta_inc)/n_H);   % refracao no TiO2
theta_L = asin(n_air*sin(theta_inc)/n_L);   % refracao no SiO2

L_H   = lambda0/(4*n_H*cos(theta_H));   % quarto de onda TiO2 (corrigido p/ angulo)
L_L   = lambda0/(4*n_L*cos(theta_L));   % quarto de onda SiO2
L_def = lambda0/(2*n_H*cos(theta_H));   % defeito: meia onda de TiO2

% --- Admitancias transversais de cada camada (TE e TM) ---
% TE: Y = n cos(theta)/eta0 ;  TM: Y = n/(cos(theta) eta0)
Y_TE_air = (n_air*cos(theta_inc))/eta0;
Y_TE_H   = (n_H*cos(theta_H))/eta0;
Y_TE_L   = (n_L*cos(theta_L))/eta0;
Y_TM_air = n_air/(cos(theta_inc)*eta0);
Y_TM_H   = n_H/(cos(theta_H)*eta0);
Y_TM_L   = n_L/(cos(theta_L)*eta0);

% --- Matrizes de interface (independentes de lambda) ---
M_TE_ini = matching_mat(Y_TE_air, Y_TE_H);   % ar -> TiO2
M_TM_ini = matching_mat(Y_TM_air, Y_TM_H);
M_TE_fin = matching_mat(Y_TE_H, Y_TE_air);   % TiO2 -> ar
M_TM_fin = matching_mat(Y_TM_H, Y_TM_air);
M_TE_H_L = matching_mat(Y_TE_H, Y_TE_L);     % TiO2 -> SiO2
M_TM_H_L = matching_mat(Y_TM_H, Y_TM_L);
M_TE_L_H = matching_mat(Y_TE_L, Y_TE_H);     % SiO2 -> TiO2
M_TM_L_H = matching_mat(Y_TM_L, Y_TM_H);

% --- Varredura em lambda: transmissividade TE e TM ---
T_TE = zeros(size(lambda_vec));
T_TM = zeros(size(lambda_vec));
for i = 1:numel(lambda_vec)
    T_TE(i) = transmissividade(lambda_vec(i), num_periodos, ...
        n_H, n_L, theta_H, theta_L, L_H, L_L, L_def, ...
        M_TE_ini, M_TE_fin, M_TE_H_L, M_TE_L_H);
    T_TM(i) = transmissividade(lambda_vec(i), num_periodos, ...
        n_H, n_L, theta_H, theta_L, L_H, L_L, L_def, ...
        M_TM_ini, M_TM_fin, M_TM_H_L, M_TM_L_H);
end

% --- Grafico 1: resposta completa em lambda ---
figure('Color','w');
plot(lambda_vec*1e9, T_TE, 'r', 'LineWidth', 1.5); hold on;
plot(lambda_vec*1e9, T_TM, 'b--', 'LineWidth', 1.5);
grid on;
title(['Resposta do Filtro DBR (N=', num2str(num_periodos), ...
       ') - Incidencia a ', num2str(theta_inc_deg), ' graus']);
xlabel('Comprimento de Onda (nm)'); ylabel('Transmissividade (T)');
legend('Polarizacao TE', 'Polarizacao TM');

% --- Grafico 2: escala semi-log (evidencia a banda proibida) ---
figure('Color','w');
semilogy(lambda_vec*1e9, T_TE, 'r', 'LineWidth', 1.5); hold on;
semilogy(lambda_vec*1e9, T_TM, 'b--', 'LineWidth', 1.5);
grid on;
title(['Resposta do Filtro (Escala Log) - N=', num2str(num_periodos)]);
xlabel('Comprimento de Onda (nm)'); ylabel('Transmissividade (escala log)');
legend('Polarizacao TE', 'Polarizacao TM');

% --- Grafico 3: resposta em frequencia (simetrica em torno da ressonancia) ---
c = 299792458;                 % velocidade da luz [m/s]
f_THz = (c ./ lambda_vec)/1e12;
figure('Color','w');
plot(f_THz, T_TE, 'r', 'LineWidth', 1.5); hold on;
plot(f_THz, T_TM, 'b--', 'LineWidth', 1.5);
grid on;
title(['Resposta do Filtro em Frequencia (N=', num2str(num_periodos), ')']);
xlabel('Frequencia (THz)'); ylabel('Transmissividade (T)');
legend('Polarizacao TE', 'Polarizacao TM');

%% 
%{
PARTE 2 - Fator de qualidade Q em funcao do numero de periodos

Q = lambda0 / FWHM. A largura a meia altura (FWHM) e obtida por
BISSECCAO das bordas onde T cruza 0.5, sem depender de uma grade fina em
lambda: para N grande o pico fica estreito demais (poucos ou nenhum
ponto de uma grade fixa cairia sobre ele), e a bisseccao resolve isso
de forma exata.
%}

max_periodos = 11;
n_vec = 1:max_periodos;
Q_TE_vec = zeros(1, max_periodos);
Q_TM_vec = zeros(1, max_periodos);

for n = 1:max_periodos
    Q_TE_vec(n) = fator_Q(lambda0, n, n_H, n_L, theta_H, theta_L, ...
        L_H, L_L, L_def, M_TE_ini, M_TE_fin, M_TE_H_L, M_TE_L_H);
    Q_TM_vec(n) = fator_Q(lambda0, n, n_H, n_L, theta_H, theta_L, ...
        L_H, L_L, L_def, M_TM_ini, M_TM_fin, M_TM_H_L, M_TM_L_H);
end

% --- Grafico: Q x N ---
figure('Color','w');
semilogy(n_vec, Q_TE_vec, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); hold on;
semilogy(n_vec, Q_TM_vec, 'b--s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on;
title('Evolucao do Fator de Qualidade com o Numero de Periodos');
xlabel('Numero de Periodos (N)'); ylabel('Fator de Qualidade (Q)');
legend('Polarizacao TE', 'Polarizacao TM', 'Location', 'best');

% --- Saida numerica ---
fprintf('Fator de qualidade Q por numero de periodos:\n');
fprintf('  N |    Q_TE     |    Q_TM\n');
for n = 1:max_periodos
    fprintf(' %2d | %10.4g | %10.4g\n', n, Q_TE_vec(n), Q_TM_vec(n));
end

%{

                        ANALISE DOS RESULTADOS

1) Resposta espectral: o filtro apresenta uma banda proibida (stopband)
   larga - onde T ~ 0, a luz e refletida - com um pico estreito de
   transmissao T = 1 no centro, em 1550 nm. Esse pico e a ressonancia da
   microcavidade formada pela camada defeito de meia onda entre os dois
   espelhos de Bragg. O pico fica bem centrado em 1550 nm porque as
   espessuras sao de quarto de onda (corrigidas pelo angulo de refracao)
   e o defeito e de meia onda.

2) Simetria: em comprimento de onda a resposta e levemente assimetrica,
   porque lambda ~ 1/f e uma relacao nao linear e o comprimento de onda
   dentro de cada camada muda ao longo do espectro; em FREQUENCIA
   (grafico 3) a resposta fica simetrica em torno da ressonancia, ja que
   a fase acumulada e linear na frequencia. Os dois efeitos sao
   consistentes com a discussao com o professor.

3) TE vs TM (incidencia obliqua a 30 graus): a polarizacao TE tem
   stopband mais larga e fator de qualidade maior que a TM. Isso ocorre
   porque, sob incidencia obliqua, a refletividade das interfaces cresce
   para TE e diminui para TM (as duas coincidem apenas na incidencia
   normal). Perto do angulo de Brewster a refletividade TM tenderia a
   zero, estreitando ainda mais a banda.

4) Fator de qualidade Q x N: Q cresce de forma aproximadamente
   exponencial com o numero de periodos - cada periodo adicional torna os
   espelhos mais refletores, estreitando a ressonancia. Para N grande o
   pico fica extremamente estreito (Q > 1e5), o que exige o calculo da
   largura por bisseccao; uma grade fixa em lambda nao resolveria o pico
   e o Q saturaria artificialmente. Em toda a faixa, Q_TE > Q_TM, pela
   mesma razao do item 3.
%}

%% 
% Funcoes locais

function T = transmissividade(lam, Nper, n_H, n_L, theta_H, theta_L, ...
                              L_H, L_L, L_def, M_ini, M_fin, M_H_L, M_L_H)
% Transmissividade do DBR com defeito central, para uma polarizacao
% (definida pelas matrizes de interface passadas) em um dado lambda.
    k0 = 2*pi/lam;
    M = M_ini;                                   % ar -> TiO2
    for p = 1:Nper                               % espelho de baixo
        M = M * prop_mat(k0*n_H*cos(theta_H), L_H) * M_H_L ...
              * prop_mat(k0*n_L*cos(theta_L), L_L) * M_L_H;
    end
    M = M * prop_mat(k0*n_H*cos(theta_H), L_def); % camada defeito (meia onda)
    for p = 1:Nper                                % espelho de cima
        M = M * M_H_L * prop_mat(k0*n_L*cos(theta_L), L_L) ...
              * M_L_H * prop_mat(k0*n_H*cos(theta_H), L_H);
    end
    M = M * M_fin;                                % TiO2 -> ar
    T = abs(1/M(1,1))^2;
end

function Q = fator_Q(lam0, Nper, n_H, n_L, theta_H, theta_L, ...
                     L_H, L_L, L_def, M_ini, M_fin, M_H_L, M_L_H)
% Fator de qualidade Q = lambda0/FWHM da ressonancia em lam0.
% A ressonancia esta centrada em lam0 (defeito de meia onda). As bordas
% de meia altura sao encontradas por bisseccao, o que resolve picos
% arbitrariamente estreitos (N grande) sem depender de grade fina.
    Tp = @(lam) transmissividade(lam, Nper, n_H, n_L, theta_H, theta_L, ...
                                 L_H, L_L, L_def, M_ini, M_fin, M_H_L, M_L_H);
    Tpico = Tp(lam0);
    if Tpico < 0.5, Q = NaN; return; end
    meia = Tpico/2;

    % borda direita: avanca ate T < meia, depois bissecciona.
    % O limite de busca precisa cobrir a banda de N=1 (FWHM ~ 130 nm em TE
    % e ~ 200 nm em TM); com um limite curto a busca parava antes da
    % borda e o Q saia superestimado.
    passo = 0.1e-9;  lim = 300e-9;  hi = lam0;
    while Tp(hi) > meia && hi < lam0 + lim, hi = hi + passo; end
    if Tp(hi) > meia, warning('fator_Q: borda direita nao encontrada (N=%d)', Nper); end
    lo = hi - passo;
    for k = 1:60
        m = (lo+hi)/2;
        if Tp(m) > meia, lo = m; else, hi = m; end
    end
    lam_R = (lo+hi)/2;

    % borda esquerda
    lo2 = lam0;
    while Tp(lo2) > meia && lo2 > lam0 - lim, lo2 = lo2 - passo; end
    if Tp(lo2) > meia, warning('fator_Q: borda esquerda nao encontrada (N=%d)', Nper); end
    hi2 = lo2 + passo;
    for k = 1:60
        m = (lo2+hi2)/2;
        if Tp(m) > meia, hi2 = m; else, lo2 = m; end
    end
    lam_L = (lo2+hi2)/2;

    FWHM = lam_R - lam_L;
    if FWHM > 0, Q = lam0/FWHM; else, Q = NaN; end
end

function M = matching_mat(Ya, Yb)
% Matriz de interface (casamento) entre dois meios de admitancias Ya, Yb.
    rho = (Ya - Yb)/(Ya + Yb);   % coeficiente de reflexao transversal
    tau = 1 + rho;               % coeficiente de transmissao transversal
    M = (1/tau) * [1, rho; rho, 1];
end

function P = prop_mat(kz, L)
% Matriz de propagacao numa camada de espessura L (kz = componente z do
% vetor de onda na camada).
    phi = kz * L;                % fase acumulada
    P = [exp(1j*phi), 0; 0, exp(-1j*phi)];
end