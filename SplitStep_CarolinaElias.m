%{
                SIMULAÇÃO 3 - MÉTODO SPLIT-STEP FOURIER
- Disciplina: Comunicações Ópticas
- Aluna: Carolina Elias

A propagação de um pulso em uma fibra é governada pela equação não
linear de Schrödinger (NLSE), que junta a dispersão cromática (parte
linear) e o efeito Kerr (parte não linear - na sílica chi(2)=0, então a
primeira não linearidade é a de terceira ordem):
  dA/dz = -i(beta2/2) d2A/dT2 + (beta3/6) d3A/dT3 + i*gamma*|A|^2*A

O método split-step, apresentado na Seção 2.4.1 do Agrawal, recebe esse
nome porque a cada passo dz separa a equação em uma parte só linear e
outra só não linear: sozinha, cada uma tem solução simples (a dispersão
é uma multiplicação no domínio da frequência, o Kerr uma multiplicação
no domínio do tempo); juntas, não. O método aplica uma de cada vez em
cada passo, com FFTs conectando os dois domínios.

O split-step "enxerga" a dispersão apenas através das derivadas da curva
beta(omega) na portadora - beta2 = d2beta/domega2 e beta3 = d3beta/domega3. 
Por isso, assim como visto em aula, o esboço de beta(omega) é importante: 
o operador linear de cada passo é:
                exp[i(beta2/2 w^2 + beta3/6 w^3) dz], 
ou seja, o polinômio de dispersão.

Fenômeno simulado: sólitons. Um sóliton aparece quando a dispersão
anômala (beta2 < 0) equilibra o Kerr. A ordem do sóliton é:
                N = sqrt(gamma P0 T0^2 / |beta2|),
que mede a razão entre as forças não linear e dispersiva (N^2 = LD/LNL).
N=1 é o sóliton fundamental (propaga sem mudar de forma - tem solução
analítica e valida o método); N inteiro > 1 são sólitons de ordem
superior, que evoluem periodicamente. Ligando o beta3 real da sílica, a
periodicidade se quebra: o sóliton fissiona e emite uma onda dispersiva
- o efeito de 3ª ordem pedido, e a semente da geração de supercontinuum
comentada em aula.
%}

clear; clc; close all;

%%
% 1. Dispersão da sílica: beta2 e % beta3 em 1550 nm (via Sellmeier)-------
c    = 299792458;           % [m/s]
lam0 = 1550e-9;             % portadora [m]

B = [0.6961663, 0.4079426, 0.8974794];
L = [0.0684043, 0.1162414, 9.896161];
nfun = @(lu) sqrt(1 + B(1)*lu.^2./(lu.^2-L(1)^2) ...
                    + B(2)*lu.^2./(lu.^2-L(2)^2) ...
                    + B(3)*lu.^2./(lu.^2-L(3)^2));

% Derivadas de n(lambda) por diferenças finitas centradas (lambda em um)
lu = lam0*1e6;  h = 1e-4;
d2n = (nfun(lu+h) - 2*nfun(lu) + nfun(lu-h))/h^2 * 1e12;                  % [1/m^2]
d3n = (nfun(lu+2*h)-2*nfun(lu+h)+2*nfun(lu-h)-nfun(lu-2*h))/(2*h^3)*1e18; % [1/m^3]

%{
Relações entre a dispersão e as derivadas do índice:
  D     = -(lambda/c) d2n/dlambda2                (dos slides da disciplina)
  beta2 =  lambda^3/(2 pi c^2) d2n/dlambda2
  beta3 = -lambda^4/(4 pi^2 c^3) (3 d2n/dlambda2 + lambda d3n/dlambda3)
%}
beta2 = (lam0^3/(2*pi*c^2) * d2n) * 1e24 * 1e3;            % [ps^2/km]
beta3 = (-lam0^4/(4*pi^2*c^3) * (3*d2n + lam0*d3n)) * 1e36 * 1e3;  % [ps^3/km]
D     = (-lam0/c * d2n) * 1e6;                             % [ps/(nm.km)]

fprintf('Dispersao da silica em 1550 nm:\n');
fprintf('  n = %.5f\n', nfun(lu));
fprintf('  beta2 = %.3f ps^2/km  (dispersao anomala, permite solitons)\n', beta2);
fprintf('  beta3 = %.4f ps^3/km\n', beta3);
fprintf('  D     = %.2f ps/(nm.km)\n', D);
lv = linspace(1.20, 1.40, 4001);
d2v = (nfun(lv+h) - 2*nfun(lv) + nfun(lv-h))/h^2;
[~,iz] = min(abs(d2v));
fprintf('  zero de dispersao material: %.0f nm (conhecido ~1273 nm)\n\n', lv(iz)*1e3);
% Numa fibra real (SMF-28) a dispersão de guia soma-se à material e leva o
% zero para ~1310 nm, com D(1550) ~ 17 ps/(nm.km).

% 1b. Ordem do polinômio de ajuste (ponto levantado em aula)
%{
Se em vez das relações acima obtivéssemos beta2 e beta3 ajustando um
polinômio a beta(omega) e derivando, a ordem do polinômio importa: para
chegar em beta3 são 3 derivadas, e um polinômio de grau 3 tem terceira
derivada CONSTANTE - não acompanha a variação de beta3 com omega e dá um
erro sistemático, maior quanto mais larga a faixa de ajuste. Grau 5-6 já
converge. Na sílica pura (curva suave) o erro do grau 3 é moderado; em
curvas com mais estrutura (ex.: beta(w) com dispersão de guia de onda)
ele fica pior - daí a recomendação de aula de usar grau > 3.
%}
lam_j  = (1.2:0.001:1.9)*1e-6;
omega  = 2*pi*c./lam_j;
beta_j = nfun(lam_j*1e6).*omega/c;
omega0 = 2*pi*c/lam0;
es = std(omega - omega0);
xj = (omega - omega0)/es;
fprintf('Ajuste de beta(omega) - beta3 [ps^3/km] vs grau do polinomio:\n');
for grau = [3 5 6]
    p  = polyfit(xj, beta_j, grau);
    b3 = polyval(polyder(polyder(polyder(p))), 0)/es^3 *1e36*1e3;
    fprintf('  grau %d: beta3 = %.4f\n', grau, b3);
end
fprintf('  valor analitico (secao 1): beta3 = %.4f\n\n', beta3);

%% 
% 2. Parâmetros do pulso e normalização (unidades de sóliton do Agrawal)
%{
  tau = T/T0,  xi = z/LD,  u = A/sqrt(P0)
  du/dxi = -i(s/2) d2u/dtau2 + delta3 d3u/dtau3 + i N^2 |u|^2 u
  s = sgn(beta2);  N = sqrt(gamma P0 T0^2/|beta2|);  delta3 = beta3/(6|beta2|T0)
%}

T0    = 0.050;                 % largura do pulso [ps]  (50 fs)
gamma = 1.3;                   % coef. não linear [1/(W.km)]

s      = sign(beta2);          % -1 => dispersão anômala (regime de sóliton)
LD     = T0^2/abs(beta2);      % comprimento de dispersão [km]
delta3 = beta3/(6*abs(beta2)*T0);
P0     = abs(beta2)/(gamma*T0^2);   % potência de pico do sóliton N=1
z0     = (pi/2)*LD;            % período do sóliton

fprintf('Pulso de %.0f fs em 1550 nm:\n', T0*1e3);
fprintf('  LD = %.2f cm | z0 = %.2f cm | P0(N=1) = %.2f kW | delta3 = %.4f\n\n',...
        LD*1e5, z0*1e5, P0/1e3, delta3);

nt   = 2048;                   % pontos da FFT (livro usa 1024; aqui 2048
                               % para resolver melhor a onda dispersiva)
Tmax = 32;                     % janela: tau em [-32, 32]  (como no livro)
dtau = 2*Tmax/nt;
tau  = (-nt/2:nt/2-1)*dtau;    % grade temporal
omega = (pi/Tmax)*[0:nt/2-1, -nt/2:-1];   % grade de frequência (do livro)

energia = @(u) sum(abs(u).^2)*dtau;
uu0 = sech(tau);               % pulso de entrada (perfil de sóliton, m=0)

%% 
% 3. Cenário A - sóliton fundamental N=1 (validação do método)
[UA, xiA] = ssfm(uu0, omega, s, 1, 0, 10, 201);
erroA = max(abs(abs(UA(end,:)) - sech(tau)));
fprintf('Cenario A (N=1, sem beta3, 10 LD):\n');
fprintf('  desvio maximo em relacao ao sech analitico: %.2e\n', erroA);
fprintf('  energia inicial/final: %.6f / %.6f\n\n', energia(uu0), energia(UA(end,:)));

figure('Name','Cenario A - soliton fundamental','Position',[90 90 900 400]);
subplot(1,2,1);
imagesc(tau, xiA, abs(UA)); axis xy; xlim([-10 10]); colorbar;
xlabel('\tau = T/T_0'); ylabel('\xi = z/L_D'); title('N = 1: |u(\xi,\tau)|');
subplot(1,2,2);
plot(tau, abs(uu0),'k-', tau, abs(UA(end,:)),'r--','LineWidth',1.3);
xlim([-10 10]); grid on; xlabel('\tau'); ylabel('|u|');
legend('entrada','saida (\xi=10)');
title('O soliton fundamental propaga sem mudar de forma');

%% 
% 4. Cenário B - sóliton de ordem N=3 (evolução periódica)
[UB, xiB] = ssfm(uu0, omega, s, 3, 0, pi, 201);
[~,kp] = min(abs(xiB - pi/2));
erroB = norm(abs(UB(kp,:)) - abs(UB(1,:)))/norm(abs(UB(1,:)));
fprintf('Cenario B (N=3): desvio relativo da forma inicial em xi=pi/2: %.3f\n\n', erroB);

figure('Name','Cenario B - soliton N=3','Position',[100 100 640 440]);
imagesc(tau, xiB, abs(UB)); axis xy; xlim([-6 6]); colorbar;
xlabel('\tau = T/T_0'); ylabel('\xi = z/L_D');
title('N = 3: evolucao periodica, periodo \xi_0 = \pi/2');
yline(pi/2,'w--','\xi_0','LineWidth',1.2); yline(pi,'w--','2\xi_0','LineWidth',1.2);

%% 
% 5. Cenário C - dispersão de 3ª ordem: fissão e onda dispersiva
%{
Usa o delta3 do beta3 real da sílica (seção 1). Precisa de N>1 e de uma
distância maior para a perturbação acumular.
%}
%{
Janela temporal propria para este cenario: a onda dispersiva viaja a
dtau/dxi ~ 16 em relacao ao soliton, ou seja, ~600 unidades de tau em 40 LD.
Na janela de +-32 ela daria varias voltas pela borda periodica da FFT e
voltaria a cruzar o soliton (interagindo por XPM), distorcendo o espectro.
Mantem-se o mesmo dtau (0.03125) e alarga-se a janela para +-1024.
%}
TmaxC  = 1024;  ntC = 65536;
tauC   = (-ntC/2:ntC/2-1)*(2*TmaxC/ntC);
omegaC = (pi/TmaxC)*[0:ntC/2-1, -ntC/2:-1];
uu0C   = sech(tauC);
idxC   = find(tauC >= -20 & tauC <= 40);   % trecho guardado para o mapa
[UC, xiC, uC_fim] = ssfm(uu0C, omegaC, s, 2, delta3, 40, 301, idxC);
ws   = fftshift(omegaC);
Sin  = abs(fftshift(ifft(uu0C))).^2;
Sout = abs(fftshift(ifft(uC_fim))).^2;

fprintf('Cenario C (N=2, delta3=%.4f, 40 LD):\n', delta3);
fprintf('  energia inicial/final: %.6f / %.6f\n', energia(uu0C), energia(uC_fim));
fprintf('  a periodicidade e quebrada pela 3a ordem -> fissao do soliton\n\n');

figure('Name','Cenario C - efeito da 3a ordem','Position',[110 110 940 400]);
subplot(1,2,1);
imagesc(tauC(idxC), xiC, abs(UC)); axis xy; xlim([-15 30]); colorbar;
xlabel('\tau'); ylabel('\xi');
title(sprintf('N = 2 com \\delta_3 = %.3f: fissao do soliton', delta3));
subplot(1,2,2);
semilogy(ws, Sin/max(Sin),'k-', ws, Sout/max(Sout),'r-','LineWidth',1.1);
xlim([-40 60]); ylim([1e-8 2]); grid on;
xlabel('(\omega-\omega_0)T_0'); ylabel('espectro (norm.)');
legend('entrada','saida (\xi=40)','Location','southwest');
title('Onda dispersiva gerada pela dispersao de 3a ordem');

%{

                        ANÁLISE DOS RESULTADOS

1) Dispersão: as relações analíticas aplicadas à Sellmeier da sílica dão
   beta2 = -27.9 ps^2/km, beta3 = +0.152 ps^3/km, D = +21.9 ps/(nm.km) e
   zero de dispersão em 1273 nm - todos coincidem com os valores
   conhecidos do material, validando a seção. O sinal negativo de beta2
   (dispersão anômala) é o que permite sólitons em 1550 nm.

2) Ordem do polinômio: ajustar beta(omega) e derivar três vezes com grau
   3 erra o beta3 sistematicamente (terceira derivada forçada a
   constante); grau 5-6 recupera o valor analítico. Na sílica pura o
   erro do grau 3 é moderado porque a curva é suave, mas o princípio
   vale - daí a recomendação de aula de usar grau maior que 3 (sem
   exagerar, para o ajuste não oscilar).

3) Cenário A: o sóliton fundamental atravessa 10 comprimentos de
   dispersão sem se deformar (desvio ~1e-4 vs sech) e a energia se
   conserva até a precisão numérica - cada passo do método é apenas uma
   multiplicação de fase, então ele é naturalmente conservativo. Como o
   sóliton só existe com o equilíbrio exato entre dispersão e Kerr, esse
   resultado valida a implementação.

4) Cenário B: o sóliton N=3 comprime, se divide e se recompõe com o
   período teórico xi0 = pi/2, reproduzindo a figura clássica do Agrawal.
   É o balanço periódico entre a automodulação de fase (que gera chirp)
   e a dispersão anômala (que comprime o pulso com chirp).

5) Cenário C: com o beta3 real da sílica (delta3 ~ 0.018), a
   periodicidade se quebra - o pulso fissiona em sólitons fundamentais
   com velocidades distintas e emite uma onda dispersiva, que aparece
   como um pico lateral no espectro de saída, no lado AZUL
   ((omega-omega0)T0 ~ +28, perto da previsão de casamento de fase
   1/(2 delta3) ~ 27.5), enquanto o sóliton se atrasa (desloca para tau>0),
   como esperado para beta3 > 0. Essa é a assinatura do
   efeito de 3ª ordem, e o mecanismo que, com potência suficiente,
   inicia a geração de supercontinuum discutida em aula.

6) delta3 = beta3/(6|beta2|T0) explica por que a 3ª ordem só importa
   para pulsos ultracurtos: com T0 = 5 ps (telecom) delta3 ~ 2e-4 e o
   efeito é imperceptível; com 50 fs ele já domina a dinâmica.

7) Limitação do modelo: para T0 = 50 fs, o espalhamento Raman
   intrapulso (T_R/T0 ~ 0.06, que desloca o sóliton para o vermelho) e o
   self-steepening são da mesma ordem que delta3 e não foram incluídos.
   A simulação isola o efeito da 3ª ordem; uma descrição completa
   exigiria a NLSE generalizada.
%}

%% ========================================================================
% Função local: split-step de Fourier
%{
Adaptada do m-file do Agrawal (Nonlinear Fiber Optics, 5ª ed., Ap. B).
Mantém o esquema do livro e acrescenta o termo de 3ª ordem no operador de 
dispersão. Aqui a função também guarda a evolução ao longo de z (matriz U) 
para montar os mapas |u(xi,tau)|.
%}

function [U, xi, ufim] = ssfm(uu, omega, s, N, delta3, dist, nsave, idx)
% idx (opcional): indices de tau guardados na matriz U (para janelas
% grandes, guardar a grade inteira em todos os passos ocuparia centenas de
% MB). ufim e o campo completo na saida.
    if nargin < 8, idx = 1:numel(uu); end
    step_num = round(40*dist*max(N,1)^2);  % livro usa 20; usei 40 p/ folga
    deltaz   = dist/step_num;

    %{
    Sinal do termo de 3a ordem: com u(tau) = sum_w u~(w) exp(-i w tau)
    (convencao das FFTs abaixo, a mesma do livro), d/dtau -> -i w, logo
    d3/dtau3 -> +i w^3, e o operador linear de
        du/dxi = -i(s/2) d2u/dtau2 + delta3 d3u/dtau3
    e exp[(i s w^2/2 + i delta3 w^3) dxi]. O termo de beta2 nao distingue
    a convencao (w^2 e par), o de beta3 sim: com o sinal trocado
    simula-se beta3 < 0 (onda dispersiva no lado vermelho).
    %}
    dispersion = exp((0.5i*s*omega.^2 + 1i*delta3*omega.^3)*deltaz);
    hhz = 1i*N^2*deltaz;

    U  = zeros(nsave, numel(idx));
    xi = linspace(0, dist, nsave);
    ksave = round(linspace(1, step_num, nsave));

    U(1,:) = uu(idx);
    temp = uu.*exp(abs(uu).^2.*hhz/2);
    isv = 1;
    for n = 1:step_num
        f_temp = ifft(temp).*dispersion;
        uu     = fft(f_temp);
        temp   = uu.*exp(abs(uu).^2.*hhz);
        while isv < nsave && n == ksave(isv+1)
            isv = isv+1;
            usv = temp.*exp(-abs(uu).^2.*hhz/2);
            U(isv,:) = usv(idx);
        end
    end
    ufim = temp.*exp(-abs(uu).^2.*hhz/2);
    U(end,:) = ufim(idx);
end