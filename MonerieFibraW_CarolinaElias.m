%{
                  SIMULAÇÃO 4 - PAPER MONERIE FIBRA W
- Disciplina: Comunicações Ópticas
- Aluna: Carolina Elias

Baseada no artigo do Monerie, eu tentei simular o PERFIL (b) da Fig. 9.
O nome "fibra W" vem do formato do perfil quando a casca interna é
DEPRIMIDA (n2 < n): o corte radial do índice lembra a letra W.

PERFIL DE ÍNDICE (três regiões, de dentro para fora):
  núcleo:        raio a, índice n1  (o mais alto)
  casca interna: raio b, índice n2
  casca externa: índice n (sílica, semi-infinita)

PARÂMETROS (Fig. 1 e Seção II do artigo):
  Dn = n1-n2 ;  Dn' = n2-n
  V = k0 a sqrt(2 n (Dn+Dn'))   -> frequência normalizada
  R = Dn'/Dn                    -> razão dos saltos de índice (R<0 = W)
  S = b/a                       -> razão dos raios
  B = (beta - k0 n)/(k0 (Dn+Dn'))  -> parâmetro de propagação, em (0,1)

PARÂMETROS DO PERFIL (b) (lidos da legenda da Fig. 9):
  b/a = 1.5           -> S = 1.5
  2a = 6.5 um         -> a = 3.25 um
  Dn+Dn' = 6.3e-3  e  Dn' = -6.3e-3  ->  Dn = 12.6e-3,  R = -0.5
O artigo afirma que essa fibra tem dispersão total nula ~1.55 um.

O QUE ESTE CODIGO FAZ:
  Parte 1: curvas universais B(V) do LP01 variando o perfil R (base do
           método; a Fig. 5 do artigo).
  Parte 2: SIMULAÇÃO DO PERFIL (b) - a dispersão cromática total dessa
           fibra em função de lambda, separada nas duas contribuições:
             D_total(lambda) = D_material(lambda) + D_guia(lambda)
           A material vem da sílica (Sellmeier, como no split-step); a de
           guia vem da derivada V d2(VB)/dV2 de B(V), avaliada ao longo
           de lambda (via V(lambda)). O resultado mostra o cancelamento
           perto de 1.55 um - o ponto central do artigo.

MÉTODO NUMÉRICO (equação de dispersão, LP01, m=0):
No guiamento fraco, o perfil radial e (Eqs. 1-2 do artigo):
  r<=a:  psi = A0 J0(u r/a)
  a<r<=b: psi = A1 P0(w r/b) + A2 Q0(w r/b)   (P,Q)=(I,K) se evanescente
                                              na casca interna, (J,Y) se
                                              oscilante
  r>b:   psi = A3 K0(v r/b)
com u=V sqrt(1-B), v=V S sqrt(B), w=V S sqrt|B-R/(1+R)|. Continuidade de
psi e dpsi/dr nas interfaces -> sistema 4x4; a equação de dispersão e
det(M(B))=0. Para cada V, a raiz de maior B e o LP01.

Nota sobre a dispersão de guia: usei a forma de Gloge,
  D_guia = -(n Delta)/(c lambda) V d2(VB)/dV2 ,  Delta = (Dn+Dn')/n.
que é uma aproximação (não inclui o termo de dispersão de perfil p da Eq. 11
do artigo), entao o cancelamento sai como um MAXIMO de D_total logo abaixo
de zero (~ -0.4 ps/nm/km) em torno de 1.5-1.55 um, em vez de um cruzamento
exato, o suficiente para demonstrar o efeito.
%}

clear; clc; close all;

%% 
% Parametros do perfil (b) da Fig. 9
c     = 299792458;      % [m/s]
S     = 1.5;            % b/a
a     = 3.25e-6;        % raio do núcleo [m]  (2a = 6.5 um)
DnTot = 6.3e-3;         % Dn + Dn'
R     = -0.5;           % Dn'/Dn (casca interna deprimida)

% Sellmeier da sílica (índice da casca externa e referência)
Bs = [0.6961663, 0.4079426, 0.8974794];
Ls = [0.0684043, 0.1162414, 9.896161];
nsil = @(lu) sqrt(1 + Bs(1)*lu.^2./(lu.^2-Ls(1)^2) ...
                    + Bs(2)*lu.^2./(lu.^2-Ls(2)^2) ...
                    + Bs(3)*lu.^2./(lu.^2-Ls(3)^2));

fprintf('Simulacao 4 - Fibra W (Monerie), perfil (b) da Fig. 9\n');
fprintf('S=%.1f, a=%.2f um, Dn+Dn''=%.1e, R=%.2f\n\n', S, a*1e6, DnTot, R);

%%
% PARTE 1 - curvas universais B(V) do LP01 variando o perfil R
Vvec = linspace(1.0, 5.0, 140);
R_list = [0.5, 0.0, -0.25, -0.5, -0.75];

figure('Name','Fibra W: B(V) universal','Position',[60 60 720 500]);
hold on; cores = lines(numel(R_list));
for i = 1:numel(R_list)
    B = arrayfun(@(V) B_LP01(V, R_list(i), S), Vvec);
    est = '-'; lw = 1.4;
    if R_list(i) == R, est = '-'; lw = 2.6; end   % destaca o perfil (b)
    plot(Vvec, B, est, 'LineWidth', lw, 'Color', cores(i,:), ...
         'DisplayName', sprintf('R = %+.2f', R_list(i)));
end
xlabel('V'); ylabel('B'); grid on; ylim([0 1]);
title(sprintf('B(V) universal do LP01 (S = %.1f). Linha grossa: perfil (b), R=%.1f', S, R));
legend('Location','southeast');

%% 
% PARTE 2 - dispersão total do perfil (b) em função de lambda
lam = linspace(1.20, 2.00, 130) * 1e-6;   % [m]  (faixa larga p/ ver a volta)
Dmat = zeros(size(lam));
Dwg  = zeros(size(lam));
for i = 1:numel(lam)
    Dmat(i) = disp_material(lam(i), c, nsil);
    Dwg(i)  = disp_guia(lam(i), c, a, S, R, DnTot, nsil);
end
Dtot = Dmat + Dwg;

% Restringi a busca para evita marcar o segundo cruzamento no zero (próximo do corte)
mask = lam*1e6 <= 1.70;
[~, ib] = min(abs(Dtot(mask)));
idxs = find(mask);
imin = idxs(ib);
lam_op = lam(imin)*1e6;   % comprimento de onda de dispersão ~zero (operação)

figure('Name','Fibra W perfil (b): dispersao total','Position',[100 100 760 520]);
plot(lam*1e6, Dmat, 'b--', 'LineWidth', 1.5); hold on;
plot(lam*1e6, Dwg,  'r--', 'LineWidth', 1.5);
plot(lam*1e6, Dtot, 'k-',  'LineWidth', 2.4);
yline(0, 'Color', [.5 .5 .5]);
xline(lam_op, 'g:', sprintf('%.3f \\mum', lam_op), 'LineWidth', 1.2);
plot(lam_op, Dtot(imin), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
xlabel('\lambda [\mum]'); ylabel('D [ps/(nm\cdotkm)]');
title(sprintf(['Perfil (b): dispersao total achatada perto de %.2f \\mum\n' ...
    'minimo |D_{tot}| = %.2f ps/(nm km) em \\lambda = %.3f \\mum'], ...
    lam_op, abs(Dtot(imin)), lam_op));
legend('material D_{mat} (silica)', 'guia de onda D_{wg}', ...
       'total D_{tot}', 'Location', 'northwest');
grid on;

%%
% Saida numerica
fprintf('Dispersao do perfil (b) em funcao de lambda:\n');
fprintf(' lam[um] |  Dmat  |  Dwg   |  Dtot  [ps/(nm km)]\n');
for lu = [1.30 1.40 1.50 1.55 1.60 1.70 1.80 1.90 2.00]
    lm = lu*1e-6;
    dm = disp_material(lm, c, nsil);
    dw = disp_guia(lm, c, a, S, R, DnTot, nsil);
    fprintf('  %.2f   | %6.2f | %6.2f | %6.2f\n', lu, dm, dw, dm+dw);
end
fprintf('\nPonto de operacao (dispersao minima na banda plana):\n');
fprintf('  D_tot = %.2f ps/(nm km) em lambda = %.3f um\n', Dtot(imin), lam(imin)*1e6);

%{

                        ANÁLISE DOS RESULTADOS

1) Curvas universais B(V) (Parte 1): reproduzem a Fig. 5 do artigo.
   Quanto mais negativo R (casca interna mais deprimida), mais baixa a
   curva - a fibra W guia mais fracamente que a comum e o LP01 passa a
   ter corte em V finito (B -> 0), ao contrário da fibra degrau. O perfil (b) que
   vamos simular (R=-0.5, linha grossa) é uma dessas curvas.

2) Simulação do perfil (b) (Parte 2): a dispersão cromática total é a
   soma da material (positiva, ~+22 ps/nm/km em 1.55 um na sílica) com a
   de guia de onda (negativa). Para o perfil (b), essas duas parcelas se
   cancelam quase exatamente em torno de 1.5-1.55 um: a dispersão total
   sobe até um MÁXIMO logo abaixo de zero (D_tot ~ -0.4 ps/(nm km)) e fica
   dentro de ~1 ps/(nm km) numa banda larga (~1.46-1.60 um). É o que o
   Monerie projeta na Fig. 9.
   Estendendo lambda, a curva total desce até um mínimo (~ -6 ps/(nm km)
   perto de 1.85 um) e depois VOLTA A SUBIR, cruzando zero perto de 2 um.
   Essa subida ocorre porque V diminui (V ~ 1/lambda) e o modo se aproxima
   do seu corte (o LP01 da fibra W tem corte em V finito), onde a dispersão
   de guia muda de comportamento e deixa de compensar a material. Por isso
   a busca do ponto de operação foi restrita a lambda <= 1.7 um.

3) Por que isso importa: numa fibra degrau comum, obter dispersão zero em
   1.55 um exige alto contraste de índice (muita dopagem de germânio),
   o que aumenta as perdas por espalhamento. A fibra W consegue o mesmo
   com Dn menor, porque a casca interna deprimida dá uma dispersão de
   guia forte e ajustável. Poder mexer em R e S (o "brincar com o
   perfil") é o que dá essa liberdade extra de projeto.

4) Sobre a aproximação: a dispersão de guia foi calculada pela forma de
   Gloge, sem o termo de dispersão de perfil p (Eq. 11 do artigo) nem o
   índice de grupo exato da casca. Por isso a curva total fica tangente a
   zero por baixo (~ -0.4 ps/nm/km) em vez de cruzar exatamente - o
   suficiente para demonstrar o cancelamento. Incluir p e o índice de grupo
   possivelmente desloca a curva o bastante para cruzar zero em 1.55 um,
   como no artigo (não verificado aqui).

5) Observação de aula: a parte oscilante do campo fica dentro do guia e a
   evanescente, fora - por isso a equação de dispersão troca de forma
   (Bessel J,Y para I,K) quando beta cruza k0*n2 (o código decide pelo
   sinal de B-R/(1+R)). Aqui B é real: não há perda de propagação; a
   evanescência é só o decaimento transversal do campo fora do guia.
%}

%% ========================================================================
% Funções locais
% =========================================================================

function D = disp_material(lam, c, nsil)
% Dispersão material da sílica [ps/(nm km)]: D = -(lambda/c) d2n/dlambda2.
    lu = lam*1e6;  h = 1e-4;
    d2n = (nsil(lu+h) - 2*nsil(lu) + nsil(lu-h))/h^2 * 1e12;   % [1/m^2]
    D = -lam/c * d2n * 1e6;
end

function D = disp_guia(lam, c, a, S, R, DnTot, nsil)
% Dispersão de guia de onda [ps/(nm km)] pela forma de Gloge:
%   D_guia = -(n Delta)/(c lambda) * V d2(VB)/dV2,  Delta=(Dn+Dn')/n.
    n  = nsil(lam*1e6);
    V  = 2*pi/lam * a * sqrt(2*n*DnTot);
    vd = Vd2VB(V, R, S);
    if ~isfinite(vd), D = NaN; return; end
    Delta = DnTot/n;
    D = -(n*Delta)/(c*lam) * vd * 1e6;
end

function y = Vd2VB(V, R, S)
% V d2(VB)/dV2 do LP01 por diferencas finitas centradas.
    h = 0.02;
    f0 = VB_of(V-h,R,S); f1 = VB_of(V,R,S); f2 = VB_of(V+h,R,S);
    if ~isfinite(f0) || ~isfinite(f1) || ~isfinite(f2)
        y = NaN;
    else
        y = V*(f2 - 2*f1 + f0)/h^2;
    end
end

function y = VB_of(V, R, S)
    b = B_LP01(V, R, S);
    if isfinite(b), y = V*b; else, y = NaN; end
end

function B = B_LP01(V, R, S)
% Parâmetro de propagação B do modo fundamental LP01 (maior raiz).
%
% Para R > 0, o argumento B - R/(1+R) troca de sinal em B = R/(1+R): ali
% o determinante passa de (I,K) para (J,Y), diverge (K0, Y0 -> inf com
% w -> 0) e troca de sinal SEM ser raiz. Esses intervalos sao ignorados,
% senao a bisseccao devolve uma raiz falsa em B = R/(1+R) (ex.: 1/3 para
% R = 0.5), que viraria o LP01 por ser a maior.
    Bs = linspace(1e-4, 1-1e-4, 2000);
    dv = arrayfun(@(b) det_disp(b, V, R, S), Bs);
    argv = Bs - R/(1+R);
    raizes = [];
    for i = 1:numel(Bs)-1
        if isfinite(dv(i)) && isfinite(dv(i+1)) && dv(i)*dv(i+1) < 0 ...
                && sign(argv(i)) == sign(argv(i+1))
            raizes(end+1) = bisseccao(@(b) det_disp(b,V,R,S), Bs(i), Bs(i+1)); %#ok<AGROW>
        end
    end
    if isempty(raizes), B = NaN; else, B = max(raizes); end
end

function d = det_disp(B, V, R, S)
%{
Determinante 4x4 da equação de dispersão do LP01 (m=0). Linhas: psi e
dpsi/dr contínuos em r=a e r=b. Derivadas radiais multiplicadas por a
(linhas de r=a) e por b (linhas de r=b) para ficarem adimensionais.
Para m=0: d/dx J0=-J1, d/dx I0=I1, d/dx K0=-K1, d/dx Y0=-Y1.
%}
    B = min(max(B,1e-12), 1-1e-12);
    u = V*sqrt(1-B);
    v = V*S*sqrt(B);
    arg = B - R/(1+R);
    aob = 1/S;

    fc  = besselj(0,u);
    dfc = -u*besselj(1,u);

    if arg > 0        % casca interna evanescente: I0, K0
        w = V*S*sqrt(arg);  wa = w*aob;
        Pa = besseli(0,wa);  Qa = besselk(0,wa);
        dPa =  wa*besseli(1,wa);  dQa = -wa*besselk(1,wa);
        Pb = besseli(0,w);   Qb = besselk(0,w);
        dPb =  w*besseli(1,w);    dQb = -w*besselk(1,w);
    else              % casca interna oscilante: J0, Y0
        w = V*S*sqrt(-arg); wa = w*aob;
        Pa = besselj(0,wa);  Qa = bessely(0,wa);
        dPa = -wa*besselj(1,wa);  dQa = -wa*bessely(1,wa);
        Pb = besselj(0,w);   Qb = bessely(0,w);
        dPb = -w*besselj(1,w);    dQb = -w*bessely(1,w);
    end

    g  = besselk(0,v);
    dg = -v*besselk(1,v);

    M = [ fc,  -Pa,  -Qa,   0 ; ...
          dfc, -dPa, -dQa,  0 ; ...
          0,    Pb,   Qb,  -g ; ...
          0,    dPb,  dQb, -dg ];
    d = det(M);
end

function m = bisseccao(f, lo, hi)
% Bissecção para a raiz de f em [lo,hi] (sinais opostos).
    flo = f(lo);
    for k = 1:80
        m = (lo+hi)/2;  fm = f(m);
        if flo*fm <= 0, hi = m; else, lo = m; flo = fm; end
    end
    m = (lo+hi)/2;
end