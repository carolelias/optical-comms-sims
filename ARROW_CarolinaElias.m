%{
                    SIMULAÇÃO 2 - GUIA ARROW
            (Antiresonant Reflecting Optical Waveguide)
- Disciplina: Comunicações Ópticas
- Aluna: Carolina Elias

ESTRUTURA SIMULADA (de cima para baixo):
  1: superstrato n=1.5 (semi-infinito)
  2: G2: núcleo n=4.0, espessura d2
  3: casca n=1.5, espessura dc1
  4: ESPELHO de silício, espessura t  <- parâmetro varrido
  5: casca n=1.5, espessura dc2
  6: substrato de silício (semi-infinito)

Como visto em aula, o guia G2 confina por reflexão interna total (núcleo
4.0 entre cascas 1.5), mas a casca inferior é finita e o substrato de Si
tem índice maior que o índice efetivo do modo (n_ef ~ 2.61 < 3.6). A
cauda evanescente do modo então "tunela" pela casca e vira onda
propagante no substrato: é um modo com fuga (leaky), com constante de
propagação complexa beta = beta_r - j*beta_i, onde a potência decai como
exp(-2*beta_i*z).

A estrutura ARROW insere na casca uma camada fina de Si (o "espelho"),
que se comporta como uma cavidade Fabry-Perot em regime de tunelamento
ressonante. Trabalhando nos pontos de ANTIRRESSONÂNCIA da camada, a
reflexão de volta ao núcleo é máxima e a perda por fuga é mínima; nas
espessuras RESSONANTES a camada fica "transparente" e a fuga dispara.
Por isso queremos operar nos pontos de antirressonância.

Assim como no exemplo de aplicação visto em aula, essa estrutura aparece
na detecção integrada: a luz é lançada em G1, acopla para G2 ao longo do
comprimento L_pi e, nesse trecho, o espelho antirressonante não acrescenta
fuga para o substrato (a perda fica praticamente igual à do guia sem
espelho, ver item 4 da análise). Onde o espelho termina, o modo volta a ser fortemente
leaky e a potência é despejada no substrato de Si, onde o detector MSM a
converte em fotocorrente. O guia G1 não precisa entrar nesta simulação: a
sobreposição evanescente entre os dois guias é desprezível para a perda
(ele só importa na análise de acoplamento), então a análise modal abaixo
considera apenas o guia G2 sobre o espelho.

MÉTODO NUMÉRICO: Para cada espessura t do espelho, a raiz complexa de n_ef 
é encontrada assim:
  1) varredura de n_ef real numa região estável (janela fixa em torno do
     modo do guia isolado), procurando os mínimos de |F(n_ef)|, onde F é
     a função de dispersão montada pelo TMM com as condições de contorno
     de modo leaky (onda de saída no substrato e campo evanescente no
     superstrato);
  2) refino de cada candidato no plano complexo (fminsearch);
  3) quando há mais de uma raiz (anticruzamento do modo de G2 com os
     modos guiados do próprio espelho), fica-se com a raiz cujo campo
     está mais confinado em G2, calculando a fração de |E|^2 dentro do
     núcleo.
A varredura em t é feita em duas passadas: uma grossa e um refinamento
automático em torno dos picos, já que as ressonâncias são estreitas
(largura de poucos nm).

São comparados dois espelhos: um sem perda (n=3.6) e o valor usado em
aula (n=3.6-j0.2), para mostrar o efeito da absorção sobre a largura dos
picos e o piso dos vales.
%}

clear; clc; close all;

%% ------------------------------------------------------------------------
% Parâmetros
lambda = 1.55;   k0 = 2*pi/lambda;      % [um], [rad/um]

n_sup = 1.5;    % superstrato
n_G2  = 4.0;    % núcleo do guia
n_cl  = 1.5;    % cascas
n_mir_A = 3.6;          % espelho sem perda
n_mir_B = 3.6 - 0.2j;   % espelho com o índice usado em aula
n_sub   = 3.6 - 0.2j;   % substrato de Si
% (em 1550 nm o c-Si real é transparente, n ~ 3.48; o valor com k=0.2 é
%  o do exemplo de aula e serve para evidenciar o papel da absorção)

d2  = 0.10;    % núcleo G2 [um]
dc1 = 0.45;    % casca entre G2 e o espelho (grossa -> ressonâncias finas)
dc2 = 0.25;    % casca entre o espelho e o substrato

t_coarse = linspace(0, 1.35, 451);   % varredura grossa de t (passo 3 nm)
w_fino   = 0.025;                    % meia-janela do refino nos picos [um]
N_fino   = 121;                      % pontos por janela (passo ~0.4 nm)

fprintf('Simulacao 2 - Guia ARROW (TMM, beta complexo)\n');
fprintf('lambda=%.0f nm | espelho A=%.2f, B=%.2f%+.2fj | dc1=%.2f dc2=%.2f um\n\n',...
    lambda*1e3, real(n_mir_A), real(n_mir_B), imag(n_mir_B), dc1, dc2);

%% 
% Modo TE0 do guia isolado (referência para a janela de busca)
%{
Equação transcendental do slab simétrico, deduzida das equações de
Maxwell + condições de contorno:
                u*tan(u) = w,   u^2 + w^2 = V^2
%}
a  = d2/2;
V  = k0*a*sqrt(n_G2^2 - n_cl^2);
u0 = fzero(@(u) u.*tan(u) - sqrt(max(V^2-u.^2,0)), [1e-6, min(V,pi/2)-1e-6]);
neff_iso = sqrt(n_G2^2 - (u0/a/k0)^2);
fprintf('Guia isolado: V = %.3f -> neff(TE0) = %.4f\n', V, neff_iso);
fprintf('(1.5 < neff < 3.6 => o modo e leaky para o substrato)\n\n');

%% 
% Varredura da espessura do espelho (dois casos)
fprintf('Varredura A (espelho sem perda)...\n');
[ttA, neff_A] = varre(t_coarse, w_fino, N_fino, neff_iso, k0, ...
    n_cl, n_mir_A, n_G2, dc1, dc2, d2, n_sup, n_sub);
fprintf('Varredura B (espelho com perda, valor da aula)...\n');
[ttB, neff_B] = varre(t_coarse, w_fino, N_fino, neff_iso, k0, ...
    n_cl, n_mir_B, n_G2, dc1, dc2, d2, n_sup, n_sub);
fprintf('Concluido.\n\n');

alpha_A = 4.343*2*k0*abs(imag(neff_A))*1e4;   % perda [dB/cm]
alpha_B = 4.343*2*k0*abs(imag(neff_B))*1e4;
base    = alpha_A(1);                          % t=0: guia sem espelho

%%
% Comparação com a teoria
%{
Condição de Fabry-Perot no espelho, incluindo a fase adquirida nas
reflexões internas contra as regiões evanescentes (análoga ao
deslocamento de fase da reflexão total):
  ressonância:      kx*t = m*pi + 2*atan(gamma/kx),  m = 0, 1, 2, ...
  antirressonância: kx*t = m*pi + 2*atan(gamma/kx) + pi/2
com kx = k0*sqrt(nSi^2 - nef^2) e gamma = k0*sqrt(nef^2 - ncl^2).
Observe que existe ressonância de ordem m=0.
%}
nmed = median(real(neff_A));
kx   = k0*sqrt(real(n_mir_A)^2 - nmed^2);
gam  = k0*sqrt(nmed^2 - n_cl^2);
phi  = atan(gam/kx);
t_res_teo = ((0:4)*pi + 2*phi)/kx;
t_AR_teo  = t_res_teo + (pi/2)/kx;

fprintf('Guia sem espelho (t=0): alpha = %.4g dB/cm\n', base);
fprintf('Teoria: kx=%.3f, gamma=%.3f, fase de interface = %.3f rad\n', kx, gam, phi);
fprintf('  t_res = %s um\n', sprintf('%.4f ', t_res_teo(t_res_teo<=ttA(end))));
fprintf('  t_AR  = %s um\n', sprintf('%.4f ', t_AR_teo(t_AR_teo<=ttA(end))));
fprintf('  espacamento entre picos: pi/kx = %.4f um\n\n', pi/kx);

lista_extremos('Espelho sem perda:', ttA, alpha_A, base);
lista_extremos('Espelho com perda (aula):', ttB, alpha_B, base);

%% ------------------------------------------------------------------------
% Gráfico: perda vs espessura do espelho
% -------------------------------------------------------------------------
figure('Name','ARROW: perda vs espessura do espelho','Position',[80 80 1000 480]);
semilogy(ttA, alpha_A, 'k-','LineWidth',1.4); hold on;
semilogy(ttB, alpha_B, 'Color',[.85 .4 0],'LineWidth',1.2);
yline(base,'b--','guia sem espelho','LineWidth',1.1);
for m = 1:numel(t_res_teo)
    if t_res_teo(m)<=ttA(end), xline(t_res_teo(m),'r--',sprintf('res m=%d',m-1)); end
    if t_AR_teo(m) <=ttA(end), xline(t_AR_teo(m),'g-.',sprintf('AR m=%d',m-1)); end
end
xlabel('Espessura do espelho de Si, t [\mum]');
ylabel('Perda de propagacao \alpha [dB/cm]');
legend('espelho sem perda (3.6)','espelho da aula (3.6-0.2j)','Location','northeast');
grid on;
title(sprintf(['Perda do modo de G2 vs espessura do espelho (n_{ef} \\approx %.3f). ' ...
    'Operar nas antirressonancias (vales)'], nmed));

%% ------------------------------------------------------------------------
% Perfis de campo: antirressonância vs ressonância (espelho sem perda)
% -------------------------------------------------------------------------
ipk = extremos(alpha_A, +1, 3*base);
imn = extremos(alpha_A, -1, 0);
if ~isempty(imn), kAR = imn(1); else, [~,kAR] = min(abs(ttA-t_AR_teo(1)));  end
if ~isempty(ipk), kRS = ipk(1); else, [~,kRS] = min(abs(ttA-t_res_teo(1))); end

figure('Name','Perfis de campo','Position',[120 120 1000 460]);
casos = [kAR kRS];
tit = {sprintf('ANTIRRESSONANCIA t=%.3f \\mum, \\alpha=%.3g dB/cm', ttA(kAR), alpha_A(kAR)),...
       sprintf('RESSONANCIA t=%.4f \\mum, \\alpha=%.3g dB/cm',     ttA(kRS), alpha_A(kRS))};
for c = 1:2
    kc = casos(c);
    [nbt,dbt] = pilha(ttA(kc), n_cl, n_mir_A, n_G2, dc1, dc2, d2);
    [xv,Ey]   = perfil_bi(neff_A(kc), k0, nbt, dbt, n_sup, n_sub);
    Ey = abs(Ey)/max(abs(Ey));
    subplot(1,2,c);
    xb = [0 cumsum(dbt)];
    patch([xb(2) xb(3) xb(3) xb(2)],[0 0 1.15 1.15],[1 .8 .8],'FaceAlpha',.5,'EdgeColor','none'); hold on;
    patch([xb(4) xb(5) xb(5) xb(4)],[0 0 1.15 1.15],[.8 1 .8],'FaceAlpha',.5,'EdgeColor','none');
    plot(xv, Ey, 'k-','LineWidth',1.5);
    text(mean(xb(2:3)),1.08,'Si','HorizontalAlignment','center','FontSize',8);
    text(mean(xb(4:5)),1.08,'G2','HorizontalAlignment','center','FontSize',8);
    text(-0.4,1.08,'substrato','HorizontalAlignment','center','FontSize',8);
    xlabel('x [\mum]'); ylabel('|E_y| normalizado'); title(tit{c},'FontSize',9);
    grid on; ylim([0 1.15]); xlim([xv(1) xv(end)]);
end
sgtitle('Na antirressonancia o campo fica confinado em G2; na ressonancia invade o espelho e vaza');

%{

                        ANÁLISE DOS RESULTADOS

1) Sem o espelho (t=0), o modo de G2 é leaky com alpha ~ 1.4 dB/cm: é a
   fuga por tunelamento através da casca de 0.70 um para o substrato de
   Si.

2) Ao inserir o espelho e varrer t, aparecem picos estreitos de perda em
   t = 0.142, 0.455, 0.767 e 1.080 um, espaçados de pi/kx = 0.313 um,
   exatamente nas espessuras ressonantes previstas por
   kx*t = m*pi + 2*atan(gamma/kx). A fase 2*atan(gamma/kx) aparece porque
   o espelho está entre regiões evanescentes; é ela que explica a
   existência do pico de ordem m=0 e o deslocamento em relação à condição
   "ingênua" de quarto de onda.

3) A largura dos picos é uma propriedade física (o Q do tunelamento
   ressonante), não uma questão de resolução numérica. Com o espelho sem
   perda os picos são agulhas de 5-10 nm de largura e contraste de até
   ~1000x sobre o guia sem espelho. Com o índice usado em aula (k = 0.2),
   a absorção do campo residual dentro do Si alarga os picos (>= 60 nm),
   "lava" as ordens altas e impõe um piso nos vales (~38-94 dB/cm). Ou
   seja, a qualidade do material do espelho controla diretamente a
   nitidez das ressonâncias.

4) Nos vales (antirressonâncias), a perda com espelho sem perda volta a
   ~0.98x a do guia sem espelho: o espelho fica "opticamente invisível".
   Este é o ponto de operação do ARROW. Na prática, é isso que permite
   ter a camada de Si próxima do guia (necessária, por exemplo, para o
   detector MSM discutido em aula) sem custo de perda - desde que a
   espessura seja antirressonante. Na espessura errada (ressonante), a
   mesma camada custa até 3 ordens de grandeza a mais.

5) Os perfis de campo confirmam o mecanismo: na antirressonância o campo
   fica confinado em G2, com penetração mínima no espelho; na ressonância
   o campo invade o espelho com amplitude comparável à do núcleo, abrindo
   o canal de fuga para o substrato.

6) Critérios de projeto que saem da simulação: usar a MENOR espessura
   antirressonante (m=0), pois é a mais fina e fácil de fabricar e, com
   espelho absorvente, a que absorve menos (os vales em si se repetem com
   período pi/kx, então a tolerância em t é parecida em todas as ordens);
   e o material do
   espelho deve ter índice bem diferente do meio em torno (aqui 3.6
   contra 1.5), como comentado em aula, para que o contraste de
   impedância produza reflexão forte.
%}

%% 
% Funções locais

function [tt, neff_c] = varre(t_coarse, w_fino, N_fino, neff_iso, k0, ...
                              n_cl, n_mir, n_G2, dc1, dc2, d2, n_sup, n_sub)
% Varredura grossa em t + refinamento automático em torno dos picos.
    Nt = numel(t_coarse);
    r  = zeros(1,Nt);  imh = -1e-4;
    for k = 1:Nt
        r(k) = resolve_um(t_coarse(k), imh, neff_iso, k0, n_cl, n_mir, ...
                          n_G2, dc1, dc2, d2, n_sup, n_sub);
        imh  = imag(r(k));
        if mod(k,150)==0, fprintf('   %d/%d\n', k, Nt); end
    end
    al  = abs(imag(r));
    ipk = extremos(al, +1, 3*al(1)+1e-30);
    tt  = t_coarse;  neff_c = r;
    for p = 1:numel(ipk)
        tp = t_coarse(ipk(p));
        tf = linspace(max(tp-w_fino,1e-6), tp+w_fino, N_fino);
        rf = zeros(1,N_fino);  imh = imag(r(ipk(p)));
        for k = 1:N_fino
            rf(k) = resolve_um(tf(k), imh, neff_iso, k0, n_cl, n_mir, ...
                               n_G2, dc1, dc2, d2, n_sup, n_sub);
            imh = imag(rf(k));
        end
        tt = [tt, tf];  neff_c = [neff_c, rf]; %#ok<AGROW>
    end
    [tt, iu] = unique(round(tt*1e7)/1e7);
    neff_c   = neff_c(iu);
end

function raiz = resolve_um(t, im_hint, neff_iso, k0, n_cl, n_mir, n_G2, ...
                           dc1, dc2, d2, n_sup, n_sub)
% Resolve o modo de G2 para um único t: varredura de n_ef real em janela
% fixa, refino complexo de cada candidato e seleção pelo confinamento.
    [nbt,dbt] = pilha(t, n_cl, n_mir, n_G2, dc1, dc2, d2);
    xb = [0 cumsum(dbt)];
    ng = linspace(neff_iso-0.10, neff_iso+0.10, 1501);
    Fa = abs(Fdisp(ng + 1j*im_hint, k0, nbt, dbt, n_sup, n_sub));
    ii = find(Fa(2:end-1)<Fa(1:end-2) & Fa(2:end-1)<Fa(3:end)) + 1;
    if isempty(ii)
        raiz = refina(neff_iso, im_hint, k0, nbt, dbt, n_sup, n_sub);
        return
    end
    rts = [];
    for c = ng(ii)
        rr = refina(c, im_hint, k0, nbt, dbt, n_sup, n_sub);
        if isempty(rts) || all(abs(rr-rts) > 1e-6)
            rts(end+1) = rr; %#ok<AGROW>
        end
    end
    if isscalar(rts)
        raiz = rts(1);
    else
        fr = zeros(size(rts));
        for q = 1:numel(rts)
            [xv,Ey] = perfil_bi(rts(q), k0, nbt, dbt, n_sup, n_sub);
            fr(q)   = confinamento(xv, Ey, xb(4), xb(5));
        end
        [~,qb] = max(fr);
        raiz   = rts(qb);
    end
end

function idx = extremos(y, sinal, limiar)
% Máximos (sinal=+1) ou mínimos (sinal=-1) locais.
    z   = sinal*y;
    idx = find(z(2:end-1)>z(1:end-2) & z(2:end-1)>z(3:end)) + 1;
    if sinal > 0, idx = idx(y(idx) > limiar); end
end

function lista_extremos(cab, tt, al, base)
    ipk = extremos(al, +1, 3*base);
    imn = extremos(al, -1, 0);
    fprintf('%s\n', cab);
    for q = 1:numel(ipk)
        fprintf('  pico t=%.4f um: alpha=%.5g dB/cm, FWHM=%.1f nm (%.0fx o guia s/ espelho)\n',...
            tt(ipk(q)), al(ipk(q)), 1e3*fwhm(tt,al,ipk(q)), al(ipk(q))/base);
    end
    for q = 1:numel(imn)
        fprintf('  vale t=%.4f um: alpha=%.4g dB/cm (%.2gx o guia s/ espelho)\n',...
            tt(imn(q)), al(imn(q)), al(imn(q))/base);
    end
    fprintf('\n');
end

function L = fwhm(tv, al, i)
% Largura a meia altura do pico em al(i).
    half = al(i)/2;  l = i;  r = i;
    while l>1 && al(l)>half,          l = l-1; end
    while r<numel(al) && al(r)>half,  r = r+1; end
    L = tv(r) - tv(l);
end

function [nbt, dbt] = pilha(t, n_cl, n_mir, n_G2, dc1, dc2, d2)
% Camadas finitas, do substrato para o superstrato.
    nbt = [n_cl, n_mir, n_cl, n_G2];
    dbt = [dc2,  t,    dc1,  d2  ];
end

function F = Fdisp(neff, k0, n_bt, d_bt, n_sup, n_sub)
%{
Função de dispersão TE via TMM (aceita vetor de n_ef). F = 0 no modo.
Condições de contorno: onda de saída no substrato (modo leaky) e
decaimento evanescente no superstrato. O vetor de estado [Ey; Ey'] é
renormalizado a cada camada para evitar estouro numérico.
%}
    b2  = (k0*neff).^2;
    kxs = sqrt((k0*n_sub)^2 - b2);
    fl  = real(kxs)<0; kxs(fl) = -kxs(fl);
    V1  = ones(size(b2));  V2 = 1j*kxs;
    for m = 1:numel(n_bt)
        d = d_bt(m);
        if d > 0
            kx = sqrt((k0*n_bt(m))^2 - b2);
            kx(abs(kx)<1e-12) = 1e-12;
            c = cos(kx*d);  s = sin(kx*d);
            W1 =  c.*V1 + (s./kx).*V2;
            W2 = -kx.*s.*V1 + c.*V2;
            nn = max(abs(W1),abs(W2)) + 1e-300;
            V1 = W1./nn;  V2 = W2./nn;
        end
    end
    gs = sqrt(b2 - (k0*n_sup)^2);
    fl = real(gs)<0; gs(fl) = -gs(fl);
    F  = (V2 + gs.*V1)./(abs(V1)+abs(V2)+1e-300);
end

function raiz = refina(cand, im0, k0, n_bt, d_bt, n_sup, n_sub)
%{
Refino local da raiz no plano complexo, parametrizado como desvio do
candidato; a penalidade mantém o refino na vizinhança (evita convergir
para a raiz de outro modo).
%}
    opts = optimset('TolX',1e-12,'TolFun',1e-12,...
                    'MaxFunEvals',3000,'MaxIter',3000,'Display','off');
    J = @(v) 1e6*max(0,abs(v(1))-0.03).^2 + ...
        log10(abs(Fdisp((cand+v(1))+1j*(im0+v(2)), k0, n_bt, d_bt, n_sup, n_sub))+1e-300);
    v = fminsearch(J, [0,0], opts);
    raiz = (cand+v(1)) + 1j*(im0+v(2));
end

function fr = confinamento(xv, Ey, x1, x2)
% Fração de |Ey|^2 contida no intervalo [x1,x2].
    I  = abs(Ey).^2;  ms = xv>=x1 & xv<=x2;
    fr = trapz(xv(ms), I(ms)) / trapz(xv, I);
end

function [xv, Ey] = perfil_bi(neff, k0, n_bt, d_bt, n_sup, n_sub)
%{
Perfil transversal Ey(x). A rotina permite integrar das duas extremidades
e casar as soluções na interface i_match (evita a instabilidade numérica
de propagar numa única direção através de regiões evanescentes espessas).
Com i_match = 4 (topo do núcleo), todas as camadas finitas são integradas
a partir do substrato e só a cauda evanescente do superstrato é anexada;
com as cascas usadas aqui (<= 0.45 um) isso é numericamente estável.
%}
    b2 = (k0*neff)^2;  Ns = 200;  ext = 0.8;
    i_match = min(4, numel(n_bt));
    kxs = sqrt((k0*n_sub)^2 - b2);
    if real(kxs)<0, kxs = -kxs; end
    xs = linspace(-ext,0,Ns);
    xv = xs;  Ey = exp(1j*kxs*xs);
    Vv = [1; 1j*kxs];  xcur = 0;
    for m = 1:i_match
        d = d_bt(m);
        if d > 0
            kx = sqrt((k0*n_bt(m))^2 - b2);
            if abs(kx)<1e-12, kx = 1e-12; end
            xi = linspace(0,d,Ns);
            Ey = [Ey, cos(kx*xi)*Vv(1) + sin(kx*xi)/kx*Vv(2)]; %#ok<AGROW>
            xv = [xv, xcur+xi]; %#ok<AGROW>
            M  = [cos(kx*d) sin(kx*d)/kx; -kx*sin(kx*d) cos(kx*d)];
            Vv = M*Vv;  xcur = xcur+d;
        end
    end
    E_low = Vv(1);
    gs = sqrt(b2 - (k0*n_sup)^2);
    if real(gs)<0, gs = -gs; end
    Vt = [1; -gs];  xtop = sum(d_bt);
    xv_up = [];  Ey_up = [];  xcu = xtop;
    for m = numel(n_bt):-1:i_match+1
        d = d_bt(m);
        if d > 0
            kx = sqrt((k0*n_bt(m))^2 - b2);
            if abs(kx)<1e-12, kx = 1e-12; end
            Mi = [cos(kx*d) -sin(kx*d)/kx; kx*sin(kx*d) cos(kx*d)];
            Vb = Mi*Vt;
            xi = linspace(0,d,Ns);
            Ey_up = [cos(kx*xi)*Vb(1) + sin(kx*xi)/kx*Vb(2), Ey_up]; %#ok<AGROW>
            xv_up = [(xcu-d)+xi, xv_up]; %#ok<AGROW>
            Vt = Vb;  xcu = xcu-d;
        end
    end
    E_up = Vt(1);
    xi  = linspace(0,ext,Ns);
    esc = E_low/(E_up+1e-300);
    xv  = [xv, xv_up, xtop+xi];
    Ey  = [Ey, esc*Ey_up, esc*exp(-gs*xi)];
end