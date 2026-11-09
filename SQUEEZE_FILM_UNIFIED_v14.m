function SQUEEZE_FILM_UNIFIED_v14(mode)
%% ========================================================================
%  SQUEEZE_FILM_UNIFIED_v14.m
%  v14 CHANGES (plotting and data only; the physics is unchanged from v13):
%    - Fig 5 normalisation bug fixed (each step now divided by the isothermal load at the same h)
%    - all figures in Times New Roman, journal widths, 300 dpi PNG + vector PDF
%    - Fig 1 widened to the full column width, no empty margins
%    - Figs 3, 4, 5, 7 show several curves (sweeps), Fig 2 shows first- and total-order indices
%    - new Part 16 one-at-a-time sweeps: roughness, slip, porosity, couple-stress length, viscosity
%      closures incl. Mooney psi, kinetic aggregation, constant-load squeeze, shear thinning (5 lambdas)
%    - every plotted array saved to CSV
%  ONE coupled model for the circular squeeze film, replacing the superposed
%  treatment of v12.  Function file written to run in MATLAB and Octave
%  (tested in Octave only).
%
%    SQUEEZE_FILM_UNIFIED_v14          full run
%    SQUEEZE_FILM_UNIFIED_v14('quick') reduced grids, for a smoke test
%
%  UNIFIED MODEL (axisymmetric; p(a)=0; Q = U r/2 per unit circumference)
%    dp/dr = -6 U r mu(p,T,phi) / ( h^3 F )                     [Reynolds]
%    F     = [ fR(h_l) * Phi_gap(h_l; ell, b1, b2) + 12 kappa_e h_p ] / h^3 ,  h_l = h + compl*p (elastic liner)
%    Phi_gap = couple stress + Navier/BJS slip solved TOGETHER (closed form,
%              derived below and checked against a direct BVP), not added
%    mu    = mu0 * rmu(phi_a) * exp(-beta*theta) * m(p)   m: Barus or bounded
%    ell^2 = eta/mu  (eta = mu0 ell0^2 (1+gc*phi))          [couple stress]
%    theta : quasi-steady energy balance with dissipation from the SAME p
%    phi   : depth-averaged particle transport, Brownian + shear-induced
%            migration (Phillips et al. 1992), advected by Q = U r/2
%    phi_a = phi * Ra^(3-D) , Ra from Smoluchowski kinetics or fixed
%  Solved by a fixed-point loop p -> theta -> phi -> (mu, ell) -> p.
%
%  CLOSURES AND SOURCES (found by web search; see SRC table printed at start)
%    - modified Krieger-Dougherty with aggregates (Chen et al. 2007):
%        mu_r = (1 - phi_a/phi_m)^(-2.5 phi_m), phi_a = phi (a_a/a)^(3-D),
%        D = 1.8, phi_m = 0.605, a_a/a = 3.34 (typical)
%    - couple-stress length ell ~ particle/aggregate size
%        (Binu et al. 2015; Awad et al. 2023): ell0 = d_agg
%    - Mooney self-crowding factor psi_M in [1.35, 1.91] (as stated in the paper)
%    - shear-induced migration constants K_c = 0.41, K_mu = 0.62
%        (Phillips et al. 1992, as quoted by later papers)
%    - roughness: Patir & Cheng (1978) isotropic: 1 - 0.9 exp(-0.56 h/sigma)
%    - slip length: ~100 nm steel/PAO/steel (Corneli, TU Darmstadt);
%        ~20 nm smooth hydrophobic (Cottin-Bizonne et al. 2005); 0 hydrophilic
%    - permeability of sintered bronze bearings 0.6-60e-3 darcy
%        (Watanabe et al., J. Japan Soc. Powder Metall.): 6e-16 .. 6e-14 m^2
%    - bounded viscosity: mu_max = glass-transition viscosity ~ 1e10 Pa s
%        (organic liquids / lubricants, Bair 2019); 1e12 for minerals
%
%  PARTS
%    0 sources | 1 viscosity ratios | 2, 2b verification | 3 coupled solve at the sourced operating point
%    4 each effect alone vs all together | 5 existence limit vs 5 MPa envelope | 6 2-D plates (v_c by Lambda and fluid)
%    7 Sobol (direct and surrogate) | 8 comparison with Lin et al. (2013) and Shiralashetti & Vatsala (2026)
%    9 TiO2/Al2O3/CuO (Abass et al. 2023) | 10 constant-load squeeze | 11 elastic liner | 12 separation + cavitation
%    13 transient thermal + radial conduction | 14 shear-thinning sensitivity | 15 figures
%
%  WHAT IS STILL NOT KNOWN (no measurements exist in the sources found)
%    colloidal stability ratio Ws, eta(phi), roughness sigma of the actual
%    plates, wall thickness h_p and saturation of the porous layer, real
%    particle properties.  They enter as ranges and are propagated by the
%    Sobol study (Part 7).  Nothing here is validated against experiment.
%
%  TESTED IN OCTAVE ONLY (not in MATLAB).  Checks are printed as PASS/FAIL.
%% ========================================================================
if nargin < 1, mode = 'full'; end
quick = strcmp(mode,'quick');
tAll = tic;
try, rng(1); catch, rand('twister',1); end %#ok<NOCOMMA>
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
outDir = fullfile(here,'results_unified_v14');
if ~exist(outDir,'dir'), [okm,~] = mkdir(outDir); if ~okm, outDir = fullfile(tempdir,'results_unified_v14'); mkdir(outDir); end, end
logf = fullfile(outDir,'console_log.txt'); if exist(logf,'file'), delete(logf); end
diary(logf); cleanupDiary = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('Output folder: %s\n', outDir);
TEST = {}; FIG = struct();

S0 = defaultParams();
if quick, S0.N = 80; end
if mod(S0.N,2)==1, S0.N = S0.N+1; end

%% ---- PART 0: parameter and source table --------------------------------
fprintf('\n=== PART 0: PARAMETERS AND SOURCES ===\n');
SRC = {
 'mu0, alpha, rho, cp, kth, beta, hw','SAE 30 base oil, Jabbar et al. 2026 / Ferron 1983 (as in v12)','sourced (v12)';
 'phi_m = 0.605, D = 1.8, a_a/a = 3.34','Chen et al. 2007 modified Krieger-Dougherty; Binu et al.; arXiv 2604.05596','literature';
 'psi_M in [1.35, 1.91]','Mooney bounds quoted in the paper (v12 used 1.1: outside)','literature';
 'ell0 = d_agg (couple-stress length)','Binu et al. 2015; Awad et al. 2023 (size as couple-stress parameter)','literature closure';
 'K_c = 0.41, K_mu = 0.62','Phillips et al. 1992 diffusive-flux model','literature; authors say hard to measure';
 'Patir-Cheng C = 0.9, r = 0.56','Patir & Cheng 1978 isotropic','literature';
 'slip b = 0 .. 100 nm','Corneli (steel-PAO ~100 nm); Cottin-Bizonne 2005 (20 nm)','literature range';
 'kappa = 6e-16 .. 6e-14 m^2','Watanabe et al., sintered bronze 0.6-60e-3 darcy','literature range';
 'mu_max = 1e10 Pa s','Bair 2019: glass-transition viscosity of organic lubricants','literature';
 'sigma (roughness rms), h_p, saturation S','typical engineering values, NOT sourced here','ASSUMED range';
 'Ws (colloidal stability ratio)','no source: unknown for this oil/particle pair','UNKNOWN -> sweep';
 'rho_p, cp_p, k_p (TiO2 default)','Abass et al. 2023 Table 2 (TiO2 4230, 710, 8.4; Al2O3, CuO in Part 9)','literature';
 'experimental squeeze-film data','none found for nanofluid squeeze films','NONE'};
for i = 1:size(SRC,1), fprintf('  %-38s | %-62s | %s\n', SRC{i,:}); end

%% ---- PART 1: viscosity ratios (limitation 1: psi_M) ---------------------
fprintf('\n=== PART 1: RELATIVE VISCOSITY, phi = 0.05 ===\n');
Sv = S0; Sv.phi = 0.05; rowsV = {};
mods = {'Einstein','Brinkman','Batchelor','Krieger','Nielsen'};
for i = 1:numel(mods), Sv.visc = mods{i}; rowsV(end+1,:) = {mods{i}, relVisc(Sv,Sv.phi)}; end %#ok<AGROW>
for psi = [1.1 1.35 1.63 1.91]
    Sv.visc = 'Mooney'; Sv.psiM = psi;
    rowsV(end+1,:) = {sprintf('Mooney psi=%.2f%s',psi,tern(psi<1.35,' (v12, outside bounds)','')), relVisc(Sv,Sv.phi)}; %#ok<AGROW>
end
for ra = [1 2 3.34]
    Sv.visc = 'KDagg'; Sv.RaEff = ra; Sv.Ra = ra;
    rowsV(end+1,:) = {sprintf('mod. KD aggregates Ra=%.2f',ra), relVisc(Sv,Sv.phi)}; %#ok<AGROW>
end
for i = 1:size(rowsV,1), fprintf('  %-34s r = %.4f\n', rowsV{i,1}, rowsV{i,2}); end
wcsvc(fullfile(outDir,'table2_viscosity_ratios.csv'),'model,ratio',rowsV);
fprintf('  NOTE: Nielsen here uses exp(phi/(1-phi_m)) as in the literature review; v12 used exp(phi/(1-phi)). Check against the original.\n');

%% ---- PART 2: VERIFICATION OF THE UNIFIED SOLVER --------------------------
fprintf('\n=== PART 2: VERIFICATION ===\n');
% unit scales: mu0=U=a=h0=1, alpha=v, ell0=k (reproduces v12 Part 1 numbers)
Su = unitParams();
% T1 isoviscous
Su.alpha = 0; Su.ellFixed = 0; o = solveCoupled(Su, 0.7);
TEST(end+1,:) = {'isoviscous: W = 1.5 pi mu U a^4/h^3 (rel<1e-6)', abs(o.W-1.5*pi/0.7^3)/(1.5*pi/0.7^3) < 1e-6};
TEST(end+1,:) = {'isoviscous: p(0) = 3 mu U a^2/h^3 (rel<1e-9)', abs(o.P0-3/0.7^3)/(3/0.7^3) < 1e-9};
% T2 Barus closed form, no couple stress
Su = unitParams(); Su.ellFixed = 0; vv = 0.1; Su.alpha = vv; h1 = 0.8;
o = solveCoupled(Su,h1); P0cf = -log(1-3*vv/h1^3)/vv;
TEST(end+1,:) = {'Barus closed form p(0) = -ln(1-3 v/h^3)/v (rel<1e-6)', abs(o.P0-P0cf)/P0cf < 1e-6};
vc = vcritAxi(Su,h1,1e-9);
TEST(end+1,:) = {'Barus existence limit v_c = h^3/3 (rel<1e-4)', abs(vc-h1^3/3)/(h1^3/3) < 1e-4};
% T3 reproduce v12 Part 1: k=0.10, H*=0.4, r=1
Su = unitParams(); Su.ellFixed = 0.10; Su.alpha = 0; o = solveCoupled(Su,0.4);
TEST(end+1,:) = {sprintf('v12 reproduction v=0: P*=%.3f (76.654), W*=%.3f (120.409)',o.P0,o.W), abs(o.P0-76.654)<2e-3 && abs(o.W-120.409)<2e-3};
vcExact = numClosed(0.4,0.10)/3;
Su.alpha = 0.5*vcExact; o = solveCoupled(Su,0.4);
TEST(end+1,:) = {sprintf('v12 reproduction v=0.5 v_c: P*=%.3f (99.880), W*=%.3f (141.850)',o.P0,o.W), abs(o.P0-99.880)<5e-3 && abs(o.W-141.850)<5e-3};
Su.visc = 'fixed'; Su.rVal = 1.5; Su.alpha = 0; o = solveCoupled(Su,0.4);
TEST(end+1,:) = {'nanofluid r multiplies the isoviscous pressure exactly (p = 3r/Phi)', abs(o.P0-1.5*76.654)/ (1.5*76.654) < 1e-4};
% T4 flow coefficient: couple stress + slip solved together
okf = true;
okf = okf && abs(flowGap(0.7,0.07,0,0) - (0.7^3-12*0.07^2*0.7+24*0.07^3*tanh(0.7/(2*0.07)))) < 1e-12;
okf = okf && abs(flowGap(0.7,0,0.2,0.3) - 0.7^2*(0.7^2+4*0.7*0.5+12*0.06)/(0.7+0.5)) < 1e-12;
TEST(end+1,:) = {'flowGap: b=0 -> couple-stress Phi; ell=0 -> two-wall slip Phi', okf};
errBVP = 0;
for cs = [0.7 0.07 0.0 0.2; 1.0 0.3 0.1 0.0; 0.5 0.2 0.05 0.15; 0.4 0.05 0.02 0.02]'
    h = cs(1); ell = cs(2); b1 = cs(3); b2 = cs(4);
    errBVP = max(errBVP, abs(flowGap(h,ell,b1,b2)-flowBVP(h,ell,b1,b2))/flowGap(h,ell,b1,b2));
end
TEST(end+1,:) = {sprintf('flowGap vs direct finite-difference BVP, couple stress + slip (max rel %.1e < 1e-4)',errBVP), errBVP < 1e-4};
TEST(end+1,:) = {'flowGap series branch continuous at large ell/h', abs(flowGap(1,50,0.1,0.1)-flowGap(1,49.999,0.1,0.1))/flowGap(1,50,0.1,0.1) < 1e-3 && flowGap(1,50,0,0) > 0};
% T5 bounded law
Su = unitParams(); Su.ellFixed = 0; Su.alpha = 0.5; h1 = 0.8; Su.law = 'barus'; ob = solveCoupled(Su,h1);
Su.law = 'bounded'; Su.muMax = 1e6; oc = solveCoupled(Su,h1);
TEST(end+1,:) = {'Barus beyond v_c is flagged non-existent (v=0.5 > h^3/3)', ~ob.exist};
TEST(end+1,:) = {'bounded law solves the same case (finite p and W)', oc.exist && isfinite(oc.W)};
Su.alpha = 1e-4; Su.law = 'bounded'; o1 = solveCoupled(Su,h1); Su.law = 'barus'; o2 = solveCoupled(Su,h1);
TEST(end+1,:) = {'bounded law -> Barus for small v (rel<1e-6)', abs(o1.W-o2.W)/o2.W < 1e-6};

%% ---- PART 2b: EXTENDED VERIFICATION (added after the line-by-line code audit) ----------
fprintf('\n=== PART 2b: EXTENDED VERIFICATION ===\n'); flushOut();
% quadrature and differentiation matrices
[Dq,xq] = chebD(12); wq = clencurt(12); [~,xq13] = chebD(13); wq13 = clencurt(13);
TEST(end+1,:) = {'Clenshaw-Curtis weights integrate x^2 and x^4 exactly (even and odd N)', abs(wq*xq.^2-2/3) < 1e-13 && abs(wq13*xq13.^4-2/5) < 1e-13 && max(abs(Dq*(xq.^3)-3*xq.^2)) < 1e-12};
% RK4 pressure integration is 4th order (Barus closed form, ell = 0)
Su = unitParams(); Su.ellFixed = 0; Su.alpha = 0.1; pex = -log(1-3*0.1/0.8^3)/0.1; er = zeros(1,3); Nl = [50 100 200];
for i = 1:3, Su.N = Nl(i); er(i) = abs(solveCoupled(Su,0.8).P0-pex); end
ordr = log2(er(1:2)./er(2:3));
TEST(end+1,:) = {sprintf('pressure ODE: RK4 convergence order %.2f, %.2f (> 3.8)',ordr(1),ordr(2)), all(ordr > 3.8)};
% elastic liner: isoviscous closed form p = ((h^4 + 12 c mu U (a^2-r^2))^(1/4) - h)/c
Sc2 = defaultParams(); Sc2.visc = 'none'; Sc2.alpha = 0; Sc2.ellSource = 'fixed'; Sc2.ellFixed = 0; Sc2.N = 200; Sc2.U = 1e-3; Sc2.compl = 1e-11; hc2 = 20e-6;
oc2 = solveCoupled(Sc2,hc2); pc2 = ((hc2^4+12*Sc2.compl*Sc2.mu0*Sc2.U*(Sc2.a^2-oc2.r.^2)).^0.25-hc2)/Sc2.compl;
TEST(end+1,:) = {sprintf('elastic liner: pressure equals the isoviscous closed form (max rel %.1e < 1e-6)',max(abs(oc2.p-pc2))/max(pc2)), max(abs(oc2.p-pc2))/max(pc2) < 1e-6};
% quasi-steady energy balance: residual of the ODE, for fast and very slow squeezes
okT = true; Sth = defaultParams(); Sth.visc = 'KDagg'; Sth.law = 'bounded'; Sth.ellSource = 'fixed'; Sth.ellFixed = 0; Sth.N = 200; Nt = Sth.N; rt = linspace(0,Sth.a,Nt+1)';
for Ut = [1e-6 1e-4 1e-3]
    Sth.U = Ut; ht = 0.4*Sth.h0; Mt = relVisc(Sth,Sth.phi)*ones(Nt+1,1); [~,auxt] = pressureODE(Sth,ht,Mt,zeros(Nt+1,1),0,0,0);
    tht = thermalField(Sth,rt,ht,auxt,Sth.phi); rct = ((1-Sth.phi)*Sth.rho*Sth.cp+Sth.phi*Sth.rhoP*Sth.cpP)*Ut; Dt = ht^3*auxt.F.*auxt.dp.^2./(12*auxt.mu);
    rest = rct*rt.*gradient(tht,rt(2)-rt(1))/2 + 2*Sth.hw*tht - Dt; okT = okT && all(isfinite(tht)) && max(abs(rest(5:end-5)))/max(Dt) < 1e-4;
end
TEST(end+1,:) = {'thermal field satisfies rho cp (U r/2) theta'' = D - 2 h_w theta for k from 0.2 to 180 (residual < 1e-4 D_max)', okT};
% couple stress + slip flux in the large-ell regime (series branch) against the independent BVP
eB = 0; for cs = [0.2 1.0 0.01 0.0; 0.05 1.0 0.002 0.0; 0.1 0.3 0.0 0.05]', eB = max(eB, abs(flowGap(cs(1),cs(2),cs(3),cs(4))/flowBVP(cs(1),cs(2),cs(3),cs(4))-1)); end
TEST(end+1,:) = {sprintf('flowGap vs BVP for ell/h up to 20 (max rel %.1e < 1e-5)',eB), eB < 1e-5};
% aggregated Krieger-Dougherty value by hand
Sk0 = defaultParams(); Sk0.RaEff = 3.34; paH = 0.05*3.34^1.2; rH = (1-paH/0.605)^(-2.5*0.605);
TEST(end+1,:) = {sprintf('modified Krieger-Dougherty with aggregates = hand value %.5f',rH), abs(relVisc(Sk0,0.05)-rH) < 1e-12};
% Sobol estimators on the Ishigami function (known S1 = 0.314 0.442 0; ST = 0.558 0.442 0.244)
try, rng(7); catch, rand('twister',7); end %#ok<NOCOMMA>
fI = @(X) sin(X(:,1))+7*sin(X(:,2)).^2+0.1*X(:,3).^4.*sin(X(:,1)); nI = 100000; AI = rand(nI,3)*2*pi-pi; BI = rand(nI,3)*2*pi-pi; YAI = fI(AI); YBI = fI(BI); VI = var([YAI;YBI]); S1I = zeros(1,3); STI = zeros(1,3);
for i = 1:3, ABi = AI; ABi(:,i) = BI(:,i); Yi = fI(ABi); S1I(i) = mean(YBI.*(Yi-YAI))/VI; STI(i) = 0.5*mean((YAI-Yi).^2)/VI; end
TEST(end+1,:) = {sprintf('Sobol estimators reproduce Ishigami (S1 %.3f %.3f %.3f, ST %.3f %.3f %.3f)',S1I,STI), all(abs(S1I-[0.314 0.442 0]) < 0.02) && all(abs(STI-[0.558 0.442 0.244]) < 0.03)};
% quadratic surrogate is exact for a quadratic function
Zq = rand(300,12)*2-1; Yqq = 1+Zq(:,1)-2*Zq(:,2).*Zq(:,3)+Zq(:,4).^2; Cqq = fitQuad(Zq,Yqq); Ztq = rand(50,12)*2-1;
TEST(end+1,:) = {'quadratic surrogate reproduces a quadratic test function (error < 1e-6)', max(abs(evalQuad(Ztq,Cqq,12)-(1+Ztq(:,1)-2*Ztq(:,2).*Ztq(:,3)+Ztq(:,4).^2))) < 1e-6};
% particle transport conserves the particle amount (diffusion only; sub-micron advection)
Sg2 = defaultParams(); Sg2.visc = 'none'; Sg2.mu0 = 1e-9; Sg2.beta = 0; Sg2.U = 1e-12; Sg2.hStart = 20e-6 + 1e-12; Ng = 200; rg = linspace(0,Sg2.a,Ng+1)';
phg0 = 0.05*(1+0.5*cos(pi*rg/Sg2.a)); phg = migrate(Sg2,rg,20e-6,phg0,1,1,zeros(Ng+1,1)); m0 = simpson(rg,phg0.*rg); m1 = simpson(rg,phg.*rg);
TEST(end+1,:) = {sprintf('particle transport: amount conserved (rel %.1e < 1e-5) and gradients relax',(m1-m0)/m0), abs(m1-m0)/m0 < 1e-5 && (max(phg)-min(phg)) < (max(phg0)-min(phg0))};
% 2-D: Clenshaw-Curtis area of the plate
sq2 = setup2D(12,2); TEST(end+1,:) = {'2-D quadrature integrates the plate area (4)', abs(sq2.w*ones(13,13)*sq2.w'-4) < 1e-12};
try, rng(1); catch, rand('twister',1); end %#ok<NOCOMMA>

%% ---- PART 3: COUPLED PHYSICS ON THE SOURCED OPERATING POINT --------------
fprintf('\n=== PART 3: COUPLED SOLVE, SOURCED OPERATING POINT ===\n');
SA = S0; SA.visc = 'KDagg'; hA = 0.4*SA.h0;
SA.thermal = true; SA.migration = true; SA.slip = 100e-9; SA.sigma = 0.5e-6; SA.porous = true;
oA = solveCoupled(SA,hA); FIG.oA = oA;
fprintf('  h = %.1f um, U = %.1e m/s, Ra = %.2f, phi_a/phi_m = %.3f, mu_r = %.3f\n', hA*1e6, SA.U, oA.Ra, oA.phi_a_ratio, oA.mu_r);
fprintf('  ell0 = %.3g m (k = ell0/h = %.2e) | iterations %d | exist %d\n', oA.ell0, oA.ell0/hA, oA.iter, oA.exist);
fprintf('  P0 = %.4g Pa | W = %.4g N | max dT = %.3g K | max |dphi|/phi = %.2e\n', oA.P0, oA.W, max(oA.th), oA.dphiRel);
% diagnostics: why effects are small or large
kB = 1.380649e-23; ap = 0.5*SA.dp*oA.Ra; DB = kB*SA.T0/(6*pi*SA.mu0*ap);
gd = 3*SA.U*SA.a/hA^2; tsq = (SA.hStart-hA)/SA.U;
fprintf('  Peclet (radial, Brownian) = %.2e | shear-migration diffusivity / Brownian = %.2e\n', SA.U*SA.a*SA.a/(2*hA*DB)/1, SA.Kc*ap^2*gd*SA.phi/DB);
fprintf('  squeeze time %.3f s vs gap-diffusion time h^2/D_B = %.2e s\n', tsq, hA^2/DB);
fprintf('  slip/h = %.2e | sigma/h = %.2e | Darcy/h^3 = %.2e | v = alpha*mu0*U*a^2/h^3 = %.2e\n', SA.slip/hA, SA.sigma/hA, 12*SA.kappa*SA.Ssat*SA.hp/hA^3, SA.alpha*SA.mu0*SA.U*SA.a^2/hA^3);
% thermal and migration tests
St = S0; St.thermal = false; St.migration = false; St.visc = 'KDagg';
o0 = solveCoupled(St,hA); St.thermal = true; o1 = solveCoupled(St,hA);
TEST(end+1,:) = {'thermal: theta>=0 and W_coupled <= W_isothermal', all(o1.th>=-1e-14) && o1.W <= o0.W*(1+1e-12)};
Sx = S0; Sx.U = 0.02; Sx.alpha = 0; Sx.visc = 'KDagg'; Sx.thermal = false; hx = 0.4*Sx.h0;
ox0 = solveCoupled(Sx,hx); Sx.thermal = true; ox1 = solveCoupled(Sx,hx);
fprintf('  fast squeeze U = 2 cm/s: max dT = %.3f K, W change = %.2f %%\n', max(ox1.th), 100*(ox1.W/ox0.W-1));
TEST(end+1,:) = {'thermal: heating lowers load at fast squeeze', ox1.W < ox0.W && max(ox1.th) > 1e-3};
Sm = S0; Sm.visc = 'KDagg'; Sm.migration = true;
om = solveCoupled(Sm,hA);
TEST(end+1,:) = {sprintf('migration with sourced constants leaves phi uniform (max rel change %.1e < 1e-5)',om.dphiRel), om.dphiRel < 1e-5};
Sm.migAmp = 1e9; om2 = solveCoupled(Sm,hA);
TEST(end+1,:) = {sprintf('stress test (migration x1e9): phi moves toward the axis (phi(0)>phi(a), change %.2e)',om2.dphiRel), om2.phi(1) > om2.phi(end) && om2.dphiRel > 1e-4 && all(om2.phi>0 & om2.phi<Sm.phim)};
% aggregation kinetics
Sg = S0; Sg.aggMode = 'kinetic'; Sg.Ws = 1;
tgrid = [0 1e-3 1e-2 1e-1 1];
Rg = arrayfun(@(t) aggRatio(Sg,t), tgrid);
TEST(end+1,:) = {'aggregation: Ra(0)=1, monotone in t, bounded by RaMax', abs(Rg(1)-1)<1e-12 && all(diff(Rg)>=0) && max(Rg)<=Sg.RaMax+1e-12};
fprintf('  aggregation ratio Ra(t) for Ws=1: '); fprintf('%.2f ', Rg); fprintf('(t = 0, 1 ms, 10 ms, 0.1 s, 1 s)\n');
rowsW = []; for Ws = 10.^(0:6), Sg.Ws = Ws; rowsW(end+1,:) = [Ws aggRatio(Sg,0.3)]; end %#ok<AGROW>
fprintf('  Ra after 0.3 s squeeze vs stability ratio Ws:\n'); fprintf('    Ws=%.0e  Ra=%.3f\n', rowsW');
wcsv(fullfile(outDir,'aggregation_vs_Ws.csv'),'Ws,Ra_after_0p3s',rowsW);
% viscosity law ordering
Sv = S0; Sv.visc = 'Mooney'; Sv.phi = 0.05; r135 = relVisc(setfield(Sv,'psiM',1.35),0.05); r191 = relVisc(setfield(Sv,'psiM',1.91),0.05); %#ok<SFLD>
TEST(end+1,:) = {'Mooney: ratio increases with psi_M', r191 > r135};

%% ---- PART 4: SUPERPOSED vs SOLVED TOGETHER (ladder) ----------------------
fprintf('\n=== PART 4: EACH EFFECT ALONE vs ALL TOGETHER ===\n');
fprintf('  (ell = 30 d_agg ~ 5 um is the size that reproduces the v12 range k = 0.1-0.25 at h* = 0.4)\n');
cases = struct('name',{'A sourced: h*=0.4, U=1e-4','B faster: h*=0.4, U=1e-3','C thin film: h*=0.3, U=2e-4'}, ...
               'U',{1e-4,1e-3,2e-4},'hs',{0.4,0.4,0.3});
rowsL = {}; supErr = zeros(numel(cases),2);
for ic = 1:numel(cases)
    Sc = S0; Sc.law = 'bounded'; Sc.U = cases(ic).U; h = cases(ic).hs*Sc.h0;
    base = Sc; base.visc = 'none'; base.ellSource = 'fixed'; base.ellFixed = 0; base.thermal = false; base.migration = false;
    W0 = solveCoupled(base,h).W;
    eff = {'viscosity (mod. KD, Ra=3.34)','couple stress ell=d_agg','couple stress ell=30 d_agg','slip 100 nm','roughness 0.5 um', ...
           'roughness 3 um','porosity (kappa=6e-15, h_p=2 mm)','thermal (on viscosity-only state)','migration (on viscosity-only state)', ...
           'kinetic aggregation Ws=1e3 (replaces Ra=3.34)'};
    Wi = zeros(1,numel(eff)); Wvis = 0;
    for ie = 1:numel(eff)
        Sb = base;
        switch ie
            case 1, Sb.visc = 'KDagg';
            case 2, Sb.ellSource = 'particle'; Sb.ellScale = 1;
            case 3, Sb.ellSource = 'particle'; Sb.ellScale = 30;
            case 4, Sb.slip = 100e-9;
            case 5, Sb.sigma = 0.5e-6;
            case 6, Sb.sigma = 3e-6;
            case 7, Sb.porous = true;
            case 8, Sb.visc = 'KDagg'; Sb.thermal = true;
            case 9, Sb.visc = 'KDagg'; Sb.migration = true; Sb.migAmp = 1;
            case 10, Sb.visc = 'KDagg'; Sb.aggMode = 'kinetic'; Sb.Ws = 1e3;
        end
        Wi(ie) = solveCoupled(Sb,h).W; if ie == 1, Wvis = Wi(1); end
    end
    rat = Wi/W0; rat(8:10) = Wi(8:10)/Wvis;
    Sall = Sc; Sall.visc = 'KDagg'; Sall.ellSource = 'particle'; Sall.ellScale = 30; Sall.slip = 100e-9; Sall.sigma = 0.5e-6;
    Sall.porous = true; Sall.thermal = true; Sall.migration = true;
    Wall = solveCoupled(Sall,h).W; used = [1 3 4 5 7 8 9];
    prodInd = prod(rat(used)); addInd = 1 + sum(rat(used)-1);
    fprintf('\n  %s   (isoviscous Newtonian W0 = %.4g N)\n', cases(ic).name, W0);
    for ie = 1:numel(eff), fprintf('    %-46s W ratio = %.5f\n', eff{ie}, rat(ie)); end
    fprintf('    ALL TOGETHER (1,3,4,5,7,8,9) = %.5f | product of singles = %.5f | additive singles = %.5f\n', Wall/W0, prodInd, addInd);
    supErr(ic,:) = [100*(prodInd/(Wall/W0)-1) 100*(addInd/(Wall/W0)-1)];
    fprintf('    superposition error: product %+.2f %%, additive %+.2f %%\n', supErr(ic,1), supErr(ic,2));
    for ie = 1:numel(eff), rowsL(end+1,:) = {cases(ic).name, eff{ie}, rat(ie)}; end %#ok<AGROW>
    rowsL(end+1,:) = {cases(ic).name, 'ALL together', Wall/W0}; %#ok<AGROW>
    rowsL(end+1,:) = {cases(ic).name, 'product of singles', prodInd}; %#ok<AGROW>
    rowsL(end+1,:) = {cases(ic).name, 'additive singles', addInd}; %#ok<AGROW>
end
wcsvc(fullfile(outDir,'ladder_single_vs_together.csv'),'case,effect,W_ratio',rowsL);
TEST(end+1,:) = {'ladder: all-together is finite and positive in every case', all(cellfun(@(x) isfinite(x) && x>0, rowsL(:,3)))};

%% ---- PART 5: EXISTENCE LIMIT vs PRESSURE-VALIDITY ENVELOPE ----------------
fprintf('\n=== PART 5: EXISTENCE LIMIT, BOUNDED VISCOSITY AND THE 5 MPa ENVELOPE ===\n');
Se = S0; Se.visc = 'KDagg'; Se.ellSource = 'particle'; Se.ellScale = 30; Se.thermal = false; Se.migration = false;
rowsE = []; hs = [0.3 0.4 0.6 1.0];
for ih = 1:numel(hs)
    h = hs(ih)*Se.h0; Se.law = 'barus'; Uc = vcritU(Se,h);
    vcBar = Se.alpha*Se.mu0*Uc*Se.a^2/h^3;
    Uenv = envU(Se,h,Se.pLimit);
    Se.law = 'bounded'; Se.U = 1.5*Uc; ob = solveCoupled(Se,h);
    Se.law = 'barus'; Se.U = 0.99*Uc; o99 = solveCoupled(Se,h); Se.U = S0.U;
    rowsE(end+1,:) = [hs(ih) Uc vcBar Uenv Uenv/Uc o99.P0 ob.exist ob.P0]; %#ok<AGROW>
    fprintf('  h*=%.1f: U_c = %.3e m/s (v_c = %.4f) | U at 5 MPa = %.3e m/s = %.3f U_c | Barus P0(0.99 U_c) = %.2e Pa | bounded P0(1.5 U_c) = %.2e Pa\n', ...
        hs(ih), Uc, vcBar, Uenv, Uenv/Uc, o99.P0, ob.P0);
end
wcsv(fullfile(outDir,'existence_limit_vs_envelope.csv'),'h_star,Uc,vc,U_5MPa,U5MPa_over_Uc,P0_barus_at_0p99Uc,bounded_exists_1p5Uc,P0_bounded_1p5Uc',rowsE);
TEST(end+1,:) = {'bounded law exists at 1.5 U_c in all cases', all(rowsE(:,7)==1)};
TEST(end+1,:) = {'5 MPa envelope is reached before the existence limit (U_5MPa < U_c)', all(rowsE(:,5) < 1)};
fprintf('  READING: the 5 MPa envelope binds at a squeeze speed well below U_c, so the existence limit is mathematical,\n');
fprintf('  not practical; the bounded law only removes the divergence. Its pressures near the 1e10 Pa s cap are not predictions.\n\n');

%% ---- PART 6: 2-D UNIFIED SOLVER ------------------------------------------
fprintf('\n=== PART 6: 2-D RECTANGULAR PLATES, UNIFIED ===\n'); flushOut();
N2 = 16; if quick, N2 = 12; end
Lams = [1 2 3]; if quick, Lams = [1 2]; end
spL = setup2D(N2,1); P2 = pars2D(0.10,0,'barus'); P2.H = 0.4*ones(spL.n,1); P2.fR = ones(spL.n,1);
[p2,ok2] = solve2Du(spL,P2,[]);
u = poissonU1(spL.X,spL.Y,1); tau0 = flowGap(0.4,0.10,0,0)/1; pex = 12*u/tau0;
TEST(end+1,:) = {sprintf('2-D linear solve = 12 u/Phi (Poisson series), rel %.1e < 2e-5 (spectral error of the N=%d grid)', max(abs(p2-pex))/max(pex), N2), ok2 && max(abs(p2-pex))/max(pex) < 2e-5};
% convergence of v_c with N (the 0.2 % bias question): Lambda = 1, k = 0.1, H = 0.4, base fluid
Ns = [12 16 24]; if quick, Ns = [10 12 16]; end
vcN = zeros(size(Ns));
for iN = 1:numel(Ns)
    sp = setup2D(Ns(iN),1); Pc = pars2D(0.10,0,'barus'); Pc.H = 0.4*ones(sp.n,1); Pc.fR = ones(sp.n,1);
    vcN(iN) = vcrit2D(sp,Pc,0.03);
end
fprintf('  v_c(Lambda=1, k=0.1, H=0.4, base fluid) for N = '); fprintf('%d ',Ns); fprintf(': '); fprintf('%.5f ',vcN); fprintf('\n');
d1 = vcN(2)-vcN(1); d2 = vcN(3)-vcN(2);
if d1*d2 > 0 && abs(d2) < abs(d1), vInf = vcN(3) - d2^2/(d2-d1); else, vInf = NaN; end
if isnan(vInf), fprintf('  spread between last two grids: %.3f %%; Aitken extrapolation not applicable (differences not geometrically decreasing)\n', 100*abs(d2)/vcN(3));
else, fprintf('  spread between last two grids: %.3f %%; Aitken-extrapolated limit %.5f (%.3f %% below the last grid)\n', 100*abs(d2)/vcN(3), vInf, 100*(vcN(3)-vInf)/vcN(3)); end
TEST(end+1,:) = {'2-D v_c converges with N (last two grids agree to 1 %)', abs(d2)/vcN(end) < 0.01};
% exact nonlinear check: for ell = 0 and Barus law, p = -ln(1 - v p_lin)/v with p_lin the linear (v = 0) solution
Pe0 = pars2D(0,0.005,'barus'); Pe0.H = 0.4*ones(spL.n,1); Pe0.fR = ones(spL.n,1); Pl0 = Pe0; Pl0.v = 0;
[pl0,okl0] = solve2Du(spL,Pl0,[]); [pn0,okn0] = contin2D(spL,Pe0); pex0 = -log(1-0.005*pl0)/0.005;
TEST(end+1,:) = {sprintf('2-D nonlinear solve = exact -ln(1-v p_lin)/v for ell=0 (max rel %.1e < tol %.0e)',max(abs(pn0-pex0))/max(pex0),tern(quick,1e-5,1e-6)), okl0 && okn0 && max(abs(pn0-pex0))/max(pex0) < tern(quick,1e-5,1e-6)};
Pe1 = pars2D(0.10,0.005,'barus'); Pe1.H = 0.4*ones(spL.n,1); Pe1.fR = ones(spL.n,1); [pn1,okn1] = contin2D(spL,Pe1); Pm = reshape(pn1,spL.N+1,spL.N+1);
TEST(end+1,:) = {'2-D nonlinear solution has the square-plate symmetries p(x,y) = p(-x,y) = p(y,x)', okn1 && max(max(abs(Pm-Pm(end:-1:1,:)))) < 1e-8*max(Pm(:)) && max(max(abs(Pm-Pm'))) < 1e-8*max(Pm(:))};
FIG.vcN = [Ns(:) vcN(:)];
% sweep: v_c(Lambda, fluid) and status of every case (none is dropped).
% Same conventions as v12 and the published tables: phi_m = 0.63, Mooney l = 1.1, Nielsen with (1-phi), and the couple-stress
% length independent of the nanofluid ratio r ('pressure' coupling). With the self-consistent 'eta' coupling (ell^2 = eta/mu)
% the limits would be somewhat different.
vSw = [0.002 0.005 0.01]; rows2 = {}; nBeyond = 0; nTot = 0; nInc = 0; vcTab = [];
names7 = {'Original','Einstein','Brinkman','Batchelor','Krieger','Mooney','Nielsen'}; vis7s = {'none','Einstein','Brinkman','Batchelor','Krieger','Mooney','NielsenPaper'};
for iL = 1:numel(Lams)
    sp = setup2D(N2,Lams(iL));
    for im = 1:7
        Sv = S0; Sv.phi = 0.05; Sv.phim = 0.63; Sv.psiM = 1.1; Sv.visc = vis7s{im}; rr = relVisc(Sv,Sv.phi);
        Pb = pars2D(0.10,0.03,'barus'); Pb.coupling = 'pressure'; Pb.M = rr; Pb.H = 0.4*ones(sp.n,1); Pb.fR = ones(sp.n,1);
        vc = vcrit2D(sp,Pb,0.03); vcTab(end+1,:) = [Lams(iL) im rr vc]; %#ok<AGROW>
        for iv = 1:numel(vSw)
            Pb.v = vSw(iv); inside = vSw(iv) < vc*(1-1e-3); Wb = NaN; okb = false;
            if inside, [pb,okb] = contin2D(sp,Pb); if okb, Wb = Lams(iL)*(sp.w*reshape(pb,sp.N+1,sp.N+1)*sp.w'); end, end
            nTot = nTot+1; nBeyond = nBeyond + ~inside; nInc = nInc + (inside && ~okb);
            rows2(end+1,:) = {Lams(iL), names7{im}, rr, vc, vSw(iv), tern(inside,'exists','beyond v_c'), Wb}; %#ok<AGROW>
        end
        flushOut();
    end
end
wcsvc(fullfile(outDir,'sweep_2D_unified.csv'),'Lambda,model,r,v_c,v,status,W',rows2);
wcsv(fullfile(outDir,'vc_2D_by_Lambda_and_fluid.csv'),'Lambda,model_index,r,v_c',vcTab);
fprintf('  %d cases: %d lie beyond the Barus existence limit v_c(Lambda, fluid) and are REPORTED with their v_c, not dropped; %d predicted-inside cases failed to converge\n', nTot, nBeyond, nInc);
for iL = 1:numel(Lams), mk = vcTab(:,1)==Lams(iL); fprintf('  Lambda = %d: v_c from %.5f (r=%.3f) to %.5f (base fluid)\n', Lams(iL), min(vcTab(mk,4)), vcTab(find(mk & vcTab(:,4)==min(vcTab(mk,4)),1),3), max(vcTab(mk,4))); end
TEST(end+1,:) = {'2-D sweep: every case is classified and every predicted-inside case converges', nInc == 0};
fprintf('  NOTE: the bounded-viscosity law could NOT be continued in 2-D beyond v_c (Newton fails as pressures jump toward the 1e10 Pa s cap); those cases are physically outside the 5 MPa envelope anyway.\n');

%% ---- PART 7: SOBOL SENSITIVITY: DIRECT, AND THROUGH A SURROGATE ---------------
fprintf('\n=== PART 7: SOBOL (coupled model, all closures uncertain) ===\n'); flushOut();
pn = {'phi','phi_m','Ra','D','d_p','ell_scale','slip','sigma','log10 kappa','h_p','saturation','log10 U'};
lo = [0.01 0.55 1.0 1.6 10e-9  1  0    0.1e-6 -15.2 1e-3 0.5 -5];
hi = [0.08 0.65 3.34 2.2 100e-9 30 100e-9 3e-6 -13.2 5e-3 1.0 -3];
Nsob = 64; if quick, Nsob = 12; end
d = numel(pn); A = rand(Nsob,d); B = rand(Nsob,d);
fy = @(Z) sobolEval(S0, lo + Z.*(hi-lo));
nFail = 0; YA = fy(A); YB = fy(B); YAB = zeros(Nsob,d); nFail = nFail + nnz(isnan(YA)) + nnz(isnan(YB)); YA = capNaN(YA); YB = capNaN(YB);
for i = 1:d, ABi = A; ABi(:,i) = B(:,i); Yi = fy(ABi); nFail = nFail + nnz(isnan(Yi)); YAB(:,i) = capNaN(Yi); flushOut(); end
V = var([YA;YB]); S1 = zeros(1,d); ST = zeros(1,d);
for i = 1:d
    S1(i) = mean(YB.*(YAB(:,i)-YA))/V; ST(i) = 0.5*mean((YA-YAB(:,i)).^2)/V;
end
fprintf('  DIRECT: %d coupled solves (N=%d, d=%d), output = ln(W/W_isoviscous) at h*=0.4, bounded law\n', Nsob*(d+2), Nsob, d);
fprintf('  %-14s  S1      ST\n', 'parameter'); for i = 1:d, fprintf('  %-14s %6.3f  %6.3f\n', pn{i}, S1(i), ST(i)); end
fprintf('  sum S1 = %.3f (1 = purely additive; <1 = interactions); statistical error of each index is large at this N\n', sum(S1));
fprintf('  failed solves among the direct samples (capped at the largest finite value): %d of %d\n', nFail, Nsob*(d+2));
wcsvc(fullfile(outDir,'sobol_coupled.csv'),'parameter,S1,ST',[pn(:) num2cell(S1(:)) num2cell(ST(:))]);
TEST(end+1,:) = {sprintf('Sobol: at most 2 %% of the direct solves failed (%d of %d)',nFail,Nsob*(d+2)), nFail <= 0.02*Nsob*(d+2)};
% surrogate: quadratic polynomial trained on random solves, then Saltelli with 20000 samples
nTr = 300; nTe = 100; Nsur = 20000; if quick, nTr = 60; nTe = 20; Nsur = 4000; end
Xtr = rand(nTr,d); Xte = rand(nTe,d); Ytr = fy(Xtr); Yte = fy(Xte); nFail2 = nnz(isnan(Ytr)) + nnz(isnan(Yte)); Ytr = capNaN(Ytr); Yte = capNaN(Yte);
fprintf('  failed solves among the surrogate training and test samples: %d of %d\n', nFail2, nTr+nTe);
Cq = fitQuad(2*Xtr-1,Ytr); Yp = evalQuad(2*Xte-1,Cq,d);
R2 = 1 - sum((Yte-Yp).^2)/sum((Yte-mean(Yte)).^2);
fprintf('  SURROGATE: quadratic polynomial, %d training solves, held-out R^2 = %.4f (%d test solves), rms error %.4f in ln W\n', nTr, R2, nTe, sqrt(mean((Yte-Yp).^2)));
As = rand(Nsur,d); Bs = rand(Nsur,d); fs = @(Z) evalQuad(2*Z-1,Cq,d); YAs = fs(As); YBs = fs(Bs); Vs = var([YAs;YBs]); S1s = zeros(1,d); STs = zeros(1,d);
for i = 1:d, ABi = As; ABi(:,i) = Bs(:,i); Yi = fs(ABi); S1s(i) = mean(YBs.*(Yi-YAs))/Vs; STs(i) = 0.5*mean((YAs-Yi).^2)/Vs; end
fprintf('  %-14s  S1(sur) ST(sur)   S1(direct) ST(direct)\n','parameter'); for i = 1:d, fprintf('  %-14s %6.3f  %6.3f     %6.3f  %6.3f\n', pn{i}, S1s(i), STs(i), S1(i), ST(i)); end
fprintf('  surrogate sum S1 = %.3f; cost: %d solves for the surrogate versus %d for the direct study at N=%d\n', sum(S1s), nTr+nTe, Nsob*(d+2), Nsob);
wcsvc(fullfile(outDir,'sobol_surrogate.csv'),'parameter,S1_surrogate,ST_surrogate,S1_direct,ST_direct',[pn(:) num2cell(S1s(:)) num2cell(STs(:)) num2cell(S1(:)) num2cell(ST(:))]);
TEST(end+1,:) = {sprintf('Sobol surrogate: held-out R^2 = %.3f (> %.1f; %d training solves%s)',R2,tern(quick,0.5,0.9),nTr,tern(quick,', quick mode: threshold relaxed','')), R2 > tern(quick,0.5,0.9)};
FIG.sobol = [S1(:) ST(:) S1s(:) STs(:)]; FIG.pn = pn;

%% ---- PART 8: COMPARISON WITH LIN ET AL. (2013) AND SHIRALASHETTI & VATSALA (2026) ----
fprintf('\n=== PART 8: COMPARISON WITH PUBLISHED RESULTS ===\n');
% --- Lin et al. 2013: Table I (V = 0, h* = 1) and the Barus check W*(h=1, V=0.02, K=0) = 4.807
Wpub = [4.712 4.843 5.213]; Kp = [0 0.05 0.10]; Wl = zeros(1,3);
for i = 1:3, Su = unitParams(); Su.ellFixed = Kp(i); Su.alpha = 0; Wl(i) = solveCoupled(Su,1).W; end
fprintf('  Lin Table I (V=0, h*=1): K = 0, 0.05, 0.10\n    published '); fprintf('%.3f ',Wpub); fprintf('\n    present   '); fprintf('%.4f ',Wl); fprintf('\n');
TEST(end+1,:) = {'Lin et al. Table I (V=0, h*=1) reproduced, rel < 2e-4', all(abs(Wl-Wpub)./Wpub < 2e-4)};
Su = unitParams(); Su.ellFixed = 0; Su.alpha = 0.02; Su.N = 400; Wb = solveCoupled(Su,1).W;
fprintf('  Lin Barus check W*(h=1, V=0.02, K=0): published 4.807, present %.4f (%.3f %% apart)\n', Wb, 100*(Wb/4.807-1));
TEST(end+1,:) = {'Lin Barus check W*(h=1,V=0.02,K=0)=4.807 within 0.1 %', abs(Wb/4.807-1) < 1e-3};
% --- Shiralashetti & Vatsala 2026: Tables 2-4, six viscosity models, phi=0.05, phi_m=0.63, l=1.1, H*=0.4
Sx = unitParams(); Sx.phim = 0.63; Sx.psiM = 1.1; Sx.N = 100;
vis7 = {'none','Einstein','Brinkman','Batchelor','Krieger','Mooney','NielsenPaper'};
names7 = {'Original','Einstein','Brinkman','Batchelor','Krieger','Mooney','Nielsen'};
r7 = zeros(1,7); for i = 1:7, Sx.visc = vis7{i}; r7(i) = relVisc(Sx,0.05); end
ks = [0.10 0.15 0.20 0.25];
pubP = [76.65 86.24 87.14 87.48 87.32 87.50 86.86; 113.31 127.48 128.82 129.32 129.08 129.34 128.40; ...
        164.60 185.17 187.12 187.85 187.50 187.88 186.51; 230.52 259.34 262.06 263.09 262.59 263.13 261.21];
pubW = [120.41 135.46 136.88 137.42 137.16 137.44 136.43; 177.99 200.24 202.35 203.14 202.75 203.17 201.68; ...
        258.55 290.87 293.93 295.07 294.52 295.12 292.96; 362.11 407.37 411.65 413.26 412.48 413.32 410.30];
pubT = [16.966 19.087 19.288 19.363 19.326 19.366 19.224; 22.588 25.411 25.678 25.778 25.730 25.782 25.594; ...
        30.442 34.248 34.608 34.742 34.677 34.748 34.494; 40.536 45.602 46.082 46.261 46.174 46.268 45.931];
[xg,wg] = gaussleg(12); hT = 0.4; hq = 0.5*(1-hT)*xg + 0.5*(1+hT);
calP = zeros(4,7); calW = calP; calT = calP;
for ik = 1:4
    for im = 1:7
        S = Sx; S.visc = 'fixed'; S.rVal = r7(im); S.alpha = 0; S.ellFixed = ks(ik);
        o = solveCoupled(S,0.4); calP(ik,im) = o.P0; calW(ik,im) = o.W;
        Wh = zeros(size(hq)); for j = 1:numel(hq), Wh(j) = solveCoupled(S,hq(j)).W; end
        calT(ik,im) = 0.5*(1-hT)*sum(wg(:).*Wh(:));
    end
end
eP = max(abs(calP(:)./pubP(:)-1)); eW = max(abs(calW(:)./pubW(:)-1)); eT = max(abs(calT(:)./pubT(:)-1));
fprintf('  v = 0 entries of Tables 2, 3, 4 (28 each), max relative deviation: P* %.1e, W* %.1e, T* %.1e\n', eP, eW, eT);
fprintf('  viscosity ratios used: '); fprintf('%.4f ', r7); fprintf('(paper: 1, 1.1250, 1.1368, 1.1413, 1.1391, 1.1416, 1.1330)\n');
TEST(end+1,:) = {sprintf('Shiralashetti Tables 2-4 at v=0 reproduced (P* %.1e, W* %.1e, T* %.1e < 2e-3)',eP,eW,eT), eP < 2e-3 && eW < 2e-3 && eT < 2e-3};
wcsv(fullfile(outDir,'compare_shiralashetti_v0.csv'),'k,model_index,P_pub,P_present,W_pub,W_present,T_pub,T_present', ...
     [repelem(ks(:),7) repmat((1:7)',4,1) reshape(pubP',[],1) reshape(calP',[],1) reshape(pubW',[],1) reshape(calW',[],1) reshape(pubT',[],1) reshape(calT',[],1)]);
% --- v = 0.04 and 0.08 entries: are they inside the exact existence region?
rat = []; 
for v = [0.04 0.08], for ik = 1:4, for im = 1:7
    vc = numClosed(0.4,ks(ik))/(3*r7(im)); rat(end+1) = v/vc; %#ok<AGROW>
end, end, end
fprintf('  v = 0.04 and 0.08 entries (56 per table): v/v_c ranges %.2f to %.2f; entries with v < v_c: %d of %d\n', min(rat), max(rat), nnz(rat<1), numel(rat));
fprintf('  (v_c(H*=0.4, base fluid) = %.4f, %.4f, %.4f, %.4f for k = 0.10, 0.15, 0.20, 0.25)\n', numClosed(0.4,0.10)/3, numClosed(0.4,0.15)/3, numClosed(0.4,0.20)/3, numClosed(0.4,0.25)/3);
TEST(end+1,:) = {'published v>0 entries (Tables 2-4) lie beyond the exact Barus existence limit (all 56 per table)', all(rat > 1)};
fprintf('  => the published v>0 numbers cannot be solutions of the stated model; v=0 entries agree with the present solver.\n');
% --- present methodology below the limit: nanofluid enhancement at the same v (k = 0.10, H* = 0.4)
fprintf('  Present results below the limit, k = 0.10, H* = 0.4, same v for all fluids (P* ratio to base fluid):\n');
fprintf('    v/v_c(base):  '); fprintf('%-10s', names7{2:end}); fprintf('\n');
rowsC = [];
for fr = [0 0.3 0.6 0.9]
    vv = fr*numClosed(0.4,0.1)/3; Pb = 0; row = zeros(1,6);
    for im = 1:7
        S = Sx; S.visc = 'fixed'; S.rVal = r7(im); S.alpha = vv; S.ellFixed = 0.1; o = solveCoupled(S,0.4);
        if im == 1, Pb = o.P0; else, row(im-1) = o.P0/Pb; end
    end
    fprintf('    %.1f           ', fr); fprintf('%-10.4f', row); fprintf('\n'); rowsC(end+1,:) = [fr row]; %#ok<AGROW>
end
fprintf('  (Inf = beyond the existence limit of that nanofluid, which is the base-fluid limit divided by r)\n');
fprintf('  (paper reports 12.5 %% to 24.7 %% enhancement; here it equals the viscosity ratio at v=0 and grows toward the limit)\n');
wcsv(fullfile(outDir,'nanofluid_enhancement_below_limit.csv'),'v_over_vc,Einstein,Brinkman,Batchelor,Krieger,Mooney,Nielsen',rowsC);

%% ---- PART 9: NANOPARTICLE TYPES (properties from Abass et al. 2023) --------
fprintf('\n=== PART 9: TiO2, Al2O3, CuO (Abass et al. 2023 Table 2 and aggregate ratios) ===\n');
DB = particleDB();
S9 = S0; S9.mu0 = 0.0788; S9.beta = 0.032; S9.rho = 875; S9.cp = 2000; S9.kth = 0.13;   % Veedol ISO VG 46, Abass Table 1
S9.ellSource = 'fixed'; S9.ellFixed = 0; S9.law = 'bounded'; S9.thermal = true; S9.migration = false; S9.N = 100;
h9 = 0.4*S9.h0; phis = [0.005 0.01 0.02 0.03];
fprintf('  oil: Veedol ISO VG46 (mu0 = %.4f Pa s, beta = %.3f 1/K); alpha not given there, kept at %.1e 1/Pa (ASSUMED)\n', S9.mu0, S9.beta, S9.alpha);
fprintf('  particle  Ra    phi    mu_r    k_nf/k_f  (rho cp)_nf/(rho cp)_f   W/W0 (U=1e-4)  W/W0 (U=1e-3)  dT_max K (U=1e-3)\n');
rows9 = {}; ratio9 = zeros(3,4);
base1 = S9; base1.visc = 'none'; base1.phi = 0; base1.thermal = true; W0a = solveCoupled(base1,h9).W; base1.U = 1e-3; W0b = solveCoupled(base1,h9).W;
for ip = 1:3
    for iphi = 1:numel(phis)
        S = S9; S.visc = 'KDagg'; S.phi = phis(iphi); S.Ra = DB.Ra(ip); S.RaEff = DB.Ra(ip); S.phim = 0.605; S.D = 1.8;
        S.rhoP = DB.rho(ip); S.cpP = DB.cp(ip); S.kP = DB.k(ip);
        oa = solveCoupled(S,h9); S.U = 1e-3; ob = solveCoupled(S,h9);
        kr = ((S.kP+2*S.kth) - 2*S.phi*(S.kth-S.kP))/((S.kP+2*S.kth) + S.phi*(S.kth-S.kP));
        hc = ((1-S.phi)*S.rho*S.cp + S.phi*S.rhoP*S.cpP)/(S.rho*S.cp);
        ratio9(ip,iphi) = oa.W/W0a;
        fprintf('  %-8s %.2f  %.3f  %.4f  %.4f    %.4f                  %.4f         %.4f         %.3f\n', DB.name{ip}, DB.Ra(ip), S.phi, oa.mu_r, kr, hc, oa.W/W0a, ob.W/W0b, max(ob.th));
        rows9(end+1,:) = {DB.name{ip}, S.phi, oa.mu_r, kr, hc, oa.W/W0a, ob.W/W0b}; %#ok<AGROW>
    end
end
wcsvc(fullfile(outDir,'nanoparticle_types_abass.csv'),'particle,phi,mu_r,k_ratio,rhocp_ratio,W_ratio_U1e-4,W_ratio_U1e-3',rows9);
fprintf('  Abass et al. (journal bearing, CFD) report for phi = 0.01 / 0.03: max-pressure rise  Al2O3 6-17.5 %%, CuO 6-18.66 %%, TiO2 8.126-19.85 %%;\n');
fprintf('  load rise Al2O3 5-13.8 %%, CuO 4.8-14 %%, TiO2 5.2-14.3 %%. Different bearing, thermal and cavitation: only the ORDER is comparable.\n');
TEST(end+1,:) = {'particle ordering at phi=0.03 (TiO2 > CuO > Al2O3) matches Abass et al.', ratio9(1,4) > ratio9(3,4) && ratio9(3,4) > ratio9(2,4)};
TEST(end+1,:) = {'load ratio increases with phi for every particle type', all(all(diff(ratio9,1,2) > 0))};

%% ---- PART 10: CONSTANT-LOAD SQUEEZE, v NOT FROZEN ---------------------------
fprintf('\n=== PART 10: CONSTANT-LOAD APPROACH TIME, TRUE TRANSIENT vs FROZEN v ===\n'); flushOut();
Sc = S0; Sc.visc = 'KDagg'; Sc.ellSource = 'fixed'; Sc.ellFixed = 0; Sc.thermal = false; Sc.migration = false; Sc.N = 100;
ha = Sc.h0; hb = 0.4*Sc.h0; nq = 12; if quick, nq = 6; end
St = Sc; St.visc = 'none'; St.alpha = 0; Wl = 2000;
tt = approachTime(St,Wl,ha,hb,nq); tex = 1.5*pi*St.mu0*St.a^4/(2*Wl)*(1/hb^2-1/ha^2);
TEST(end+1,:) = {sprintf('constant-load approach time = analytic isoviscous result (rel %.1e < 1e-3)',abs(tt/tex-1)), abs(tt/tex-1) < 1e-3};
rowsT = []; FIG.sq = [];
for Ua = [1e-4 1e-3]
    Sb = Sc; Sb.law = 'barus'; Sb.U = Ua; Wl = solveCoupled(Sb,ha).W;
    [tT,tF,hq,Uq] = approachTime(Sb,Wl,ha,hb,nq);
    tN = 1.5*pi*Sb.mu0*relVisc(Sb,Sb.phi)*Sb.a^4/(2*Wl)*(1/hb^2-1/ha^2);
    fprintf('  initial U = %.0e m/s, load %.1f N: t(true) = %.4f s | t(frozen v) = %.4f s | t(no piezo-viscosity) = %.4f s | frozen/true = %.3f, Newtonian/true = %.3f\n', Ua, Wl, tT, tF, tN, tF/tT, tN/tT);
    rowsT(end+1,:) = [Ua Wl tT tF tN]; %#ok<AGROW>
    if Ua == 1e-3, FIG.sq = [hq(:)/Sc.h0 Uq(:)]; end
end
wcsv(fullfile(outDir,'constant_load_approach.csv'),'U_initial,load_N,t_true_s,t_frozen_s,t_newtonian_s',rowsT);
TEST(end+1,:) = {'frozen-v approach time differs from the true transient (piezoviscous case)', abs(rowsT(2,4)/rowsT(2,3)-1) > 1e-4};

%% ---- PART 11: ELASTIC LINER (Winkler foundation) ----------------------------
fprintf('\n=== PART 11: ELASTIC LINER (local Winkler model, quasi-static) ===\n'); flushOut();
lin = {'rigid',Inf,0.33; 'bronze liner E=100 GPa, 2 mm',100e9,0.33; 'polymer liner E=1 GPa, 2 mm',1e9,0.40}; tl = 2e-3;
Se2 = S0; Se2.visc = 'KDagg'; Se2.ellSource = 'fixed'; Se2.ellFixed = 0; Se2.law = 'bounded'; Se2.N = 100; h = 0.4*Se2.h0; Pq = zeros(3,2); iu = 0;
for Ua = [1e-4 1e-3]
    iu = iu+1; Se2.U = Ua;
    for il = 1:3
        S = Se2; S.compl = tl*(1+lin{il,3})*(1-2*lin{il,3})/((1-lin{il,3})*lin{il,2}); o = solveCoupled(S,h); Pq(il,iu) = o.P0;
        if il == 1, Wr = o.W; end
        if S.compl*o.P0/h > 0.10, warn = ' <-- deflection > 10 % of h: Winkler/quasi-static model no longer valid'; else, warn = ''; end
        fprintf('  U = %.0e: %-30s P0 = %.4e Pa | W/W_rigid = %.4f | max deflection = %.3f um (%.2f %% of h)%s\n', Ua, lin{il,1}, o.P0, o.W/Wr, 1e6*S.compl*o.P0, 100*S.compl*o.P0/h, warn);
    end
end
TEST(end+1,:) = {'elastic liner lowers the pressure: rigid >= bronze >= polymer', all(Pq(1,:) >= Pq(2,:)) && all(Pq(2,:) > Pq(3,:))};
fprintf('  Liner modulus, thickness and Poisson ratio are typical values, NOT sourced here (ASSUMED). Abass et al. also report lower pressure with deformation.\n');

%% ---- PART 12: FILM SEPARATION WITH CAVITATION ---------------------------------
fprintf('\n=== PART 12: SEPARATION (plates pulled apart) WITH CAVITATION ===\n'); flushOut();
% analytic isoviscous check: p(r) = -(6 mu Us/h^3)[(a^2-r^2)/2 - rc^2 ln(a/r)], p(rc) = p_cav
Ss = S0; Ss.alpha = 0; Ss.visc = 'none'; Ss.ellSource = 'fixed'; Ss.ellFixed = 0; Ss.N = 400; hs = 20e-6; Us = 1e-3; pcav = -(101325-29185);
os = sepSolve(Ss,hs,Us);
fr = @(rc) -(6*Ss.mu0*Us/hs^3)*((Ss.a^2-rc^2)/2 - rc^2*log(Ss.a/rc)) - pcav; rcA = fzero(fr,[1e-6 Ss.a*(1-1e-9)]);
TEST(end+1,:) = {sprintf('separation: cavity radius vs analytic (%.4f vs %.4f of a, rel %.1e < 1e-3)',os.rc/Ss.a,rcA/Ss.a,abs(os.rc/rcA-1)), abs(os.rc/rcA-1) < 1e-3};
Sn = S0; Sn.visc = 'KDagg'; Sn.ellSource = 'fixed'; Sn.ellFixed = 0; Sn.law = 'barus'; Sn.N = 150; rowsS = []; FIG.sep = [];
fprintf('  separation speed | cavity radius/a | min p (Pa)  | W with cavitation (N) | W Gumbel full film (N) | ratio\n');
Uslist = [1e-5 1e-4 1e-3]; FIG.sep = cell(1,3);
for ik2 = 1:3
    Us2 = Uslist(ik2);
    o = sepSolve(Sn,hs,Us2);
    fprintf('  %.0e m/s        | %.4f          | %.3e | %.4e          | %.4e            | %.4f\n', Us2, o.rc/Sn.a, o.pmin, o.W, o.Wfull, o.W/o.Wfull);
    rowsS(end+1,:) = [Us2 o.rc/Sn.a o.pmin o.W o.Wfull]; %#ok<AGROW>
    FIG.sep{ik2} = [o.r(:)/Sn.a o.p(:)];
end
wcsv(fullfile(outDir,'separation_cavitation.csv'),'Us,rc_over_a,pmin_Pa,W_cavitating_N,W_gumbel_N',rowsS);
TEST(end+1,:) = {'cavitation limits the suction force (|W| below the Gumbel full-film value)', all(abs(rowsS(:,4)) < abs(rowsS(:,5)))};
fprintf('  Cavitation pressure = ambient minus vapour pressure of Abass Table 1 (29,185 Pa); gas-cavity (Swift-Stieber) condition, no oil-film reformation model.\n');

%% ---- PART 13: THERMAL: TRANSIENT ENERGY BALANCE AND RADIAL CONDUCTION -------
fprintf('\n=== PART 13: THERMAL, TRANSIENT vs QUASI-STEADY, AND RADIAL CONDUCTION ===\n'); flushOut();
Sh = S0; Sh.visc = 'KDagg'; Sh.ellSource = 'fixed'; Sh.ellFixed = 0; Sh.law = 'bounded'; Sh.N = 100; Sh.migration = false;
nT = 40; if quick, nT = 20; end; ha = Sh.h0; hb = 0.4*Sh.h0; rhoCp = (1-Sh.phi)*Sh.rho*Sh.cp + Sh.phi*Sh.rhoP*Sh.cpP;
fprintf('  film thermal time constant rho cp h/(2 h_w) = %.3f s at h = 20 um; squeeze times below\n', rhoCp*hb/(2*Sh.hw));
% quasi-steady limit test: very slow squeeze must agree with the quasi-steady balance
Sl = Sh; Sl.U = 1e-6; ot = thermalTransient(Sl,ha,hb,nT); Sq = Sl; Sq.thermal = true; oq = solveCoupled(Sq,ot.h(end));
TEST(end+1,:) = {sprintf('transient energy balance -> quasi-steady for a very slow squeeze (theta_max ratio %.3f within 10 %%)',ot.thm(end)/max(oq.th)), abs(ot.thm(end)/max(oq.th)-1) < 0.10};
rowsH = []; FIG.thU = cell(1,2); iu = 0;
for Ua = [1e-4 1e-3]
    Sx = Sh; Sx.U = Ua; ot = thermalTransient(Sx,ha,hb,nT); hl = ot.h(end); Sx.thermal = true; oq = solveCoupled(Sx,hl);
    Sy = Sx; Sy.thermal = false; o0 = solveCoupled(Sy,hl);
    fprintf('  U = %.0e m/s (squeeze time %.3f s): at h*=%.3f  theta_max transient %.4f K, quasi-steady %.4f K | W/W_isothermal: transient %.4f, quasi-steady %.4f\n', ...
        Ua,(ha-hb)/Ua,hl/Sh.h0,ot.thm(end),max(oq.th),ot.W(end)/o0.W,oq.W/o0.W);
    rowsH(end+1,:) = [Ua ot.thm(end) max(oq.th) ot.W(end)/o0.W oq.W/o0.W]; %#ok<AGROW>
    iu = iu + 1; hg = ot.h(:); Wiso_t = zeros(size(hg)); Wqs_t = zeros(size(hg)); thq_t = zeros(size(hg));
    Sq2 = Sx; Sq2.thermal = true; Sn2 = Sx; Sn2.thermal = false;   % v14: isothermal load at the SAME film thickness
    for nn = 1:numel(hg)
        Wiso_t(nn) = solveCoupled(Sn2,hg(nn)).W; oqn = solveCoupled(Sq2,hg(nn)); Wqs_t(nn) = oqn.W; thq_t(nn) = max(oqn.th);
    end
    FIG.thU{iu} = [hg/Sh.h0, ot.W(:)./Wiso_t, Wqs_t./Wiso_t, ot.thm(:), thq_t];
end
wcsv(fullfile(outDir,'thermal_transient_vs_quasisteady.csv'),'U,theta_transient,theta_quasi,Wratio_transient,Wratio_quasi',rowsH);
% radial conduction
Sx = Sh; Sx.U = 1e-3; Sx.thermal = true; oq = solveCoupled(Sx,hb);
M = relVisc(Sx,Sx.phi)*ones(Sx.N+1,1).*exp(-Sx.beta*oq.th); [~,aux] = pressureODE(Sx,hb,M,zeros(Sx.N+1,1),0,0,0);
[th0,thc] = thermalWithConduction(Sx,oq.r,hb,aux,Sx.phi,[]); [~,thk] = thermalWithConduction(Sx,oq.r,hb,aux,Sx.phi,1e-9);
TEST(end+1,:) = {sprintf('conduction solver with k -> 0 reproduces the advective solution (max diff %.1f %% < 4 %%)',100*max(abs(thk-th0))/max(th0)), max(abs(thk-th0))/max(th0) < 0.04};
Pe = rhoCp*Sx.U*Sx.a^2/(2*hb*Sx.kth);
fprintf('  radial conduction: radial Peclet number rho cp U a^2/(2 h k) = %.1e; theta_max with conduction %.4f K, without %.4f K (%.2f %% apart)\n', Pe, max(thc), max(th0), 100*(max(thc)/max(th0)-1));
fprintf('  Newtonian, depth-averaged temperature, constant wall-loss coefficient h_w: still NOT validated.\n');

%% ---- PART 14: SHEAR-THINNING (Carreau) SENSITIVITY ----------------------------
fprintf('\n=== PART 14: SHEAR-THINNING SENSITIVITY (parameters ASSUMED, no data found) ===\n'); flushOut();
Sk = S0; Sk.visc = 'KDagg'; Sk.ellSource = 'fixed'; Sk.ellFixed = 0; Sk.law = 'bounded'; Sk.N = 100; h = 0.4*Sk.h0; W0 = solveCoupled(Sk,h).W; rowsK = [];
fprintf('  max wall shear rate at the rim: %.2e 1/s\n', 3*Sk.U*Sk.a/h^2);
fprintf('  lambda (s)   n=0.95   n=0.80   n=0.60   (W / W_Newtonian)\n');
for lam = [1e-5 3e-5 1e-4 3e-4 1e-3]
    row = zeros(1,3); jj = 0;
    for nn = [0.95 0.8 0.6], jj = jj+1; S = Sk; S.lam = lam; S.nCar = nn; row(jj) = solveCoupled(S,h).W/W0; end
    fprintf('  %.0e     %.4f   %.4f   %.4f\n', lam, row); rowsK(end+1,:) = [lam row]; %#ok<AGROW>
end
S = Sk; S.lam = 1e-4; S.nCar = 1; TEST(end+1,:) = {'Carreau with n = 1 reproduces the Newtonian load', abs(solveCoupled(S,h).W/W0-1) < 1e-12};
TEST(end+1,:) = {'shear-thinning lowers the load', all(all(rowsK(:,2:4) < 1))};
wcsv(fullfile(outDir,'shear_thinning_sensitivity.csv'),'lambda_s,n0p95,n0p80,n0p60',rowsK);
fprintf('  Film-averaged wall shear rate 3 U r/h^2 used in the Carreau factor (couple-stress and thermal feedback on the rate neglected).\n');

%% ---- PART 16: ONE-AT-A-TIME PARAMETER SWEEPS (v14: data for Figs 3, 8, 9, 11) ----
fprintf('\n=== PART 16: PARAMETER SWEEPS (one parameter at a time, others at base values) ===\n'); flushOut();
Sb0 = S0; Sb0.law = 'bounded'; Sb0.visc = 'KDagg'; Sb0.ellSource = 'fixed'; Sb0.ellFixed = 0;
Sb0.thermal = false; Sb0.migration = false; Sb0.U = 1e-4;
hsw = [0.3 0.4 0.6];                        % film heights h*
hdf = 0.4*S0.h0;                            % reference film for single-height sweeps
Wref = solveCoupled(Sb0,hdf).W;
% (a) roughness, Patir-Cheng isotropic: W / W(smooth)
sig = logspace(log10(0.1e-6),log10(5e-6),9)';
FIG.rough = [sig zeros(numel(sig),numel(hsw))];
for ih = 1:numel(hsw)
    h = hsw(ih)*S0.h0; Wr0 = solveCoupled(Sb0,h).W;
    for is = 1:numel(sig), Sr = Sb0; Sr.sigma = sig(is); FIG.rough(is,ih+1) = solveCoupled(Sr,h).W/Wr0; end
end
fprintf('  roughness: W/W(smooth), sigma in um | h*=0.3 | h*=0.4 | h*=0.6\n');
fprintf('    %8.3f  %.5f  %.5f  %.5f\n', [FIG.rough(:,1)*1e6 FIG.rough(:,2:end)]');
% (b) slip length, both walls (100 nm = sourced upper value; 200-1000 nm = extension beyond the sources)
bsl = [0 10 20 50 100 200 500 1000]'*1e-9;
FIG.slip = [bsl zeros(numel(bsl),numel(hsw))];
for ih = 1:numel(hsw)
    h = hsw(ih)*S0.h0; Ws0 = solveCoupled(Sb0,h).W;
    for ib = 1:numel(bsl), Sr = Sb0; Sr.slip = bsl(ib); FIG.slip(ib,ih+1) = solveCoupled(Sr,h).W/Ws0; end
end
fprintf('  slip: W/W(no slip), b in nm | h*=0.3 | h*=0.4 | h*=0.6\n');
fprintf('    %8.0f  %.5f  %.5f  %.5f\n', [FIG.slip(:,1)*1e9 FIG.slip(:,2:end)]');
% (c) porosity: bronze-type permeability kappa and wall thickness h_p (h* = 0.4); Darcy term + BJS slip sqrt(kappa)
kap = [6e-16 2e-15 6e-15 2e-14 6e-14 2e-13]';
hpl = [1 2 5]*1e-3;
FIG.por = [kap zeros(numel(kap),numel(hpl))];
for ik3 = 1:numel(kap)
    for ip3 = 1:numel(hpl), Sr = Sb0; Sr.porous = true; Sr.kappa = kap(ik3); Sr.hp = hpl(ip3); FIG.por(ik3,ip3+1) = solveCoupled(Sr,hdf).W/Wref; end
end
fprintf('  porosity: W/W(solid), kappa (m^2) | hp=1 mm | hp=2 mm | hp=5 mm\n');
fprintf('    %9.2e  %.5f  %.5f  %.5f\n', FIG.por');
% (d) couple-stress length ell (absolute; UNSOURCED closure: d_agg = 0.17 um): W / W(ell = 0)
ellv = [0 0.1 0.2 0.5 1 2 5 10]'*1e-6;
FIG.ell = [ellv zeros(numel(ellv),numel(hsw))];
for ih = 1:numel(hsw)
    h = hsw(ih)*S0.h0; We0 = solveCoupled(Sb0,h).W;
    for ie = 1:numel(ellv), Sr = Sb0; Sr.ellSource = 'fixed'; Sr.ellFixed = ellv(ie); FIG.ell(ie,ih+1) = solveCoupled(Sr,h).W/We0; end
end
fprintf('  couple-stress length: W/W(ell=0), ell in um | h*=0.3 | h*=0.4 | h*=0.6\n');
fprintf('    %8.2f  %.5f  %.5f  %.5f\n', [FIG.ell(:,1)*1e6 FIG.ell(:,2:end)]');
% (e) viscosity closures: mu_r(phi) and load ratio W/W_Newtonian(phi) at h* = 0.4 (Mooney psi = 1.35, 1.63, 1.91)
vmod = {'Einstein',1.63; 'Brinkman',1.63; 'Batchelor',1.63; 'Krieger',1.63; 'Mooney',1.35; 'Mooney',1.63; 'Mooney',1.91; 'KDagg',1.63};
vnam = {'Einstein','Brinkman','Batchelor','Krieger','Mooney psi=1.35','Mooney psi=1.63','Mooney psi=1.91','mod. KD Ra=3.34'};
phg = linspace(0.001,0.07,35)'; phs = [0.01 0.02 0.03 0.05 0.07]';
FIG.visc_r = [phg zeros(numel(phg),size(vmod,1))];
FIG.visc_W = [phs zeros(numel(phs),size(vmod,1))];
Sn0 = Sb0; Sn0.visc = 'none'; Wnw = solveCoupled(Sn0,hdf).W;
for im = 1:size(vmod,1)
    Sv = Sb0; Sv.visc = vmod{im,1}; Sv.psiM = vmod{im,2}; Sv.phim = 0.605;
    FIG.visc_r(:,im+1) = relVisc(Sv,phg);
    for j = 1:numel(phs), Sw = Sv; Sw.phi = phs(j); FIG.visc_W(j,im+1) = solveCoupled(Sw,hdf).W/Wnw; end
end
fprintf('  load ratio W/W_Newtonian at h*=0.4, U=0.1 mm/s; rows = phi\n    phi   ');
fprintf('%16s', vnam{:}); fprintf('\n');
for j = 1:numel(phs), fprintf('    %.3f ', FIG.visc_W(j,1)); fprintf('%16.5f', FIG.visc_W(j,2:end)); fprintf('\n'); end
% (f) kinetic aggregation: Ra(t) for several stability ratios Ws, and the load ratio at the squeeze time (0.3 s)
Wsl = [1 10 100 1e3 1e4 1e6]';
tg = logspace(-4,0,60)';
FIG.Rat = [tg zeros(numel(tg),numel(Wsl))];
Sg = Sb0; Sg.aggMode = 'kinetic';
for iw = 1:numel(Wsl)
    Sg.Ws = Wsl(iw);
    for it = 1:numel(tg), FIG.Rat(it,iw+1) = aggRatio(Sg,tg(it)); end
end
FIG.agg = zeros(numel(Wsl),2);
for iw = 1:numel(Wsl), Sg.Ws = Wsl(iw); FIG.agg(iw,:) = [Wsl(iw) solveCoupled(Sg,hdf).W/Wref]; end
fprintf('  kinetic aggregation: W / W(Ra=3.34) at h*=0.4 (squeeze time 0.3 s)\n');
fprintf('    Ws=%8.0e  W ratio=%.5f\n', FIG.agg');
% (g) constant-load squeeze: U(h) for two initial speeds, aggregated nanofluid (Barus) vs Newtonian at the same load
hcl = linspace(0.4,1,13)'*S0.h0;
FIG.cl = [hcl/S0.h0 zeros(numel(hcl),4)];
Uinit = [1e-4 1e-3];
for iu = 1:2
    Ua0 = Uinit(iu);
    Sk1 = Sc; Sk1.visc = 'KDagg'; Sk1.U = Ua0; Wl0 = solveCoupled(Sk1,S0.h0).W;
    Sn1 = Sk1; Sn1.visc = 'none';
    for ih2 = 1:numel(hcl)
        FIG.cl(ih2,2*iu)   = speedForLoad(Sk1,hcl(ih2),Wl0);
        FIG.cl(ih2,2*iu+1) = speedForLoad(Sn1,hcl(ih2),Wl0);
    end
    fprintf('  constant load, U0 = %.0e m/s: U(h0) = %.4e, U(0.4 h0) KD = %.4e, Newtonian = %.4e m/s\n', ...
        Ua0, FIG.cl(end,2*iu), FIG.cl(1,2*iu), FIG.cl(1,2*iu+1));
end
flushOut();

%% ---- PART 15: FIGURES (v14) ------------------------------------------------
% Times New Roman; journal widths (single column 3.5 in, double column 7.2 in); PNG 300 dpi + vector PDF.
% Every plotted array is also written to CSV so each figure can be rebuilt without rerunning the model.
fprintf('\n=== PART 15: FIGURES ===\n'); flushOut();
try, set(0,'defaultAxesFontName','Times New Roman'); set(0,'defaultTextFontName','Times New Roman'); set(0,'defaultLegendFontName','Times New Roman'); catch, end %#ok<NOCOMMA>
nOK = 0;
% ---- Fig 1: coupled fields, three panels filling the double-column width
try
    o = FIG.oA; r = o.r/S0.a;
    wd = (1-0.085-0.02-2*0.085)/3; xs = 0.085 + (0:2)*(wd+0.085);
    f = figure('visible','off','Units','inches','Position',[1 1 7.2 2.7]);
    axes('Position',[xs(1) 0.17 wd 0.74]); plot(r,o.p/1e6,'-','LineWidth',1.2); xlabel('r/a'); ylabel('p (MPa)'); title('(a) pressure'); grid on; xlim([0 1]);
    axes('Position',[xs(2) 0.17 wd 0.74]); plot(r,o.th,'-','LineWidth',1.2); xlabel('r/a'); ylabel('\theta (K)'); title('(b) temperature rise'); grid on; xlim([0 1]);
    axes('Position',[xs(3) 0.17 wd 0.74]); plot(r,(o.phi/SA.phi-1)*1e6,'-','LineWidth',1.2); xlabel('r/a'); ylabel('\Delta\phi/\phi_0 (ppm)'); title('(c) particle fraction'); grid on; xlim([0 1]);
    wcsv(fullfile(outDir,'fig1_coupled_fields.csv'),'r_over_a,p_Pa,theta_K,dphi_over_phi0_ppm',[r o.p(:) o.th(:) (o.phi(:)/SA.phi-1)*1e6]);
    saveFig(f,fullfile(outDir,'fig1_coupled_fields')); nOK = nOK+1;
catch ME, fprintf('  Fig 1 skipped: %s\n', ME.message); end
% ---- Fig 2: Sobol indices, direct and surrogate, first-order and total-order
try
    snames = {'phi','phi_m','Ra','D','d_p','ell_s','slip','sigma','log kappa','h_p','S_sat','log U'};
    f = figure('visible','off','Units','inches','Position',[1 1 7.2 3.0]);
    subplot(1,2,1); bar(FIG.sobol(:,[1 3])); set(gca,'XTick',1:12,'XTickLabel',snames); ylabel('first-order S_1'); legend('direct','surrogate','Location','northeast'); grid on; title('(a)');
    subplot(1,2,2); bar(FIG.sobol(:,[2 4])); set(gca,'XTick',1:12,'XTickLabel',snames); ylabel('total-order S_T'); legend('direct','surrogate','Location','northeast'); grid on; title('(b)');
    try, xtickangle(subplot(1,2,1),45); xtickangle(subplot(1,2,2),45); catch, end %#ok<NOCOMMA>
    set(subplot(1,2,1),'Position',[0.07 0.30 0.41 0.64]); set(subplot(1,2,2),'Position',[0.56 0.30 0.41 0.64]);
    saveFig(f,fullfile(outDir,'fig2_sobol')); nOK = nOK+1;
catch ME, fprintf('  Fig 2 skipped: %s\n', ME.message); end
% ---- Fig 3: constant-load squeeze speed vs film thickness, two initial speeds, KD vs Newtonian
try
    wcsv(fullfile(outDir,'constant_load_speed_curves.csv'),'h_over_h0,U_KD_U0_1e-4,U_Newt_U0_1e-4,U_KD_U0_1e-3,U_Newt_U0_1e-3',FIG.cl);
    f = figure('visible','off','Units','inches','Position',[1 1 3.5 3.1]);
    hh = FIG.cl(:,1);
    plot(hh,FIG.cl(:,2)*1e3,'-o','MarkerSize',3,'LineWidth',1.1); hold on;
    plot(hh,FIG.cl(:,3)*1e3,'--','LineWidth',1.1);
    plot(hh,FIG.cl(:,4)*1e3,'-s','MarkerSize',3,'LineWidth',1.1);
    plot(hh,FIG.cl(:,5)*1e3,'--','LineWidth',1.1); hold off;
    xlabel('h/h_0'); ylabel('squeeze speed U (mm/s)'); grid on; xlim([0.4 1]);
    legend('U_0 = 0.1 mm/s, aggregated nanofluid','U_0 = 0.1 mm/s, Newtonian (same load)','U_0 = 1 mm/s, aggregated nanofluid','U_0 = 1 mm/s, Newtonian (same load)','Location','northwest','FontSize',7);
    saveFig(f,fullfile(outDir,'fig3_constant_load_speed')); nOK = nOK+1;
catch ME, fprintf('  Fig 3 skipped: %s\n', ME.message); end
% ---- Fig 4: separation pressure profiles for three separation speeds (cavitation at p = -72 kPa)
try
    M4 = []; for ik2 = 1:3, P = FIG.sep{ik2}; M4 = [M4; repmat(Uslist(ik2),size(P,1),1) P]; end
    wcsv(fullfile(outDir,'separation_profiles.csv'),'Us_m_per_s,r_over_a,p_Pa',M4);
    f = figure('visible','off','Units','inches','Position',[1 1 3.5 3.1]); hold on; hs4 = zeros(1,3); lg4 = cell(1,3);
    for ik2 = 1:3, P = FIG.sep{ik2}; hs4(ik2) = plot(P(:,1),P(:,2)/1e3,'-','LineWidth',1.1); lg4{ik2} = sprintf('U_s = %.0e m/s',Uslist(ik2)); end
    hold off; xlabel('r/a'); ylabel('p (kPa)'); legend(hs4,lg4,'Location','northwest','FontSize',8); grid on; xlim([0 1]);
    saveFig(f,fullfile(outDir,'fig4_separation')); nOK = nOK+1;
catch ME, fprintf('  Fig 4 skipped: %s\n', ME.message); end
% ---- Fig 5: transient thermal effect. v14 FIX: W is divided by the isothermal load at the SAME film thickness
%      (v13 divided every time step by the isothermal load at the final thickness, so the curve was wrong)
try
    T1 = FIG.thU{1}; T2 = FIG.thU{2};
    TEST(end+1,:) = {sprintf('Fig 5 data: transient W/W_iso = 1 at the first step (rel %.1e < 2 %%)', abs(T2(1,2)-1)), abs(T2(1,2)-1) < 0.02};
    wcsv(fullfile(outDir,'thermal_transient_curves.csv'),'h_over_h0,Wratio_tr_1e-4,Wratio_qs_1e-4,theta_tr_1e-4,theta_qs_1e-4,Wratio_tr_1e-3,Wratio_qs_1e-3,theta_tr_1e-3,theta_qs_1e-3',[T1 T2(:,2:5)]);
    f = figure('visible','off','Units','inches','Position',[1 1 7.2 2.9]);
    subplot(1,2,1); hold on;
    hA1 = plot(T1(:,1),T1(:,2),'-','LineWidth',1.2); hA2 = plot(T1(:,1),T1(:,3),'--','LineWidth',1.2);
    hB1 = plot(T2(:,1),T2(:,2),'-','LineWidth',1.2); hB2 = plot(T2(:,1),T2(:,3),'--','LineWidth',1.2); hold off;
    xlabel('h/h_0'); ylabel('W / W_{isothermal}'); grid on; xlim([0.4 1]);
    legend([hA1 hA2 hB1 hB2],{'U = 0.1 mm/s, transient','U = 0.1 mm/s, quasi-steady','U = 1 mm/s, transient','U = 1 mm/s, quasi-steady'},'Location','southeast','FontSize',7); title('(a)');
    subplot(1,2,2); hold on;
    hC1 = plot(T1(:,1),T1(:,4),'-','LineWidth',1.2); hC2 = plot(T1(:,1),T1(:,5),'--','LineWidth',1.2);
    hD1 = plot(T2(:,1),T2(:,4),'-','LineWidth',1.2); hD2 = plot(T2(:,1),T2(:,5),'--','LineWidth',1.2); hold off;
    xlabel('h/h_0'); ylabel('\theta_{max} (K)'); grid on; xlim([0.4 1]); title('(b)');
    legend([hC1 hC2 hD1 hD2],{'U = 0.1 mm/s, transient','U = 0.1 mm/s, quasi-steady','U = 1 mm/s, transient','U = 1 mm/s, quasi-steady'},'Location','northeast','FontSize',6);
    saveFig(f,fullfile(outDir,'fig5_thermal_transient')); nOK = nOK+1;
catch ME, fprintf('  Fig 5 skipped: %s\n', ME.message); end
% ---- Fig 6: particle types (Abass properties), load ratio vs volume fraction
try
    f = figure('visible','off','Units','inches','Position',[1 1 3.5 3.0]);
    plot(phis,ratio9','o-','LineWidth',1.1,'MarkerSize',4); xlabel('\phi'); ylabel('W / W_0'); legend(DB.name,'Location','northwest'); grid on;
    wcsv(fullfile(outDir,'fig6_particle_types.csv'),'phi,W_ratio_TiO2,W_ratio_Al2O3,W_ratio_CuO',[phis(:) ratio9']);
    saveFig(f,fullfile(outDir,'fig6_particle_types')); nOK = nOK+1;
catch ME, fprintf('  Fig 6 skipped: %s\n', ME.message); end
% ---- Fig 7: 2-D existence limit: grid convergence (left) and v_c by fluid model and Lambda (right)
try
    wcsv(fullfile(outDir,'vc_convergence.csv'),'N,v_c_Lambda1',FIG.vcN);
    nL = numel(Lams); vcM = zeros(7,nL);
    for iL = 1:nL, for im = 1:7, mk = vcTab(:,1)==Lams(iL) & vcTab(:,2)==im; vcM(im,iL) = vcTab(mk,4); end, end
    f = figure('visible','off','Units','inches','Position',[1 1 7.2 2.9]);
    subplot(1,2,1); plot(FIG.vcN(:,1),(FIG.vcN(:,2)/FIG.vcN(end,2)-1)*100,'o-','LineWidth',1.1); xlabel('grid N'); ylabel('v_c deviation from N = 24 (%)'); grid on; title('(a)');
    subplot(1,2,2); bar(vcM); set(gca,'XTick',1:7,'XTickLabel',names7); ylabel('v_c (Barus limit)');
    lgL = cell(1,nL); for iL = 1:nL, lgL{iL} = sprintf('\\Lambda = %d',Lams(iL)); end
    ylim([0 1.3*max(vcM(:))]); legend(lgL,'Location','north','Orientation','horizontal','FontSize',7); grid on; title('(b)');
    try, xtickangle(subplot(1,2,2),45); catch, end %#ok<NOCOMMA>
    saveFig(f,fullfile(outDir,'fig7_vc_convergence')); nOK = nOK+1;
catch ME, fprintf('  Fig 7 skipped: %s\n', ME.message); end
% ---- Fig 8: viscosity closures: mu_r(phi) (left) and load ratio vs phi (right), Mooney psi family included
try
    wcsv(fullfile(outDir,'viscosity_ratio_vs_phi.csv'),strjoin([{'phi'} vnam],','),FIG.visc_r);
    wcsv(fullfile(outDir,'viscosity_load_vs_phi.csv'),strjoin([{'phi'} vnam],','),FIG.visc_W);
    f = figure('visible','off','Units','inches','Position',[1 1 7.2 3.0]);
    subplot(1,2,1); plot(FIG.visc_r(:,1),FIG.visc_r(:,2:end),'LineWidth',1.0); xlabel('\phi'); ylabel('\mu_r'); legend(vnam,'Location','northwest','FontSize',6); grid on; title('(a)');
    subplot(1,2,2); plot(FIG.visc_W(:,1),FIG.visc_W(:,2:end),'o-','MarkerSize',3,'LineWidth',1.0); xlabel('\phi'); ylabel('W / W_{Newtonian}'); legend(vnam,'Location','northwest','FontSize',6); grid on; title('(b)');
    saveFig(f,fullfile(outDir,'fig8_viscosity_models')); nOK = nOK+1;
catch ME, fprintf('  Fig 8 skipped: %s\n', ME.message); end
% ---- Fig 9: one-at-a-time sweeps of the closures, 2 x 3 panels
try
    wcsv(fullfile(outDir,'roughness_sweep.csv'),'sigma_m,ratio_h0p3,ratio_h0p4,ratio_h0p6',FIG.rough);
    wcsv(fullfile(outDir,'slip_sweep.csv'),'slip_m,ratio_h0p3,ratio_h0p4,ratio_h0p6',FIG.slip);
    wcsv(fullfile(outDir,'porosity_sweep.csv'),'kappa_m2,ratio_hp1mm,ratio_hp2mm,ratio_hp5mm',FIG.por);
    wcsv(fullfile(outDir,'couple_stress_length_sweep.csv'),'ell_m,ratio_h0p3,ratio_h0p4,ratio_h0p6',FIG.ell);
    wcsv(fullfile(outDir,'aggregation_load_vs_Ws.csv'),'Ws,W_ratio_vs_Ra3p34',FIG.agg);
    f = figure('visible','off','Units','inches','Position',[1 1 7.2 5.8]);
    subplot(2,3,1); semilogx(rowsK(:,1),rowsK(:,2:4),'o-','MarkerSize',3); xlabel('\lambda (s)'); ylabel('W / W_{Newt}'); legend('n=0.95','n=0.80','n=0.60','Location','southwest','FontSize',7); grid on; title('(a) shear thinning');
    subplot(2,3,2); plot(FIG.rough(:,1)*1e6,FIG.rough(:,2:end),'-o','MarkerSize',3); xlabel('\sigma (\mum)'); ylabel('W / W_{smooth}'); legend('h*=0.3','h*=0.4','h*=0.6','Location','northwest','FontSize',7); grid on; title('(b) roughness');
    subplot(2,3,3); plot(FIG.slip(:,1)*1e9,FIG.slip(:,2:end),'-o','MarkerSize',3); hold on; Y9 = FIG.slip(:,2:end); yl = [min(Y9(:))*0.99 max(Y9(:))*1.01]; plot([100 100],yl,'k:'); hold off; ylim(yl); xlabel('slip length b (nm)'); ylabel('W / W_{no slip}'); grid on; title('(c) slip (dotted: 100 nm)');
    subplot(2,3,4); semilogx(FIG.por(:,1),FIG.por(:,2:end),'-o','MarkerSize',3); xlabel('\kappa (m^2)'); ylabel('W / W_{solid}'); legend('h_p=1 mm','h_p=2 mm','h_p=5 mm','Location','southwest','FontSize',7); grid on; title('(d) porosity');
    subplot(2,3,5); plot(FIG.ell(:,1)*1e6,FIG.ell(:,2:end),'-o','MarkerSize',3); xlabel('couple-stress length l (\mum)'); ylabel('W / W_{l=0}'); legend('h*=0.3','h*=0.4','h*=0.6','Location','northwest','FontSize',7); grid on; title('(e) couple stress');
    subplot(2,3,6); semilogx(FIG.agg(:,1),FIG.agg(:,2),'o-','MarkerSize',3); xlabel('stability ratio W_s'); ylabel('W / W_{Ra=3.34}'); grid on; title('(f) aggregation');
    cx = [0.08 0.41 0.74]; ry = [0.56 0.12];
    for k9 = 1:6, set(subplot(2,3,k9),'Position',[cx(mod(k9-1,3)+1) ry(ceil(k9/3)) 0.23 0.29]); end
    saveFig(f,fullfile(outDir,'fig9_parametric_sweeps')); nOK = nOK+1;
catch ME, fprintf('  Fig 9 skipped: %s\n', ME.message); end
% ---- Fig 10: single effects vs all together, from the Part 4 ladder (three cases)
try
    effs = {'viscosity (KD)','couple stress d_{agg}','couple stress 30 d_{agg}','slip 100 nm','roughness 0.5 \mum','roughness 3 \mum','porosity','thermal (qs)','migration','kinetic agg. W_s=10^3'};
    ttl = {'(a) A: h*=0.4, U=0.1 mm/s','(b) B: h*=0.4, U=1 mm/s','(c) C: h*=0.3, U=0.2 mm/s'};
    f = figure('visible','off','Units','inches','Position',[1 1 7.2 3.4]);
    px10 = [0.25 0.50 0.76]; pw10 = [0.20 0.22 0.22];
    for ic = 1:3
        b0 = (ic-1)*13;
        axes('Position',[px10(ic) 0.17 pw10(ic) 0.68]); barh(1:10,cell2mat(rowsL(b0+(1:10),3)));
        set(gca,'YDir','reverse'); xlim([0 3.2]); grid on; hold on;
        allv = rowsL{b0+11,3}; pv = rowsL{b0+12,3};
        line([allv allv],[0.5 10.5],'Color','r','LineStyle','--'); line([pv pv],[0.5 10.5],'Color','b','LineStyle',':'); hold off;
        xlabel('load ratio'); title(ttl{ic});
        if ic == 1, set(gca,'YTick',1:10,'YTickLabel',effs,'FontSize',7); else, set(gca,'YTick',1:10,'YTickLabel',repmat({' '},10,1)); end
    end
    saveFig(f,fullfile(outDir,'fig10_single_vs_together')); nOK = nOK+1;
catch ME, fprintf('  Fig 10 skipped: %s\n', ME.message); end
% ---- Fig 11: kinetic aggregation: stability-limited growth of Ra(t) (left), load ratio at 0.3 s (right)
try
    wcsv(fullfile(outDir,'aggregation_Ra_vs_time.csv'),'time_s,Ra_Ws1,Ra_Ws10,Ra_Ws100,Ra_Ws1e3,Ra_Ws1e4,Ra_Ws1e6',FIG.Rat);
    f = figure('visible','off','Units','inches','Position',[1 1 7.2 2.9]);
    subplot(1,2,1); semilogx(FIG.Rat(:,1),FIG.Rat(:,2:end),'LineWidth',1.0); xlabel('time (s)'); ylabel('R_a'); grid on; title('(a)');
    legend('W_s=1','10','10^2','10^3','10^4','10^6','Location','northwest','FontSize',7);
    subplot(1,2,2); semilogx(FIG.agg(:,1),FIG.agg(:,2),'o-','MarkerSize',4); xlabel('stability ratio W_s'); ylabel('W / W_{R_a=3.34}'); grid on; title('(b)');
    saveFig(f,fullfile(outDir,'fig11_aggregation_kinetics')); nOK = nOK+1;
catch ME, fprintf('  Fig 11 skipped: %s\n', ME.message); end
fprintf('  %d of 11 figures written (PNG 300 dpi and PDF), with their data in CSV files\n', nOK);


%% ---- SUMMARY --------------------------------------------------------------
fprintf('\n==================== SELF-TEST SUMMARY ====================\n');
for q = 1:size(TEST,1), fprintf('  [%s] %s\n', tern(TEST{q,2},'PASS','FAIL'), TEST{q,1}); end
fprintf('  %d / %d passed. Run time %.1f min.\n', nnz([TEST{:,2}]), size(TEST,1), toc(tAll)/60);
fprintf('\nREMAINING LIMITS: closures marked ASSUMED/UNKNOWN in Part 0; no experiment; thermal is quasi-steady Newtonian;\n  bounded-law pressures exceed the 5 MPa envelope where elastic deformation matters; no experimental validation.\n');
diary off
end

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================
function s = tern(c,a,b), if c, s = a; else, s = b; end, end

function S = defaultParams()
S.mu0 = 0.065; S.alpha = 1.8e-8; S.a = 20e-3; S.h0 = 50e-6; S.U = 1e-4;       % Jabbar et al. 2026 / v12
S.beta = 0.0285; S.rho = 890; S.cp = 1900; S.kth = 0.145; S.hw = 80;
S.law = 'barus'; S.muMax = 1e10; S.pLimit = 5e6;                               % Bair 2019 glass viscosity
S.visc = 'KDagg'; S.phi = 0.05; S.phim = 0.605; S.psiM = 1.63; S.D = 1.8;     % Chen et al. 2007
S.Ra = 3.34; S.RaEff = 3.34; S.RaMax = 3.34; S.aggMode = 'fixed'; S.Ws = 1e3; S.rVal = 1;   % RaMax caps kinetic growth only (literature range 3.33-4)
S.dp = 50e-9; S.ellSource = 'particle'; S.ellScale = 1; S.ellFixed = 0; S.ellCoupling = 'eta'; S.gc = 0;
S.slip = 0; S.sigma = 0; S.C_PC = 0.9; S.r_PC = 0.56;                           % Patir-Cheng
S.porous = false; S.kappa = 6e-15; S.hp = 2e-3; S.Ssat = 1; S.alphaBJS = 1;     % kappa: Watanabe range, hp ASSUMED
S.thermal = false; S.migration = false; S.migAmp = 1;
S.compl = 0; S.lam = 0; S.nCar = 1;                                          % elastic liner compliance (m/Pa); Carreau lambda (s), n
S.rhoP = 4230; S.cpP = 710; S.kP = 8.4; S.T0 = 313.15; S.Kc = 0.41; S.Kmu = 0.62;   % TiO2: Abass et al. 2023 Table 2; Kc, Kmu: Phillips et al. 1992
S.N = 200; S.relax = 0.7; S.tol = 1e-7; S.maxIt = 80; S.hStart = S.h0;
end

function S = unitParams()
S = defaultParams();
S.mu0 = 1; S.a = 1; S.h0 = 1; S.U = 1; S.alpha = 0; S.beta = 0; S.hStart = 1;
S.visc = 'none'; S.ellSource = 'fixed'; S.ellFixed = 0; S.ellCoupling = 'pressure'; S.law = 'barus'; S.muMax = 1e6;
S.N = 400; S.porous = false; S.slip = 0; S.sigma = 0; S.thermal = false; S.migration = false;
end

function wcsv(fn,hdr,M)
fid = fopen(fn,'w'); fprintf(fid,'%s\n',hdr);
for i = 1:size(M,1), fprintf(fid,[repmat('%.10g,',1,size(M,2)-1) '%.10g\n'],M(i,:)); end
fclose(fid);
end

function wcsvc(fn,hdr,C)
fid = fopen(fn,'w'); fprintf(fid,'%s\n',hdr);
for i = 1:size(C,1)
    for j = 1:size(C,2)
        x = C{i,j};
        if ischar(x), fprintf(fid,'"%s"',x); else, fprintf(fid,'%.10g',x); end
        if j < size(C,2), fprintf(fid,','); else, fprintf(fid,'\n'); end
    end
end
fclose(fid);
end

% ---------------------------------------------------------------- closures
function r = relVisc(S,phi)
switch S.visc
    case 'none',      r = ones(size(phi));
    case 'fixed',     r = S.rVal*ones(size(phi));
    case 'Einstein',  r = 1+2.5*phi;
    case 'Brinkman',  r = (1-phi).^(-2.5);
    case 'Batchelor', r = 1+2.5*phi+6.5*phi.^2;
    case 'Krieger',   r = (1-phi/S.phim).^(-2.5*S.phim);
    case 'Mooney',    r = exp(2.5*phi./(1-S.psiM*phi));
    case 'Nielsen',   r = (1+1.5*phi).*exp(phi/(1-S.phim));      % original Nielsen (Mishra et al. 2014 review)
    case 'NielsenPaper', r = (1+1.5*phi).*exp(phi./(1-phi));    % form printed in Shiralashetti & Vatsala 2026
    case 'KDagg'
        pa = phi*S.RaEff^(3-S.D);
        pa = min(pa,0.999*S.phim);
        r = (1-pa/S.phim).^(-2.5*S.phim);
    otherwise, error('unknown viscosity model');
end
end

function Ra = aggRatio(S,t)
% Smoluchowski perikinetic kinetics with stability ratio Ws, fractal clusters:
%   N_p(t) = 1 + t/(tau_c Ws),  tau_c = 3 mu0/(4 kB T n0),  Ra = N_p^(1/D)
if strcmp(S.aggMode,'fixed'), Ra = S.Ra; return; end
kB = 1.380649e-23; n0 = S.phi/(pi*S.dp^3/6); tc = 3*S.mu0/(4*kB*S.T0*n0);
Ra = min((1+t/(tc*S.Ws))^(1/S.D), S.RaMax);
end

function l0 = ellBase(S,Ra)
if strcmp(S.ellSource,'particle'), l0 = S.ellScale*S.dp*Ra; else, l0 = S.ellFixed; end
end

% ------------------------------------------------- couple stress + slip flow
function g = gfun(y)
A = 3*[2/15,-17/315,62/2835,-1382/155925,21844/6081075,-929569/638512875,6404582/10854718875,-443861162/1856156927625];
g = zeros(size(y)); sm = y < 0.3;
if any(sm(:))
    ys = y(sm); y2 = ys.^2; gs = zeros(size(ys));
    for n = 8:-1:1, gs = gs.*y2 + A(n); end
    g(sm) = gs.*y2;
end
if any(~sm(:)), yl = y(~sm); g(~sm) = 1 - 3./yl.^2 + 3*tanh(yl)./yl.^3; end
end

function w = wfun(c)
% w = 1 - tanh(c)/c, series for small c (avoids cancellation when ell >> h)
w = zeros(size(c)); sm = c < 0.3;
if any(sm(:))
    x = c(sm); x2 = x.^2;
    w(sm) = x2.*(1/3 + x2.*(-2/15 + x2.*(17/315 + x2.*(-62/2835 + x2.*1382/155925))));
end
if any(~sm(:)), xl = c(~sm); w(~sm) = 1 - tanh(xl)./xl; end
end

function Phi = flowGap(h,ell,b1,b2)
% Flux of a couple-stress fluid with Navier slip b1 (y=0), b2 (y=h) and vanishing
% couple stress at both walls:  Q = -(G/(12 mu)) Phi, derived in closed form:
%   Phi = h^3 g(c) + 6 h s (h(b1+b2) + 4 b1 b2)/(h+b1+b2),  c = h/(2 ell),
%   s = (h/2)(1 - tanh(c)/c),  g(c) = 1 - 3/c^2 + 3 tanh(c)/c^3
% ell -> 0: Navier two-wall slip; b -> 0: Lin et al. (2013) couple-stress Phi.
c = h./(2*max(ell,1e-300));
s = 0.5*h.*wfun(c);
Phi = h.^3.*gfun(c) + 6*h.*s.*(h.*(b1+b2)+4*b1.*b2)./(h+b1+b2);
end

function Phi = flowBVP(h,ell,b1,b2)
% independent check: solve u'' - ell^2 u'''' = -1 (G = mu = 1) as two second-order BVPs:
%   w = u'':  ell^2 w'' - w = 1 , w(0) = w(h) = 0 (vanishing couple stress);
%   u'' = w with Navier slip u(0) = b1 u'(0), u(h) = -b2 u'(h);  Phi = 12 int u dy
n = 1500; y = linspace(0,h,n+1)'; dy = y(2)-y(1);
if ell < 1e-9
    Phi = h^2*(h^2+4*h*(b1+b2)+12*b1*b2)/(h+b1+b2); return;
end
e = ones(n+1,1);
Aw = full(spdiags([e*ell^2/dy^2, -2*e*ell^2/dy^2 - e, e*ell^2/dy^2],-1:1,n+1,n+1));
Aw(1,:) = 0; Aw(1,1) = 1; Aw(n+1,:) = 0; Aw(n+1,n+1) = 1;
rw = ones(n+1,1); rw(1) = 0; rw(n+1) = 0; w = Aw\rw;
Au = full(spdiags([e/dy^2, -2*e/dy^2, e/dy^2],-1:1,n+1,n+1));
Au(1,:) = 0; Au(1,1) = 1 + 3*b1/(2*dy); Au(1,2) = -4*b1/(2*dy); Au(1,3) = b1/(2*dy);
Au(n+1,:) = 0; Au(n+1,n+1) = 1 + 3*b2/(2*dy); Au(n+1,n) = -4*b2/(2*dy); Au(n+1,n-1) = b2/(2*dy);
ru = w; ru(1) = 0; ru(n+1) = 0; u = Au\ru;
Phi = 12*trapz(y,u);
end

function v = numClosed(H,k)
I = integral(@(u) tanh(u*H/(2*k))./u.^6,1,Inf,'AbsTol',1e-16,'RelTol',1e-13);
v = (H^3-6*k^2*H) + 48*k^3*I;
end

% --------------------------------------------------------- viscosity laws
function mp = mlawS(S,p)
if strcmp(S.law,'barus')
    mp = exp(min(S.alpha*p,700));
else
    L = log(S.muMax/S.mu0); pp = max(p,0);
    mp = exp(L*(1-exp(min(-S.alpha*pp/L,50)))).*exp(S.alpha*min(p,0));   % saturates for p>0, Barus for suction
end
end

function [dpdr,F,mu] = rhsP(S,h,p,r,M,ell0,b1,b2,darcy)
mp = mlawS(S,p); m = M.*mp;
if strcmp(S.ellCoupling,'eta'), ell = ell0./sqrt(m); else, ell = ell0./sqrt(mp); end
hl = h + S.compl*p;                                   % local film thickness with elastic liner (quasi-static)
Phi = flowGap(hl,ell,b1,b2);
if S.sigma > 0, Phi = Phi.*(1 - S.C_PC*exp(-S.r_PC*hl/S.sigma)); end
Phi = Phi + darcy;
F = Phi/h^3; mu = S.mu0*m;
dpdr = -6*S.U*r.*mu./Phi;
dpdr(Phi <= 0) = NaN;
end

function [b1,b2,darcy] = wallPars(S)
if S.porous
    ke = S.Ssat*S.kappa; b1 = sqrt(ke)/S.alphaBJS; b2 = S.slip; darcy = 12*ke*S.hp;
else
    b1 = S.slip; b2 = S.slip; darcy = 0;
end
end

function [p,aux,ok] = pressureODE(S,h,M,ell0,b1,b2,darcy)
N = S.N; r = linspace(0,S.a,N+1)'; dr = S.a/N; p = zeros(N+1,1); ok = true;
f = @(pp,rr,MM,LL) rhsP(S,h,pp,rr,MM,LL,b1,b2,darcy);
for i = N+1:-1:2
    r1 = r(i); r0 = r(i-1); rm = 0.5*(r1+r0); Mm = 0.5*(M(i)+M(i-1)); Lm = 0.5*(ell0(i)+ell0(i-1));
    k1 = f(p(i),r1,M(i),ell0(i));
    k2 = f(p(i)-0.5*dr*k1,rm,Mm,Lm);
    k3 = f(p(i)-0.5*dr*k2,rm,Mm,Lm);
    k4 = f(p(i)-dr*k3,r0,M(i-1),ell0(i-1));
    p(i-1) = p(i) - dr/6*(k1+2*k2+2*k3+k4);
    if ~isfinite(p(i-1)) || (strcmp(S.law,'barus') && S.alpha*p(i-1) > 400)
        ok = false; aux = []; return;
    end
end
[dp,F,mu] = f(p,r,M,ell0);
aux.dp = dp; aux.F = F; aux.mu = mu;
if ~all(isfinite(dp)), ok = false; end
end

% ----------------------------------------------------------- energy balance
function th = thermalField(S,r,h,aux,phmean)
% quasi-steady energy balance  rho cp (U r/2) th' = D - 2 h_w th  (D = dissipation per unit area).
% Exact exponential integrator for D linear between nodes; stable for any k = 4 h_w/(rho cp U) (ratio <= 1).
rhoCp = (1-phmean)*S.rho*S.cp + phmean*S.rhoP*S.cpP;
rc = rhoCp*S.U; k = 4*S.hw/rc;
D = h^3*aux.F.*aux.dp.^2./(12*aux.mu);
th = zeros(size(r)); dr = r(2)-r(1);
for j = 2:numel(r)
    r0 = r(j-1); r1 = r(j); q = r0/r1; qk = q^k; qk1 = q^(k+1);
    beta = (D(j)-D(j-1))/dr; alp = D(j-1) - beta*r0;
    th(j) = qk*th(j-1) + (2/rc)*(alp*(1-qk)/k + beta*r1*(1-qk1)/(k+1));
end
end

% -------------------------------------------------------- particle transport
function ph = migrate(S,r,h,ph0,tsq,Ra,th)
% h phi_t + Q phi_r = -(1/r)(r h J)_r ,  Q = U r/2
% J = -[D_B + a^2 gd (Kc phi + Kmu phi^2 dln mu/dphi)] phi_r - Kc a^2 phi^2 gd_r  (Phillips et al. 1992)
kB = 1.380649e-23; N = numel(r)-1; dr = r(2)-r(1); ap = 0.5*S.dp*Ra; ph = ph0(:);
if tsq <= 0, return; end
rmid = 0.5*(r(1:end-1)+r(2:end));
t = 0; nst = 0; hmin = max(h,1e-9);
while t < tsq && nst < 40000
    hh = max(S.hStart - S.U*t, hmin); gd = 3*S.U*r/hh^2; gdr = 3*S.U/hh^2;
    mub = S.mu0*exp(-S.beta*th); DB = kB*S.T0./(6*pi*mub*ap);
    dl = 1e-5; Sv = S; 
    dlnmu = (log(relVisc(Sv,ph+dl)) - log(relVisc(Sv,max(ph-dl,0)))) ./ ((ph+dl)-max(ph-dl,0));
    Dtot = DB + S.migAmp*ap^2*gd.*(S.Kc*ph + S.Kmu*ph.^2.*dlnmu);
    Dm = 0.5*(Dtot(1:end-1)+Dtot(2:end)); phm = 0.5*(ph(1:end-1)+ph(2:end));
    J = -Dm.*diff(ph)/dr - S.migAmp*S.Kc*ap^2*phm.^2*gdr;           % at midpoints
    rJ = rmid.*J;
    div = zeros(N+1,1);
    div(1) = 4*J(1)/dr;                                               % axis cell
    div(2:N) = (rJ(2:N)-rJ(1:N-1))./(r(2:N)*dr);
    div(N+1) = (0 - rJ(N))/(0.5*(r(N+1)^2-rmid(N)^2));                % no diffusive flux through the rim
    adv = zeros(N+1,1); adv(2:end) = (S.U*r(2:end)/(2*hh)).*diff(ph)/dr;
    dt = 0.4*min(dr/max(S.U*S.a/(2*hh),eps), dr^2/max(max(Dtot),eps));
    dt = min(dt, tsq-t);
    ph = ph + dt*(-adv - div);
    ph = min(max(ph,1e-12),0.999*S.phim);
    t = t + dt; nst = nst + 1;
end
end

% -------------------------------------------------------- the coupled solver
function out = solveCoupled(S,h)
N = S.N; r = linspace(0,S.a,N+1)'; [b1,b2,darcy] = wallPars(S);
tsq = max(S.hStart-h,0)/S.U; Ra = aggRatio(S,tsq); S.RaEff = Ra;
ell0s = ellBase(S,Ra);
ph = S.phi*ones(N+1,1); th = zeros(N+1,1); out.exist = false; out.iter = 0;
coupled = S.thermal || S.migration; dphiRel = 0; p = []; aux = [];
for it = 1:S.maxIt
    M = relVisc(S,ph).*exp(-S.beta*th).*carreauF(S,r,h);
    ell0 = ell0s*sqrt(1+S.gc*ph);
    [p,aux,ok] = pressureODE(S,h,M,ell0,b1,b2,darcy);
    if ~ok, out.W = Inf; out.P0 = Inf; out.p = p; out.th = th; out.phi = ph; out.r = r; out.Ra = Ra; out.ell0 = ell0s; ...
        out.phi_a_ratio = NaN; out.mu_r = NaN; out.dphiRel = 0; out.iter = it; return; end
    if ~coupled, break; end
    thn = th; phn = ph;
    if S.thermal, thn = thermalField(S,r,h,aux,mean(ph)); end
    if S.migration, phn = migrate(S,r,h,S.phi*ones(N+1,1),tsq,Ra,th); end
    err = max(abs(thn-th))/(1+max(abs(thn))) + max(abs(phn-ph))/S.phi;
    th = (1-S.relax)*th + S.relax*thn; ph = (1-S.relax)*ph + S.relax*phn;
    if err < S.tol, out.iter = it; break; end
    out.iter = it;
end
if ~coupled, out.iter = 1; end
if S.migration, dphiRel = max(abs(ph-S.phi))/S.phi; end
out.exist = true; out.p = p; out.r = r; out.th = th; out.phi = ph; out.Ra = Ra; out.ell0 = ell0s;
out.W = 2*pi*simpson(r,p.*r); out.P0 = p(1); out.dphiRel = dphiRel;
out.phi_a_ratio = S.phi*Ra^(3-S.D)/S.phim; out.mu_r = relVisc(S,S.phi);
end

function q = simpson(x,f)
% composite Simpson (uniform grid, even number of intervals)
n = numel(x)-1; h = x(2)-x(1);
if mod(n,2) == 1, q = trapz(x,f); return; end
q = h/3*(f(1) + f(end) + 4*sum(f(2:2:end-1)) + 2*sum(f(3:2:end-2)));
end

function a = vcritAxi(S,h,tol)
% largest Barus number v = alpha (unit scales) for which a solution exists
S.law = 'barus'; S.thermal = false; S.migration = false;
ex = @(al) solveCoupled(setfield(S,'alpha',al),h).exist; %#ok<SFLD>
lo = 0; hi = 1e-3; it = 0;
while ex(hi) && it < 60, lo = hi; hi = 2*hi; it = it+1; end
while (hi-lo) > tol*hi
    mid = 0.5*(lo+hi); if ex(mid), lo = mid; else, hi = mid; end
end
a = 0.5*(lo+hi);
end

function Uc = vcritU(S,h)
S.law = 'barus'; S.thermal = false; S.migration = false;
ex = @(U) solveCoupled(setfield(S,'U',U),h).exist; %#ok<SFLD>
lo = 1e-7; hi = 1e-6; it = 0;
while ex(hi) && it < 80, lo = hi; hi = 2*hi; it = it+1; end
while (hi-lo) > 1e-5*hi
    mid = 0.5*(lo+hi); if ex(mid), lo = mid; else, hi = mid; end
end
Uc = 0.5*(lo+hi);
end

function Ue = envU(S,h,pLim)
% squeeze speed at which the Barus-law centre pressure reaches pLim
S.law = 'barus'; S.thermal = false; S.migration = false;
f = @(lu) pcen(S,h,10^lu) - pLim;
Ue = 10^fzero(f,[-7 log10(vcritU(S,h))-1e-6]);
end
function P = pcen(S,h,U)
S.U = U; o = solveCoupled(S,h); P = o.P0; if ~isfinite(P), P = 1e30; end
end

function Y = capNaN(Y)
% replace failed solves (NaN) by the largest finite value so that variance estimates stay defined (count reported by caller)
if all(isnan(Y)), error('all Sobol solves failed'); end
Y(isnan(Y)) = max(Y(~isnan(Y)));
end

function C = fitQuad(Z,Y)
B = quadBasis(Z); C = (B'*B + 1e-8*eye(size(B,2)))\(B'*Y);
end
function Y = evalQuad(Z,C,~)
Y = quadBasis(Z)*C;
end
function B = quadBasis(Z)
[n,d] = size(Z); nb = 1+d+d*(d+1)/2; B = ones(n,nb); B(:,2:d+1) = Z; k = d+2;
for i = 1:d, for j = i:d, B(:,k) = Z(:,i).*Z(:,j); k = k+1; end, end
end

function Y = sobolEval(S0,X)
n = size(X,1); Y = zeros(n,1);
for i = 1:n
    S = S0; S.law = 'bounded'; S.visc = 'KDagg'; S.thermal = true; S.migration = true; S.aggMode = 'fixed';
    S.N = min(S0.N,80);
    S.phi = X(i,1); S.phim = X(i,2); S.Ra = X(i,3); S.D = X(i,4); S.dp = X(i,5); S.ellScale = X(i,6);
    S.slip = X(i,7); S.sigma = X(i,8); S.kappa = 10^X(i,9); S.hp = X(i,10); S.Ssat = X(i,11); S.U = 10^X(i,12);
    S.porous = true; h = 0.4*S.h0;
    o = solveCoupled(S,h);
    Wiso = 1.5*pi*S.mu0*S.U*S.a^4/h^3;
    if o.exist && isfinite(o.W) && o.W > 0, Y(i) = log(o.W/Wiso); else, Y(i) = NaN; end
end
% failed solves are returned as NaN; the caller counts them and caps them explicitly (capNaN)
end

function [x,w] = gaussleg(n)
b = (1:n-1)./sqrt(4*(1:n-1).^2-1); J = diag(b,1)+diag(b,-1); [V,D] = eig(J); [x,i] = sort(diag(D)); w = 2*V(1,i)'.^2;
end

function DB = particleDB()
% Abass, Ahmed, Kadhim (2023) Table 2 and aggregate ratios (from their ref. [16], Tushar et al. 2021)
DB.name = {'TiO2','Al2O3','CuO'}; DB.rho = [4230 3970 6310]; DB.cp = [710 765 565]; DB.k = [8.4 36 18]; DB.Ra = [4 3.33 3.5];
end

function flushOut()
try, fflush(stdout); catch, end %#ok<NOCOMMA>
end

function F = carreauF(S,r,h)
if S.lam > 0, gd = 3*S.U*abs(r)/h^2; F = (1+(S.lam*gd).^2).^((S.nCar-1)/2); else, F = ones(size(r)); end
end

function U = speedForLoad(S,h,Wload)
% squeeze speed at which the film carries Wload at thickness h
Ui = Wload*h^3/(1.5*pi*S.mu0*relVisc(S,S.phi)*S.a^4); lo = log(Ui*1e-3); hi = log(Ui*10);
for it = 1:28
    mid = 0.5*(lo+hi); Sj = S; Sj.U = exp(mid); o = solveCoupled(Sj,h);
    if o.exist && o.W < Wload, lo = mid; else, hi = mid; end
end
U = exp(0.5*(lo+hi));
end

function [t,tf,hq,Uq] = approachTime(S,Wload,ha,hb,nq)
% t: true constant-load approach time; tf: time with v frozen at its initial value (v12 convention)
[x,w] = gaussleg(nq); hq = 0.5*(ha-hb)*x + 0.5*(ha+hb); Uq = zeros(size(hq));
for j = 1:nq, Uq(j) = speedForLoad(S,hq(j),Wload); end
t = 0.5*(ha-hb)*sum(w(:)./Uq(:));
if nargout < 2, tf = []; return; end
Ua = speedForLoad(S,ha,Wload); Wper = zeros(size(hq));
for j = 1:nq, Sj = S; Sj.U = Ua; Wper(j) = solveCoupled(Sj,hq(j)).W/Ua; end
tf = 0.5*(ha-hb)*sum(w(:).*Wper(:))/Wload;
end

function [p,r] = sepInt(S,h,rc,M,l0,b1,b2,dar)
% separation with a gas cavity r < rc: inward flux 2 pi r Q = pi (r^2-rc^2) Us  (S.U = -Us)
N = S.N; r = linspace(rc,S.a,N+1)'; dr = (S.a-rc)/N; p = zeros(N+1,1);
rw = @(x) (x.^2-rc^2)./max(x,1e-12);
f = @(pp,x) rhsP(S,h,pp,rw(x),M,l0,b1,b2,dar);
for i = N+1:-1:2
    r1 = r(i); r0 = r(i-1); rm = 0.5*(r1+r0);
    k1 = f(p(i),r1); k2 = f(p(i)-0.5*dr*k1,rm); k3 = f(p(i)-0.5*dr*k2,rm); k4 = f(p(i)-dr*k3,r0);
    p(i-1) = p(i) - dr/6*(k1+2*k2+2*k3+k4);
end
end

function out = sepSolve(S,h,Us)
pcav = -(101325-29185); S.RaEff = aggRatio(S,0); M = relVisc(S,S.phi); l0 = ellBase(S,S.RaEff);
[b1,b2,dar] = wallPars(S); S.U = -Us;
[pf,rf] = sepInt(S,h,0,M,l0,b1,b2,dar); out.Wfull = 2*pi*simpson(rf,pf.*rf);
if pf(1) >= pcav
    out.rc = 0; out.p = pf; out.r = rf; out.pmin = pf(1); out.W = out.Wfull; return;
end
lo = 0; hi = S.a*(1-1e-6);
for it = 1:45
    mid = 0.5*(lo+hi); pm = sepInt(S,h,mid,M,l0,b1,b2,dar);
    if pm(1) < pcav, lo = mid; else, hi = mid; end
end
rc = 0.5*(lo+hi); [p,r] = sepInt(S,h,rc,M,l0,b1,b2,dar);
out.rc = rc; out.p = p; out.r = r; out.pmin = pcav; out.W = 2*pi*(pcav*rc^2/2 + simpson(r,p.*r));
end

function out = thermalTransient(S,ha,hb,nT)
% constant squeeze speed U from ha to hb; implicit upwind energy balance, explicit p-theta coupling
N = S.N; r = linspace(0,S.a,N+1)'; dr = S.a/N; ts = (ha-hb)/S.U; dt = ts/nT; th = zeros(N+1,1);
[b1,b2,dar] = wallPars(S); S.RaEff = aggRatio(S,0); ell0 = ellBase(S,S.RaEff)*ones(N+1,1); ph = S.phi*ones(N+1,1);
rhoCp = (1-S.phi)*S.rho*S.cp + S.phi*S.rhoP*S.cpP; out.W = zeros(nT,1); out.thm = zeros(nT,1); out.h = zeros(nT,1);
for n = 1:nT
    h = ha - S.U*(n-0.5)*dt; M = relVisc(S,ph).*exp(-S.beta*th).*carreauF(S,r,h);
    [p,aux,ok] = pressureODE(S,h,M,ell0,b1,b2,dar);
    if ~ok, out.W(n:end) = Inf; out.thm(n:end) = Inf; out.h(n:end) = h; break; end
    D = h^3*aux.F.*aux.dp.^2./(12*aux.mu); u = S.U*r/(2*h); s = 2*S.hw/(rhoCp*h); src = D/(rhoCp*h); thn = th;
    for j = 1:N+1
        if j == 1, thn(j) = (th(j)/dt + src(j))/(1/dt + s);
        else, thn(j) = (th(j)/dt + u(j)*thn(j-1)/dr + src(j))/(1/dt + u(j)/dr + s); end
    end
    th = thn; out.W(n) = 2*pi*simpson(r,p.*r); out.thm(n) = max(th); out.h(n) = h;
end
end

function [th,thc] = thermalWithConduction(S,r,h,aux,phmean,kover)
% steady energy balance with radial conduction (upwind advection): h k (th''+th'/r) - rho cp Q th' - 2 h_w th + D = 0
N = numel(r)-1; dr = r(2)-r(1); rhoCp = (1-phmean)*S.rho*S.cp + phmean*S.rhoP*S.cpP;
kf = S.kth; kp = S.kP; knf = kf*((kp+2*kf)-2*phmean*(kf-kp))/((kp+2*kf)+phmean*(kf-kp));
if ~isempty(kover), knf = kover; end
D = h^3*aux.F.*aux.dp.^2./(12*aux.mu); Q = S.U*r/2; th = thermalField(S,r,h,aux,phmean);
A = zeros(N+1); b = -D;
for j = 2:N
    c = h*knf/dr^2; a1 = h*knf/(2*r(j)*dr); ad = rhoCp*Q(j)/dr;
    A(j,j-1) = c - a1 + ad; A(j,j) = -2*c - ad - 2*S.hw; A(j,j+1) = c + a1;
end
A(1,1) = -4*h*knf/dr^2 - 2*S.hw; A(1,2) = 4*h*knf/dr^2; A(N+1,N+1) = 1; A(N+1,N) = -1; b(N+1) = 0;
thc = A\b;
end

% ================================================================== 2-D solver
function [D,x] = chebD(N)
x = cos(pi*(0:N)'/N); c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1); dX = X-X';
D = (c*(1./c)')./(dX+eye(N+1)); D = D - diag(sum(D,2));
end

function w = clencurt(N)
theta = pi*(0:N)'/N; w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
if mod(N,2)==0
    w(1) = 1/(N^2-1); w(N+1) = w(1);
    for k = 1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    v = v - cos(N*theta(ii))/(N^2-1);
else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k = 1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
end
w(ii) = 2*v/N;
end

function sp = setup2D(N,Lam)
[D,x] = chebD(N); I = eye(N+1); [X,Y] = ndgrid(x,x);
sp.N = N; sp.Lam = Lam; sp.x = x; sp.X = X(:); sp.Y = Y(:); sp.n = (N+1)^2;
sp.Dx = kron(I,D); sp.Dy = kron(D,I)/Lam;
sp.bnd = abs(sp.X)>1-1e-12 | abs(sp.Y)>1-1e-12; sp.ib = find(sp.bnd);
sp.w = clencurt(N);
end

function u = poissonU1(x,y,Lam)
u = (8*poisU(x,y,Lam,4001) - poisU(x,y,Lam,2001))/7;
end
function u = poisU(x,y,Lam,M)
m = 1:2:M-1; n = m; cm = (-1).^((m-1)/2);
lam = (m'*pi/2).^2 + (n*pi/(2*Lam)).^2;
coef = 16*(cm'*cm)./(pi^2*(m'*n).*lam);
Cx = cos(x(:)*m*pi/2); Cy = cos(y(:)*n*pi/(2*Lam));
u = sum((Cx*coef).*Cy,2); u = reshape(u,size(x));
end

function P2 = pars2D(k,v,law)
% dimensionless: lengths by h0, pressure by mu0 U a^2/h0^3, v = alpha mu0 U a^2/h0^3
P2.ell0 = k; P2.v = v; P2.law = law; P2.L = log(1e10/0.065); P2.M = 1;   % L = ln(mu_max/mu0), used only by the bounded law
P2.coupling = 'eta'; P2.b1 = 0; P2.b2 = 0; P2.darcy = 0; P2.sigma = 0; P2.C_PC = 0.9; P2.r_PC = 0.56;
P2.H = 1; P2.fR = 1;
end

function t = tau2v(p,P2)
if strcmp(P2.law,'barus'), mp = exp(min(P2.v*p,700)); else, mp = exp(P2.L*(1-exp(min(-P2.v*p/P2.L,50)))); end
m = P2.M.*mp;
if strcmp(P2.coupling,'eta'), ell = P2.ell0./sqrt(m); else, ell = P2.ell0./sqrt(mp); end
Phi = flowGap(P2.H,ell,P2.b1,P2.b2).*P2.fR + P2.darcy;
t = Phi./m;
end

function [tau,tp] = tau2(p,P2)
tau = tau2v(p,P2); d = 1e-6*max(1,abs(p));
tp = (tau2v(p+d,P2)-tau2v(p-d,P2))./(2*d);
end

function [p,ok] = solve2Du(sp,P2,p0)
n = sp.n; bnd = sp.bnd; ib = sp.ib;
if isempty(p0), p = zeros(n,1); else, p = p0; end
p(bnd) = 0; ok = false;
for it = 1:30
    [tau,tp] = tau2(p,P2); px = sp.Dx*p; py = sp.Dy*p;
    F = sp.Dx*(tau.*px) + sp.Dy*(tau.*py) + 12; F(bnd) = p(bnd);
    if ~all(isfinite(F)), return; end
    if norm(F(~bnd),inf) < 1e-10*12, ok = true; return; end
    J = sp.Dx*(diag(tp.*px) + diag(tau)*sp.Dx) + sp.Dy*(diag(tp.*py) + diag(tau)*sp.Dy);
    J(bnd,:) = 0; J(sub2ind(size(J),ib,ib)) = 1;
    dp = -(J\F); if ~all(isfinite(dp)), return; end
    lam = 1; n0 = norm(F); acc = false;
    while lam > 1e-4
        pn = p+lam*dp; tn = tau2v(pn,P2);
        Fn = sp.Dx*(tn.*(sp.Dx*pn)) + sp.Dy*(tn.*(sp.Dy*pn)) + 12; Fn(bnd) = pn(bnd);
        if all(isfinite(Fn)) && norm(Fn) < (1-1e-4*lam)*n0, acc = true; break; end
        lam = lam/2;
    end
    if ~acc, return; end
    p = p+lam*dp;
end
end

function [p,ok,vreach] = contin2D(sp,P2)
% natural continuation in v from the linear solution; vreach = furthest v reached
vT = P2.v; P0 = P2; P0.v = 0; [p,ok0] = solve2Du(sp,P0,[]); vreach = 0; ok = false;
if ~ok0, return; end
if vT == 0, ok = true; return; end
vc = 0; dv = vT/4; cnt = 0;
while vc < vT && cnt < 300
    cnt = cnt+1; vt = min(vT,vc+dv); Pt = P2; Pt.v = vt;
    [p2,ok2] = solve2Du(sp,Pt,p);
    if ok2, p = p2; vc = vt; dv = dv*1.5; else, dv = dv/2; if dv < 1e-6*vT, break; end, end
end
vreach = vc; ok = (vc >= vT*(1-1e-12));
end

function vc = vcrit2D(sp,P2,vmax)
% largest v with a Newton solution: continuation then bisection from the last good field
P2.v = vmax; [p,~,vr] = contin2D(sp,P2);
lo = vr; hi = vmax; pLo = p; if vr >= vmax, vc = vmax; return; end
[pLo,okLo] = solve2Du(sp,setfield(P2,'v',lo),p); %#ok<SFLD>
if ~okLo, pLo = p; end
for q = 1:30
    mid = 0.5*(lo+hi); [pm,okm] = solve2Du(sp,setfield(P2,'v',mid),pLo); %#ok<SFLD>
    if okm, lo = mid; pLo = pm; else, hi = mid; end
end
vc = lo;
end

function saveFig(f,base)
% v14: Times New Roman on every text object, 10 pt axes, PNG at 300 dpi and vector PDF, then close
set(findall(f,'-property','FontName'),'FontName','Times New Roman');
ax = findall(f,'type','axes');
for k = 1:numel(ax), set(ax(k),'FontSize',10,'LineWidth',0.7,'Box','on'); end
lg = findall(f,'type','legend'); for k = 1:numel(lg), set(lg(k),'Box','off'); end
set(f,'PaperPositionMode','auto');
print(f,[base '.png'],'-dpng','-r300');
print(f,[base '.pdf'],'-dpdf');
close(f);
end
