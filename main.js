// main.js - Refactored version (same behavior, cleaner structure)

// =======================
// Common DOM / number utils
// =======================
function num(id) {
  const el = document.getElementById(id);
  if (!el) return null;
  const v = parseFloat(el.value);
  return Number.isFinite(v) ? v : null;
}

function setText(id, value, digits = 3) {
  const el = document.getElementById(id);
  if (!el) return;
  if (!Number.isFinite(value)) {
    el.textContent = '-';
  } else {
    el.textContent = value.toFixed(digits);
  }
}

function getVal(id) {
  const el = document.getElementById(id);
  return el ? el.value : '';
}

function makeTsv(headers, row) {
  return headers.join('\t') + '\n' + row.join('\t');
}

function bindInputs(ids, handler, event = 'input') {
  ids.forEach(id => {
    const el = document.getElementById(id);
    if (el) el.addEventListener(event, handler);
  });
}

// Temperature helper
function cToK(tC) {
  return tC + 273.15;
}

// =======================
// 1. Mass / Volume
// =======================
function calcBlock1() {
  const mass = num('mass1');
  const mw = num('mw1');

  if (mass === null || mw === null || mw === 0) {
    setText('vol1', NaN);
    return;
  }
  const vol = (mass * 22.414) / mw;
  setText('vol1', vol);
}

function calcBlock2() {
  const vol = num('vol2_in');
  const mw = num('mw2');

  if (vol === null || mw === null || mw === 0) {
    setText('mass2', NaN);
    return;
  }
  const mass = (vol * mw) / 22.414;
  setText('mass2', mass);
}

// TSV
function getBlock1Tsv() {
  const header = [
    'Mass Flowrate (kg/h)',
    'MW (kg/kmol)',
    'Volume Flowrate (Nm3/h)'
  ];
  const row = [
    getVal('mass1'),
    getVal('mw1'),
    '=A2*22.414/B2'
  ];
  return makeTsv(header, row);
}

function getBlock2Tsv() {
  const header = [
    'Volume Flowrate (Nm3/h)',
    'MW (kg/kmol)',
    'Mass Flowrate (kg/h)'
  ];
  const row = [
    getVal('vol2_in'),
    getVal('mw2'),
    '=A2*B2/22.414'
  ];
  return makeTsv(header, row);
}

// =======================
// 2. Static Head
// =======================
function calcStaticHeadFromPressure() {
  const p = num('p_kgcm2');
  const sg = num('sg1');

  if (p === null || sg === null || sg === 0) {
    setText('head_m', NaN);
    return;
  }
  const head = (p * 10) / sg;
  setText('head_m', head);
}

function calcPressureFromStaticHead() {
  const head = num('head_m_in');
  const sg = num('sg2');

  if (head === null || sg === null) {
    setText('p_kgcm2_out', NaN);
    return;
  }
  const p = (head * sg) / 10;
  setText('p_kgcm2_out', p);
}

// TSV
function getBlock3Tsv() {
  const header = [
    'D. Pressure (kg/cm2)',
    'S.G (-)',
    'Head (m)'
  ];
  const row = [
    getVal('p_kgcm2'),
    getVal('sg1'),
    '=A2*10/B2'
  ];
  return makeTsv(header, row);
}

function getBlock4Tsv() {
  const header = [
    'Head (m)',
    'S.G (-)',
    'D. Pressure (kg/cm2)'
  ];
  const row = [
    getVal('head_m_in'),
    getVal('sg2'),
    '=A2*B2/10'
  ];
  return makeTsv(header, row);
}

// =======================
// 3. Vessel Volume
// =======================
// Vertical
function calcVerticalVessel() {
  const headType = getVal('v_head_type');
  const diam_mm = num('v_diam');
  const level_mm = num('v_level');

  if (
    diam_mm === null ||
    level_mm === null ||
    diam_mm <= 0 ||
    level_mm <= 0
  ) {
    setText('v_volume', NaN);
    return;
  }

  const D = diam_mm / 1000;
  const Level = level_mm / 1000;
  const R = D / 2;

  let volume;
  if (headType === '2:1E') {
    volume =
      ((2 / 3) / 2) * Math.PI * Math.pow(R, 3) +
      Math.PI * (Math.pow(D, 2) / 4) * Level;
  } else if (headType === 'HS') {
    volume = (Math.PI / 12) * Math.pow(D, 3) +
      Math.PI * (Math.pow(D, 2) / 4) * Level;
  } else {
    // FLAT
    volume = Math.PI * (Math.pow(D, 2) / 4) * Level;
  }

  setText('v_volume', volume);
}

// Horizontal
function calcHorizontalVessel() {
  const headType = getVal('h_head_type'); // "2:1E" / "HS"
  const diam_mm = num('h_diam');
  const len_mm = num('h_length');
  const level_mm = num('h_level');

  if (
    diam_mm === null ||
    len_mm === null ||
    level_mm === null ||
    diam_mm <= 0 ||
    len_mm <= 0 ||
    level_mm < 0
  ) {
    setText('h_volume', NaN);
    return;
  }

  const B = diam_mm;
  const C = len_mm;
  const D = level_mm;

  const acos_term = Math.acos(1 - 2 * (D / B));
  const factor =
    (Math.PI / 8 +
      Math.pow(0.5, 2) *
        ((acos_term - Math.acos(0)) -
          0.5 *
            (Math.sin(2 * acos_term) - Math.sin(2 * Math.acos(0))))) /
    (Math.PI / 4);

  const cyl_full = (Math.PI / 4 * B * B * C) / 1e9;

  let volume;
  if (headType === '2:1E') {
    const head_liq =
      ((Math.PI / 4 / 2 * B * D * D * (1 - (2 * D) / (3 * B))) / 1e9) * 2;
    volume = head_liq + factor * cyl_full;
  } else {
    // HS
    const head_liq = (2 * Math.PI / 24 * (D * D) * (3 * B - 2 * D) * 2) / 1e9;
    volume = head_liq + factor * cyl_full;
  }

  setText('h_volume', volume);
}

// TSV
function getBlock5Tsv() {
  const type = getVal('v_head_type');
  const diam = getVal('v_diam');
  const len = getVal('v_length');
  const level = getVal('v_level');

  const header = [
    'Head Type',
    'Diameter (mm)',
    'Length (mm)',
    'Liquid Level (mm)',
    'Volume (m3)'
  ];

  let formula;
  if (type === '2:1E') {
    formula =
      '=((2/3)/2*PI()*(B2/2/1000)^3+PI()*((B2/1000)^2)/4*(D2/1000))';
  } else if (type === 'HS') {
    formula = '=(PI()/12*(B2/1000)^3)+PI()*((B2/1000)^2)/4*(D2/1000)';
  } else {
    formula = '=PI()*((B2/1000)^2)/4*(D2/1000)';
  }

  const row = [type, diam, len, level, formula];
  return makeTsv(header, row);
}

function getBlock6Tsv() {
  const type = getVal('h_head_type');
  const diam = getVal('h_diam');
  const len = getVal('h_length');
  const level = getVal('h_level');

  const header = [
    'Type',
    'Diameter (mm)',
    'Length (mm)',
    'Liquid Level (mm)',
    'Volume (m3)'
  ];

  let formula;
  if (type === '2:1E') {
    formula =
      '=((PI()/4/2*B2*D2^2*(1-2*D2/(3*B2)))/10^9)*2+' +
      '((PI()/8+0.5^2*((ACOS(1-2*(D2/B2))-ACOS(0))-' +
      '0.5*(SIN(2*ACOS(1-2*(D2/B2)))-SIN(2*ACOS(0)))))/(PI()/4))' +
      '*(PI()/4*B2*B2*C2/10^9)';
  } else {
    // HS
    formula =
      '=(2*PI()/24*(D2^2)*(3*B2-2*D2)*2)/10^9+' +
      '((PI()/8+0.5^2*((ACOS(1-2*(D2/B2))-ACOS(0))-' +
      '0.5*(SIN(2*ACOS(1-2*(D2/B2)))-SIN(2*ACOS(0)))))/(PI()/4))' +
      '*(PI()/4*B2*B2*C2/10^9)';
  }

  const row = [type, diam, len, level, formula];
  return makeTsv(header, row);
}

// =======================
// 4. Combined Gas Law
// =======================
// V2
function calcCombinedV2() {
  const atm = num('atm_v2');
  const p1 = num('p1_v2');
  const t1 = num('t1_v2');
  const v1 = num('v1_v2');
  const p2 = num('p2_v2');
  const t2 = num('t2_v2');

  if (
    [atm, p1, t1, v1, p2, t2].some(v => v === null) ||
    p2 + atm === 0
  ) {
    setText('v2_out', NaN);
    return;
  }

  const v2 = ((p1 + atm) * v1 * cToK(t2)) / ((p2 + atm) * cToK(t1));
  setText('v2_out', v2);
}

// P2
function calcCombinedP2() {
  const atm = num('atm_p2');
  const p1 = num('p1_p2');
  const t1 = num('t1_p2');
  const v1 = num('v1_p2');
  const t2 = num('t2_p2');
  const v2 = num('v2_p2');

  if ([atm, p1, t1, v1, t2, v2].some(v => v === null) || v2 === 0) {
    setText('p2_out', NaN);
    return;
  }

  const p2_abs = ((p1 + atm) * v1 * cToK(t2)) / (v2 * cToK(t1));
  const p2_g = p2_abs - atm;
  setText('p2_out', p2_g);
}

// T2
function calcCombinedT2() {
  const atm = num('atm_t2');
  const p1 = num('p1_t2');
  const t1 = num('t1_t2');
  const v1 = num('v1_t2');
  const p2 = num('p2_t2');
  const v2 = num('v2_t2');

  if (
    [atm, p1, t1, v1, p2, v2].some(v => v === null) ||
    p1 + atm === 0
  ) {
    setText('t2_out', NaN);
    return;
  }

  const t2K = ((p2 + atm) * v2 * cToK(t1)) / ((p1 + atm) * v1);
  const t2C = t2K - 273.15;
  setText('t2_out', t2C);
}

// TSV
function getBlock7Tsv() {
  const header = [
    '1 ATM (kg/cm2g)',
    'P1 (kg/cm2g)',
    'T1 (degC)',
    'V1 (m3)',
    'P2 (kg/cm2g)',
    'T2 (degC)',
    'V2 (m3)'
  ];
  const row = [
    getVal('atm_v2'),
    getVal('p1_v2'),
    getVal('t1_v2'),
    getVal('v1_v2'),
    getVal('p2_v2'),
    getVal('t2_v2'),
    '=((B2+A2)*D2*(F2+273.15))/((E2+A2)*(C2+273.15))'
  ];
  return makeTsv(header, row);
}

function getBlock8Tsv() {
  const header = [
    '1 ATM (kg/cm2g)',
    'P1 (kg/cm2g)',
    'T1 (degC)',
    'V1 (m3)',
    'P2 (kg/cm2g)',
    'T2 (degC)',
    'V2 (m3)'
  ];
  const row = [
    getVal('atm_p2'),
    getVal('p1_p2'),
    getVal('t1_p2'),
    getVal('v1_p2'),
    '=((B2+A2)*D2*(F2+273.15))/(G2*(C2+273.15))-A2',
    getVal('t2_p2'),
    getVal('v2_p2')
  ];
  return makeTsv(header, row);
}

function getBlock9Tsv() {
  const header = [
    '1 ATM (kg/cm2g)',
    'P1 (kg/cm2g)',
    'T1 (degC)',
    'V1 (m3)',
    'P2 (kg/cm2g)',
    'T2 (degC)',
    'V2 (m3)'
  ];
  const row = [
    getVal('atm_t2'),
    getVal('p1_t2'),
    getVal('t1_t2'),
    getVal('v1_t2'),
    getVal('p2_t2'),
    '=((E2+A2)*G2*(C2+273.15))/((B2+A2)*D2)-273.15',
    getVal('v2_t2')
  ];
  return makeTsv(header, row);
}

// =======================
// 5. Gas Density
// =======================
function calcGasDensity() {
  const mw = num('gd_mw');
  const z = num('gd_z');
  const p = num('gd_p');
  const tC = num('gd_t');

  if (
    [mw, z, p, tC].some(v => v === null) ||
    z === 0
  ) {
    setText('gd_rho', NaN);
    return;
  }

  // Excel: =((C3+1.03323)*98.0665*A3)/(B3*8.314*(D3+273.15))
  const rho =
    ((p + 1.03323) * 98.0665 * mw) /
    (z * 8.314 * (tC + 273.15));
  setText('gd_rho', rho);
}

function getBlock10Tsv() {
  const header = [
    'MW (-)',
    'Comp. Factor Z (-)',
    'Pressure (kg/cm2g)',
    'Temperature (degC)',
    'Gas Density (kg/m3)'
  ];
  const row = [
    getVal('gd_mw'),
    getVal('gd_z'),
    getVal('gd_p'),
    getVal('gd_t'),
    '=((C2+1.03323)*98.0665*A2)/(B2*8.314*(D2+273.15))'
  ];
  return makeTsv(header, row);
}

// =======================
// 6. Line Sizing
// =======================
// common incompressible DP (100 m 기준)
function ls_incompressible_dp(massFlow_kg_h, rho, mu_cP, d_inch, rough_ft) {
  const mdot = massFlow_kg_h / 3600.0;
  const D = d_inch * 0.0254;
  const eps = rough_ft * 0.3048;
  const area = (Math.PI * D * D) / 4.0;
  const v = mdot / (rho * area);
  const mu = mu_cP * 1e-3; // Pa·s
  const Re = (rho * v * D) / mu;

  if (!isFinite(Re) || Re <= 0 || !isFinite(v)) {
    return { v: NaN, dp100: NaN };
  }

  const rel = eps / D;
  const f =
    0.25 /
    Math.pow(
      Math.log10(rel / 3.7 + 5.74 / Math.pow(Re, 0.9)),
      2
    );

  const dP_per_m_Pa = (f * rho * v * v) / (2 * D);
  const dP_100m_kgcm2 = (dP_per_m_Pa * 100.0) / 98066.5;
  return { v, dp100: dP_100m_kgcm2 };
}

// LS1: Liquid from Mass Flow
function calcLS1() {
  const m = num('ls1_mass');
  const rho = num('ls1_rho');
  const mu = num('ls1_mu');
  const d = num('ls1_d');
  const eps = num('ls1_eps');

  if (
    [m, rho, mu, d, eps].some(v => v === null) ||
    rho <= 0 ||
    d <= 0
  ) {
    ['ls1_q', 'ls1_dp100', 'ls1_v'].forEach(id =>
      setText(id, NaN)
    );
    return;
  }

  const q = m / rho; // m3/h
  const res = ls_incompressible_dp(m, rho, mu, d, eps);

  setText('ls1_q', q);
  setText('ls1_dp100', res.dp100);
  setText('ls1_v', res.v);
}

// LS2: Liquid from Vol Flow
function calcLS2() {
  const q = num('ls2_q');
  const rho = num('ls2_rho');
  const mu = num('ls2_mu');
  const d = num('ls2_d');
  const eps = num('ls2_eps');

  if (
    [q, rho, mu, d, eps].some(v => v === null) ||
    rho <= 0 ||
    d <= 0
  ) {
    ['ls2_mass', 'ls2_dp100', 'ls2_v'].forEach(id =>
      setText(id, NaN)
    );
    return;
  }

  const m = q * rho; // kg/h
  const res = ls_incompressible_dp(m, rho, mu, d, eps);

  setText('ls2_mass', m);
  setText('ls2_dp100', res.dp100);
  setText('ls2_v', res.v);
}

// Gas helpers
function ls_gas_density(mw, z, p_g_kgcm2, tC) {
  const p_abs = (p_g_kgcm2 + 1.03323) * 98066.5; // Pa
  const T = cToK(tC);
  return (p_abs * mw) / (z * 8.314 * 1000.0 * T);
}

function ls_sonic_speed(k, mw, tC) {
  const T = cToK(tC);
  const Rspec = (8.314 * 1000.0) / mw; // J/kg/K
  return Math.sqrt(k * Rspec * T);
}

// LS3: Vapor Darcy (incompressible approximation)
function calcLS3() {
  const m = num('ls3_mass');
  const mw = num('ls3_mw');
  const mu = num('ls3_mu');
  const z = num('ls3_z');
  const k = num('ls3_k');
  const side = num('ls3_side');
  const pKnown = num('ls3_pknown');
  const tC = num('ls3_t');
  const d = num('ls3_d');
  const eps = num('ls3_eps');
  const L = num('ls3_leq');

  const ids = [
    'ls3_pin',
    'ls3_pout',
    'ls3_dp',
    'ls3_dp100',
    'ls3_rho_in',
    'ls3_rho_out',
    'ls3_v_in',
    'ls3_v_out',
    'ls3_mach_in',
    'ls3_mach_out',
    'ls3_a'
  ];

  if (
    [m, mw, mu, z, k, side, pKnown, tC, d, eps, L].some(
      v => v === null
    ) ||
    mw <= 0 ||
    z <= 0 ||
    d <= 0 ||
    L <= 0
  ) {
    ids.forEach(id => setText(id, NaN));
    return;
  }

  const rhoKnown = ls_gas_density(mw, z, pKnown, tC);
  const res = ls_incompressible_dp(m, rhoKnown, mu, d, eps);
  const dp_total = res.dp100 * (L / 100.0);

  let p_in = pKnown;
  let p_out = pKnown;
  if (side === 1) {
    p_in = pKnown;
    p_out = pKnown - dp_total;
  } else {
    p_out = pKnown;
    p_in = pKnown + dp_total;
  }

  const rho_in = ls_gas_density(mw, z, p_in, tC);
  const rho_out = ls_gas_density(mw, z, p_out, tC);

  const D = d * 0.0254;
  const area = (Math.PI * D * D) / 4.0;
  const mdot = m / 3600.0;
  const v_in = mdot / (rho_in * area);
  const v_out = mdot / (rho_out * area);
  const a = ls_sonic_speed(k, mw, tC);
  const mach_in = v_in / a;
  const mach_out = v_out / a;

  setText('ls3_pin', p_in);
  setText('ls3_pout', p_out);
  setText('ls3_dp', dp_total);
  setText('ls3_dp100', res.dp100);
  setText('ls3_rho_in', rho_in, 2);
  setText('ls3_rho_out', rho_out, 2);
  setText('ls3_v_in', v_in, 2);
  setText('ls3_v_out', v_out, 2);
  setText('ls3_mach_in', mach_in, 3);
  setText('ls3_mach_out', mach_out, 3);
  setText('ls3_a', a, 1);
}

// TSV for Line Sizing
function getLS1Tsv() {
  const header = [
    'Mass Flow (kg/h)',
    'Density (kg/m3)',
    'Viscosity (cP)',
    'Press (kg/cm2g)',
    'Temp (degC)',
    'Pipe ID (inch)',
    'Rough (ft)',
    'Vol Flow (m3/h)',
    'DP/100m (kg/cm2)',
    'Velocity (m/s)'
  ];
  const row = [
    getVal('ls1_mass'),
    getVal('ls1_rho'),
    getVal('ls1_mu'),
    getVal('ls1_p'),
    getVal('ls1_t'),
    getVal('ls1_d'),
    getVal('ls1_eps'),
    '=A2/B2',
    '=0.25/(LOG10((G2*0.3048)/(3.7*(F2*0.0254))+5.74/(((B2*J2*(F2*0.0254))/(C2/1000))^0.9))^2)*B2*J2^2/(2*(F2*0.0254))*100/98066.5',
    '=4*(A2/(3600*B2))/(PI()*(F2*0.0254)^2)'
  ];
  return makeTsv(header, row);
}

function getLS2Tsv() {
  const header = [
    'Vol Flow (m3/h)',
    'Density (kg/m3)',
    'Viscosity (cP)',
    'Press (kg/cm2g)',
    'Temp (degC)',
    'Pipe ID (inch)',
    'Rough (ft)',
    'Mass Flow (kg/h)',
    'DP/100m (kg/cm2)',
    'Velocity (m/s)'
  ];
  const row = [
    getVal('ls2_q'),
    getVal('ls2_rho'),
    getVal('ls2_mu'),
    getVal('ls2_p'),
    getVal('ls2_t'),
    getVal('ls2_d'),
    getVal('ls2_eps'),
    '=A2*B2',
    '=0.25/(LOG10((G2*0.3048)/(3.7*(F2*0.0254))+5.74/(((B2*J2*(F2*0.0254))/(C2/1000))^0.9))^2)*B2*J2^2/(2*(F2*0.0254))*100/98066.5',
    '=4*(A2/3600)/(PI()*(F2*0.0254)^2)'
  ];
  return makeTsv(header, row);
}

function getLS3Tsv() {
  const header = [
    'Mass Flow (kg/h)',
    'MW (kg/kmol)',
    'Viscosity (cP)',
    'Z-factor',
    'k (Cp/Cv)',
    'Known Side',
    'Known P (kg/cm2g)',
    'Temp (degC)',
    'Pipe ID (inch)',
    'Rough (ft)',
    'EQ Length (m)',
    'P_in (kg/cm2g)',
    'P_out (kg/cm2g)',
    'Total dP (kg/cm2)',
    'DP/100m (kg/cm2)',
    'rho_in (kg/m3)',
    'rho_out (kg/m3)',
    'v_in (m/s)',
    'v_out (m/s)',
    'Mach_in',
    'Mach_out',
    'Sonic Vel (m/s)'
  ];

  const fPin = '=IF(F2=1,G2,G2+N2)';
  const fPout = '=IF(F2=1,G2-N2,G2)';
  const fDP100 =
    '=LET(' +
    'rhoKnown,((G2+1.03323)*98.0665*B2)/(D2*8.314*(H2+273.15)),' +
    'Dm,I2*0.0254,' +
    'area,PI()*Dm^2/4,' +
    'v,(A2/3600)/(rhoKnown*area),' +
    'Re,rhoKnown*v*Dm/(C2/1000),' +
    'rel,(J2*0.3048)/Dm,' +
    'f,0.25/(LOG10(rel/3.7+5.74/(Re^0.9))^2),' +
    'dPperM,f*rhoKnown*v^2/(2*Dm),' +
    'dPperM*100/98066.5)';
  const fDPtot = '=O2*K2/100';

  const fRhoIn =
    '=((L2+1.03323)*98.0665*B2)/(D2*8.314*(H2+273.15))';
  const fRhoOut =
    '=((M2+1.03323)*98.0665*B2)/(D2*8.314*(H2+273.15))';

  const fVIn = '=(A2/3600)/(P2*PI()*(I2*0.0254)^2/4)';
  const fVOut = '=(A2/3600)/(Q2*PI()*(I2*0.0254)^2/4)';

  const fA = '=SQRT(E2*8.314*1000/B2*(H2+273.15))';
  const fMachIn = '=R2/V2';
  const fMachOut = '=S2/V2';

  const row = [
    getVal('ls3_mass'),
    getVal('ls3_mw'),
    getVal('ls3_mu'),
    getVal('ls3_z'),
    getVal('ls3_k'),
    getVal('ls3_side'),
    getVal('ls3_pknown'),
    getVal('ls3_t'),
    getVal('ls3_d'),
    getVal('ls3_eps'),
    getVal('ls3_leq'),
    fPin,
    fPout,
    fDPtot,
    fDP100,
    fRhoIn,
    fRhoOut,
    fVIn,
    fVOut,
    fMachIn,
    fMachOut,
    fA
  ];

  // 여기서는 탭을 직접 써도 무방하지만 위와 통일
  return header.join('\t') + '\n' + row.join('\t');
}

// =======================
// 7. Clipboard
// =======================
function copyTsvToClipboard(tsv, label) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard
      .writeText(tsv)
      .then(() => {
        if (label) {
          alert(label + ' In Excel, select cell A1 and paste.');
        }
      })
      .catch(() => alert('Failed to copy to clipboard.'));
  } else {
    alert('This browser does not support the Clipboard API.');
  }
}

// =======================
// 8. Mode switching
// =======================
const sectionsByMode = {
  'mass-volume': 'section-mass-volume',
  'static-head': 'section-static-head',
  'vessel-volume': 'section-vessel-volume',
  'combined-gas': 'section-combined-gas',
  'gas-density': 'section-gas-density',
  'line-sizing': 'section-line-sizing'
};

function applyMode(mode) {
  Object.entries(sectionsByMode).forEach(([key, id]) => {
    const el = document.getElementById(id);
    if (!el) return;
    el.style.display = key === mode ? '' : 'none';
  });
}

// =======================
// 9. Init
// =======================
window.addEventListener('DOMContentLoaded', () => {
  const modeSelect = document.getElementById('calc-mode');
  if (modeSelect) {
    applyMode(modeSelect.value);
    modeSelect.addEventListener('change', function () {
      applyMode(this.value);
    });
  }

  // Mass/Volume
  bindInputs(['mass1', 'mw1'], calcBlock1);
  bindInputs(['vol2_in', 'mw2'], calcBlock2);

  // Static Head
  bindInputs(['p_kgcm2', 'sg1'], calcStaticHeadFromPressure);
  bindInputs(['head_m_in', 'sg2'], calcPressureFromStaticHead);

  // Vertical / Horizontal Vessel
  bindInputs(['v_head_type'], calcVerticalVessel, 'change');
  bindInputs(['v_diam', 'v_level'], calcVerticalVessel);

  bindInputs(['h_head_type'], calcHorizontalVessel, 'change');
  bindInputs(['h_diam', 'h_length', 'h_level'], calcHorizontalVessel);

  // Gas Density
  bindInputs(['gd_mw', 'gd_z', 'gd_p', 'gd_t'], calcGasDensity);

  // Line Sizing
  bindInputs(
    ['ls1_mass', 'ls1_rho', 'ls1_mu', 'ls1_p', 'ls1_t', 'ls1_d', 'ls1_eps'],
    calcLS1
  );
  bindInputs(
    ['ls2_q', 'ls2_rho', 'ls2_mu', 'ls2_p', 'ls2_t', 'ls2_d', 'ls2_eps'],
    calcLS2
  );
  bindInputs(
    [
      'ls3_mass',
      'ls3_mw',
      'ls3_mu',
      'ls3_z',
      'ls3_k',
      'ls3_side',
      'ls3_pknown',
      'ls3_t',
      'ls3_d',
      'ls3_eps',
      'ls3_leq'
    ],
    calcLS3
  );

  // Combined Gas
  bindInputs(
    ['atm_v2', 'p1_v2', 't1_v2', 'v1_v2', 'p2_v2', 't2_v2'],
    calcCombinedV2
  );
  bindInputs(
    ['atm_p2', 'p1_p2', 't1_p2', 'v1_p2', 't2_p2', 'v2_p2'],
    calcCombinedP2
  );
  bindInputs(
    ['atm_t2', 'p1_t2', 't1_t2', 'v1_t2', 'p2_t2', 'v2_t2'],
    calcCombinedT2
  );

  // Copy buttons
  const copyButtons = document.querySelectorAll('.btn-copy');
  copyButtons.forEach(btn => {
    btn.addEventListener('click', function () {
      const block = this.getAttribute('data-block');
      let tsv = null;
      switch (block) {
        case '1':
          tsv = getBlock1Tsv();
          break;
        case '2':
          tsv = getBlock2Tsv();
          break;
        case '3':
          tsv = getBlock3Tsv();
          break;
        case '4':
          tsv = getBlock4Tsv();
          break;
        case '5':
          tsv = getBlock5Tsv();
          break;
        case '6':
          tsv = getBlock6Tsv();
          break;
        case '7':
          tsv = getBlock7Tsv();
          break;
        case '8':
          tsv = getBlock8Tsv();
          break;
        case '9':
          tsv = getBlock9Tsv();
          break;
        case '10':
          tsv = getBlock10Tsv();
          break;
        case 'ls1':
          tsv = getLS1Tsv();
          break;
        case 'ls2':
          tsv = getLS2Tsv();
          break;
        case 'ls3':
          tsv = getLS3Tsv();
          break;
        default:
          break;
      }
      if (tsv) {
        copyTsvToClipboard(tsv, 'Current block table has been copied.');
      }
    });
  });

  // Initial calculations
  if (document.getElementById('mass1') && document.getElementById('mw1')) calcBlock1();
  if (document.getElementById('vol2_in') && document.getElementById('mw2')) calcBlock2();
  if (document.getElementById('p_kgcm2') && document.getElementById('sg1')) calcStaticHeadFromPressure();
  if (document.getElementById('head_m_in') && document.getElementById('sg2')) calcPressureFromStaticHead();
  if (document.getElementById('v_head_type')) calcVerticalVessel();
  if (document.getElementById('h_head_type')) calcHorizontalVessel();
  if (document.getElementById('ls1_mass')) calcLS1();
  if (document.getElementById('ls2_q')) calcLS2();
  if (document.getElementById('ls3_mass')) calcLS3();
  if (document.getElementById('atm_v2')) calcCombinedV2();
  if (document.getElementById('atm_p2')) calcCombinedP2();
  if (document.getElementById('atm_t2')) calcCombinedT2();
  if (document.getElementById('gd_mw')) calcGasDensity();
});
