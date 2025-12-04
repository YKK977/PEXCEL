// main.js - All calculations + Excel copy + mode switching in one file

// --- Mass / Volume ---
function calcBlock1() {
  const mass = parseFloat(document.getElementById('mass1').value);
  const mw   = parseFloat(document.getElementById('mw1').value);
  const volSpan = document.getElementById('vol1');

  if (!isFinite(mass) || !isFinite(mw) || mw === 0) {
    volSpan.textContent = '-';
    return;
  }
  const vol = mass * 22.414 / mw;
  volSpan.textContent = vol.toFixed(3);
}

function calcBlock2() {
  const vol = parseFloat(document.getElementById('vol2_in').value);
  const mw  = parseFloat(document.getElementById('mw2').value);
  const massSpan = document.getElementById('mass2');

  if (!isFinite(vol) || !isFinite(mw) || mw === 0) {
    massSpan.textContent = '-';
    return;
  }
  const mass = vol * mw / 22.414;
  massSpan.textContent = mass.toFixed(3);
}

// --- Static Head ---
function calcStaticHeadFromPressure() {
  const p = parseFloat(document.getElementById('p_kgcm2').value);
  const sg = parseFloat(document.getElementById('sg1').value);
  const headSpan = document.getElementById('head_m');

  if (!isFinite(p) || !isFinite(sg) || sg === 0) {
    headSpan.textContent = '-';
    return;
  }
  const head = p * 10 / sg;
  headSpan.textContent = head.toFixed(3);
}

function calcPressureFromStaticHead() {
  const head = parseFloat(document.getElementById('head_m_in').value);
  const sg = parseFloat(document.getElementById('sg2').value);
  const pSpan = document.getElementById('p_kgcm2_out');

  if (!isFinite(head) || !isFinite(sg)) {
    pSpan.textContent = '-';
    return;
  }
  const p = head * sg / 10;
  pSpan.textContent = p.toFixed(3);
}

// --- Vessel: Vertical (using the provided Excel logic as-is) ---

function calcVerticalVessel() {
  const headType = document.getElementById('v_head_type').value;
  const diam_mm  = parseFloat(document.getElementById('v_diam').value);
  const len_mm   = parseFloat(document.getElementById('v_length').value);
  const level_mm = parseFloat(document.getElementById('v_level').value);
  const volSpan  = document.getElementById('v_volume');

  if (!isFinite(diam_mm) || !isFinite(len_mm) || !isFinite(level_mm) ||
      diam_mm <= 0 || len_mm <= 0 || level_mm < 0) {
    volSpan.textContent = '-';
    updateVerticalVesselVisual(diam_mm || 0, len_mm || 0, 0, headType);
    return;
  }

  const D = diam_mm / 1000;      // m
  const L = len_mm / 1000;       // m (straight length)
  const Level = Math.max(0, Math.min(level_mm, len_mm)) / 1000; // m, clamp 0~Length
  const R = D / 2;

  let volume = 0;
  if (headType === "2:1E") {
    volume = ((2/3)/2) * Math.PI * Math.pow(R, 3) + Math.PI * (Math.pow(D, 2)/4) * Level;
  } else if (headType === "HS") {
    volume = (Math.PI/12) * Math.pow(D, 3) + Math.PI * (Math.pow(D, 2)/4) * Level;
  } else { // FLAT
    volume = Math.PI * (Math.pow(D, 2)/4) * Level;
  }

  volSpan.textContent = volume.toFixed(3);

  // SVG 업데이트 (길이에 비례해서 채움)
  updateVerticalVesselVisual(diam_mm, len_mm, Level * 1000, headType);
}

function updateVerticalHeadShape(headType) {
  const head2 = document.getElementById('vv-head-2_1E');
  const headHS = document.getElementById('vv-head-HS');
  if (!head2 || !headHS) return;

  if (headType === '2:1E') {
    head2.style.display = '';
    headHS.style.display = 'none';
  } else if (headType === 'HS') {
    head2.style.display = 'none';
    headHS.style.display = '';
  } else { // FLAT: head 없음, 박스만
    head2.style.display = 'none';
    headHS.style.display = 'none';
  }
}

function updateVerticalHeadLiquidShape(headType) {
  const liq2 = document.getElementById('vv-head-2_1E-liquid');
  const liqHS = document.getElementById('vv-head-HS-liquid');

  if (!liq2 || !liqHS) return;

  if (headType === '2:1E') {
    // 2:1E head일 때: 2:1E liquid만 표시
    liq2.style.display = '';
    liqHS.style.display = 'none';
  } else if (headType === 'HS') {
    // HS head일 때: HS liquid만 표시
    liq2.style.display = 'none';
    liqHS.style.display = '';
  } else {
    // FLAT일 때: head liquid는 둘 다 숨김
    liq2.style.display = 'none';
    liqHS.style.display = 'none';
  }
}

function updateVerticalVesselVisual(diam_mm, len_mm, level_mm, headType) {
  const liquidRect = document.getElementById('vessel-vertical-liquid');
  if (!liquidRect) return;

  updateVerticalHeadShape(headType);
  updateVerticalHeadLiquidShape(headType);

  const L = len_mm;                                 // mm (straight length)
  const Level = Math.max(0, Math.min(level_mm, L)); // 0~L
  const frac = L === 0 ? 0 : Level / L;             // 0~1

  // SVG 기준: 탱크 내부 높이 = 120 (y=40~160)
  const tankTop = 40;
  const tankBottom = 160;
  const tankHeight = tankBottom - tankTop;

  const liquidHeight = tankHeight * frac;
  const liquidY = tankBottom - liquidHeight;

  liquidRect.setAttribute('y', liquidY);
  liquidRect.setAttribute('height', liquidHeight);
}


// --- Vessel: Horizontal (implementing the provided Excel formulas in JS) ---
function calcHorizontalVessel() {
  const headType = document.getElementById('h_head_type').value; // "2:1E" / "HS"
  const diam_mm  = parseFloat(document.getElementById('h_diam').value);
  const len_mm   = parseFloat(document.getElementById('h_length').value);
  const level_mm = parseFloat(document.getElementById('h_level').value);
  const volSpan  = document.getElementById('h_volume');

  if (!isFinite(diam_mm) || !isFinite(len_mm) || !isFinite(level_mm) ||
      diam_mm <= 0 || len_mm <= 0 || level_mm < 0) {
    volSpan.textContent = '-';
    updateHorizontalVesselVisual(diam_mm || 0, 0, headType);
    return;
  }

  const B = diam_mm;
  const C = len_mm;
  const D = level_mm;

  const acos_term = Math.acos(1 - 2 * (D / B));
  const factor = (Math.PI/8 + Math.pow(0.5, 2) * (
                    (acos_term - Math.acos(0)) -
                    0.5 * (Math.sin(2 * acos_term) - Math.sin(2 * Math.acos(0)))
                  )) / (Math.PI/4);

  const cyl_full = (Math.PI/4 * B * B * C) / 1e9;

  let volume = 0;
  if (headType === "2:1E") {
    const head_liq = ((Math.PI/4/2 * B * D * D * (1 - 2*D/(3*B))) / 1e9) * 2;
    volume = head_liq + factor * cyl_full;
  } else { // HS
    const head_liq = (2*Math.PI/24 * (D*D) * (3*B - 2*D) * 2) / 1e9;
    volume = head_liq + factor * cyl_full;
  }

  volSpan.textContent = volume.toFixed(3);

  // SVG 업데이트
  updateHorizontalVesselVisual(diam_mm, level_mm, headType);
}

function updateHorizontalHeadShape(headType) {
  const pairs = [
    ['2:1E', 'hv-head-2_1E-left', 'hv-head-2_1E-right'],
    ['HS',   'hv-head-HS-left',   'hv-head-HS-right']
  ];

  pairs.forEach(([type, leftId, rightId]) => {
    const show = (type === headType);
    const left = document.getElementById(leftId);
    const right = document.getElementById(rightId);
    if (left) left.style.display = show ? '' : 'none';
    if (right) right.style.display = show ? '' : 'none';
  });

  const liq2 = document.getElementById('vessel-horizontal-liquid-2_1E');
  const liqHS = document.getElementById('vessel-horizontal-liquid-HS');
  if (liq2 && liqHS) {
    if (headType === '2:1E') {
      liq2.style.display = '';
      liqHS.style.display = 'none';
    } else if (headType === 'HS') {
      liq2.style.display = 'none';
      liqHS.style.display = '';
    } else { // FLAT: head 없음, body만
      liq2.style.display = 'none';
      liqHS.style.display = 'none';
    }
  }
}

function updateHorizontalVesselVisual(diam_mm, level_mm, headType) {
  const liq2 = document.getElementById('vessel-horizontal-liquid-2_1E');
  const liqHS = document.getElementById('vessel-horizontal-liquid-HS');
  if (!liq2 || !liqHS) return;

  updateHorizontalHeadShape(headType);

  const D = diam_mm;
  const Level = Math.max(0, Math.min(level_mm, D));   // 0~D
  const frac = D === 0 ? 0 : Level / D;               // 0~1

  // SVG 기준: 탱크 내부 높이 = 80 (y=40~120)
  const tankTop = 40;
  const tankBottom = 120;
  const tankHeight = tankBottom - tankTop;

  const liquidHeight = tankHeight * frac;
  const liquidY = tankBottom - liquidHeight;

  [liq2, liqHS].forEach(rect => {
    rect.setAttribute('y', liquidY);
    rect.setAttribute('height', liquidHeight);
  });
}


// --- Combined Gas Law ---
// helper
function cToK(tC) {
  return tC + 273.15;
}

// V2
function calcCombinedV2() {
  const atm = parseFloat(document.getElementById('atm_v2').value);
  const p1  = parseFloat(document.getElementById('p1_v2').value);
  const t1  = parseFloat(document.getElementById('t1_v2').value);
  const v1  = parseFloat(document.getElementById('v1_v2').value);
  const p2  = parseFloat(document.getElementById('p2_v2').value);
  const t2  = parseFloat(document.getElementById('t2_v2').value);
  const out = document.getElementById('v2_out');

  if (![atm,p1,t1,v1,p2,t2].every(Number.isFinite) || (p2+atm) === 0) {
    out.textContent = '-';
    return;
  }
  const v2 = ((p1 + atm) * v1 * cToK(t2)) / ((p2 + atm) * cToK(t1));
  out.textContent = v2.toFixed(3);
}

// P2
function calcCombinedP2() {
  const atm = parseFloat(document.getElementById('atm_p2').value);
  const p1  = parseFloat(document.getElementById('p1_p2').value);
  const t1  = parseFloat(document.getElementById('t1_p2').value);
  const v1  = parseFloat(document.getElementById('v1_p2').value);
  const t2  = parseFloat(document.getElementById('t2_p2').value);
  const v2  = parseFloat(document.getElementById('v2_p2').value);
  const out = document.getElementById('p2_out');

  if (![atm,p1,t1,v1,t2,v2].every(Number.isFinite) || v2 === 0) {
    out.textContent = '-';
    return;
  }
  const p2_abs = ((p1 + atm) * v1 * cToK(t2)) / (v2 * cToK(t1));
  const p2_g = p2_abs - atm;
  out.textContent = p2_g.toFixed(3);
}

// T2
function calcCombinedT2() {
  const atm = parseFloat(document.getElementById('atm_t2').value);
  const p1  = parseFloat(document.getElementById('p1_t2').value);
  const t1  = parseFloat(document.getElementById('t1_t2').value);
  const v1  = parseFloat(document.getElementById('v1_t2').value);
  const p2  = parseFloat(document.getElementById('p2_t2').value);
  const v2  = parseFloat(document.getElementById('v2_t2').value);
  const out = document.getElementById('t2_out');

  if (![atm,p1,t1,v1,p2,v2].every(Number.isFinite) || (p1+atm) === 0) {
    out.textContent = '-';
    return;
  }
  const t2K = ((p2 + atm) * v2 * cToK(t1)) / ((p1 + atm) * v1);
  const t2C = t2K - 273.15;
  out.textContent = t2C.toFixed(3);
}

// TSV generators for Combined Gas Law
function getBlock7Tsv() {
  const atm = document.getElementById('atm_v2').value || '';
  const p1  = document.getElementById('p1_v2').value || '';
  const t1  = document.getElementById('t1_v2').value || '';
  const v1  = document.getElementById('v1_v2').value || '';
  const p2  = document.getElementById('p2_v2').value || '';
  const t2  = document.getElementById('t2_v2').value || '';
  const header = [
    '1 ATM (kg/cm2g)',
    'P1 (kg/cm2g)',
    'T1 (degC)',
    'V1 (m3)',
    'P2 (kg/cm2g)',
    'T2 (degC)',
    'V2 (m3)'
  ].join('\t');
  const row = [
    atm,
    p1,
    t1,
    v1,
    p2,
    t2,
    '=((B2+A2)*D2*(F2+273.15))/((E2+A2)*(C2+273.15))'
  ].join('\t');
  return header + '\n' + row;
}

function getBlock8Tsv() {
  const atm = document.getElementById('atm_p2').value || '';
  const p1  = document.getElementById('p1_p2').value || '';
  const t1  = document.getElementById('t1_p2').value || '';
  const v1  = document.getElementById('v1_p2').value || '';
  const t2  = document.getElementById('t2_p2').value || '';
  const v2  = document.getElementById('v2_p2').value || '';
  const header = [
    '1 ATM (kg/cm2g)',
    'P1 (kg/cm2g)',
    'T1 (degC)',
    'V1 (m3)',
    'P2 (kg/cm2g)',
    'T2 (degC)',
    'V2 (m3)'
  ].join('\t');
  const row = [
    atm,
    p1,
    t1,
    v1,
    '=((B2+A2)*D2*(F2+273.15))/(G2*(C2+273.15))-A2',
    t2,
    v2
  ].join('\t');
  return header + '\n' + row;
}

function getBlock9Tsv() {
  const atm = document.getElementById('atm_t2').value || '';
  const p1  = document.getElementById('p1_t2').value || '';
  const t1  = document.getElementById('t1_t2').value || '';
  const v1  = document.getElementById('v1_t2').value || '';
  const p2  = document.getElementById('p2_t2').value || '';
  const v2  = document.getElementById('v2_t2').value || '';
  const header = [
    '1 ATM (kg/cm2g)',
    'P1 (kg/cm2g)',
    'T1 (degC)',
    'V1 (m3)',
    'P2 (kg/cm2g)',
    'T2 (degC)',
    'V2 (m3)'
  ].join('\t');
  const row = [
    atm,
    p1,
    t1,
    v1,
    p2,
    '=((E2+A2)*G2*(C2+273.15))/((B2+A2)*D2)-273.15',
    v2
  ].join('\t');
  return header + '\n' + row;
}


// --- Gas Density ---
function calcGasDensity() {
  const mw = parseFloat(document.getElementById('gd_mw')?.value);
  const z  = parseFloat(document.getElementById('gd_z')?.value);
  const p  = parseFloat(document.getElementById('gd_p')?.value);
  const tC = parseFloat(document.getElementById('gd_t')?.value);
  const rhoSpan = document.getElementById('gd_rho');

  if (!rhoSpan || ![mw,z,p,tC].every(Number.isFinite) || z === 0) {
    if (rhoSpan) rhoSpan.textContent = '-';
    return;
  }
  // Excel: =((C3+1.03323)*98.0665*A3)/(B3*8.314*(D3+273.15))
  const rho = ((p + 1.03323) * 98.0665 * mw) / (z * 8.314 * (tC + 273.15));
  rhoSpan.textContent = rho.toFixed(3);
}

function getBlock10Tsv() {
  const mw = document.getElementById('gd_mw')?.value || '';
  const z  = document.getElementById('gd_z')?.value || '';
  const p  = document.getElementById('gd_p')?.value || '';
  const tC = document.getElementById('gd_t')?.value || '';
  const header = [
    'MW (-)',
    'Comp. Factor Z (-)',
    'Pressure (kg/cm2g)',
    'Temperature (degC)',
    'Gas Density (kg/m3)'
  ].join('\t');
  const row = [
    mw,
    z,
    p,
    tC,
    '=((C2+1.03323)*98.0665*A2)/(B2*8.314*(D2+273.15))'
  ].join('\t');
  return header + '\n' + row;
}

// --- Clipboard common ---
function copyTsvToClipboard(tsv, label) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(tsv)
      .then(() => {
        if (label) alert(label + ' In Excel, select cell A1 and paste.');
      })
      .catch(() => alert('Failed to copy to clipboard.'));
  } else {
    alert('This browser does not support the Clipboard API.');
  }
}

// --- TSV generators ---
function getBlock1Tsv() {
  const mass = document.getElementById('mass1').value || '';
  const mw   = document.getElementById('mw1').value || '';
  const header = [
    'Mass Flowrate (kg/h)',
    'MW (kg/kmol)',
    'Volume Flowrate (Nm3/h)'
  ].join('\t');
  const row = [
    mass,
    mw,
    '=A2*22.414/B2'
  ].join('\t');
  return header + '\n' + row;
}

function getBlock2Tsv() {
  const vol = document.getElementById('vol2_in').value || '';
  const mw  = document.getElementById('mw2').value || '';
  const header = [
    'Volume Flowrate (Nm3/h)',
    'MW (kg/kmol)',
    'Mass Flowrate (kg/h)'
  ].join('\t');
  const row = [
    vol,
    mw,
    '=A2*B2/22.414'
  ].join('\t');
  return header + '\n' + row;
}

function getBlock3Tsv() {
  const p  = document.getElementById('p_kgcm2').value || '';
  const sg = document.getElementById('sg1').value || '';
  const header = [
    'D. Pressure (kg/cm2)',
    'S.G (-)',
    'Head (m)'
  ].join('\t');
  const row = [
    p,
    sg,
    '=A2*10/B2'
  ].join('\t');
  return header + '\n' + row;
}

function getBlock4Tsv() {
  const head = document.getElementById('head_m_in').value || '';
  const sg   = document.getElementById('sg2').value || '';
  const header = [
    'Head (m)',
    'S.G (-)',
    'D. Pressure (kg/cm2)'
  ].join('\t');
  const row = [
    head,
    sg,
    '=A2*B2/10'
  ].join('\t');
  return header + '\n' + row;
}

function getBlock5Tsv() {
  const type  = document.getElementById('v_head_type').value || '';
  const diam  = document.getElementById('v_diam').value || '';
  const len   = document.getElementById('v_length').value || '';
  const level = document.getElementById('v_level').value || '';

  const header = [
    'Head Type',
    'Diameter (mm)',
    'Length (mm)',
    'Liquid Level (mm)',
    'Volume (m3)'
  ].join('\t');

  let formula = '';
  if (type === "2:1E") {
    formula = '=((2/3)/2*PI()*(B2/2/1000)^3+PI()*((B2/1000)^2)/4*(D2/1000))';
  } else if (type === "HS") {
    formula = '=(PI()/12*(B2/1000)^3)+PI()*((B2/1000)^2)/4*(D2/1000)';
  } else {
    // FLAT or others
    formula = '=PI()*((B2/1000)^2)/4*(D2/1000)';
  }

  const row = [
    type,
    diam,
    len,
    level,
    formula
  ].join('\t');

  return header + '\n' + row;
}


function getBlock6Tsv() {
  const type  = document.getElementById('h_head_type').value || '';
  const diam  = document.getElementById('h_diam').value || '';
  const len   = document.getElementById('h_length').value || '';
  const level = document.getElementById('h_level').value || '';

  const header = [
    'Type',
    'Diameter (mm)',
    'Length (mm)',
    'Liquid Level (mm)',
    'Volume (m3)'
  ].join('\t');

  let formula = '';
  if (type === "2:1E") {
    formula = '=((PI()/4/2*B2*D2^2*(1-2*D2/(3*B2)))/10^9)*2+((PI()/8+0.5^2*((ACOS(1-2*(D2/B2))-ACOS(0))-0.5*(SIN(2*ACOS(1-2*(D2/B2)))-SIN(2*ACOS(0)))))/(PI()/4))*(PI()/4*B2*B2*C2/10^9)';
  } else {
    // HS
    formula = '=(2*PI()/24*(D2^2)*(3*B2-2*D2)*2)/10^9+((PI()/8+0.5^2*((ACOS(1-2*(D2/B2))-ACOS(0))-0.5*(SIN(2*ACOS(1-2*(D2/B2)))-SIN(2*ACOS(0)))))/(PI()/4))*(PI()/4*B2*B2*C2/10^9)';
  }

  const row = [
    type,
    diam,
    len,
    level,
    formula
  ].join('\t');

  return header + '\n' + row;
}





// --- Line Sizing common helpers ---
function ls_incompressible_dp(massFlow_kg_h, rho, mu_cP, d_inch, rough_ft, length_m) {
  const mdot = massFlow_kg_h / 3600.0;
  const D = d_inch * 0.0254;
  const eps = rough_ft * 0.3048;
  const area = Math.PI * D * D / 4.0;
  const v = mdot / (rho * area);
  const mu = mu_cP * 1e-3; // Pa·s
  const Re = (rho * v * D) / mu;
  if (!isFinite(Re) || Re <= 0 || !isFinite(v)) {
    return { v: NaN, dp100: NaN };
  }
  const rel = eps / D;
  const f = 0.25 / Math.pow(Math.log10(rel / 3.7 + 5.74 / Math.pow(Re, 0.9)), 2);
  const dP_per_m_Pa = f * rho * v * v / (2 * D);
  const dP_100m_kgcm2 = dP_per_m_Pa * 100.0 / 98066.5;
  return { v, dp100: dP_100m_kgcm2 };
}

// Liquid from Mass Flow
function calcLS1() {
  const m = parseFloat(document.getElementById('ls1_mass')?.value);
  const rho = parseFloat(document.getElementById('ls1_rho')?.value);
  const mu = parseFloat(document.getElementById('ls1_mu')?.value);
  const d = parseFloat(document.getElementById('ls1_d')?.value);
  const eps = parseFloat(document.getElementById('ls1_eps')?.value);
  const qSpan = document.getElementById('ls1_q');
  const dpSpan = document.getElementById('ls1_dp100');
  const vSpan = document.getElementById('ls1_v');

  if (![m,rho,mu,d,eps].every(Number.isFinite) || rho <= 0 || d <= 0) {
    if (qSpan) qSpan.textContent = '-';
    if (dpSpan) dpSpan.textContent = '-';
    if (vSpan) vSpan.textContent = '-';
    return;
  }

  const q = m / rho; // m3/h
  const res = ls_incompressible_dp(m, rho, mu, d, eps, 100.0);
  if (qSpan) qSpan.textContent = q.toFixed(3);
  if (dpSpan) dpSpan.textContent = isFinite(res.dp100) ? res.dp100.toFixed(3) : '-';
  if (vSpan) vSpan.textContent = isFinite(res.v) ? res.v.toFixed(3) : '-';
}

// Liquid from Volumetric Flow
function calcLS2() {
  const q = parseFloat(document.getElementById('ls2_q')?.value);
  const rho = parseFloat(document.getElementById('ls2_rho')?.value);
  const mu = parseFloat(document.getElementById('ls2_mu')?.value);
  const d = parseFloat(document.getElementById('ls2_d')?.value);
  const eps = parseFloat(document.getElementById('ls2_eps')?.value);
  const mSpan = document.getElementById('ls2_mass');
  const dpSpan = document.getElementById('ls2_dp100');
  const vSpan = document.getElementById('ls2_v');

  if (![q,rho,mu,d,eps].every(Number.isFinite) || rho <= 0 || d <= 0) {
    if (mSpan) mSpan.textContent = '-';
    if (dpSpan) dpSpan.textContent = '-';
    if (vSpan) vSpan.textContent = '-';
    return;
  }

  const m = q * rho; // kg/h
  const res = ls_incompressible_dp(m, rho, mu, d, eps, 100.0);
  if (mSpan) mSpan.textContent = m.toFixed(3);
  if (dpSpan) dpSpan.textContent = isFinite(res.dp100) ? res.dp100.toFixed(3) : '-';
  if (vSpan) vSpan.textContent = isFinite(res.v) ? res.v.toFixed(3) : '-';
}

// Vapor Darcy (incompressible approximation)
function ls_gas_density(mw, z, p_g_kgcm2, tC) {
  const p_abs = (p_g_kgcm2 + 1.03323) * 98066.5; // Pa
  const T = tC + 273.15;
  const Rspec = 8.314 * 1000.0 / mw; // J/kg/K
  return p_abs * mw / (z * 8.314 * 1000.0 * T); // same as p/(z*Rspec*T)
}

function ls_sonic_speed(k, mw, tC) {
  const T = tC + 273.15;
  const Rspec = 8.314 * 1000.0 / mw;
  return Math.sqrt(k * Rspec * T);
}

function calcLS3() {
  const m = parseFloat(document.getElementById('ls3_mass')?.value);
  const mw = parseFloat(document.getElementById('ls3_mw')?.value);
  const mu = parseFloat(document.getElementById('ls3_mu')?.value);
  const z = parseFloat(document.getElementById('ls3_z')?.value);
  const k = parseFloat(document.getElementById('ls3_k')?.value);
  const side = parseFloat(document.getElementById('ls3_side')?.value);
  const pKnown = parseFloat(document.getElementById('ls3_pknown')?.value);
  const tC = parseFloat(document.getElementById('ls3_t')?.value);
  const d = parseFloat(document.getElementById('ls3_d')?.value);
  const eps = parseFloat(document.getElementById('ls3_eps')?.value);
  const L = parseFloat(document.getElementById('ls3_leq')?.value);

  const pinSpan = document.getElementById('ls3_pin');
  const poutSpan = document.getElementById('ls3_pout');
  const dpSpan = document.getElementById('ls3_dp');
  const dp100Span = document.getElementById('ls3_dp100');
  const rhoInSpan = document.getElementById('ls3_rho_in');
  const rhoOutSpan = document.getElementById('ls3_rho_out');
  const vInSpan = document.getElementById('ls3_v_in');
  const vOutSpan = document.getElementById('ls3_v_out');
  const machInSpan = document.getElementById('ls3_mach_in');
  const machOutSpan = document.getElementById('ls3_mach_out');
  const aSpan = document.getElementById('ls3_a');
  const veSpan = document.getElementById('ls3_ve');
  const vratioSpan = document.getElementById('ls3_vratio');

  if (![m,mw,mu,z,k,side,pKnown,tC,d,eps,L].every(Number.isFinite) || mw<=0 || z<=0 || d<=0 || L<=0) {
    [pinSpan,poutSpan,dpSpan,dp100Span,rhoInSpan,rhoOutSpan,vInSpan,vOutSpan,machInSpan,machOutSpan,aSpan,veSpan,vratioSpan].forEach(el=>{
      if (el) el.textContent='-';
    });
    return;
  }

  const rhoKnown = ls_gas_density(mw, z, pKnown, tC);
  const res = ls_incompressible_dp(m, rhoKnown, mu, d, eps, L);
  const dp_total = res.dp100 * (L/100.0);

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
  const area = Math.PI * D * D / 4.0;

  const mdot = m / 3600.0;
  const v_in = mdot / (rho_in * area);
  const v_out = mdot / (rho_out * area);
  const a = ls_sonic_speed(k, mw, tC);
  const mach_in = v_in / a;
  const mach_out = v_out / a;

  // erosional velocity (API) using C=100 ft/s
  let ve = NaN;
  let vratio = NaN;
  if (rho_out > 0) {
    const C_ft_s = 100.0;
    const C_ms = C_ft_s * 0.3048;
    ve = C_ms / Math.sqrt(rho_out);
    vratio = v_out / ve;
  }

  if (pinSpan) pinSpan.textContent = p_in.toFixed(3);
  if (poutSpan) poutSpan.textContent = p_out.toFixed(3);
  if (dpSpan) dpSpan.textContent = dp_total.toFixed(3);
  if (dp100Span) dp100Span.textContent = res.dp100.toFixed(3);
  if (rhoInSpan) rhoInSpan.textContent = rho_in.toFixed(2);
  if (rhoOutSpan) rhoOutSpan.textContent = rho_out.toFixed(2);
  if (vInSpan) vInSpan.textContent = v_in.toFixed(2);
  if (vOutSpan) vOutSpan.textContent = v_out.toFixed(2);
  if (machInSpan) machInSpan.textContent = mach_in.toFixed(3);
  if (machOutSpan) machOutSpan.textContent = mach_out.toFixed(3);
  if (aSpan) aSpan.textContent = a.toFixed(1);
  if (veSpan) veSpan.textContent = Number.isFinite(ve) ? ve.toFixed(2) : '-';
  if (vratioSpan) vratioSpan.textContent = Number.isFinite(vratio) ? vratio.toFixed(2) : '-';
}


// --- Residence Time ---// --- Residence Time ---
function calcResidence() {
  const V = parseFloat(document.getElementById('res_vol')?.value);
  const Q = parseFloat(document.getElementById('res_q')?.value);
  const thSpan = document.getElementById('res_time_h');
  const tmSpan = document.getElementById('res_time_min');
  if (![V,Q].every(Number.isFinite) || Q <= 0) {
    if (thSpan) thSpan.textContent = '-';
    if (tmSpan) tmSpan.textContent = '-';
    return;
  }
  const t_h = V / Q;
  const t_min = t_h * 60.0;
  if (thSpan) thSpan.textContent = t_h.toFixed(3);
  if (tmSpan) tmSpan.textContent = t_min.toFixed(1);
}

// --- NPSHa ---
function calcNPSH() {
  const ps_g = parseFloat(document.getElementById('npsh_ps')?.value);
  const hs   = parseFloat(document.getElementById('npsh_hs')?.value);
  const hf   = parseFloat(document.getElementById('npsh_hf')?.value);
  const sg   = parseFloat(document.getElementById('npsh_sg')?.value);
  const pv_a = parseFloat(document.getElementById('npsh_pv')?.value);
  const out  = document.getElementById('npsh_avail');
  if (![ps_g,hs,hf,sg,pv_a].every(Number.isFinite) || sg <= 0) {
    if (out) out.textContent = '-';
    return;
  }
  const ATM = 1.03323; // kg/cm2 abs
  const p_surface_abs = ps_g + ATM;
  const head_surface = p_surface_abs * 10.0 / sg;
  const head_vapor   = pv_a * 10.0 / sg;
  const npsh = head_surface + hs - hf - head_vapor;
  if (out) out.textContent = npsh.toFixed(3);
}

// --- Pump Power ---
function calcPumpPower() {
  const Q = parseFloat(document.getElementById('hp_q')?.value);
  const dp = parseFloat(document.getElementById('hp_dp')?.value);
  const eff = parseFloat(document.getElementById('hp_eff')?.value);
  const phSpan = document.getElementById('hp_ph');
  const psSpan = document.getElementById('hp_ps');
  if (![Q,dp,eff].every(Number.isFinite) || Q < 0 || dp < 0 || eff <= 0) {
    if (phSpan) phSpan.textContent = '-';
    if (psSpan) psSpan.textContent = '-';
    return;
  }
  const Q_m3s = Q / 3600.0;
  const dp_Pa = dp * 98066.5;
  const P_h_kw = dp_Pa * Q_m3s / 1000.0;
  const P_s_kw = P_h_kw / (eff / 100.0);
  if (phSpan) phSpan.textContent = P_h_kw.toFixed(3);
  if (psSpan) psSpan.textContent = P_s_kw.toFixed(3);
}

// --- Heat Duty / Area ---
function calcHeatDuty() {
  const m  = parseFloat(document.getElementById('hd_m')?.value);
  const cp = parseFloat(document.getElementById('hd_cp')?.value);
  const tin = parseFloat(document.getElementById('hd_tin')?.value);
  const tout = parseFloat(document.getElementById('hd_tout')?.value);
  const U = parseFloat(document.getElementById('hd_u')?.value);
  const lmtd = parseFloat(document.getElementById('hd_lmtd')?.value);
  const qSpan = document.getElementById('hd_q');
  const aSpan = document.getElementById('hd_a');
  if (![m,cp,tin,tout].every(Number.isFinite)) {
    if (qSpan) qSpan.textContent = '-';
    if (aSpan) aSpan.textContent = '-';
    return;
  }
  const Q_kw = m * cp * (tout - tin) / 3600.0;
  let A = NaN;
  if (Number.isFinite(U) && Number.isFinite(lmtd) && U > 0 && lmtd !== 0) {
    const Q_W = Q_kw * 1000.0;
    A = Q_W / (U * lmtd);
  }
  if (qSpan) qSpan.textContent = Q_kw.toFixed(3);
  if (aSpan) aSpan.textContent = Number.isFinite(A) ? A.toFixed(3) : '-';
}

// --- Steam Duty ---
function calcSteamDuty() {
  const Q = parseFloat(document.getElementById('sd_q')?.value);
  const lambda = parseFloat(document.getElementById('sd_lambda')?.value);
  const mSteamSpan = document.getElementById('sd_msteam');
  const mCondSpan = document.getElementById('sd_mcond');
  if (![Q,lambda].every(Number.isFinite) || lambda <= 0) {
    if (mSteamSpan) mSteamSpan.textContent = '-';
    if (mCondSpan) mCondSpan.textContent = '-';
    return;
  }
  const m_kgph = Q * 3600.0 / lambda;
  if (mSteamSpan) mSteamSpan.textContent = m_kgph.toFixed(3);
  if (mCondSpan) mCondSpan.textContent = m_kgph.toFixed(3);
}

// --- Control Valve Cv Inversion ---
function calcCvInvert() {
  const q_in = parseFloat(document.getElementById('cv_q_in')?.value);
  const dp_in = parseFloat(document.getElementById('cv_dp_in')?.value);
  const sg = parseFloat(document.getElementById('cv_sg')?.value);
  const cv = parseFloat(document.getElementById('cv_cv')?.value);
  const qSpan = document.getElementById('cv_q_calc');
  const dpSpan = document.getElementById('cv_dp_calc');
  if (![sg,cv].every(Number.isFinite) || sg <= 0 || cv <= 0) {
    if (qSpan) qSpan.textContent = '-';
    if (dpSpan) dpSpan.textContent = '-';
    return;
  }
  const GPM_PER_M3H = 4.402867;
  const KGCM2_PER_PSI = 0.07030697;
  // From Cv & ΔP -> Q
  let q_from = NaN;
  if (Number.isFinite(dp_in) && dp_in > 0) {
    const dp_psi = dp_in / KGCM2_PER_PSI;
    const q_gpm = cv * Math.sqrt(dp_psi / sg);
    q_from = q_gpm / GPM_PER_M3H;
  }
  // From Cv & Q -> ΔP
  let dp_from = NaN;
  if (Number.isFinite(q_in) && q_in > 0) {
    const q_gpm2 = q_in * GPM_PER_M3H;
    const dp_psi2 = (q_gpm2 / cv) ** 2 * sg;
    dp_from = dp_psi2 * KGCM2_PER_PSI;
  }
  if (qSpan) qSpan.textContent = Number.isFinite(q_from) ? q_from.toFixed(3) : '-';
  if (dpSpan) dpSpan.textContent = Number.isFinite(dp_from) ? dp_from.toFixed(3) : '-';
}

// --- Restriction Orifice ---
function calcOrifice() {
  const dp = parseFloat(document.getElementById('ro_dp')?.value);
  const d_mm = parseFloat(document.getElementById('ro_d')?.value);
  const C = parseFloat(document.getElementById('ro_c')?.value);
  const rho = parseFloat(document.getElementById('ro_rho')?.value);
  const qSpan = document.getElementById('ro_q');
  const vSpan = document.getElementById('ro_v');
  if (![dp,d_mm,C,rho].every(Number.isFinite) || dp <= 0 || d_mm <= 0 || C <= 0 || rho <= 0) {
    if (qSpan) qSpan.textContent = '-';
    if (vSpan) vSpan.textContent = '-';
    return;
  }
  const dp_Pa = dp * 98066.5;
  const d = d_mm / 1000.0;
  const A = Math.PI * d * d / 4.0;
  const Q_m3s = C * A * Math.sqrt(2.0 * dp_Pa / rho);
  const Q_m3h = Q_m3s * 3600.0;
  const v = Q_m3s / A;
  if (qSpan) qSpan.textContent = Q_m3h.toFixed(3);
  if (vSpan) vSpan.textContent = v.toFixed(3);
}

// --- Mixing Temperature ---
function calcMixing() {
  const m1 = parseFloat(document.getElementById('mix_m1')?.value);
  const cp1 = parseFloat(document.getElementById('mix_cp1')?.value);
  const t1 = parseFloat(document.getElementById('mix_t1')?.value);
  const m2 = parseFloat(document.getElementById('mix_m2')?.value);
  const cp2 = parseFloat(document.getElementById('mix_cp2')?.value);
  const t2 = parseFloat(document.getElementById('mix_t2')?.value);
  const out = document.getElementById('mix_tmix');
  if (![m1,cp1,t1,m2,cp2,t2].every(Number.isFinite) || m1<0 || m2<0 || cp1<=0 || cp2<=0) {
    if (out) out.textContent = '-';
    return;
  }
  const num = m1*cp1*t1 + m2*cp2*t2;
  const den = m1*cp1 + m2*cp2;
  if (den === 0) {
    if (out) out.textContent = '-';
    return;
  }
  const tmix = num / den;
  if (out) out.textContent = tmix.toFixed(2);
}
// TSV generators for Line Sizing

// TSV generators for Line Sizing

function getResidenceTsv() {
  const header = ['Holdup Volume (m3)','Flowrate (m3/h)','Residence Time (h)','Residence Time (min)'].join('\t');
  const row = [
    document.getElementById('res_vol')?.value || '',
    document.getElementById('res_q')?.value || '',
    '=A2/B2',
    '=(A2/B2)*60'
  ].join('\t');
  return header + '\n' + row;
}

function getNPSHTsv() {
  const header = ['Surface P (kg/cm2g)','Static Head (m)','Line Loss (m)','SG','Vapor P (kg/cm2 abs)','NPSHa (m)'].join('\t');
  const row = [
    document.getElementById('npsh_ps')?.value || '',
    document.getElementById('npsh_hs')?.value || '',
    document.getElementById('npsh_hf')?.value || '',
    document.getElementById('npsh_sg')?.value || '',
    document.getElementById('npsh_pv')?.value || '',
    '=((A2+1.03323)*10/D2)+B2-C2-(E2*10/D2)'
  ].join('\t');
  return header + '\n' + row;
}

function getPumpPowerTsv() {
  const header = ['Flow (m3/h)','DP (kg/cm2)','P_hyd (kW)','Eff (%)','P_shaft (kW)'].join('\t');
  const row = [
    document.getElementById('hp_q')?.value || '',
    document.getElementById('hp_dp')?.value || '',
    '=B2*98066.5*(A2/3600)/1000',
    document.getElementById('hp_eff')?.value || '',
    '=(B2*98066.5*(A2/3600)/1000)/(D2/100)'
  ].join('\t');
  return header + '\n' + row;
}

function getHeatDutyTsv() {
  const header = ['Mass Flow (kg/h)','Cp (kJ/kgK)','Tin (C)','Tout (C)','U (W/m2K)','LMTD (K)','Q (kW)','Area (m2)'].join('\t');
  const row = [
    document.getElementById('hd_m')?.value || '',
    document.getElementById('hd_cp')?.value || '',
    document.getElementById('hd_tin')?.value || '',
    document.getElementById('hd_tout')?.value || '',
    document.getElementById('hd_u')?.value || '',
    document.getElementById('hd_lmtd')?.value || '',
    '=A2*B2*(D2-C2)/3600',
    '=IF(AND(E2>0,F2<>0),A2*B2*(D2-C2)/3600*1000/(E2*F2),"")'
  ].join('\t');
  return header + '\n' + row;
}

function getSteamDutyTsv() {
  const header = ['Heat Duty (kW)','Latent Heat (kJ/kg)','Steam Flow (kg/h)','Condensate Flow (kg/h)'].join('\t');
  const row = [
    document.getElementById('sd_q')?.value || '',
    document.getElementById('sd_lambda')?.value || '',
    '=A2*3600/B2',
    '=A2*3600/B2'
  ].join('\t');
  return header + '\n' + row;
}

function getCvInvertTsv() {
  const header = ['Q_in (m3/h)','DP_in (kg/cm2)','SG','Cv','Q_from_Cv_DP (m3/h)','DP_from_Cv_Q (kg/cm2)'].join('\t');
  const row = [
    document.getElementById('cv_q_in')?.value || '',
    document.getElementById('cv_dp_in')?.value || '',
    document.getElementById('cv_sg')?.value || '',
    document.getElementById('cv_cv')?.value || '',
    '', // will be calculated in JS
    ''  // will be calculated in JS
  ].join('\t');
  return header + '\n' + row;
}

function getOrificeTsv() {
  const header = ['DP (kg/cm2)','d (mm)','C','rho (kg/m3)','Q (m3/h)','v (m/s)'].join('\t');
  const row = [
    document.getElementById('ro_dp')?.value || '',
    document.getElementById('ro_d')?.value || '',
    document.getElementById('ro_c')?.value || '',
    document.getElementById('ro_rho')?.value || '',
    '', // JS calc
    ''
  ].join('\t');
  return header + '\n' + row;
}

function getMixingTsv() {
  const header = ['m1 (kg/h)','Cp1 (kJ/kgK)','T1 (C)','m2 (kg/h)','Cp2 (kJ/kgK)','T2 (C)','Tmix (C)'].join('\t');
  const row = [
    document.getElementById('mix_m1')?.value || '',
    document.getElementById('mix_cp1')?.value || '',
    document.getElementById('mix_t1')?.value || '',
    document.getElementById('mix_m2')?.value || '',
    document.getElementById('mix_cp2')?.value || '',
    document.getElementById('mix_t2')?.value || '',
    '=(A2*B2*C2+D2*E2*F2)/(A2*B2+D2*E2)'
  ].join('\t');
  return header + '\n' + row;
}

function getLS1Tsv() {
  const header = [
    'Mass Flow (kg/h)','Density (kg/m3)','Viscosity (cP)','Press (kg/cm2g)','Temp (degC)',
    'Pipe ID (inch)','Rough (ft)','Vol Flow (m3/h)','DP/100m (kg/cm2)','Velocity (m/s)'
  ].join('\t');
  const row = [
    document.getElementById('ls1_mass')?.value || '',
    document.getElementById('ls1_rho')?.value || '',
    document.getElementById('ls1_mu')?.value || '',
    document.getElementById('ls1_p')?.value || '',
    document.getElementById('ls1_t')?.value || '',
    document.getElementById('ls1_d')?.value || '',
    document.getElementById('ls1_eps')?.value || '',
    '=A2/B2',
    '=0.25/(LOG10((G2*0.3048)/(3.7*(F2*0.0254))+5.74/(((B2*J2*(F2*0.0254))/(C2/1000))^0.9))^2)*B2*J2^2/(2*(F2*0.0254))*100/98066.5',
    '=4*(A2/(3600*B2))/(PI()*(F2*0.0254)^2)'
  ].join('\t');
  return header + '\n' + row;
}

function getLS2Tsv() {
  const header = [
    'Vol Flow (m3/h)','Density (kg/m3)','Viscosity (cP)','Press (kg/cm2g)','Temp (degC)',
    'Pipe ID (inch)','Rough (ft)','Mass Flow (kg/h)','DP/100m (kg/cm2)','Velocity (m/s)'
  ].join('\t');
  const row = [
    document.getElementById('ls2_q')?.value || '',
    document.getElementById('ls2_rho')?.value || '',
    document.getElementById('ls2_mu')?.value || '',
    document.getElementById('ls2_p')?.value || '',
    document.getElementById('ls2_t')?.value || '',
    document.getElementById('ls2_d')?.value || '',
    document.getElementById('ls2_eps')?.value || '',
    '=A2*B2',
    '=0.25/(LOG10((G2*0.3048)/(3.7*(F2*0.0254))+5.74/(((B2*J2*(F2*0.0254))/(C2/1000))^0.9))^2)*B2*J2^2/(2*(F2*0.0254))*100/98066.5',
    '=4*(A2/3600)/(PI()*(F2*0.0254)^2)'
  ].join('\t');
  return header + '\n' + row;
}

// LS3 TSV: input + output in one row for Excel
function getLS3Tsv() {
  const header = [
    'Mass Flow (kg/h)','MW (kg/kmol)','Viscosity (cP)','Z-factor','k (Cp/Cv)','Known Side',
    'Known P (kg/cm2g)','Temp (degC)','Pipe ID (inch)','Rough (ft)','EQ Length (m)',
    'P_in (kg/cm2g)','P_out (kg/cm2g)','Total dP (kg/cm2)','DP/100m (kg/cm2)',
    'rho_in (kg/m3)','rho_out (kg/m3)','v_in (m/s)','v_out (m/s)','Mach_in','Mach_out','Sonic Vel (m/s)'
  ].join('\t');
  const row = [
    document.getElementById('ls3_mass')?.value || '',
    document.getElementById('ls3_mw')?.value || '',
    document.getElementById('ls3_mu')?.value || '',
    document.getElementById('ls3_z')?.value || '',
    document.getElementById('ls3_k')?.value || '',
    document.getElementById('ls3_side')?.value || '',
    document.getElementById('ls3_pknown')?.value || '',
    document.getElementById('ls3_t')?.value || '',
    document.getElementById('ls3_d')?.value || '',
    document.getElementById('ls3_eps')?.value || '',
    document.getElementById('ls3_leq')?.value || '',
    document.getElementById('ls3_pin')?.textContent || '',
    document.getElementById('ls3_pout')?.textContent || '',
    document.getElementById('ls3_dp')?.textContent || '',
    document.getElementById('ls3_dp100')?.textContent || '',
    document.getElementById('ls3_rho_in')?.textContent || '',
    document.getElementById('ls3_rho_out')?.textContent || '',
    document.getElementById('ls3_v_in')?.textContent || '',
    document.getElementById('ls3_v_out')?.textContent || '',
    document.getElementById('ls3_mach_in')?.textContent || '',
    document.getElementById('ls3_mach_out')?.textContent || '',
    document.getElementById('ls3_a')?.textContent || ''
  ].join('\t');
  return header + '\n' + row;
}

// --- Init / events ---
window.addEventListener('DOMContentLoaded', function () {
  const categorySelect = document.getElementById('category-select');

  const sections = {
    'section-mass-volume': document.getElementById('section-mass-volume'),
    'section-static-head': document.getElementById('section-static-head'),
    'section-vessel-volume': document.getElementById('section-vessel-volume'),
    'section-combined-gas': document.getElementById('section-combined-gas'),
    'section-gas-density': document.getElementById('section-gas-density'),
    'section-line-sizing': document.getElementById('section-line-sizing'),
    'section-residence': document.getElementById('section-residence'),
    'section-npsh': document.getElementById('section-npsh'),
    'section-pump-power': document.getElementById('section-pump-power'),
    'section-heat-duty': document.getElementById('section-heat-duty'),
    'section-steam-duty': document.getElementById('section-steam-duty'),
    'section-cv-invert': document.getElementById('section-cv-invert'),
    'section-orifice': document.getElementById('section-orifice'),
    'section-mixing-temp': document.getElementById('section-mixing-temp')
  };

  const indexes = {
    'vessel': document.getElementById('index-vessel'),
    'pump': document.getElementById('index-pump'),
    'heat-exchanger': document.getElementById('index-heat-exchanger'),
    'instrumentation': document.getElementById('index-instrumentation'),
    'line-sizing': document.getElementById('index-line-sizing'),
    'physical-properties': document.getElementById('index-physical-properties')
  };

  const CATEGORY_SECTIONS = {
    'vessel': ['section-vessel-volume', 'section-residence'],
    'pump': ['section-static-head', 'section-npsh', 'section-pump-power'],
    'heat-exchanger': ['section-heat-duty', 'section-steam-duty'],
    'instrumentation': ['section-cv-invert', 'section-orifice'],
    'line-sizing': ['section-line-sizing'],
    'physical-properties': ['section-mass-volume', 'section-combined-gas', 'section-gas-density', 'section-mixing-temp']
  };

  function applyCategory(category) {
    // hide all sections
    Object.keys(sections).forEach(function (id) {
      const sec = sections[id];
      if (sec) sec.style.display = 'none';
    });
    // hide all indexes
    Object.keys(indexes).forEach(function (key) {
      const nav = indexes[key];
      if (nav) nav.style.display = 'none';
    });
    const list = CATEGORY_SECTIONS[category] || [];
    list.forEach(function (id) {
      const sec = sections[id];
      if (sec) sec.style.display = '';
    });
    const nav = indexes[category];
    if (nav) nav.style.display = '';
  }

  if (categorySelect) {
    applyCategory(categorySelect.value);
    categorySelect.addEventListener('change', function () {
      applyCategory(this.value);
    });
  }

// hook up input events
  const mass1 = document.getElementById('mass1');
  const mw1 = document.getElementById('mw1');
  const vol2_in = document.getElementById('vol2_in');
  const mw2 = document.getElementById('mw2');
  if (mass1) mass1.addEventListener('input', calcBlock1);
  if (mw1) mw1.addEventListener('input', calcBlock1);
  if (vol2_in) vol2_in.addEventListener('input', calcBlock2);
  if (mw2) mw2.addEventListener('input', calcBlock2);

  const p_kgcm2 = document.getElementById('p_kgcm2');
  const sg1 = document.getElementById('sg1');
  const head_m_in = document.getElementById('head_m_in');
  const sg2 = document.getElementById('sg2');
  if (p_kgcm2) p_kgcm2.addEventListener('input', calcStaticHeadFromPressure);
  if (sg1) sg1.addEventListener('input', calcStaticHeadFromPressure);
  if (head_m_in) head_m_in.addEventListener('input', calcPressureFromStaticHead);
  if (sg2) sg2.addEventListener('input', calcPressureFromStaticHead);

  const v_head_type = document.getElementById('v_head_type');
  const v_diam = document.getElementById('v_diam');
  const v_length = document.getElementById('v_length');
  const v_level = document.getElementById('v_level');
  if (v_head_type) v_head_type.addEventListener('change', calcVerticalVessel);
  if (v_diam) v_diam.addEventListener('input', calcVerticalVessel);
  if (v_length) v_length.addEventListener('input', calcVerticalVessel);
  if (v_level) v_level.addEventListener('input', calcVerticalVessel);

  const h_head_type = document.getElementById('h_head_type');
  const h_diam = document.getElementById('h_diam');
  const h_length = document.getElementById('h_length');
  const h_level = document.getElementById('h_level');
  if (h_head_type) h_head_type.addEventListener('change', calcHorizontalVessel);
  if (h_diam) h_diam.addEventListener('input', calcHorizontalVessel);
  if (h_length) h_length.addEventListener('input', calcHorizontalVessel);
  if (h_level) h_level.addEventListener('input', calcHorizontalVessel);

  // Gas Density listeners
  const gd_mw = document.getElementById('gd_mw');
  const gd_z  = document.getElementById('gd_z');
  const gd_p  = document.getElementById('gd_p');
  const gd_t  = document.getElementById('gd_t');
  if (gd_mw) gd_mw.addEventListener('input', calcGasDensity);
  if (gd_z)  gd_z.addEventListener('input', calcGasDensity);
  if (gd_p)  gd_p.addEventListener('input', calcGasDensity);
  if (gd_t)  gd_t.addEventListener('input', calcGasDensity);

  // Line sizing listeners
  const ls1_ids = ['ls1_mass','ls1_rho','ls1_mu','ls1_p','ls1_t','ls1_d','ls1_eps'];
  ls1_ids.forEach(id=>{ const el=document.getElementById(id); if(el) el.addEventListener('input', calcLS1); });
  const ls2_ids = ['ls2_q','ls2_rho','ls2_mu','ls2_p','ls2_t','ls2_d','ls2_eps'];
  ls2_ids.forEach(id=>{ const el=document.getElementById(id); if(el) el.addEventListener('input', calcLS2); });
  const ls3_ids = ['ls3_mass','ls3_mw','ls3_mu','ls3_z','ls3_k','ls3_side','ls3_pknown','ls3_t','ls3_d','ls3_eps','ls3_leq'];
  ls3_ids.forEach(id=>{ const el=document.getElementById(id); if(el) el.addEventListener('input', calcLS3); });

  // Residence time listeners
  const res_vol = document.getElementById('res_vol');
  const res_q = document.getElementById('res_q');
  if (res_vol) res_vol.addEventListener('input', calcResidence);
  if (res_q) res_q.addEventListener('input', calcResidence);

  // NPSH listeners
  ['npsh_ps','npsh_hs','npsh_hf','npsh_sg','npsh_pv'].forEach(id=>{
    const el = document.getElementById(id);
    if (el) el.addEventListener('input', calcNPSH);
  });

  // Pump power listeners
  ['hp_q','hp_dp','hp_eff'].forEach(id=>{
    const el = document.getElementById(id);
    if (el) el.addEventListener('input', calcPumpPower);
  });

  // Heat duty listeners
  ['hd_m','hd_cp','hd_tin','hd_tout','hd_u','hd_lmtd'].forEach(id=>{
    const el = document.getElementById(id);
    if (el) el.addEventListener('input', calcHeatDuty);
  });

  // Steam duty listeners
  ['sd_q','sd_lambda'].forEach(id=>{
    const el = document.getElementById(id);
    if (el) el.addEventListener('input', calcSteamDuty);
  });

  // Cv inversion listeners
  ['cv_q_in','cv_dp_in','cv_sg','cv_cv'].forEach(id=>{
    const el = document.getElementById(id);
    if (el) el.addEventListener('input', calcCvInvert);
  });

  // Orifice listeners
  ['ro_dp','ro_d','ro_c','ro_rho'].forEach(id=>{
    const el = document.getElementById(id);
    if (el) el.addEventListener('input', calcOrifice);
  });

  // Mixing temperature listeners
  ['mix_m1','mix_cp1','mix_t1','mix_m2','mix_cp2','mix_t2'].forEach(id=>{
    const el = document.getElementById(id);
    if (el) el.addEventListener('input', calcMixing);
  });


  // Combined Gas Law listeners
  const atm_v2 = document.getElementById('atm_v2');
  const p1_v2 = document.getElementById('p1_v2');
  const t1_v2 = document.getElementById('t1_v2');
  const v1_v2 = document.getElementById('v1_v2');
  const p2_v2 = document.getElementById('p2_v2');
  const t2_v2 = document.getElementById('t2_v2');
  if (atm_v2) atm_v2.addEventListener('input', calcCombinedV2);
  if (p1_v2) p1_v2.addEventListener('input', calcCombinedV2);
  if (t1_v2) t1_v2.addEventListener('input', calcCombinedV2);
  if (v1_v2) v1_v2.addEventListener('input', calcCombinedV2);
  if (p2_v2) p2_v2.addEventListener('input', calcCombinedV2);
  if (t2_v2) t2_v2.addEventListener('input', calcCombinedV2);

  const atm_p2 = document.getElementById('atm_p2');
  const p1_p2 = document.getElementById('p1_p2');
  const t1_p2 = document.getElementById('t1_p2');
  const v1_p2 = document.getElementById('v1_p2');
  const t2_p2 = document.getElementById('t2_p2');
  const v2_p2 = document.getElementById('v2_p2');
  if (atm_p2) atm_p2.addEventListener('input', calcCombinedP2);
  if (p1_p2) p1_p2.addEventListener('input', calcCombinedP2);
  if (t1_p2) t1_p2.addEventListener('input', calcCombinedP2);
  if (v1_p2) v1_p2.addEventListener('input', calcCombinedP2);
  if (t2_p2) t2_p2.addEventListener('input', calcCombinedP2);
  if (v2_p2) v2_p2.addEventListener('input', calcCombinedP2);

  const atm_t2 = document.getElementById('atm_t2');
  const p1_t2 = document.getElementById('p1_t2');
  const t1_t2 = document.getElementById('t1_t2');
  const v1_t2 = document.getElementById('v1_t2');
  const p2_t2 = document.getElementById('p2_t2');
  const v2_t2 = document.getElementById('v2_t2');
  if (atm_t2) atm_t2.addEventListener('input', calcCombinedT2);
  if (p1_t2) p1_t2.addEventListener('input', calcCombinedT2);
  if (t1_t2) t1_t2.addEventListener('input', calcCombinedT2);
  if (v1_t2) v1_t2.addEventListener('input', calcCombinedT2);
  if (p2_t2) p2_t2.addEventListener('input', calcCombinedT2);
  if (v2_t2) v2_t2.addEventListener('input', calcCombinedT2);

  // copy buttons
  const copyButtons = document.querySelectorAll('.btn-copy');
  copyButtons.forEach(btn => {
    btn.addEventListener('click', function () {
      const block = this.getAttribute('data-block');
      let tsv = null;
      if (block === '1') tsv = getBlock1Tsv();
      else if (block === '2') tsv = getBlock2Tsv();
      else if (block === '3') tsv = getBlock3Tsv();
      else if (block === '4') tsv = getBlock4Tsv();
      else if (block === '5') tsv = getBlock5Tsv();
      else if (block === '6') tsv = getBlock6Tsv();
      else if (block === '7') tsv = getBlock7Tsv();
      else if (block === '8') tsv = getBlock8Tsv();
      else if (block === '9') tsv = getBlock9Tsv();
      else if (block === '10') tsv = getBlock10Tsv();
      else if (block === 'residence') tsv = getResidenceTsv();
      else if (block === 'npsh') tsv = getNPSHTsv();
      else if (block === 'pump-power') tsv = getPumpPowerTsv();
      else if (block === 'heat-duty') tsv = getHeatDutyTsv();
      else if (block === 'steam-duty') tsv = getSteamDutyTsv();
      else if (block === 'cv-invert') tsv = getCvInvertTsv();
      else if (block === 'orifice') tsv = getOrificeTsv();
      else if (block === 'mixing-temp') tsv = getMixingTsv();
      else if (block === 'ls1') tsv = getLS1Tsv();
      else if (block === 'ls2') tsv = getLS2Tsv();
      else if (block === 'ls3') tsv = getLS3Tsv();
      if (tsv) copyTsvToClipboard(tsv, 'Current block table has been copied.');
    });
  });


  // initial one-time calculations
  if (mass1 && mw1) calcBlock1();
  if (vol2_in && mw2) calcBlock2();
  if (p_kgcm2 && sg1) calcStaticHeadFromPressure();
  if (head_m_in && sg2) calcPressureFromStaticHead();
  if (v_head_type && v_diam && document.getElementById('v_length') && v_level) calcVerticalVessel();
  if (h_head_type && h_diam && h_length && h_level) calcHorizontalVessel();
  if (document.getElementById('ls1_mass')) calcLS1();
  if (document.getElementById('ls2_q')) calcLS2();
  if (document.getElementById('ls3_mass')) calcLS3();

  if (document.getElementById('res_vol')) calcResidence();
  if (document.getElementById('npsh_ps')) calcNPSH();
  if (document.getElementById('hp_q')) calcPumpPower();
  if (document.getElementById('hd_m')) calcHeatDuty();
  if (document.getElementById('sd_q')) calcSteamDuty();
  if (document.getElementById('cv_q_in')) calcCvInvert();
  if (document.getElementById('ro_dp')) calcOrifice();
  if (document.getElementById('mix_m1')) calcMixing();

  // Combined Gas initial
  if (document.getElementById('atm_v2')) calcCombinedV2();
  if (document.getElementById('atm_p2')) calcCombinedP2();
  if (document.getElementById('atm_t2')) calcCombinedT2();

  // Gas Density initial
  if (document.getElementById('gd_mw')) calcGasDensity();
});
