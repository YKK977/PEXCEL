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
  const level_mm = parseFloat(document.getElementById('v_level').value);
  const volSpan  = document.getElementById('v_volume');

  if (!isFinite(diam_mm) || !isFinite(level_mm) || diam_mm <= 0 || level_mm <= 0) {
    volSpan.textContent = '-';
    return;
  }

  const D = diam_mm / 1000;      // m
  const Level = level_mm / 1000; // m
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
    return;
  }

  const B = diam_mm;
  const C = len_mm;
  const D = level_mm;

  // factor = (PI()/8+0.5^2*((ACOS(1-2*(D/B))-ACOS(0))-0.5*(SIN(2*ACOS(1-2*(D/B)))-SIN(2*ACOS(0)))))/(PI()/4)
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



// --- Init / events ---
window.addEventListener('DOMContentLoaded', function () {
  const modeSelect = document.getElementById('calc-mode');
  const sectionMassVol = document.getElementById('section-mass-volume');
  const sectionStaticHead = document.getElementById('section-static-head');
  const sectionVessel = document.getElementById('section-vessel-volume');

  const sectionCombined = document.getElementById('section-combined-gas');
  const sectionGasDensity = document.getElementById('section-gas-density');

  function applyMode(mode) {
    sectionMassVol.style.display   = (mode === 'mass-volume')  ? ''    : 'none';
    sectionStaticHead.style.display = (mode === 'static-head')  ? ''    : 'none';
    sectionVessel.style.display    = (mode === 'vessel-volume') ? ''    : 'none';
    if (sectionCombined) {
      sectionCombined.style.display = (mode === 'combined-gas') ? '' : 'none';
    }
    if (sectionGasDensity) {
      sectionGasDensity.style.display = (mode === 'gas-density') ? '' : 'none';
    }
  }

  if (modeSelect) {
    applyMode(modeSelect.value);
    modeSelect.addEventListener('change', function () {
      applyMode(this.value);
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
  const v_level = document.getElementById('v_level');
  if (v_head_type) v_head_type.addEventListener('change', calcVerticalVessel);
  if (v_diam) v_diam.addEventListener('input', calcVerticalVessel);
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
      if (tsv) copyTsvToClipboard(tsv, 'Current block table has been copied.');
    });
  });


  // initial one-time calculations
  if (mass1 && mw1) calcBlock1();
  if (vol2_in && mw2) calcBlock2();
  if (p_kgcm2 && sg1) calcStaticHeadFromPressure();
  if (head_m_in && sg2) calcPressureFromStaticHead();
  if (v_head_type && v_diam && v_level) calcVerticalVessel();
  if (h_head_type && h_diam && h_length && h_level) calcHorizontalVessel();

  // Combined Gas initial
  if (document.getElementById('atm_v2')) calcCombinedV2();
  if (document.getElementById('atm_p2')) calcCombinedP2();
  if (document.getElementById('atm_t2')) calcCombinedT2();

  // Gas Density initial
  if (document.getElementById('gd_mw')) calcGasDensity();
});
