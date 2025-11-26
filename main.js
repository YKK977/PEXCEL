// main.js - 모든 계산 + 엑셀복사 + 모드 전환까지 한 번에 정리본

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

// --- Vessel: Vertical (주신 Excel 로직 그대로) ---
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

// --- Vessel: Horizontal (주신 Excel 수식을 JS로 구현) ---
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

// --- 클립보드 공통 ---
function copyTsvToClipboard(tsv, label) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(tsv)
      .then(() => {
        if (label) alert(label + ' Excel에서 A1 셀을 선택하고 붙여넣기 하세요.');
      })
      .catch(() => alert('클립보드 복사에 실패했습니다.'));
  } else {
    alert('브라우저가 클립보드 API를 지원하지 않습니다.');
  }
}

// --- TSV 생성들 ---
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
  const formula = '=IF(A2="2:1E",((2/3)/2*PI()*(B2/2/1000)^3+PI()*((B2/1000)^2)/4*(D2/1000)),IF(A2="HS",(PI()/12*(B2/1000)^3)+PI()*((B2/1000)^2)/4*(D2/1000),PI()*((B2/1000)^2)/4*(D2/1000)))';
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
  const formula = '=IF(A2="2:1E",(((PI()/4/2*B2*D2^2*(1-2*D2/(3*B2)))/10^9)*2+((PI()/8+0.5^2*((ACOS(1-2*(D2/B2))-ACOS(0))-0.5*(SIN(2*ACOS(1-2*(D2/B2)))-SIN(2*ACOS(0)))))/(PI()/4))*(PI()/4*B2*B2*C2/10^9)),(2*PI()/24*(D2^2)*(3*B2-2*D2)*2)/10^9+((PI()/8+0.5^2*((ACOS(1-2*(D2/B2))-ACOS(0))-0.5*(SIN(2*ACOS(1-2*(D2/B2)))-SIN(2*ACOS(0)))))/(PI()/4))*(PI()/4*B2*B2*C2/10^9))';
  const row = [
    type,
    diam,
    len,
    level,
    formula
  ].join('\t');
  return header + '\n' + row;
}

// --- 초기화 / 이벤트 ---
window.addEventListener('DOMContentLoaded', function () {
  const modeSelect = document.getElementById('calc-mode');
  const sectionMassVol = document.getElementById('section-mass-volume');
  const sectionStaticHead = document.getElementById('section-static-head');
  const sectionVessel = document.getElementById('section-vessel-volume');

  function applyMode(mode) {
    sectionMassVol.style.display   = (mode === 'mass-volume')  ? ''    : 'none';
    sectionStaticHead.style.display = (mode === 'static-head')  ? ''    : 'none';
    sectionVessel.style.display    = (mode === 'vessel-volume') ? ''    : 'none';
  }

  if (modeSelect) {
    applyMode(modeSelect.value);
    modeSelect.addEventListener('change', function () {
      applyMode(this.value);
    });
  }

  // 입력 이벤트 연결
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

  // 복사 버튼
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
      if (tsv) copyTsvToClipboard(tsv, '현재 블록 표가 복사되었습니다.');
    });
  });

  // Ctrl+C 가로채기
  document.addEventListener('copy', function (e) {
    const active = document.activeElement;
    if (!active) return;
    const blocks = [
      {el: document.getElementById('block-1'), fn: getBlock1Tsv},
      {el: document.getElementById('block-2'), fn: getBlock2Tsv},
      {el: document.getElementById('block-3'), fn: getBlock3Tsv},
      {el: document.getElementById('block-4'), fn: getBlock4Tsv},
      {el: document.getElementById('block-5'), fn: getBlock5Tsv},
      {el: document.getElementById('block-6'), fn: getBlock6Tsv}
    ];
    let tsv = null;
    blocks.forEach(info => {
      if (info.el && info.el.contains(active)) {
        tsv = info.fn();
      }
    });
    if (tsv) {
      e.preventDefault();
      e.clipboardData.setData('text/plain', tsv);
      alert('현재 블록 내용(수식 포함)이 Excel용으로 복사되었습니다. A1 셀에 붙여넣으세요.');
    }
  });

  // 초기 계산 한 번씩
  if (mass1 && mw1) calcBlock1();
  if (vol2_in && mw2) calcBlock2();
  if (p_kgcm2 && sg1) calcStaticHeadFromPressure();
  if (head_m_in && sg2) calcPressureFromStaticHead();
  if (v_head_type && v_diam && v_level) calcVerticalVessel();
  if (h_head_type && h_diam && h_length && h_level) calcHorizontalVessel();
});
