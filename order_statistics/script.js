"use strict";

const nInput = document.getElementById("sampleSize");
const bInput = document.getElementById("simulationCount");
const seedInput = document.getElementById("seedInput");

const statusEl = document.getElementById("status");
const latestEl = document.getElementById("latest");

const runButton = document.getElementById("runButton");
const onceButton = document.getElementById("onceButton");
const resetButton = document.getElementById("resetButton");

const TYPES = ["min", "median", "max", "range"];
const DISTRIBUTIONS = ["u", "n"];
const SQRT_2PI = Math.sqrt(2 * Math.PI);

const results = {
  u: { min: [], median: [], max: [], range: [] },
  n: { min: [], median: [], max: [], range: [] }
};

let currentN = null;
let currentSeed = null;
let random = null;
let busy = false;

/* ==================== 带种子的随机数 ==================== */

function mulberry32(seed) {
  let state = seed >>> 0;

  return function () {
    state = (state + 0x6D2B79F5) | 0;

    let t = Math.imul(state ^ (state >>> 15), 1 | state);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;

    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function makeRandom(seed) {
  const uniform = mulberry32(seed);
  let spareNormal = null;

  function normal() {
    if (spareNormal !== null) {
      const value = spareNormal;
      spareNormal = null;
      return value;
    }

    const u1 = 1 - uniform();
    const u2 = uniform();
    const radius = Math.sqrt(-2 * Math.log(u1));
    const angle = 2 * Math.PI * u2;

    spareNormal = radius * Math.sin(angle);
    return radius * Math.cos(angle);
  }

  // Marsaglia–Tsang Gamma(shape, 1) 抽样。
  function gamma(shape) {
    const d = shape - 1 / 3;
    const c = 1 / Math.sqrt(9 * d);

    while (true) {
      const z = normal();
      const v0 = 1 + c * z;

      if (v0 <= 0) continue;

      const v = v0 ** 3;
      const u = uniform();

      if (u < 1 - 0.0331 * z ** 4) return d * v;

      if (
        Math.log(u) <
        0.5 * z * z + d * (1 - v + Math.log(v))
      ) {
        return d * v;
      }
    }
  }

  function beta(a, b) {
    const x = gamma(a);
    const y = gamma(b);
    return x / (x + y);
  }

  return { uniform, beta };
}

/* ==================== 正态分布函数 ==================== */

function normalPDF(x) {
  return Math.exp(-x * x / 2) / SQRT_2PI;
}

// 标准正态 CDF 数值近似。
function normalCDF(x) {
  const a = Math.abs(x);
  const t = 1 / (1 + 0.2316419 * a);

  const polynomial =
    t * (
      0.319381530 +
      t * (
        -0.356563782 +
        t * (
          1.781477937 +
          t * (-1.821255978 + t * 1.330274429)
        )
      )
    );

  const upper = normalPDF(a) * polynomial;
  return x >= 0 ? 1 - upper : upper;
}

// Acklam 标准正态分位数近似。
function normalQuantile(p) {
  p = Math.max(1e-15, Math.min(1 - 1e-15, p));

  const a = [
    -3.969683028665376e+01,
     2.209460984245205e+02,
    -2.759285104469687e+02,
     1.383577518672690e+02,
    -3.066479806614716e+01,
     2.506628277459239e+00
  ];

  const b = [
    -5.447609879822406e+01,
     1.615858368580409e+02,
    -1.556989798598866e+02,
     6.680131188771972e+01,
    -1.328068155288572e+01
  ];

  const c = [
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
     4.374664141464968e+00,
     2.938163982698783e+00
  ];

  const d = [
     7.784695709041462e-03,
     3.224671290700398e-01,
     2.445134137142996e+00,
     3.754408661907416e+00
  ];

  if (p < 0.02425 || p > 0.97575) {
    const q = Math.sqrt(
      -2 * Math.log(p < 0.02425 ? p : 1 - p)
    );

    const numerator =
      (((((c[0] * q + c[1]) * q + c[2]) * q +
          c[3]) * q + c[4]) * q + c[5]);

    const denominator =
      ((((d[0] * q + d[1]) * q + d[2]) * q +
          d[3]) * q + 1);

    const value = numerator / denominator;
    return p < 0.02425 ? value : -value;
  }

  const q = p - 0.5;
  const r = q * q;

  return (
    (((((a[0] * r + a[1]) * r + a[2]) * r +
        a[3]) * r + a[4]) * r + a[5]) * q
  ) / (
    (((((b[0] * r + b[1]) * r + b[2]) * r +
        b[3]) * r + b[4]) * r + 1)
  );
}

/* ==================== 参数管理 ==================== */

function readInteger(input, min, max, name) {
  const raw = input.value.trim();
  const value = Number(raw);

  if (
    !/^\d+$/.test(raw) ||
    !Number.isSafeInteger(value) ||
    value < min ||
    value > max
  ) {
    alert(`${name}必须是 ${min} 到 ${max} 之间的整数。`);
    input.focus();
    return null;
  }

  return value;
}

function readParameters() {
  const n = readInteger(nInput, 3, 10000, "样本量 n");
  if (n === null) return null;

  const seed = readInteger(
    seedInput, 0, 4294967295, "随机种子"
  );
  if (seed === null) return null;

  return { n, seed };
}

function setControlsDisabled(value) {
  for (const element of [
    nInput, bInput, seedInput,
    runButton, onceButton, resetButton
  ]) {
    element.disabled = value;
  }
}

function resetWith({ n, seed }) {
  currentN = n;
  currentSeed = seed;
  random = makeRandom(seed);

  for (const dist of DISTRIBUTIONS) {
    for (const type of TYPES) {
      results[dist][type].length = 0;
    }
  }

  latestEl.textContent = "尚未运行实验。";
  render();
}

function resetFromInputs() {
  if (busy) return;

  const parameters = readParameters();
  if (parameters !== null) resetWith(parameters);
}

/* ==================== 两栏配对模拟 ==================== */

/*
 * n <= 30：实际生成 n 个均匀数并排序。
 * 将同一组排序结果通过 Φ⁻¹ 变换，得到正态样本的排序结果。
 */
function simulateBySorting(n) {
  const sorted = Array.from(
    { length: n },
    () => random.uniform()
  ).sort((a, b) => a - b);

  const k = Math.floor((n + 1) / 2);

  return {
    uMin: sorted[0],
    uMedian: sorted[k - 1],
    uMax: sorted[n - 1],
    sorted
  };
}

/*
 * n > 30：模拟三个均匀次序统计量的精确联合分布。
 * 再经 Φ⁻¹ 变换得到对应的正态次序统计量。
 */
function simulateFast(n) {
  const k = Math.floor((n + 1) / 2);

  const uMin =
    -Math.expm1(Math.log1p(-random.uniform()) / n);

  const middleFraction = random.beta(
    k - 1,
    n - k + 1
  );

  const uMedian =
    uMin + (1 - uMin) * middleFraction;

  const upperFraction =
    Math.pow(random.uniform(), 1 / (n - k));

  const uMax =
    uMedian + (1 - uMedian) * upperFraction;

  return { uMin, uMedian, uMax, sorted: null };
}

function appendResult(sample) {
  const { uMin, uMedian, uMax } = sample;

  const nMin = normalQuantile(uMin);
  const nMedian = normalQuantile(uMedian);
  const nMax = normalQuantile(uMax);

  results.u.min.push(uMin);
  results.u.median.push(uMedian);
  results.u.max.push(uMax);
  results.u.range.push(uMax - uMin);

  results.n.min.push(nMin);
  results.n.median.push(nMedian);
  results.n.max.push(nMax);
  results.n.range.push(nMax - nMin);

  return { nMin, nMedian, nMax };
}

function run(count) {
  if (busy) return;

  const parameters = readParameters();
  if (parameters === null) return;

  if (
    parameters.n !== currentN ||
    parameters.seed !== currentSeed
  ) {
    resetWith(parameters);
  }

  busy = true;
  setControlsDisabled(true);

  let completed = 0;
  let lastSample = null;
  let lastNormal = null;

  function processBatch() {
    const take = Math.min(500, count - completed);

    for (let i = 0; i < take; i++) {
      lastSample = currentN <= 30
        ? simulateBySorting(currentN)
        : simulateFast(currentN);

      lastNormal = appendResult(lastSample);
    }

    completed += take;

    if (completed < count) {
      statusEl.textContent =
        `正在模拟：本次已完成 ${completed.toLocaleString()} / ` +
        `${count.toLocaleString()} 轮……`;

      setTimeout(processBatch, 0);
      return;
    }

    if (lastSample.sorted !== null) {
      latestEl.textContent =
        "最近一轮的均匀样本（已排序）：" +
        lastSample.sorted.map(x => x.toFixed(4)).join("，") +
        "；正态栏由这组样本经 Φ⁻¹ 变换得到。";
    } else {
      latestEl.textContent =
        "最近一轮——均匀分布：" +
        `最小值 ${lastSample.uMin.toFixed(5)}，` +
        `中间值 ${lastSample.uMedian.toFixed(5)}，` +
        `最大值 ${lastSample.uMax.toFixed(5)}；` +
        "标准正态：" +
        `最小值 ${lastNormal.nMin.toFixed(5)}，` +
        `中间值 ${lastNormal.nMedian.toFixed(5)}，` +
        `最大值 ${lastNormal.nMax.toFixed(5)}。`;
    }

    busy = false;
    setControlsDisabled(false);
    render();
  }

  setTimeout(processBatch, 0);
}

function runCustom() {
  const count = readInteger(
    bInput, 1, 20000, "模拟次数 B"
  );

  if (count !== null) run(count);
}

/* ==================== 理论密度 ==================== */

const logFactorial = [0];

for (let i = 1; i <= 10000; i++) {
  logFactorial[i] =
    logFactorial[i - 1] + Math.log(i);
}

function normalRangeDensity(r, n) {
  if (r <= 0) return 0;

  // f_R(r) = n(n−1)∫φ(t)φ(t+r)
  //          [Φ(t+r)−Φ(t)]^(n−2) dt
  // 用复合 Simpson 法进行数值积分。
  const left = -9;
  const right = 9;
  const segments = 120;
  const step = (right - left) / segments;

  function integrand(t) {
    const probability = Math.max(
      0,
      normalCDF(t + r) - normalCDF(t)
    );

    return (
      normalPDF(t) *
      normalPDF(t + r) *
      Math.pow(probability, n - 2)
    );
  }

  let sum = integrand(left) + integrand(right);

  for (let i = 1; i < segments; i++) {
    sum += (i % 2 === 0 ? 2 : 4) *
      integrand(left + i * step);
  }

  return n * (n - 1) * sum * step / 3;
}

function theoreticalDensity(dist, type, x, n) {
  const k = Math.floor((n + 1) / 2);

  if (dist === "u") {
    if (x < 0 || x > 1) return 0;

    if (type === "min") {
      return n * Math.pow(1 - x, n - 1);
    }

    if (type === "max") {
      return n * Math.pow(x, n - 1);
    }

    if (type === "range") {
      return n * (n - 1) *
        Math.pow(x, n - 2) * (1 - x);
    }
  }

  if (dist === "n" && type === "range") {
    return normalRangeDensity(x, n);
  }

  const F = dist === "u" ? x : normalCDF(x);
  const f = dist === "u" ? 1 : normalPDF(x);

  if (type === "min") {
    return n * Math.pow(1 - F, n - 1) * f;
  }

  if (type === "max") {
    return n * Math.pow(F, n - 1) * f;
  }

  if (F <= 0 || F >= 1) return 0;

  const logCoefficient =
    logFactorial[n] -
    logFactorial[k - 1] -
    logFactorial[n - k];

  return Math.exp(
    logCoefficient +
    (k - 1) * Math.log(F) +
    (n - k) * Math.log1p(-F)
  ) * f;
}

/* ==================== 绘图 ==================== */

function chartDomain(dist, type, n) {
  if (dist === "u") {
    if (n <= 30) return [0, 1];

    if (type === "min") {
      return [0, Math.min(1, 8 / n)];
    }

    if (type === "max" || type === "range") {
      return [Math.max(0, 1 - 10 / n), 1];
    }

    const k = Math.floor((n + 1) / 2);
    const center = k / (n + 1);
    const halfWidth = Math.min(0.5, 4 / Math.sqrt(n));

    return [
      Math.max(0, center - halfWidth),
      Math.min(1, center + halfWidth)
    ];
  }

  if (n <= 30) {
    return type === "range" ? [0, 8] : [-5, 5];
  }

  const q = normalQuantile(n / (n + 1));
  const extremeSD = 1 / (n * normalPDF(q));

  if (type === "min") {
    return [
      -q - 5 * extremeSD,
      -q + 5 * extremeSD
    ];
  }

  if (type === "max") {
    return [
      q - 5 * extremeSD,
      q + 5 * extremeSD
    ];
  }

  if (type === "range") {
    const center = 2 * q;
    const halfWidth = 6 * Math.SQRT2 * extremeSD;

    return [
      Math.max(0, center - halfWidth),
      center + halfWidth
    ];
  }

  const k = Math.floor((n + 1) / 2);
  const p = k / (n + 1);
  const center = normalQuantile(p);

  const variance =
    k * (n - k + 1) /
    ((n + 1) ** 2 * (n + 2));

  const sd = Math.sqrt(variance) / normalPDF(center);

  return [
    center - 5 * sd,
    center + 5 * sd
  ];
}

function formatTick(x, span) {
  if (Math.abs(x) < 1e-12) x = 0;

  if (span < 0.02) return x.toFixed(4);
  if (span < 0.2) return x.toFixed(3);
  return x.toFixed(2);
}

function drawChart(canvas, values, dist, type, n) {
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;

  const plot = {
    left: 62,
    right: 16,
    top: 16,
    bottom: 40
  };

  const plotWidth = width - plot.left - plot.right;
  const plotHeight = height - plot.top - plot.bottom;

  const [lo, hi] = chartDomain(dist, type, n);
  const span = hi - lo;
  const bins = 35;
  const binWidth = span / bins;

  const counts = Array(bins).fill(0);

  for (const value of values) {
    if (value < lo || value > hi) continue;

    const index = Math.min(
      bins - 1,
      Math.floor((value - lo) / binWidth)
    );

    counts[index]++;
  }

  // 即使放大横轴，分母仍是全部模拟次数。
  const histogram = counts.map(count =>
    values.length
      ? count / (values.length * binWidth)
      : 0
  );

  const theoretical = [];
  const pointCount = dist === "n" && type === "range"
    ? 110
    : 240;

  let maxY = Math.max(1, ...histogram);

  for (let i = 0; i <= pointCount; i++) {
    const x = lo + span * i / pointCount;
    const y = theoreticalDensity(dist, type, x, n);

    theoretical.push({ x, y });
    maxY = Math.max(maxY, y);
  }

  maxY *= 1.1;

  const xPixel = x =>
    plot.left + (x - lo) / span * plotWidth;

  const yPixel = y =>
    plot.top + plotHeight - y / maxY * plotHeight;

  ctx.clearRect(0, 0, width, height);
  ctx.font = "12px system-ui, sans-serif";

  // 网格及纵轴刻度。
  for (let i = 0; i <= 4; i++) {
    const yValue = maxY * i / 4;
    const py = yPixel(yValue);

    ctx.strokeStyle = "#e2e8f0";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(plot.left, py);
    ctx.lineTo(width - plot.right, py);
    ctx.stroke();

    ctx.fillStyle = "#64748b";
    ctx.fillText(yValue.toFixed(1), 4, py + 4);
  }

  // 经验密度直方图。
  ctx.fillStyle = "rgba(37, 99, 235, 0.48)";

  for (let i = 0; i < bins; i++) {
    const barWidth = plotWidth / bins;
    const barHeight =
      histogram[i] / maxY * plotHeight;

    ctx.fillRect(
      plot.left + i * barWidth + 1,
      plot.top + plotHeight - barHeight,
      Math.max(0, barWidth - 2),
      barHeight
    );
  }

  // 理论密度曲线。
  ctx.beginPath();
  ctx.strokeStyle = "#ea580c";
  ctx.lineWidth = 3;

  theoretical.forEach(({ x, y }, index) => {
    if (index === 0) {
      ctx.moveTo(xPixel(x), yPixel(y));
    } else {
      ctx.lineTo(xPixel(x), yPixel(y));
    }
  });

  ctx.stroke();

  // 横轴及刻度。
  ctx.strokeStyle = "#94a3b8";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(plot.left, plot.top + plotHeight);
  ctx.lineTo(width - plot.right, plot.top + plotHeight);
  ctx.stroke();

  ctx.fillStyle = "#64748b";
  ctx.font = "11px system-ui, sans-serif";

  for (let i = 0; i <= 4; i++) {
    const x = lo + span * i / 4;
    const label = formatTick(x, span);
    const textWidth = ctx.measureText(label).width;

    ctx.fillText(
      label,
      Math.min(
        width - plot.right - textWidth,
        Math.max(plot.left, xPixel(x) - textWidth / 2)
      ),
      height - 12
    );
  }
}

/* ==================== 更新理论表格与页面 ==================== */

function setText(id, text) {
  document.getElementById(id).textContent = text;
}

function updateTheory(n) {
  const k = Math.floor((n + 1) / 2);

  const medianLabel = n % 2 === 0
    ? `较小的中间观测值 X(${k})`
    : `中位数 X(${k})`;

  setText("theory-median-title", medianLabel);
  setText("chart-u-median-title", medianLabel);
  setText("chart-n-median-title", medianLabel);

  setText(
    "theory-u-min",
    `f(x) = ${n}(1−x)^${n - 1}`
  );

  setText(
    "theory-n-min",
    `f(x) = ${n}[1−Φ(x)]^${n - 1} φ(x)`
  );

  setText(
    "theory-u-median",
    n === 5
      ? "f(x) = 30x²(1−x)²"
      : `f(x) = n!/[(${k - 1})!(${n - k})!] · ` +
        `x^${k - 1}(1−x)^${n - k}`
  );

  setText(
    "theory-n-median",
    `f(x) = n!/[(${k - 1})!(${n - k})!] · ` +
    `Φ(x)^${k - 1}[1−Φ(x)]^${n - k} φ(x)`
  );

  setText(
    "theory-u-max",
    `f(x) = ${n}x^${n - 1}`
  );

  setText(
    "theory-n-max",
    `f(x) = ${n}Φ(x)^${n - 1} φ(x)`
  );

  setText(
    "theory-u-range",
    `f(r) = ${n * (n - 1)}r^${n - 2}(1−r)`
  );

  setText(
    "theory-n-range",
    "f(r) = n(n−1)∫φ(t)φ(t+r)" +
    "[Φ(t+r)−Φ(t)]^(n−2)dt"
  );
}

function render() {
  const completed = results.u.min.length;

  statusEl.textContent =
    `n = ${currentN}；随机种子 = ${currentSeed}；` +
    `两栏各累计完成 ${completed.toLocaleString()} 轮实验。` +
    (
      currentN > 30
        ? " 大 n 采用精确联合分布抽样，横轴按需放大。"
        : ""
    );

  updateTheory(currentN);

  for (const dist of DISTRIBUTIONS) {
    for (const type of TYPES) {
      drawChart(
        document.getElementById(`chart-${dist}-${type}`),
        results[dist][type],
        dist,
        type,
        currentN
      );
    }
  }
}

runButton.addEventListener("click", runCustom);
onceButton.addEventListener("click", () => run(1));
resetButton.addEventListener("click", resetFromInputs);

nInput.addEventListener("change", resetFromInputs);
seedInput.addEventListener("change", resetFromInputs);

resetFromInputs();
