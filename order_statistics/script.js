"use strict";

const distributionInput = document.getElementById("distribution");
const nInput = document.getElementById("sampleSize");
const bInput = document.getElementById("simulationCount");
const seedInput = document.getElementById("seedInput");

const statusEl = document.getElementById("status");
const latestEl = document.getElementById("latest");

const runButton = document.getElementById("runButton");
const onceButton = document.getElementById("onceButton");
const resetButton = document.getElementById("resetButton");

const data = {
  min: [],
  median: [],
  max: [],
  range: []
};

let currentParameters = null;
let random = null;
let busy = false;

const SQRT_2PI = Math.sqrt(2 * Math.PI);

/* ---------- 带种子的伪随机数 ---------- */

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

    const u1 = 1 - uniform(); // 避免 log(0)
    const u2 = uniform();
    const radius = Math.sqrt(-2 * Math.log(u1));
    const angle = 2 * Math.PI * u2;

    spareNormal = radius * Math.sin(angle);
    return radius * Math.cos(angle);
  }

  // Marsaglia–Tsang Gamma(shape, 1)；本程序中 shape >= 1。
  function gamma(shape) {
    const d = shape - 1 / 3;
    const c = 1 / Math.sqrt(9 * d);

    while (true) {
      const z = normal();
      const v0 = 1 + c * z;

      if (v0 <= 0) continue;

      const v = v0 ** 3;
      const u = uniform();

      if (u < 1 - 0.0331 * z ** 4) {
        return d * v;
      }

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

  return { uniform, normal, beta };
}

/* ---------- 标准正态 PDF、CDF 和分位数 ---------- */

function normalPDF(x) {
  return Math.exp(-x * x / 2) / SQRT_2PI;
}

// 标准正态 CDF 的数值近似。
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

/*
 * 标准正态分位数的 Acklam 近似。
 * 用于将 Uniform(0,1) 次序统计量映射为正态次序统计量。
 */
function normalQuantile(p) {
  const eps = 1e-15;
  p = Math.max(eps, Math.min(1 - eps, p));

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

  const lower = 0.02425;
  const upper = 1 - lower;

  if (p < lower || p > upper) {
    const q = Math.sqrt(
      -2 * Math.log(p < lower ? p : 1 - p)
    );

    const numerator =
      (((((c[0] * q + c[1]) * q + c[2]) * q +
          c[3]) * q + c[4]) * q + c[5]);

    const denominator =
      ((((d[0] * q + d[1]) * q + d[2]) * q +
          d[3]) * q + 1);

    const result = numerator / denominator;
    return p < lower ? result : -result;
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

/* ---------- 输入与状态管理 ---------- */

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

  return {
    distribution: distributionInput.value,
    n,
    seed
  };
}

function sameParameters(a, b) {
  return (
    a !== null &&
    a.distribution === b.distribution &&
    a.n === b.n &&
    a.seed === b.seed
  );
}

function setControlsDisabled(disabled) {
  for (const control of [
    distributionInput,
    nInput,
    bInput,
    seedInput,
    runButton,
    onceButton,
    resetButton
  ]) {
    control.disabled = disabled;
  }
}

function resetWith(parameters) {
  currentParameters = parameters;
  random = makeRandom(parameters.seed);

  for (const values of Object.values(data)) {
    values.length = 0;
  }

  latestEl.textContent = "尚未运行实验。";
  render();
}

function resetFromInputs() {
  if (busy) return;

  const parameters = readParameters();
  if (parameters !== null) resetWith(parameters);
}

/* ---------- 模拟 ---------- */

/*
 * n <= 30：生成完整样本并排序。
 */
function simulateBySorting(n, distribution) {
  const draw = distribution === "normal"
    ? random.normal
    : random.uniform;

  const values = Array.from(
    { length: n },
    () => draw()
  ).sort((a, b) => a - b);

  const k = Math.floor((n + 1) / 2);

  return {
    min: values[0],
    median: values[k - 1],
    max: values[n - 1],
    sorted: values
  };
}

/*
 * n > 30：先模拟三个 Uniform 次序统计量的精确联合分布，
 * 再用目标分布的分位数函数逐个映射。
 *
 * 单调变换保留次序和联合关系，因此极差也正确。
 */
function simulateFast(n, distribution) {
  const k = Math.floor((n + 1) / 2);

  const uMin =
    -Math.expm1(Math.log1p(-random.uniform()) / n);

  const middleFraction = random.beta(
    k - 1,
    n - k + 1
  );

  const uMiddle =
    uMin + (1 - uMin) * middleFraction;

  const upperFraction =
    Math.pow(random.uniform(), 1 / (n - k));

  const uMax =
    uMiddle + (1 - uMiddle) * upperFraction;

  const transform = distribution === "normal"
    ? normalQuantile
    : x => x;

  return {
    min: transform(uMin),
    median: transform(uMiddle),
    max: transform(uMax),
    sorted: null
  };
}

function run(count) {
  if (busy) return;

  const parameters = readParameters();
  if (parameters === null) return;

  if (!sameParameters(currentParameters, parameters)) {
    resetWith(parameters);
  }

  busy = true;
  setControlsDisabled(true);

  let completedInThisRun = 0;
  let last = null;

  function processBatch() {
    const batchSize = 500;
    const take = Math.min(
      batchSize,
      count - completedInThisRun
    );

    for (let i = 0; i < take; i++) {
      const result = parameters.n <= 30
        ? simulateBySorting(
            parameters.n,
            parameters.distribution
          )
        : simulateFast(
            parameters.n,
            parameters.distribution
          );

      data.min.push(result.min);
      data.median.push(result.median);
      data.max.push(result.max);
      data.range.push(result.max - result.min);

      last = result;
    }

    completedInThisRun += take;

    if (completedInThisRun < count) {
      statusEl.textContent =
        `正在模拟：本次已完成 ` +
        `${completedInThisRun.toLocaleString()} / ` +
        `${count.toLocaleString()} 轮……`;

      setTimeout(processBatch, 0);
      return;
    }

    if (last.sorted !== null) {
      latestEl.textContent =
        "最近一轮，排序后的样本：" +
        last.sorted.map(x => x.toFixed(4)).join("，");
    } else {
      latestEl.textContent =
        "最近一轮：" +
        `最小值 ${last.min.toFixed(5)}；` +
        `中间次序统计量 ${last.median.toFixed(5)}；` +
        `最大值 ${last.max.toFixed(5)}；` +
        `极差 ${(last.max - last.min).toFixed(5)}。` +
        "大 n 时使用精确联合分布抽样，不生成完整样本。";
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

/* ---------- 理论密度 ---------- */

const logFactorial = [0];

for (let j = 1; j <= 10000; j++) {
  logFactorial[j] =
    logFactorial[j - 1] + Math.log(j);
}

/*
 * 标准正态极差的理论密度：
 *
 * f_R(r) = n(n-1) ∫ φ(t)φ(t+r)
 *          [Φ(t+r)-Φ(t)]^(n-2) dt, r > 0。
 *
 * 这里对 t 使用复合 Simpson 法数值积分。
 */
function normalRangeDensity(r, n) {
  if (r <= 0) return 0;

  const left = -9;
  const right = 9;
  const segments = 160; // Simpson 法要求偶数
  const step = (right - left) / segments;

  function integrand(t) {
    const probability =
      normalCDF(t + r) - normalCDF(t);

    if (probability <= 0) return 0;

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

function density(type, x, n, distribution) {
  const k = Math.floor((n + 1) / 2);

  if (distribution === "uniform") {
    if (x < 0 || x > 1) return 0;

    if (type === "min") {
      return n * Math.pow(1 - x, n - 1);
    }

    if (type === "max") {
      return n * Math.pow(x, n - 1);
    }

    if (type === "range") {
      return x <= 0 || x >= 1
        ? 0
        : n * (n - 1) *
          Math.pow(x, n - 2) * (1 - x);
    }
  }

  if (distribution === "normal" && type === "range") {
    return normalRangeDensity(x, n);
  }

  const cdf = distribution === "normal"
    ? normalCDF(x)
    : x;

  const pdf = distribution === "normal"
    ? normalPDF(x)
    : 1;

  if (type === "min") {
    return n * Math.pow(1 - cdf, n - 1) * pdf;
  }

  if (type === "max") {
    return n * Math.pow(cdf, n - 1) * pdf;
  }

  if (cdf <= 0 || cdf >= 1) return 0;

  const logCoefficient =
    logFactorial[n] -
    logFactorial[k - 1] -
    logFactorial[n - k];

  return Math.exp(
    logCoefficient +
    (k - 1) * Math.log(cdf) +
    (n - k) * Math.log1p(-cdf)
  ) * pdf;
}

/* ---------- 图表 ---------- */

function chartDomain(type, n, distribution) {
  if (distribution === "uniform") {
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
    return [-q - 5 * extremeSD,
            -q + 5 * extremeSD];
  }

  if (type === "max") {
    return [q - 5 * extremeSD,
            q + 5 * extremeSD];
  }

  if (type === "range") {
    const center = 2 * q;
    const halfWidth =
      6 * Math.SQRT2 * extremeSD;

    return [
      Math.max(0, center - halfWidth),
      center + halfWidth
    ];
  }

  const k = Math.floor((n + 1) / 2);
  const p = k / (n + 1);
  const center = normalQuantile(p);

  const uniformVariance =
    k * (n - k + 1) /
    ((n + 1) ** 2 * (n + 2));

  const sd =
    Math.sqrt(uniformVariance) /
    normalPDF(center);

  return [center - 5 * sd, center + 5 * sd];
}

function formatTick(x, domainWidth) {
  if (Math.abs(x) < 1e-12) x = 0;

  if (domainWidth < 0.02) return x.toFixed(4);
  if (domainWidth < 0.2) return x.toFixed(3);
  return x.toFixed(2);
}

function drawChart(canvas, values, type, parameters) {
  const { n, distribution } = parameters;
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

  const [lo, hi] = chartDomain(
    type, n, distribution
  );

  const domainWidth = hi - lo;
  const bins = 35;
  const binWidth = domainWidth / bins;
  const counts = Array(bins).fill(0);

  for (const value of values) {
    if (value < lo || value > hi) continue;

    const index = Math.min(
      bins - 1,
      Math.floor((value - lo) / binWidth)
    );

    counts[index]++;
  }

  // 放大横轴时，仍以全部实验次数为分母。
  const histogram = counts.map(count =>
    values.length
      ? count / (values.length * binWidth)
      : 0
  );

  /*
   * 先计算理论曲线，避免标准正态极差的数值积分
   * 在求纵轴范围和画线时重复执行。
   */
  const curvePoints = 180;
  const theoretical = [];

  let maxY = Math.max(1, ...histogram);

  for (let i = 0; i <= curvePoints; i++) {
    const x = lo + domainWidth * i / curvePoints;
    const y = density(type, x, n, distribution);

    theoretical.push({ x, y });
    maxY = Math.max(maxY, y);
  }

  maxY *= 1.1;

  const xPixel = x =>
    plot.left +
    (x - lo) / domainWidth * plotWidth;

  const yPixel = y =>
    plot.top + plotHeight -
    y / maxY * plotHeight;

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

  // 经验直方图。
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

  // 理论曲线。
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

  // 横轴。
  ctx.strokeStyle = "#94a3b8";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(plot.left, plot.top + plotHeight);
  ctx.lineTo(width - plot.right, plot.top + plotHeight);
  ctx.stroke();

  ctx.fillStyle = "#64748b";
  ctx.font = "11px system-ui, sans-serif";

  for (let i = 0; i <= 4; i++) {
    const x = lo + domainWidth * i / 4;
    const label = formatTick(x, domainWidth);
    const textWidth = ctx.measureText(label).width;

    ctx.fillText(
      label,
      Math.min(
        width - plot.right - textWidth,
        Math.max(
          plot.left,
          xPixel(x) - textWidth / 2
        )
      ),
      height - 12
    );
  }
}

/* ---------- 文字和事件 ---------- */

function updateLabels(parameters) {
  const { n, distribution } = parameters;
  const k = Math.floor((n + 1) / 2);

  document.getElementById("median-title").innerHTML =
    n % 2 === 0
      ? `较小的中间观测值 X<sub>(${k})</sub>`
      : `中位数 X<sub>(${k})</sub>`;

  if (distribution === "normal") {
    document.getElementById("formula-min").textContent =
      "f(x) = n[1 − Φ(x)]^(n−1) φ(x)";

    document.getElementById("formula-median").textContent =
      `f(x) = n!/[(${k - 1})!(${n - k})!] · ` +
      `Φ(x)^${k - 1}[1−Φ(x)]^${n - k} φ(x)`;

    document.getElementById("formula-max").textContent =
      "f(x) = nΦ(x)^(n−1) φ(x)";

    document.getElementById("formula-range").textContent =
      "f(r) = n(n−1)∫φ(t)φ(t+r)[Φ(t+r)−Φ(t)]^(n−2)dt，r > 0";

    return;
  }

  document.getElementById("formula-min").textContent =
    `f(x) = ${n}(1−x)^${n - 1}`;

  document.getElementById("formula-median").textContent =
    n === 5
      ? "f(x) = 30x²(1−x)²"
      : `f(x) = n!/[(${k - 1})!(${n - k})!] · ` +
        `x^${k - 1}(1−x)^${n - k}`;

  document.getElementById("formula-max").textContent =
    `f(x) = ${n}x^${n - 1}`;

  document.getElementById("formula-range").textContent =
    `f(r) = ${n * (n - 1)}r^${n - 2}(1−r)`;
}

function render() {
  const parameters = currentParameters;
  const completed = data.min.length;

  const distributionName =
    parameters.distribution === "normal"
      ? "标准正态 N(0,1)"
      : "Uniform(0,1)";

  statusEl.textContent =
    `分布：${distributionName}；` +
    `n = ${parameters.n}；` +
    `随机种子 = ${parameters.seed}；` +
    `累计完成 ${completed.toLocaleString()} 轮实验。` +
    (
      parameters.n > 30
        ? " 大 n 使用精确联合分布抽样；图表横轴按需放大。"
        : ""
    );

  updateLabels(parameters);

  for (const type of ["min", "median", "max", "range"]) {
    drawChart(
      document.getElementById(`chart-${type}`),
      data[type],
      type,
      parameters
    );
  }
}

runButton.addEventListener("click", runCustom);

onceButton.addEventListener("click", () => {
  run(1);
});

resetButton.addEventListener(
  "click",
  resetFromInputs
);

distributionInput.addEventListener(
  "change",
  resetFromInputs
);

nInput.addEventListener(
  "change",
  resetFromInputs
);

seedInput.addEventListener(
  "change",
  resetFromInputs
);

resetFromInputs();
