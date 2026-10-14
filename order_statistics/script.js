"use strict";

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

let currentN = null;
let currentSeed = null;
let random = null;
let busy = false;

/*
 * Mulberry32：可重复的教学用伪随机数生成器。
 * 不适用于密码学。
 */
function mulberry32(seed) {
  let state = seed >>> 0;

  return function () {
    state = (state + 0x6D2B79F5) | 0;

    let t = Math.imul(state ^ (state >>> 15), 1 | state);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;

    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/*
 * 所有抽样共用一条由种子确定的随机数序列。
 * n 较大时，用 Gamma 抽样构造 Beta 随机变量。
 */
function makeRandom(seed) {
  const uniform = mulberry32(seed);
  let spareNormal = null;

  function normal() {
    if (spareNormal !== null) {
      const result = spareNormal;
      spareNormal = null;
      return result;
    }

    const u1 = 1 - uniform();
    const u2 = uniform();
    const radius = Math.sqrt(-2 * Math.log(u1));
    const angle = 2 * Math.PI * u2;

    spareNormal = radius * Math.sin(angle);
    return radius * Math.cos(angle);
  }

  // Marsaglia–Tsang Gamma(shape, 1) 抽样；此处 shape 均 >= 1。
  function gamma(shape) {
    const d = shape - 1 / 3;
    const c = 1 / Math.sqrt(9 * d);

    while (true) {
      const z = normal();
      const v0 = 1 + c * z;

      if (v0 <= 0) continue;

      const v = v0 * v0 * v0;
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

  return { uniform, beta };
}

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
    seedInput,
    0,
    4294967295,
    "随机种子"
  );
  if (seed === null) return null;

  return { n, seed };
}

function setControlsDisabled(disabled) {
  nInput.disabled = disabled;
  bInput.disabled = disabled;
  seedInput.disabled = disabled;
  runButton.disabled = disabled;
  onceButton.disabled = disabled;
  resetButton.disabled = disabled;
}

function resetWith(parameters) {
  currentN = parameters.n;
  currentSeed = parameters.seed;
  random = makeRandom(currentSeed);

  for (const values of Object.values(data)) {
    values.length = 0;
  }

  latestEl.textContent = "尚未运行实验。";
  render();
}

function resetFromInputs() {
  if (busy) return;

  const parameters = readParameters();

  if (parameters !== null) {
    resetWith(parameters);
  }
}

/*
 * n <= 30：按 slides 的步骤生成完整样本并排序。
 */
function simulateBySorting(n) {
  const values = Array.from(
    { length: n },
    () => random.uniform()
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
 * n > 30：直接模拟三个次序统计量的联合分布。
 *
 * M = X_(1) ~ Beta(1, n)
 * (T - M) / (1 - M) ~ Beta(k - 1, n - k + 1)
 * (H - T) / (1 - T) ~ Beta(n - k, 1)
 *
 * T = X_(k)，H = X_(n)。
 * 保留正确的联合分布，因此 H - M 也是正确分布的极差。
 */
function simulateFast(n) {
  const k = Math.floor((n + 1) / 2);

  const minimum =
    -Math.expm1(Math.log1p(-random.uniform()) / n);

  const middleFraction = random.beta(
    k - 1,
    n - k + 1
  );

  const middle =
    minimum + (1 - minimum) * middleFraction;

  const upperFraction =
    Math.pow(random.uniform(), 1 / (n - k));

  const maximum =
    middle + (1 - middle) * upperFraction;

  return {
    min: minimum,
    median: middle,
    max: maximum,
    sorted: null
  };
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

  let completedInThisRun = 0;
  let last = null;

  // 分批计算，避免一次运行时长时间阻塞页面。
  function processBatch() {
    const batchSize = 500;
    const take = Math.min(
      batchSize,
      count - completedInThisRun
    );

    for (let i = 0; i < take; i++) {
      const result = currentN <= 30
        ? simulateBySorting(currentN)
        : simulateFast(currentN);

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
        `最小值 ${last.min.toFixed(6)}；` +
        `中间次序统计量 ${last.median.toFixed(6)}；` +
        `最大值 ${last.max.toFixed(6)}；` +
        `极差 ${(last.max - last.min).toFixed(6)}。` +
        "大 n 时采用等价的联合分布抽样，不生成完整排序样本。";
    }

    busy = false;
    setControlsDisabled(false);
    render();
  }

  setTimeout(processBatch, 0);
}

function runCustom() {
  const count = readInteger(
    bInput,
    1,
    20000,
    "模拟次数 B"
  );

  if (count !== null) {
    run(count);
  }
}

/*
 * 预先计算 log(j!)，避免 n=10000 时阶乘溢出。
 */
const logFactorial = [0];

for (let j = 1; j <= 10000; j++) {
  logFactorial[j] =
    logFactorial[j - 1] + Math.log(j);
}

function density(type, x, n) {
  if (x < 0 || x > 1) return 0;

  const k = Math.floor((n + 1) / 2);

  if (type === "min") {
    return n * Math.pow(1 - x, n - 1);
  }

  if (type === "max") {
    return n * Math.pow(x, n - 1);
  }

  if (type === "median") {
    if (x === 0 || x === 1) return 0;

    const logCoefficient =
      logFactorial[n] -
      logFactorial[k - 1] -
      logFactorial[n - k];

    return Math.exp(
      logCoefficient +
      (k - 1) * Math.log(x) +
      (n - k) * Math.log1p(-x)
    );
  }

  if (type === "range") {
    if (x === 0 || x === 1) return 0;

    return (
      n * (n - 1) *
      Math.pow(x, n - 2) *
      (1 - x)
    );
  }

  return 0;
}

/*
 * n 较大时放大密度集中的横轴区域。
 * 坐标刻度仍显示实际 x 值。
 */
function chartDomain(type, n) {
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

function formatTick(x, domainWidth) {
  return domainWidth < 0.02
    ? x.toFixed(4)
    : x.toFixed(2);
}

function drawChart(canvas, values, type, n) {
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;

  const plot = {
    left: 58,
    right: 16,
    top: 16,
    bottom: 40
  };

  const plotWidth = width - plot.left - plot.right;
  const plotHeight = height - plot.top - plot.bottom;

  const [lo, hi] = chartDomain(type, n);
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

  // 使用总模拟次数作为分母，保持密度尺度正确。
  const histogram = counts.map(count =>
    values.length
      ? count / (values.length * binWidth)
      : 0
  );

  let maxY = Math.max(1, ...histogram);

  for (let i = 0; i <= 400; i++) {
    const x = lo + domainWidth * i / 400;
    maxY = Math.max(maxY, density(type, x, n));
  }

  maxY *= 1.1;

  const xPixel = x =>
    plot.left + (x - lo) / domainWidth * plotWidth;

  const yPixel = y =>
    plot.top + plotHeight - y / maxY * plotHeight;

  ctx.clearRect(0, 0, width, height);
  ctx.font = "12px system-ui, sans-serif";

  // 网格线和纵轴刻度。
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

  // 模拟直方图。
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

  for (let i = 0; i <= 600; i++) {
    const x = lo + domainWidth * i / 600;
    const px = xPixel(x);
    const py = yPixel(density(type, x, n));

    if (i === 0) {
      ctx.moveTo(px, py);
    } else {
      ctx.lineTo(px, py);
    }
  }

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
    const x = lo + domainWidth * i / 4;
    const label = formatTick(x, domainWidth);
    const px = xPixel(x);
    const textWidth = ctx.measureText(label).width;

    ctx.fillText(
      label,
      Math.min(
        width - plot.right - textWidth,
        Math.max(plot.left, px - textWidth / 2)
      ),
      height - 12
    );
  }
}

function updateLabels(n) {
  const k = Math.floor((n + 1) / 2);

  document.getElementById("median-title").innerHTML =
    n % 2 === 0
      ? `较小的中间观测值 X<sub>(${k})</sub>`
      : `中位数 X<sub>(${k})</sub>`;

  document.getElementById("formula-min").textContent =
    `f(x) = ${n}(1 − x)^${n - 1}`;

  document.getElementById("formula-median").textContent =
    n === 5
      ? "f(x) = 30x²(1 − x)²"
      : `f(x) = n! / [(${k - 1})!(${n - k})!] ` +
        `· x^${k - 1}(1 − x)^${n - k}`;

  document.getElementById("formula-max").textContent =
    `f(x) = ${n}x^${n - 1}`;

  document.getElementById("formula-range").textContent =
    `f(r) = ${n * (n - 1)}r^${n - 2}(1 − r)`;
}

function render() {
  const completed = data.min.length;

  statusEl.textContent =
    `样本量 n = ${currentN}；` +
    `随机种子 = ${currentSeed}；` +
    `累计完成 ${completed.toLocaleString()} 轮实验。` +
    (
      currentN > 30
        ? " 当前采用精确联合分布抽样；图表横轴按需放大。"
        : ""
    );

  updateLabels(currentN);

  for (const type of ["min", "median", "max", "range"]) {
    drawChart(
      document.getElementById(`chart-${type}`),
      data[type],
      type,
      currentN
    );
  }
}

runButton.addEventListener("click", runCustom);

onceButton.addEventListener("click", () => {
  run(1);
});

resetButton.addEventListener("click", resetFromInputs);

nInput.addEventListener("change", resetFromInputs);
seedInput.addEventListener("change", resetFromInputs);

resetFromInputs();
