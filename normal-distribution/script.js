"use strict";

const $ = id => document.getElementById(id);
let current = null;
let timer = null;

function fmt(x) {
  if (!Number.isFinite(x)) return "—";
  if (x !== 0 && (Math.abs(x) < 0.001 || Math.abs(x) >= 100000)) {
    return x.toExponential(3);
  }
  return Number(x.toFixed(4)).toString();
}

/* 可复现的伪随机数生成器 */
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t ^= t + Math.imul(t ^ (t >>> 7), 61 | t);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* Box–Muller：生成标准正态随机数 */
function normalGenerator(random) {
  let spare = null;

  return function () {
    if (spare !== null) {
      const value = spare;
      spare = null;
      return value;
    }

    const u = 1 - random();
    const v = random();
    const radius = Math.sqrt(-2 * Math.log(u));
    const angle = 2 * Math.PI * v;

    spare = radius * Math.sin(angle);
    return radius * Math.cos(angle);
  };
}

/* Lanczos 近似：log Gamma */
function logGamma(z) {
  const c = [
    676.5203681218851,
    -1259.1392167224028,
    771.3234287776531,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.984369578019572e-6,
    1.5056327351493116e-7
  ];

  if (z < 0.5) {
    return Math.log(Math.PI)
      - Math.log(Math.sin(Math.PI * z))
      - logGamma(1 - z);
  }

  z -= 1;
  let a = 0.9999999999998099;
  for (let i = 0; i < c.length; i++) {
    a += c[i] / (z + i + 1);
  }

  const t = z + c.length - 0.5;
  return 0.5 * Math.log(2 * Math.PI)
    + (z + 0.5) * Math.log(t) - t + Math.log(a);
}

function normalDensity(x, mu, sd) {
  const z = (x - mu) / sd;
  return Math.exp(-0.5 * z * z) / (sd * Math.sqrt(2 * Math.PI));
}

function chiSquareDensity(x, df) {
  if (x < 0) return 0;
  if (x === 0) {
    if (df < 2) return Infinity;
    return df === 2 ? 0.5 : 0;
  }

  const a = df / 2;
  return Math.exp(
    (a - 1) * Math.log(x)
    - x / 2 - a * Math.log(2) - logGamma(a)
  );
}

/* 经验均值与无偏经验方差 */
function summarize(values) {
  let mean = 0;
  let m2 = 0;

  for (let i = 0; i < values.length; i++) {
    const d = values[i] - mean;
    mean += d / (i + 1);
    m2 += d * (values[i] - mean);
  }

  return { mean, variance: m2 / (values.length - 1) };
}

function correlation(x, y, mx, my) {
  let xx = 0, yy = 0, xy = 0;

  for (let i = 0; i < x.length; i++) {
    const dx = x[i] - mx;
    const dy = y[i] - my;
    xx += dx * dx;
    yy += dy * dy;
    xy += dx * dy;
  }

  return xy / Math.sqrt(xx * yy);
}

/* Canvas 密度直方图与理论密度曲线 */
function drawPlot(canvas, values, density, title, xLabel, bounds) {
  const dpr = window.devicePixelRatio || 1;
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;

  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);

  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);

  const margin = { left: 58, right: 15, top: 33, bottom: 41 };
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;

  let min = bounds[0], max = bounds[1];

  for (const value of values) {
    min = Math.min(min, value);
    max = Math.max(max, value);
  }

  const bins = 36;
  const binW = (max - min) / bins;
  const counts = new Array(bins).fill(0);

  for (const value of values) {
    const i = Math.max(
      0, Math.min(bins - 1, Math.floor((value - min) / binW))
    );
    counts[i]++;
  }

  const heights = counts.map(count => count / (values.length * binW));

  let ymax = Math.max(...heights);
  for (let i = 0; i < bins; i++) {
    ymax = Math.max(ymax, density(min + (i + 0.5) * binW));
  }
  ymax *= 1.2;

  const X = x => margin.left + (x - min) / (max - min) * plotW;
  const Y = y => margin.top + plotH - y / ymax * plotH;

  ctx.clearRect(0, 0, width, height);
  ctx.font = "bold 13px system-ui";
  ctx.fillStyle = "#1e293b";
  ctx.fillText(title, margin.left, 18);

  ctx.font = "11px system-ui";

  for (let i = 0; i <= 4; i++) {
    const y = ymax * i / 4;
    const py = Y(y);

    ctx.strokeStyle = "#e2e8f0";
    ctx.beginPath();
    ctx.moveTo(margin.left, py);
    ctx.lineTo(width - margin.right, py);
    ctx.stroke();

    ctx.fillStyle = "#64748b";
    ctx.textAlign = "right";
    ctx.fillText(fmt(y), margin.left - 7, py + 4);
  }

  ctx.fillStyle = "#93c5fd";
  for (let i = 0; i < bins; i++) {
    const left = X(min + i * binW);
    const top = Y(heights[i]);
    ctx.fillRect(
      left + 0.5, top,
      Math.max(0.5, plotW / bins - 1),
      margin.top + plotH - top
    );
  }

  /* 裁剪曲线：df=1 时，卡方密度在原点趋于无穷 */
  ctx.save();
  ctx.beginPath();
  ctx.rect(margin.left, margin.top, plotW, plotH);
  ctx.clip();

  ctx.strokeStyle = "#dc2626";
  ctx.lineWidth = 2;
  ctx.beginPath();

  const points = 600;
  for (let i = 0; i <= points; i++) {
    const x = min + (max - min) * i / points;
    const y = Math.min(density(x), ymax * 2);
    if (i === 0) ctx.moveTo(X(x), Y(y));
    else ctx.lineTo(X(x), Y(y));
  }

  ctx.stroke();
  ctx.restore();

  ctx.strokeStyle = "#94a3b8";
  ctx.beginPath();
  ctx.moveTo(margin.left, margin.top);
  ctx.lineTo(margin.left, margin.top + plotH);
  ctx.lineTo(width - margin.right, margin.top + plotH);
  ctx.stroke();

  ctx.textAlign = "center";
  ctx.fillStyle = "#64748b";

  for (let i = 0; i <= 4; i++) {
    const x = min + (max - min) * i / 4;
    ctx.fillText(fmt(x), X(x), height - 23);
  }

  ctx.fillText(xLabel, margin.left + plotW / 2, height - 5);

  ctx.save();
  ctx.translate(12, margin.top + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Density", 0, 0);
  ctx.restore();
}

function renderPlots() {
  if (!current) return;

  const { xbar, qvar, n, mu, sigma } = current;
  const sd = sigma / Math.sqrt(n);
  const df = n - 1;

  drawPlot(
    $("meanPlot"), xbar,
    x => normalDensity(x, mu, sd),
    "样本均值的抽样分布",
    "X̄",
    [mu - 4 * sd, mu + 4 * sd]
  );

  drawPlot(
    $("variancePlot"), qvar,
    x => chiSquareDensity(x, df),
    "标准化样本方差的抽样分布",
    "Q = (n − 1)S² / σ²",
    [0, df + 5 * Math.sqrt(2 * df)]
  );
}

async function runSimulation() {
  const form = $("controls");
  if (!form.reportValidity()) return;

  const n = Number($("n").value);
  const B = Number($("B").value);
  const mu = Number($("mu").value);
  const sigma = Number($("sigma").value);
  const seed = Number($("seed").value);

  $("run").disabled = true;
  $("status").textContent = "正在模拟……";

  /* 先让浏览器显示状态，再执行计算 */
  await new Promise(resolve => requestAnimationFrame(resolve));

  const start = performance.now();
  const normal = normalGenerator(mulberry32(seed));
  const xbar = new Float64Array(B);
  const qvar = new Float64Array(B);

  for (let b = 0; b < B; b++) {
    let meanZ = 0;
    let m2Z = 0;

    for (let i = 1; i <= n; i++) {
      const z = normal();
      const delta = z - meanZ;
      meanZ += delta / i;
      m2Z += delta * (z - meanZ);
    }

    /*
      X_i = mu + sigma * Z_i
      (n-1)S_X²/sigma² = sum (Z_i - Zbar)²
    */
    xbar[b] = mu + sigma * meanZ;
    qvar[b] = m2Z;
  }

  const xs = summarize(xbar);
  const qs = summarize(qvar);
  const corr = correlation(xbar, qvar, xs.mean, qs.mean);

  current = { xbar, qvar, n, B, mu, sigma };

  $("meanTheory").textContent =
    `X̄ ∼ N(${fmt(mu)}, ${fmt(sigma * sigma / n)})`;

  $("varianceTheory").textContent =
    `Q ∼ χ²(${n - 1})`;

  $("rCode").textContent =
`set.seed(${seed})

B <- ${B}
n <- ${n}
mu <- ${mu}
sigma <- ${sigma}

xbar <- numeric(B)
qvar <- numeric(B)

for (b in 1:B) {
  x <- rnorm(n, mu, sigma)
  xbar[b] <- mean(x)
  qvar[b] <- (n - 1) *
    var(x) / sigma^2
}

hist(xbar, probability = TRUE)
curve(
  dnorm(x, mu, sigma / sqrt(n)),
  add = TRUE, col = "red"
)

hist(qvar, probability = TRUE)
curve(
  dchisq(x, df = n - 1),
  add = TRUE, col = "red"
)

cor(xbar, qvar)`;

  const rows = [
    ["Mean(X̄)", xs.mean, mu],
    ["Var(X̄)", xs.variance, sigma * sigma / n],
    ["Mean(Q)", qs.mean, n - 1],
    ["Var(Q)", qs.variance, 2 * (n - 1)]
  ];

  $("stats").innerHTML = rows.map(([name, empirical, theory]) =>
    `<tr>
       <td>${name}</td>
       <td>${fmt(empirical)}</td>
       <td>${fmt(theory)}</td>
     </tr>`
  ).join("");

  $("corr").textContent = corr.toFixed(5);
  renderPlots();

  $("status").textContent =
    `完成 ${B.toLocaleString()} 次模拟，每次抽取 ${n} 个观测值。` +
    `耗时 ${Math.round(performance.now() - start)} ms。`;

  $("run").disabled = false;
}

$("controls").addEventListener("submit", event => {
  event.preventDefault();
  clearTimeout(timer);
  runSimulation();
});

/* 修改参数后自动更新，防抖避免连续输入时反复计算 */
$("controls").addEventListener("input", () => {
  clearTimeout(timer);
  timer = setTimeout(() => {
    if ($("controls").checkValidity()) runSimulation();
  }, 450);
});

$("newSeed").addEventListener("click", () => {
  clearTimeout(timer);
  const seed = new Uint32Array(1);
  crypto.getRandomValues(seed);
  $("seed").value = seed[0];
  runSimulation();
});

let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(renderPlots, 100);
});

runSimulation();
