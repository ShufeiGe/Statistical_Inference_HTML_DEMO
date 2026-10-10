"use strict";

const $ = id => document.getElementById(id);
let state = null;
let timer;
let generation = 0;

const pause = () => new Promise(resolve => setTimeout(resolve, 0));

function fmt(x) {
  if (!Number.isFinite(x)) return "—";
  if (x !== 0 && (Math.abs(x) < 0.001 || Math.abs(x) >= 100000)) {
    return x.toExponential(3);
  }
  return Number(x.toFixed(4)).toString();
}

// 固定种子的伪随机数生成器
function mulberry32(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t ^= t + Math.imul(t ^ (t >>> 7), 61 | t);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Box–Muller 标准正态生成器
function normalGenerator(random) {
  let spare = null;
  return () => {
    if (spare !== null) {
      const value = spare;
      spare = null;
      return value;
    }
    const radius = Math.sqrt(-2 * Math.log(1 - random()));
    const angle = 2 * Math.PI * random();
    spare = radius * Math.sin(angle);
    return radius * Math.cos(angle);
  };
}

// Lanczos log Gamma 近似
function logGamma(z) {
  const c = [
    676.5203681218851, -1259.1392167224028,
    771.3234287776531, -176.6150291621406,
    12.507343278686905, -0.13857109526572012,
    9.984369578019572e-6, 1.5056327351493116e-7
  ];
  if (z < 0.5) {
    return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z))
      - logGamma(1 - z);
  }
  z -= 1;
  let a = 0.9999999999998099;
  for (let i = 0; i < c.length; i++) a += c[i] / (z + i + 1);
  const t = z + c.length - 0.5;
  return 0.5 * Math.log(2 * Math.PI)
    + (z + 0.5) * Math.log(t) - t + Math.log(a);
}

function normalPDF(x) {
  return Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI);
}

function tPDF(x, df) {
  return Math.exp(
    logGamma((df + 1) / 2) - logGamma(df / 2)
    - 0.5 * Math.log(df * Math.PI)
    - ((df + 1) / 2) * Math.log1p(x * x / df)
  );
}

function summarize(a) {
  let mean = 0, m2 = 0;
  for (let i = 0; i < a.length; i++) {
    const d = a[i] - mean;
    mean += d / (i + 1);
    m2 += d * (a[i] - mean);
  }
  return { mean, variance: m2 / (a.length - 1) };
}

function corr(a, b, meanA, meanB) {
  let aa = 0, bb = 0, ab = 0;
  for (let i = 0; i < a.length; i++) {
    const x = a[i] - meanA, y = b[i] - meanB;
    aa += x * x;
    bb += y * y;
    ab += x * y;
  }
  return ab / Math.sqrt(aa * bb);
}

// 通用密度直方图：支持叠加多个样本、密度曲线和参考竖线
function draw(canvas, series, curves, min, max, xLabel, reference) {
  const dpr = window.devicePixelRatio || 1;
  const width = canvas.clientWidth, height = canvas.clientHeight;
  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);

  const m = { left: 58, right: 18, top: 30, bottom: 46 };
  const w = width - m.left - m.right;
  const h = height - m.top - m.bottom;
  const bins = 48, binWidth = (max - min) / bins;

  const histograms = series.map(s => {
    const counts = new Array(bins).fill(0);
    for (const x of s.values) {
      if (x < min || x > max) continue;
      const index = Math.min(bins - 1, Math.floor((x - min) / binWidth));
      counts[index]++;
    }
    return counts.map(c => c / (s.values.length * binWidth));
  });

  let ymax = 0;
  for (const hist of histograms) {
    for (const y of hist) ymax = Math.max(ymax, y);
  }
  for (const curve of curves) {
    for (let i = 0; i <= 500; i++) {
      ymax = Math.max(ymax, curve.pdf(min + (max - min) * i / 500));
    }
  }
  ymax = (ymax || 1) * 1.18;

  const X = x => m.left + (x - min) / (max - min) * w;
  const Y = y => m.top + h - y / ymax * h;

  ctx.font = "11px system-ui";
  for (let i = 0; i <= 4; i++) {
    const y = ymax * i / 4;
    ctx.strokeStyle = "#e2e8f0";
    ctx.beginPath();
    ctx.moveTo(m.left, Y(y));
    ctx.lineTo(width - m.right, Y(y));
    ctx.stroke();
    ctx.fillStyle = "#64748b";
    ctx.textAlign = "right";
    ctx.fillText(fmt(y), m.left - 7, Y(y) + 4);
  }

  histograms.forEach((hist, k) => {
    ctx.fillStyle = series[k].color;
    hist.forEach((y, i) => {
      ctx.fillRect(
        X(min + i * binWidth) + 0.5, Y(y),
        Math.max(0.5, w / bins - 1), Y(0) - Y(y)
      );
    });
  });

  for (const curve of curves) {
    ctx.strokeStyle = curve.color;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    for (let i = 0; i <= 500; i++) {
      const x = min + (max - min) * i / 500;
      if (i === 0) ctx.moveTo(X(x), Y(curve.pdf(x)));
      else ctx.lineTo(X(x), Y(curve.pdf(x)));
    }
    ctx.stroke();
  }

  if (reference) {
    ctx.strokeStyle = "#dc2626";
    ctx.lineWidth = 2;
    ctx.setLineDash([6, 4]);
    ctx.beginPath();
    ctx.moveTo(X(reference.value), m.top);
    ctx.lineTo(X(reference.value), Y(0));
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#dc2626";
    ctx.textAlign = "center";
    ctx.font = "bold 12px system-ui";
    ctx.fillText(reference.label, X(reference.value), m.top - 9);
  }

  ctx.lineWidth = 1;
  ctx.strokeStyle = "#94a3b8";
  ctx.beginPath();
  ctx.moveTo(m.left, m.top);
  ctx.lineTo(m.left, Y(0));
  ctx.lineTo(width - m.right, Y(0));
  ctx.stroke();

  ctx.font = "11px system-ui";
  ctx.fillStyle = "#64748b";
  ctx.textAlign = "center";
  for (let i = 0; i <= 4; i++) {
    const x = min + (max - min) * i / 4;
    ctx.fillText(fmt(x), X(x), height - 25);
  }
  ctx.fillText(xLabel, m.left + w / 2, height - 6);
  ctx.save();
  ctx.translate(13, m.top + h / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Density", 0, 0);
  ctx.restore();
}

function renderPlots() {
  if (!state) return;
  const { z, t, s, df, sigma, range } = state;

  draw(
    $("ztPlot"),
    [
      { values: z, color: "rgba(37,99,235,.36)" },
      { values: t, color: "rgba(234,88,12,.36)" }
    ],
    [
      { pdf: normalPDF, color: "#2563eb" },
      { pdf: x => tPDF(x, df), color: "#ea580c" }
    ],
    -range, range, "Z / T"
  );

  let maxS = sigma * 1.25;
  for (const value of s) maxS = Math.max(maxS, value);
  draw(
    $("sPlot"),
    [{ values: s, color: "rgba(16,185,129,.5)" }],
    [], 0, maxS * 1.05, "S",
    { value: sigma, label: `σ = ${fmt(sigma)}` }
  );
}

async function runSimulation() {
  if (!$("controls").reportValidity()) return;
  const id = ++generation;
  const n = Number($("n").value);
  const B = Number($("B").value);
  const mu = Number($("mu").value);
  const sigma = Number($("sigma").value);
  const seed = Number($("seed").value);
  const range = Number($("range").value);
  const df = n - 1;
  const normal = normalGenerator(mulberry32(seed));
  const z = new Float64Array(B);
  const t = new Float64Array(B);
  const s = new Float64Array(B);

  $("status").textContent = "正在模拟……";
  const start = performance.now();

  // 分块计算，保持页面响应；新任务会取消旧任务。
  for (let begin = 0; begin < B; begin += 500) {
    await pause();
    if (id !== generation) return;
    const end = Math.min(B, begin + 500);

    for (let b = begin; b < end; b++) {
      let mean = 0, m2 = 0;
      for (let i = 1; i <= n; i++) {
        const x = normal();
        const delta = x - mean;
        mean += delta / i;
        m2 += delta * (x - mean);
      }

      /*
        原始样本为 X_i = mu + sigma * Y_i，Y_i ~ N(0,1)。
        在标准化尺度上计算可避免 mu 很大时的相减误差：
        Z = sqrt(n) * Ybar
        S = sigma * sqrt(sum(Y_i-Ybar)^2 / (n-1))
        T = Z / (S/sigma)
      */
      const relativeS = Math.sqrt(m2 / df);
      z[b] = Math.sqrt(n) * mean;
      s[b] = sigma * relativeS;
      t[b] = z[b] / relativeS;
    }
    $("status").textContent =
      `正在模拟：${end.toLocaleString()} / ${B.toLocaleString()}`;
  }

  if (id !== generation) return;
  state = { z, t, s, n, B, mu, sigma, df, range };
  const zs = summarize(z), ts = summarize(t), ss = summarize(s);

  // E(S) = sigma*c4；Var(S) = sigma²*(1-c4²)
  const c4 = Math.exp(
    0.5 * Math.log(2 / df)
    + logGamma((df + 1) / 2) - logGamma(df / 2)
  );

  const rows = [
    ["Z", zs.mean, "0", zs.variance, "1"],
    ["T", ts.mean, df > 1 ? "0" : "不存在",
      ts.variance, df > 2 ? fmt(df / (df - 2))
        : (df === 2 ? "不有限（∞）" : "不存在")],
    ["S", ss.mean, fmt(sigma * c4), ss.variance,
      fmt(sigma * sigma * Math.max(0, 1 - c4 * c4))]
  ];

  $("stats").innerHTML = rows.map(row => `
    <tr>
      <td>${row[0]}</td>
      <td>${fmt(row[1])}</td>
      <td>${row[2]}</td>
      <td>${fmt(row[3])}</td>
      <td>${row[4]}</td>
    </tr>`).join("");

  const correlations = [
    ["Corr(Z, T)", corr(z, t, zs.mean, ts.mean)],
    ["Corr(Z, S)", corr(z, s, zs.mean, ss.mean)],
    ["Corr(T, S)", corr(t, s, ts.mean, ss.mean)]
  ];
  $("correlations").innerHTML = correlations.map(([name, value]) => `
    <div class="badge">${name}<br><strong>${fmt(value)}</strong></div>
  `).join("");

  $("momentNote").textContent = df <= 2
    ? "注意：当前 T 没有有限方差，涉及 T 的理论 Pearson 相关系数未定义。" +
      "页面显示的有限样本均值、方差和相关系数仍可计算，但可能极不稳定。"
    : `当前 T 的理论均值为 0，理论方差为 ${fmt(df / (df - 2))}。` +
      "增大 B 通常改善估计精度，但重尾可能使收敛较慢。";

  let outZ = 0, outT = 0, small = 0, half = 0;
  for (let i = 0; i < B; i++) {
    if (Math.abs(z[i]) > range) outZ++;
    if (Math.abs(t[i]) > range) outT++;
    if (s[i] < sigma) small++;
    if (s[i] < sigma / 2) half++;
  }

  $("outside").textContent =
    `图外比例（|值| > ${range}）：` +
    `Z ${(100 * outZ / B).toFixed(2)}%；` +
    `T ${(100 * outT / B).toFixed(2)}%。`;

  $("smallS").textContent =
    `经验比例：S < σ 为 ${(100 * small / B).toFixed(2)}%；` +
    `S < σ/2 为 ${(100 * half / B).toFixed(2)}%。`;

  $("theory").textContent =
    `当前：μ = ${fmt(mu)}，σ = ${fmt(sigma)}，ν = n − 1 = ${df}`;

  renderPlots();
  $("status").textContent =
    `完成 ${B.toLocaleString()} 次模拟，每次样本量 n = ${n}；` +
    `耗时 ${Math.round(performance.now() - start)} ms。`;
}

$("controls").addEventListener("submit", event => {
  event.preventDefault();
  clearTimeout(timer);
  runSimulation();
});

$("controls").addEventListener("input", () => {
  ++generation;
  clearTimeout(timer);
  $("status").textContent = "参数已修改，等待更新……";
  timer = setTimeout(() => {
    if ($("controls").checkValidity()) runSimulation();
    else $("status").textContent = "请输入有效参数。";
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
