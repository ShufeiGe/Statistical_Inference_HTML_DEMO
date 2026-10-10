"use strict";

const $ = id => document.getElementById(id);

let state = null;
let generation = 0;
let inputTimer = null;
let resizeTimer = null;
let scrollFrame = null;
let plotCache = null;
let followAll = false;

const pause = () => new Promise(resolve => setTimeout(resolve, 0));

function fmt(x) {
  if (!Number.isFinite(x)) return "—";
  if (x !== 0 && (Math.abs(x) < 0.001 || Math.abs(x) >= 100000)) {
    return x.toExponential(3);
  }
  return Number(x.toFixed(4)).toString();
}

/* 固定种子的伪随机数生成器 */
function mulberry32(seed) {
  let a = seed >>> 0;

  return () => {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t ^= t + Math.imul(t ^ (t >>> 7), 61 | t);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* Box–Muller：生成标准正态随机数 */
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

/* 不完全 Beta 函数的连分式 */
function betaCF(a, b, x) {
  const tiny = 1e-300;
  const eps = 3e-14;
  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;

  let c = 1;
  let d = 1 - qab * x / qap;

  if (Math.abs(d) < tiny) d = tiny;
  d = 1 / d;

  let h = d;

  for (let m = 1; m <= 300; m++) {
    const m2 = 2 * m;

    let aa = m * (b - m) * x /
      ((qam + m2) * (a + m2));

    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;

    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;

    d = 1 / d;
    h *= d * c;

    aa = -(a + m) * (qab + m) * x /
      ((a + m2) * (qap + m2));

    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;

    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;

    d = 1 / d;
    const delta = d * c;
    h *= delta;

    if (Math.abs(delta - 1) < eps) break;
  }

  return h;
}

function regularizedBeta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;

  const factor = Math.exp(
    logGamma(a + b) - logGamma(a) - logGamma(b)
    + a * Math.log(x) + b * Math.log1p(-x)
  );

  if (x < (a + 1) / (a + b + 2)) {
    return factor * betaCF(a, b, x) / a;
  }

  return 1 - factor * betaCF(b, a, 1 - x) / b;
}

/* Student's t 累积分布函数 */
function tCDF(t, df) {
  if (t === 0) return 0.5;

  const x = df / (df + t * t);
  const halfTail = 0.5 * regularizedBeta(x, df / 2, 0.5);

  return t > 0 ? 1 - halfTail : halfTail;
}

/* 二分法计算 t 分位数 */
function tQuantile(p, df) {
  if (p === 0.5) return 0;
  if (p < 0.5) return -tQuantile(1 - p, df);

  let low = 0;
  let high = 1;

  while (tCDF(high, df) < p) high *= 2;

  for (let i = 0; i < 80; i++) {
    const middle = (low + high) / 2;

    if (tCDF(middle, df) < p) low = middle;
    else high = middle;
  }

  return (low + high) / 2;
}

/* 所选区间统一横轴的范围缓存 */
function getPlotCache(count) {
  if (
    plotCache &&
    plotCache.state === state &&
    plotCache.count === count
  ) {
    return plotCache;
  }

  let min = 0;
  let max = 0;
  let covered = 0;

  for (let i = 0; i < count; i++) {
    const ci = state.intervals[i];

    min = Math.min(min, ci.low);
    max = Math.max(max, ci.high);

    if (ci.covered) covered++;
  }

  const padding = Math.max((max - min) * 0.08, 1e-12);

  plotCache = {
    state,
    count,
    min: min - padding,
    max: max + padding,
    covered
  };

  return plotCache;
}

/*
  虚拟滚动绘图：
  提供全部区间所需的滚动高度，
  但 Canvas 仅保持一个可视窗口的尺寸。
*/
function drawPlot() {
  if (!state) return;

  const { intervals, mu, sigma, B } = state;

  const requested = Number($("shown").value);

  const count = followAll
    ? B
    : Math.min(
        B,
        Math.max(1, Number.isFinite(requested) ? Math.floor(requested) : 1)
      );

  $("shown").max = B;
  $("shown").value = count;

  const viewport = $("plotScroll");
  const space = $("plotSpace");
  const canvas = $("plot");

  const rowHeight = 12;
  const margin = { left: 65, right: 25, top: 40, bottom: 50 };

  const width = viewport.clientWidth;
  const height = viewport.clientHeight;

  if (width <= 0 || height <= 0) return;

  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  const totalHeight = Math.max(
    height,
    margin.top + count * rowHeight + margin.bottom
  );

  space.style.height = totalHeight + "px";

  // 若展示数量减少，确保滚动位置仍在有效范围内
  const maxScroll = Math.max(0, totalHeight - height);

  if (viewport.scrollTop > maxScroll) {
    viewport.scrollTop = maxScroll;
  }

  const scrollTop = viewport.scrollTop;

  canvas.style.height = height + "px";
  canvas.style.transform = `translateY(${scrollTop}px)`;

  const dpr = window.devicePixelRatio || 1;

  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);

  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);

  const cache = getPlotCache(count);
  const { min, max, covered } = cache;

  /*
    内部用标准化坐标 (x − μ)/σ，
    刻度再换算回原始尺度，避免大 μ 下的数值相减误差。
  */
  const X = x =>
    margin.left + (x - min) / (max - min) * plotWidth;

  const Y = i =>
    margin.top + (i + 0.5) * rowHeight - scrollTop;

  const first = Math.max(0, Math.floor(scrollTop / rowHeight));
  const last = Math.min(
    count,
    Math.ceil((scrollTop + plotHeight) / rowHeight)
  );

  ctx.fillStyle = "white";
  ctx.fillRect(0, 0, width, height);

  const ticks = width < 600 ? 4 : 6;

  ctx.save();
  ctx.beginPath();
  ctx.rect(margin.left, margin.top, plotWidth, plotHeight);
  ctx.clip();

  // 横轴网格
  ctx.strokeStyle = "#e2e8f0";
  ctx.lineWidth = 1;

  for (let i = 0; i <= ticks; i++) {
    const value = min + (max - min) * i / ticks;

    ctx.beginPath();
    ctx.moveTo(X(value), margin.top);
    ctx.lineTo(X(value), margin.top + plotHeight);
    ctx.stroke();
  }

  // 当前可见的置信区间
  for (let i = first; i < last; i++) {
    const ci = intervals[i];
    const y = Y(i);

    ctx.strokeStyle = ci.covered ? "#2563eb" : "#dc2626";
    ctx.fillStyle = ctx.strokeStyle;
    ctx.lineWidth = 2;

    ctx.beginPath();
    ctx.moveTo(X(ci.low), y);
    ctx.lineTo(X(ci.high), y);
    ctx.stroke();

    // 端点
    ctx.beginPath();
    ctx.moveTo(X(ci.low), y - 3);
    ctx.lineTo(X(ci.low), y + 3);
    ctx.moveTo(X(ci.high), y - 3);
    ctx.lineTo(X(ci.high), y + 3);
    ctx.stroke();

    // 样本均值
    ctx.beginPath();
    ctx.arc(X(ci.mean), y, 2.6, 0, 2 * Math.PI);
    ctx.fill();
  }

  // 真实均值 μ 的参考竖线
  ctx.strokeStyle = "#059669";
  ctx.lineWidth = 2;
  ctx.setLineDash([6, 4]);

  ctx.beginPath();
  ctx.moveTo(X(0), margin.top);
  ctx.lineTo(X(0), margin.top + plotHeight);
  ctx.stroke();

  ctx.setLineDash([]);
  ctx.restore();

  // 模拟编号
  ctx.fillStyle = "#64748b";
  ctx.font = "11px system-ui";
  ctx.textAlign = "right";

  for (let i = first; i < last; i++) {
    const y = Y(i);

    if (
      y >= margin.top &&
      y <= margin.top + plotHeight &&
      (i === 0 || (i + 1) % 5 === 0)
    ) {
      ctx.fillText(i + 1, margin.left - 10, y + 4);
    }
  }

  // μ 标签固定在窗口顶部
  ctx.fillStyle = "#059669";
  ctx.font = "bold 13px system-ui";
  ctx.textAlign = "center";
  ctx.fillText(`μ = ${fmt(mu)}`, X(0), 23);

  // 坐标轴
  ctx.strokeStyle = "#94a3b8";
  ctx.lineWidth = 1;

  ctx.beginPath();
  ctx.moveTo(margin.left, margin.top);
  ctx.lineTo(margin.left, margin.top + plotHeight);
  ctx.lineTo(width - margin.right, margin.top + plotHeight);
  ctx.stroke();

  // 横轴刻度：原始尺度
  ctx.fillStyle = "#64748b";
  ctx.font = "11px system-ui";
  ctx.textAlign = "center";

  for (let i = 0; i <= ticks; i++) {
    const value = min + (max - min) * i / ticks;

    ctx.fillText(
      fmt(mu + sigma * value),
      X(value),
      height - 27
    );
  }

  ctx.font = "12px system-ui";
  ctx.fillText(
    "均值的取值",
    margin.left + plotWidth / 2,
    height - 7
  );

  ctx.save();
  ctx.translate(15, margin.top + plotHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("模拟编号", 0, 0);
  ctx.restore();

  $("displayMode").textContent = followAll
    ? `全部展示模式：共 ${B.toLocaleString()} 条。修改 B 后自动跟随新的 B。`
    : `自定义展示模式：显示前 ${count.toLocaleString()} 条，可点击按钮展示全部。`;

  $("plotNote").textContent =
    `选中的 ${count.toLocaleString()} 条 CI 中：` +
    `${covered.toLocaleString()} 条包含 μ，` +
    `${(count - covered).toLocaleString()} 条未包含 μ。` +
    `上方统计使用全部 ${B.toLocaleString()} 次模拟。`;

  $("visibleNote").textContent =
    `当前窗口显示模拟编号约 ${first + 1}–${last}。`;
}

/* 执行重复抽样 */
async function runSimulation() {
  if (!$("controls").reportValidity()) return;

  const id = ++generation;

  const n = Number($("n").value);
  const B = Number($("B").value);
  const mu = Number($("mu").value);
  const sigma = Number($("sigma").value);
  const seed = Number($("seed").value);

  $("status").textContent = "正在模拟……旧图将在完成后更新。";

  const start = performance.now();
  const critical = tQuantile(0.975, n - 1);
  const normal = normalGenerator(mulberry32(seed));

  const intervals = new Array(B);
  let covered = 0;
  let meanLengthRelative = 0;

  // 每一块约生成不超过 10 万个观测值，保持页面响应
  const chunkSize = Math.max(1, Math.floor(100000 / n));

  for (let begin = 0; begin < B; begin += chunkSize) {
    await pause();

    // 参数已更新时，取消旧任务
    if (id !== generation) return;

    const end = Math.min(B, begin + chunkSize);

    for (let b = begin; b < end; b++) {
      let mean = 0;
      let m2 = 0;

      for (let i = 1; i <= n; i++) {
        const z = normal();

        // Welford 算法：均值与离差平方和
        const delta = z - mean;
        mean += delta / i;
        m2 += delta * (z - mean);
      }

      /*
        原始样本 X_i = μ + σ Z_i。

        X̄ = μ + σ mean
        S = σ sqrt(m2/(n−1))

        以下区间端点保存在 (x−μ)/σ 坐标中。
      */
      const relativeSD = Math.sqrt(m2 / (n - 1));
      const halfWidth = critical * relativeSD / Math.sqrt(n);

      const low = mean - halfWidth;
      const high = mean + halfWidth;
      const containsMu = low <= 0 && high >= 0;

      intervals[b] = {
        mean,
        low,
        high,
        covered: containsMu
      };

      if (containsMu) covered++;

      meanLengthRelative +=
        (2 * halfWidth - meanLengthRelative) / (b + 1);
    }

    $("status").textContent =
      `正在模拟：${end.toLocaleString()} / ${B.toLocaleString()}。`;
  }

  if (id !== generation) return;

  state = { intervals, n, B, mu, sigma };
  plotCache = null;

  $("critical").textContent =
    `自由度 ν = ${n - 1}；` +
    `t₀.₉₇₅,ν = ${critical.toFixed(5)}。` +
    `样本方差 S² 使用 n − 1 作为分母。`;

  $("coverage").textContent = (100 * covered / B).toFixed(2) + "%";
  $("covered").textContent = covered.toLocaleString();
  $("missed").textContent = (B - covered).toLocaleString();
  $("length").textContent = fmt(sigma * meanLengthRelative);

  const mcSE = Math.sqrt(0.95 * 0.05 / B);

  $("mcError").textContent =
    `在理论覆盖率 95% 下，经验覆盖率的蒙特卡洛标准误约为 ` +
    `${(100 * mcSE).toFixed(2)} 个百分点。` +
    `有限次模拟不必恰好得到 95% 的覆盖率。`;

  $("shown").max = B;

  if (followAll) {
    $("shown").value = B;
  } else {
    $("shown").value = Math.min(
      B,
      Math.max(1, Math.floor(Number($("shown").value) || 1))
    );
  }

  $("plotScroll").scrollTop = 0;
  drawPlot();

  $("status").textContent =
    `完成 ${B.toLocaleString()} 次模拟，每次样本量 n = ${n}。` +
    `耗时 ${Math.round(performance.now() - start)} ms。`;
}

/* 提交参数 */
$("controls").addEventListener("submit", event => {
  event.preventDefault();
  clearTimeout(inputTimer);
  runSimulation();
});

/* 参数修改后自动模拟，防抖 */
$("controls").addEventListener("input", () => {
  ++generation;
  clearTimeout(inputTimer);

  $("status").textContent = "参数已修改，等待更新……";

  inputTimer = setTimeout(() => {
    if ($("controls").checkValidity()) {
      runSimulation();
    } else {
      $("status").textContent = "请输入有效参数；图表仍显示上一次完成的模拟。";
    }
  }, 450);
});

/* 更换种子 */
$("newSeed").addEventListener("click", () => {
  clearTimeout(inputTimer);

  const seed = new Uint32Array(1);
  crypto.getRandomValues(seed);

  $("seed").value = seed[0];
  runSimulation();
});

/* 展示数量只影响绘图，不重新模拟 */
$("shown").addEventListener("input", () => {
  if (!state || !$("shown").checkValidity()) return;

  followAll = false;
  $("plotScroll").scrollTop = 0;
  drawPlot();
});

$("shown").addEventListener("change", () => {
  if (!state) return;

  followAll = false;

  const value = Number($("shown").value);

  $("shown").value = Math.min(
    state.B,
    Math.max(1, Number.isFinite(value) ? Math.floor(value) : 1)
  );

  $("plotScroll").scrollTop = 0;
  drawPlot();
});

/* 展示全部，并在后续模拟中自动跟随 B */
$("showAll").addEventListener("click", () => {
  if (!state) return;

  followAll = true;
  $("shown").value = state.B;
  $("plotScroll").scrollTop = 0;

  drawPlot();
});

$("backTop").addEventListener("click", () => {
  $("plotScroll").scrollTop = 0;
  drawPlot();
});

/* 滚动时仅重绘当前窗口 */
$("plotScroll").addEventListener("scroll", () => {
  if (scrollFrame !== null) return;

  scrollFrame = requestAnimationFrame(() => {
    scrollFrame = null;
    drawPlot();
  });
});

/* 自适应窗口尺寸 */
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(drawPlot, 100);
});

/* 首次打开页面自动模拟 */
runSimulation();
