const canvas = document.getElementById("canvas")
const canvas2 = document.getElementById("canvas2")
const ctx = canvas.getContext("2d")
const ctx2 = canvas2.getContext("2d")
const optimMethod = document.getElementById("optimMethod")
const learningRate = document.getElementById("learningRate")
const learningRateValue = document.getElementById("learningRateValue")

const momentum = document.getElementById("momentum")
const momentumValue = document.getElementById("momentumValue")
const momentumGroup = document.getElementById("momentumGroup")

const noise = document.getElementById("noise")
const noiseValue = document.getElementById("noiseValue")
const noiseGroup = document.getElementById("noiseGroup")

const betaInput = document.getElementById("beta");
const betaValue = document.getElementById("betaValue");
const betaGroup = document.getElementById("betaGroup");

const beta1Input = document.getElementById("beta1");
const beta1Value = document.getElementById("beta1Value");
const beta1Group = document.getElementById("beta1Group");

const beta2Input = document.getElementById("beta2");
const beta2Value = document.getElementById("beta2Value");
const beta2Group = document.getElementById("beta2Group");

const SIZE = 500
canvas.height = SIZE; canvas.width = SIZE
canvas2.height = SIZE; canvas2.width = SIZE

const width = SIZE, height = SIZE;

const to_px = x => (x + 2) * 0.25 * SIZE;
const to_real = px => (px / SIZE) * 4 - 2;

function f_default(x, y) { return x * x * x * x - 1 * x * x + 0.5 * x + y * y; }

let z = new Float32Array(width * height);

let zmin = +Infinity, zmax = -Infinity;
function computeGrid(f) {
    zmin = +Infinity; zmax = -Infinity;
    for (let j = 0; j < height; j++) {
        const y = to_real(j + 0.5);
        for (let i = 0; i < width; i++) {
            const x = to_real(i + 0.5);
            const v = f(x, y);
            z[j * width + i] = v;
            if (v < zmin) zmin = v;
            if (v > zmax) zmax = v;
        }
    }
}

function drawSurface(f) {
    computeGrid(f);
    const img = ctx.createImageData(width, height);
    const rng = (zmax - zmin) || 1e-12;
    for (let j = 0; j < height; j++) {
        for (let i = 0; i < width; i++) {
            const v = z[j * width + i];
            const t = (v - zmin) / rng;
            const gb = Math.round(255 * (1 - t));
            const idx = 4 * (j * width + i);
            img.data[idx + 0] = 255;
            img.data[idx + 1] = gb * 1.2;
            img.data[idx + 2] = gb * 1.2;
            img.data[idx + 3] = 255;
        }
    }
    ctx.putImageData(img, 0, 0);
}

function drawLevelCurves(f, levels) {
    ctx.save();
    ctx.lineWidth = 0.5;
    ctx.strokeStyle = "red";

    const lerp = (a, b, t) => a + (b - a) * t;
    const eps = 1e-12;

    for (const level of levels) {
        for (let j = 0; j < height - 1; j++) {
            for (let i = 0; i < width - 1; i++) {
                const v0 = z[j * width + i];
                const v1 = z[j * width + (i + 1)];
                const v2 = z[(j + 1) * width + (i + 1)];
                const v3 = z[(j + 1) * width + i];

                const b0 = v0 > level ? 1 : 0;
                const b1 = v1 > level ? 2 : 0;
                const b2 = v2 > level ? 4 : 0;
                const b3 = v3 > level ? 8 : 0;
                const bits = b0 | b1 | b2 | b3;
                if (bits === 0 || bits === 15) continue;

                const c0 = { x: i + 0.5, y: j + 0.5 };
                const c1 = { x: i + 1.5, y: j + 0.5 };
                const c2 = { x: i + 1.5, y: j + 1.5 };
                const c3 = { x: i + 0.5, y: j + 1.5 };

                const edgeInterp = (A, B, vA, vB) => {
                    const t = (level - vA) / (vB - vA + eps);
                    return { x: lerp(A.x, B.x, t), y: lerp(A.y, B.y, t) };
                };

                const segsByCase = {
                    1: [[3, 0]],
                    2: [[0, 1]],
                    3: [[3, 1]],
                    4: [[1, 2]],
                    5: [[3, 0], [1, 2]],
                    6: [[0, 2]],
                    7: [[3, 2]],
                    8: [[2, 3]],
                    9: [[2, 0]],
                    10: [[0, 1], [2, 3]],
                    11: [[1, 3]],
                    12: [[1, 3]],
                    13: [[0, 1]],
                    14: [[0, 3]]
                };

                const segs = segsByCase[bits];
                if (!segs) continue;

                for (const seg of segs) {
                    const eA = seg[0], eB = seg[1];
                    const P = (e) => {
                        if (e === 0) return edgeInterp(c0, c1, v0, v1);
                        if (e === 1) return edgeInterp(c1, c2, v1, v2);
                        if (e === 2) return edgeInterp(c2, c3, v2, v3);
                        return edgeInterp(c3, c0, v3, v0);
                    };
                    const pA = P(eA), pB = P(eB);
                    ctx.beginPath();
                    ctx.moveTo(pA.x, pA.y);
                    ctx.lineTo(pB.x, pB.y);
                    ctx.stroke();
                }
            }
        }
    }
    ctx.restore();
}

const to_px_z = z => (2 - z) * 0.25 * SIZE;

function drawAxes() {
    ctx2.clearRect(0, 0, SIZE, SIZE);
    ctx2.save();
    ctx2.strokeStyle = "black";
    ctx2.fillStyle = "black";
    ctx2.lineWidth = 0.5;
    ctx2.font = "12px monospace";

    const x0 = to_px(0);
    const y0 = to_px_z(0);

    ctx2.beginPath();
    ctx2.moveTo(0, y0); ctx2.lineTo(SIZE, y0);
    ctx2.moveTo(x0, 0); ctx2.lineTo(x0, SIZE);
    ctx2.stroke();

    const ticks = [-2, -1, 1, 2];
    for (const t of ticks) {
        const tx = to_px(t);
        const tz = to_px_z(t);

        ctx2.beginPath();
        ctx2.moveTo(tx, y0 - 6); ctx2.lineTo(tx, y0 + 6);
        ctx2.stroke();
        ctx2.textAlign = "center"; ctx2.textBaseline = "top";
        ctx2.fillText(String(t), tx, y0 + 8);

        ctx2.beginPath();
        ctx2.moveTo(x0 - 6, tz); ctx2.lineTo(x0 + 6, tz);
        ctx2.stroke();
        ctx2.textAlign = "left"; ctx2.textBaseline = "middle";
        ctx2.fillText(String(t), x0 + 8, tz);
    }

    ctx2.textAlign = "right"; ctx2.textBaseline = "top";
    ctx2.fillText("x", SIZE - 6, y0 + 6);
    ctx2.textAlign = "left"; ctx2.textBaseline = "top";
    ctx2.fillText("z", x0 + 6, 4);

    ctx2.restore();
}

function drawCurve(f) {
    ctx2.save();
    ctx2.strokeStyle = "red";
    ctx2.lineWidth = 1.5;
    ctx2.beginPath();

    for (let i = 0; i < width; i++) {
        const x = to_real(i + 0.5);
        const zv = f(x, 0);
        const px = to_px(x);
        const py = to_px_z(zv);

        if (i === 0) ctx2.moveTo(px, py);
        else ctx2.lineTo(px, py);
    }
    ctx2.stroke();
    ctx2.restore();
}

function drawMarkerOnSurface(pt) {
    const px = to_px(pt.x);
    const py = to_px(pt.y);
    ctx.save();
    ctx.fillStyle = "white";
    ctx.strokeStyle = "black";
    ctx.beginPath();
    ctx.arc(px, py, 4, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke()
    ctx.restore();
}

function drawMarkerOnPlot(pt, f) {
    const x = pt.x;
    const zval = f(x, 0);
    const px = to_px(x);
    const py = to_px_z(zval);
    ctx2.save();
    ctx2.fillStyle = "black";
    ctx2.beginPath();
    ctx2.arc(px, py, 4, 0, Math.PI * 2);
    ctx2.fill();

    ctx2.fillStyle = "black";
    ctx2.font = "11px monospace";
    ctx2.textAlign = "left";
    ctx2.textBaseline = "top";
    const label = `x=${pt.x.toFixed(2)}, z=${zval.toFixed(2)}`;
    ctx2.fillText(label, px + 8, py + 6);
    ctx2.restore();
}

let surfaceImage = null;
let rafPending = false;
let isDragging = false;
let draggable = { x: 1.55, y: -1 };
const dragRadiusPx = 8;

function renderStatic(f) {
    drawSurface(f);

    const nLevels = 20;
    const levels = [];
    for (let k = 1; k < nLevels; k++) levels.push(zmin + (zmax - zmin) * k / nLevels);
    drawLevelCurves(f, levels);

    surfaceImage = ctx.getImageData(0, 0, SIZE, SIZE);
}

function updateDuringDrag(f) {
    if (surfaceImage) ctx.putImageData(surfaceImage, 0, 0);
    drawMarkerOnSurface(draggable);

    drawAxes();
    drawCurve(f);
    drawMarkerOnPlot(draggable, f);
    renderCanvases(f);
    rafPending = false;
}

function renderCanvases(f = f_default) {
    renderStatic(f);

    const alpha = parseFloat(learningRate.value);
    const beta = parseFloat(momentum.value);
    const beta1 = parseFloat(beta1Input.value);
    const beta2 = parseFloat(beta2Input.value);
    const noiseLevel = parseFloat(noise.value);
    const maxIter = 200;
    const tol = 1e-4;

    let path;
    if (optimMethod.value === "sgd") {
        path = SGDincrement(f, { alpha, noiseLevel, maxIter, tol });
    } else if (optimMethod.value === "momgd") {
        path = MomGDincrement(f, { alpha, beta, maxIter, tol });
    } else if (optimMethod.value === "nesterov") {
        path = NesterovGDincrement(f, { alpha, beta, maxIter, tol });
    } else if (optimMethod.value === "rmsprop") {
        let beta = parseFloat(betaInput.value)
        path = RMSpropIncrement(f, { alpha, beta, maxIter, tol });
    } else if (optimMethod.value === "adam") {
        path = AdamIncrement(f, { alpha, beta1, beta2, maxIter, tol });
    } else if (optimMethod.value === "adagrad") {
        path = AdagradIncrement(f, { alpha, maxIter, tol });
    }
    else {
        path = GDincrement(f, { alpha, maxIter, tol });
    }

    drawIncrements(path, f);

    drawMarkerOnSurface(draggable);
    drawMarkerOnPlot(draggable, f);
}

function getMousePosOnCanvas(evt) {
    const rect = canvas.getBoundingClientRect();
    const cx = evt.clientX - rect.left;
    const cy = evt.clientY - rect.top;
    return {
        px: Math.max(0, Math.min(SIZE, cx)),
        py: Math.max(0, Math.min(SIZE, cy))
    };
}

canvas.addEventListener("pointerdown", (e) => {
    canvas.setPointerCapture(e.pointerId);
    const { px, py } = getMousePosOnCanvas(e);
    const dx = px - to_px(draggable.x);
    const dy = py - to_px(draggable.y);
    if (dx * dx + dy * dy <= dragRadiusPx * dragRadiusPx) {
        isDragging = true;
    } else {
        draggable.x = Math.max(-2, Math.min(2, to_real(px)));
        draggable.y = Math.max(-2, Math.min(2, to_real(py)));
        if (surfaceImage) {
            ctx.putImageData(surfaceImage, 0, 0);
            drawMarkerOnSurface(draggable);
        } else {
            renderCanvases(f_default);
        }
        drawAxes(); drawCurve(f_default); drawMarkerOnPlot(draggable, f_default);
    }
});

canvas.addEventListener("pointermove", (e) => {
    if (!isDragging) return;
    const { px, py } = getMousePosOnCanvas(e);
    draggable.x = Math.max(-2, Math.min(2, to_real(px)));
    draggable.y = Math.max(-2, Math.min(2, to_real(py)));
    if (!rafPending) {
        rafPending = true;
        requestAnimationFrame(() => updateDuringDrag(f_default));
    }
});

canvas.addEventListener("pointerup", (e) => {
    isDragging = false;
    try { canvas.releasePointerCapture(e.pointerId); } catch (_) { }
    if (!rafPending) updateDuringDrag(f_default);
});

function GDincrement(f, opts = {}) {
    const alpha = opts.alpha ?? 0.05;
    const maxIter = opts.maxIter ?? 200;
    const tol = opts.tol ?? 1e-4;
    const epsGrad = opts.epsGrad ?? 1e-6;

    let x = draggable.x;
    let y = draggable.y;
    const path = [{ x, y, z: f(x, y) }];

    for (let it = 0; it < maxIter; it++) {
        const h = epsGrad;
        const dfdx = (f(x + h, y) - f(x - h, y)) / (2 * h);
        const dfdy = (f(x, y + h) - f(x, y - h)) / (2 * h);

        const step_x = -alpha * dfdx;
        const step_y = -alpha * dfdy;
        const stepNorm = Math.hypot(step_x, step_y);

        x = Math.max(-2, Math.min(2, x + step_x));
        y = Math.max(-2, Math.min(2, y + step_y));

        path.push({ x, y, z: f(x, y) });

        if (stepNorm < tol) break;
    }

    return path;
}

function SGDincrement(f, opts = {}) {
    const alpha = opts.alpha ?? 0.05;
    const maxIter = opts.maxIter ?? 200;
    const tol = opts.tol ?? 1e-4;
    const epsGrad = opts.epsGrad ?? 1e-6;
    const noiseLevel = opts.noiseLevel ?? 0.4;

    let x = draggable.x;
    let y = draggable.y;
    const path = [{ x, y, z: f(x, y) }];

    for (let it = 0; it < maxIter; it++) {
        const h = epsGrad;
        const dfdx = (f(x + h, y) - f(x - h, y)) / (2 * h);
        const dfdy = (f(x, y + h) - f(x, y - h)) / (2 * h);

        const noisy_dfdx = dfdx + (Math.random() - 0.5) * noiseLevel;
        const noisy_dfdy = dfdy + (Math.random() - 0.5) * noiseLevel;

        const step_x = -alpha * noisy_dfdx;
        const step_y = -alpha * noisy_dfdy;
        const stepNorm = Math.hypot(step_x, step_y);

        x = Math.max(-2, Math.min(2, x + step_x));
        y = Math.max(-2, Math.min(2, y + step_y));

        path.push({ x, y, z: f(x, y) });

        if (stepNorm < tol) break;
    }

    return path;
}

function MomGDincrement(f, opts = {}) {
    const alpha = opts.alpha ?? 0.05;
    const beta = opts.beta ?? 0.9;
    const maxIter = opts.maxIter ?? 200;
    const tol = opts.tol ?? 1e-4;
    const epsGrad = opts.epsGrad ?? 1e-6;

    let x = draggable.x;
    let y = draggable.y;
    let vx = 0, vy = 0;
    const path = [{ x, y, z: f(x, y) }];

    for (let it = 0; it < maxIter; it++) {
        const h = epsGrad;
        const dfdx = (f(x + h, y) - f(x - h, y)) / (2 * h);
        const dfdy = (f(x, y + h) - f(x, y - h)) / (2 * h);

        vx = beta * vx - alpha * dfdx;
        vy = beta * vy - alpha * dfdy;

        x = Math.max(-2, Math.min(2, x + vx));
        y = Math.max(-2, Math.min(2, y + vy));

        path.push({ x, y, z: f(x, y) });

        if (Math.hypot(vx, vy) < tol) break;
    }
    return path;
}

function NesterovGDincrement(f, opts = {}) {
    const alpha = opts.alpha ?? 0.05;
    const beta = opts.beta ?? 0.9;
    const maxIter = opts.maxIter ?? 200;
    const tol = opts.tol ?? 1e-4;
    const epsGrad = opts.epsGrad ?? 1e-6;

    let x = draggable.x;
    let y = draggable.y;
    let vx = 0, vy = 0;
    const path = [{ x, y, z: f(x, y) }];

    for (let it = 0; it < maxIter; it++) {
        const h = epsGrad;

        const lookahead_x = x + beta * vx;
        const lookahead_y = y + beta * vy;

        const dfdx = (f(lookahead_x + h, lookahead_y) - f(lookahead_x - h, lookahead_y)) / (2 * h);
        const dfdy = (f(lookahead_x, lookahead_y + h) - f(lookahead_x, lookahead_y - h)) / (2 * h);

        vx = beta * vx - alpha * dfdx;
        vy = beta * vy - alpha * dfdy;

        x = Math.max(-2, Math.min(2, x + vx));
        y = Math.max(-2, Math.min(2, y + vy));

        path.push({ x, y, z: f(x, y) });

        if (Math.hypot(vx, vy) < tol) break;
    }

    return path;
}

function RMSpropIncrement(f, opts = {}) {
    const alpha = opts.alpha ?? 0.05;
    const beta = opts.beta ?? 0.9;
    const epsilon = opts.epsilon ?? 1e-8;
    const maxIter = opts.maxIter ?? 200;
    const tol = opts.tol ?? 1e-4;
    const epsGrad = opts.epsGrad ?? 1e-6;

    let x = draggable.x;
    let y = draggable.y;
    let s_x = 0, s_y = 0;
    const path = [{ x, y, z: f(x, y) }];

    for (let it = 0; it < maxIter; it++) {
        const h = epsGrad;
        const dfdx = (f(x + h, y) - f(x - h, y)) / (2 * h);
        const dfdy = (f(x, y + h) - f(x, y - h)) / (2 * h);

        s_x = beta * s_x + (1 - beta) * dfdx * dfdx;
        s_y = beta * s_y + (1 - beta) * dfdy * dfdy;

        const step_x = -alpha * dfdx / (Math.sqrt(s_x) + epsilon);
        const step_y = -alpha * dfdy / (Math.sqrt(s_y) + epsilon);
        const stepNorm = Math.hypot(step_x, step_y);

        x = Math.max(-2, Math.min(2, x + step_x));
        y = Math.max(-2, Math.min(2, y + step_y));

        path.push({ x, y, z: f(x, y) });

        if (stepNorm < tol) break;
    }

    return path;
}

function AdamIncrement(f, opts = {}) {
    const alpha = opts.alpha ?? 0.05;
    const beta1 = opts.beta1 ?? 0.9;
    const beta2 = opts.beta2 ?? 0.999;
    const epsilon = opts.epsilon ?? 1e-8;
    const maxIter = opts.maxIter ?? 200;
    const tol = opts.tol ?? 1e-4;
    const epsGrad = opts.epsGrad ?? 1e-6;

    let x = draggable.x;
    let y = draggable.y;
    let m_x = 0, m_y = 0;
    let v_x = 0, v_y = 0;
    const path = [{ x, y, z: f(x, y) }];

    for (let it = 1; it <= maxIter; it++) {
        const h = epsGrad;
        const dfdx = (f(x + h, y) - f(x - h, y)) / (2 * h);
        const dfdy = (f(x, y + h) - f(x, y - h)) / (2 * h);

        m_x = beta1 * m_x + (1 - beta1) * dfdx;
        m_y = beta1 * m_y + (1 - beta1) * dfdy;

        v_x = beta2 * v_x + (1 - beta2) * dfdx * dfdx;
        v_y = beta2 * v_y + (1 - beta2) * dfdy * dfdy;

        const m_hat_x = m_x / (1 - Math.pow(beta1, it));
        const m_hat_y = m_y / (1 - Math.pow(beta1, it));

        const v_hat_x = v_x / (1 - Math.pow(beta2, it));
        const v_hat_y = v_y / (1 - Math.pow(beta2, it));

        const step_x = -alpha * m_hat_x / (Math.sqrt(v_hat_x) + epsilon);
        const step_y = -alpha * m_hat_y / (Math.sqrt(v_hat_y) + epsilon);
        const stepNorm = Math.hypot(step_x, step_y);

        x = Math.max(-2, Math.min(2, x + step_x));
        y = Math.max(-2, Math.min(2, y + step_y));

        path.push({ x, y, z: f(x, y) });

        if (stepNorm < tol) break;
    }

    return path;
}

function AdagradIncrement(f, opts = {}) {
    const alpha = opts.alpha ?? 0.05;
    const epsilon = opts.epsilon ?? 1e-8;
    const maxIter = opts.maxIter ?? 200;
    const tol = opts.tol ?? 1e-4;
    const epsGrad = opts.epsGrad ?? 1e-6;

    let x = draggable.x;
    let y = draggable.y;
    let s_x = 0, s_y = 0;
    const path = [{ x, y, z: f(x, y) }];

    for (let it = 0; it < maxIter; it++) {
        const h = epsGrad;
        const dfdx = (f(x + h, y) - f(x - h, y)) / (2 * h);
        const dfdy = (f(x, y + h) - f(x, y - h)) / (2 * h);

        s_x += dfdx * dfdx;
        s_y += dfdy * dfdy;

        const step_x = -alpha * dfdx / (Math.sqrt(s_x) + epsilon);
        const step_y = -alpha * dfdy / (Math.sqrt(s_y) + epsilon);
        const stepNorm = Math.hypot(step_x, step_y);

        x = Math.max(-2, Math.min(2, x + step_x));
        y = Math.max(-2, Math.min(2, y + step_y));

        path.push({ x, y, z: f(x, y) });

        if (stepNorm < tol) break;
    }

    return path;
}

function drawIncrements(path, f, opts = {}) {
    if (!path || path.length === 0) return;

    if (!surfaceImage) renderStatic(f);
    ctx.putImageData(surfaceImage, 0, 0);

    ctx.save();
    ctx.strokeStyle = "black";
    ctx.lineWidth = 1;

    ctx.beginPath();
    for (let i = 0; i < path.length; i++) {
        const p = path[i];
        const px = to_px(p.x);
        const py = to_px(p.y);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
    }
    ctx.stroke();

    for (let i = 0; i < path.length; i++) {
        const p = path[i];
        const px = to_px(p.x);
        const py = to_px(p.y);

        if (i < path.length - 1) {
            ctx.beginPath();
            ctx.arc(px, py, 1, 0, Math.PI * 2);
            ctx.fillStyle = "black";
            ctx.fill();
        } else {
            ctx.beginPath();
            ctx.moveTo(px - 3, py - 3);
            ctx.lineTo(px + 3, py + 3);
            ctx.moveTo(px - 3, py + 3);
            ctx.lineTo(px + 3, py - 3);
            ctx.stroke();
        }
    }
    ctx.restore();

    drawAxes();
    drawCurve(f);

    ctx2.save();
    ctx2.strokeStyle = "black";
    ctx2.lineWidth = 1;

    ctx2.beginPath();
    for (let i = 0; i < path.length; i++) {
        const p = path[i];
        const px = to_px(p.x);
        const pz = to_px_z(f(p.x, 0));
        if (i === 0) ctx2.moveTo(px, pz);
        else ctx2.lineTo(px, pz);
    }
    ctx2.stroke();

    for (let i = 0; i < path.length; i++) {
        const p = path[i];
        const px = to_px(p.x);
        const pz = to_px_z(f(p.x, 0));

        if (i < path.length - 1) {
            ctx2.beginPath();
            ctx2.arc(px, pz, 1.5, 0, Math.PI * 2);
            ctx2.fillStyle = "black";
            ctx2.fill();
        } else {
            ctx2.beginPath();
            ctx2.moveTo(px - 3, pz - 3);
            ctx2.lineTo(px + 3, pz + 3);
            ctx2.moveTo(px - 3, pz + 3);
            ctx2.lineTo(px + 3, pz - 3);
            ctx2.stroke();
        }
    }
    ctx2.restore();
}

optimMethod.addEventListener("change", function () {
    const method = this.value;
    momentumGroup.style.display = (method === "momgd" || method === "nesterov") ? "flex" : "none";
    noiseGroup.style.display = (method === "sgd") ? "flex" : "none";
    renderCanvases(f_default);
    betaGroup.style.display = (method === "rmsprop") ? "flex" : "none";
    beta1Group.style.display = (method === "adam") ? "flex" : "none";
    beta2Group.style.display = (method === "adam") ? "flex" : "none";
});

learningRate.addEventListener("input", function () {
    learningRateValue.textContent = this.value;
    renderCanvases(f_default);
});

momentum.addEventListener("input", function () {
    momentumValue.textContent = this.value;
    if (optimMethod.value === "momgd") {
        renderCanvases(f_default);
    }
});

noise.addEventListener("input", function () {
    noiseValue.textContent = this.value;
    if (optimMethod.value === "sgd") {
        renderCanvases(f_default);
    }
});

betaInput.addEventListener("input", function () {
    betaValue.textContent = this.value;
    if (optimMethod.value === "rmsprop") {
        renderCanvases(f_default);
    }
});

beta1Input.addEventListener("input", function () {
    beta1Value.textContent = this.value;
    if (optimMethod.value === "adam") {
        renderCanvases(f_default);
    }
});

beta2Input.addEventListener("input", function () {
    beta2Value.textContent = this.value;
    if (optimMethod.value === "adam") {
        renderCanvases(f_default);
    }
});

renderCanvases(f_default);