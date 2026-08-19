/* Exo-Daisy World — interactive applet for the pyEDW documentation.
 *
 * The physics is a direct port of pyEDW/model.py: the same dimensionless
 * 4-dof state x = [f_B, f_W, T, L], the same equations of motion (Eqs. B8-B12),
 * and the same strong-order-1 stochastic Runge-Kutta step (Roberts 2012) that
 * updateExoDaisyWorld.m uses. Drag the planet to spin it; move the luminosity
 * slider and watch the daisies rein the surface temperature back toward T_opt.
 */
(function () {
    'use strict';

    /* ── Parameters (Table 1) ───────────────────────────────────────────── */
    const P = {
        f: 0.88, A_G: 0.3, A_B: 0.1, A_W: 0.6,
        T_opt: 300, Q: 0.1, gG: 1.0, gD: 0.2,
        tauS: 3.0, tauE: 5.0, dT: 30.0, delta: 0.05, lam: 0.0
    };

    const DT = 0.1;          // integration step, as in the MATLAB
    const SUBSTEPS = 8;      // steps per animation frame
    const N_SITES = 3000;    // daisy sites scattered over the sphere

    /* Seed bank. The growth terms are multiplicative in f_B and f_W, so zero is
     * an absorbing state: a biome that collapses can never return, and the
     * luminosity slider would go dead. The paper's ensembles sidestep this by
     * re-seeding every run. Here we instead hold a dormant spore population, as
     * Watson & Lovelock's original model does, so the planet can be
     * recolonized when it becomes habitable again. At 1e-4 the bank is far too
     * small to draw (it rounds to zero flowers) and does not perturb the
     * dynamics of a living biome. */
    const SEED = 1e-4;

    /* theta layout mirrors model.py exactly. */
    function theta() {
        return [
            P.f,                       // 0  f
            1 - P.A_G,                 // 1  1 - A_G
            P.A_B - P.A_G,             // 2  dA_B
            P.A_W - P.A_G,             // 3  dA_W
            P.Q,                       // 4  Q
            P.gD / P.gG,               // 5  gamma_D / gamma_G
            1 / (P.gG * P.tauE),       // 6  1 / (gamma_G tau_E)
            1 / (P.gG * P.tauS),       // 7  1 / (gamma_G tau_S)
            8 * Math.pow(P.T_opt / P.dT, 4),  // 8  alpha
            P.delta,                   // 9  delta
            P.lam                      // 10 lambda
        ];
    }

    /* Growth window w(T_a - 1) = exp(-alpha (T_a - 1)^4), where the scaled
     * daisy temperature T_a = base^(1/4). When base < 0 the fourth root is
     * complex; MATLAB carries it and takes the real part of the increment, so
     * we return Re(w) from the principal branch. The window is clamped to
     * [0, 1] — a growth rate cannot exceed gamma_G, and the complex branch can
     * otherwise blow up in a live simulation. */
    function growthRe(base, alpha) {
        let w;
        if (base >= 0) {
            const t = Math.pow(base, 0.25) - 1;
            w = Math.exp(-alpha * t * t * t * t);
        } else {
            const r = Math.pow(-base, 0.25), c = Math.SQRT1_2;
            const zr = r * c - 1, zi = r * c;             // z = T_a - 1
            const z2r = zr * zr - zi * zi, z2i = 2 * zr * zi;
            const z4r = z2r * z2r - z2i * z2i, z4i = 2 * z2r * z2i;
            w = Math.exp(-alpha * z4r) * Math.cos(-alpha * z4i);
        }
        if (!isFinite(w) || w < 0) return 0;
        return w > 1 ? 1 : w;
    }

    /* Drift of the four dof (ports the EoM subfunction). */
    function eom(x, th) {
        const fB = x[0], fW = x[1], T = x[2], L = x[3];
        const dAf = th[2] * fB + th[3] * fW;
        const Tg4 = T * T * T * T + th[4] * dAf;
        const df = th[0] - fB - fW;
        const wB = growthRe(Tg4 - th[4] * th[2], th[8]);
        const wW = growthRe(Tg4 - th[4] * th[3], th[8]);
        return [
            df * wB * fB - th[5] * fB,
            df * wW * fW - th[5] * fW,
            th[6] * ((1 - dAf / th[1]) * L - T * T * T * T),
            th[7] * (1 + th[10] - L)
        ];
    }

    let spare = null;
    function randn() {                        // Box-Muller
        if (spare !== null) { const s = spare; spare = null; return s; }
        let u = 0, v = 0, s = 0;
        do {
            u = Math.random() * 2 - 1;
            v = Math.random() * 2 - 1;
            s = u * u + v * v;
        } while (s >= 1 || s === 0);
        const m = Math.sqrt(-2 * Math.log(s) / s);
        spare = v * m;
        return u * m;
    }

    /* One SRK1 step. Noise enters only through L (additive), so the driving
     * Wiener process is scalar and the strong-order-1 scheme applies. */
    function srk1(x, th) {
        const bL = Math.sqrt(2 * th[7]) * th[9] * (1 + th[10]);
        const z = randn(), s = Math.random() < 0.5 ? -1 : 1;
        const sq = Math.sqrt(DT);

        const a = eom(x, th);
        const k1 = [a[0] * DT, a[1] * DT, a[2] * DT, a[3] * DT + (z - s) * bL * sq];
        const x1 = [x[0] + k1[0], x[1] + k1[1], x[2] + k1[2], x[3] + k1[3]];
        const b = eom(x1, th);
        const k2 = [b[0] * DT, b[1] * DT, b[2] * DT, b[3] * DT + (z + s) * bL * sq];

        const out = [
            x[0] + 0.5 * (k1[0] + k2[0]),
            x[1] + 0.5 * (k1[1] + k2[1]),
            x[2] + 0.5 * (k1[2] + k2[2]),
            x[3] + 0.5 * (k1[3] + k2[3])
        ];
        if (out[0] < SEED) out[0] = SEED;     // non-negative, with a seed bank
        if (out[1] < SEED) out[1] = SEED;
        if (!(out[2] > 0.05)) out[2] = 0.05;  // keep T, L physical for the view
        if (!(out[3] > 0.02)) out[3] = 0.02;
        return out;
    }

    /* ── State ──────────────────────────────────────────────────────────── */
    let X = null, running = true, frameId = null;
    const trail = [];                          // (L, T) history for the phase plot

    function seed() {
        const L0 = 1 + P.lam;
        X = [0.1 * P.f, 0.1 * P.f, Math.pow(L0, 0.25), L0];
        trail.length = 0;
    }

    function albedo() { return P.A_G + (P.A_B - P.A_G) * X[0] + (P.A_W - P.A_G) * X[1]; }

    /* ── Three.js scene ─────────────────────────────────────────────────── */
    let renderer, scene, camera, planet, blackMesh, whiteMesh, starLight, starGlow;
    let STAR_POS;
    let sites = [], perm = [];
    let dragging = false, lastX = 0, lastY = 0, spinY = 0, spinX = 0.25;

    /* Fibonacci sphere: N_SITES near-uniform points on the unit sphere. */
    function buildSites() {
        sites = [];
        const ga = Math.PI * (3 - Math.sqrt(5));
        for (let i = 0; i < N_SITES; i++) {
            const y = 1 - (i / (N_SITES - 1)) * 2;
            const r = Math.sqrt(Math.max(0, 1 - y * y));
            const th = ga * i;
            sites.push(new THREE.Vector3(Math.cos(th) * r, y, Math.sin(th) * r));
        }
        // A fixed random ranking. Black daisies fill sites from the front of
        // this permutation, white from the back, so populations grow and
        // shrink smoothly instead of flickering to new sites each frame.
        perm = Array.from({ length: N_SITES }, (_, i) => i);
        for (let i = N_SITES - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [perm[i], perm[j]] = [perm[j], perm[i]];
        }
    }

    /* A daisy: a golden centre disc ringed by eight petals, lying in the local
     * XZ plane with +Y as its outward axis. Built as one indexed geometry with
     * baked vertex colours, so black and white daisies differ only in their
     * petal colour and each still draws as a single InstancedMesh. */
    function makeFlower(petalHex) {
        const NP = 8;          // petals
        const NC = 10;         // centre-disc segments
        const cr = 0.0060;     // centre radius
        const pl = 0.0135;     // petal length
        const pw = 0.0042;     // petal half-width
        const ph = 0.0022;     // petal lift, so petals cup upward

        const pos = [], col = [], idx = [];
        const petal = new THREE.Color(petalHex);
        const heart = new THREE.Color(0xe8b93b);

        // Centre: a triangle fan.
        pos.push(0, 0.0016, 0);
        col.push(heart.r, heart.g, heart.b);
        for (let i = 0; i < NC; i++) {
            const a = i / NC * Math.PI * 2;
            pos.push(Math.cos(a) * cr, 0.0008, Math.sin(a) * cr);
            col.push(heart.r, heart.g, heart.b);
        }
        for (let i = 0; i < NC; i++) idx.push(0, 1 + i, 1 + (i + 1) % NC);

        // Petals: a lifted diamond quad each.
        for (let p = 0; p < NP; p++) {
            const a = p / NP * Math.PI * 2;
            const ca = Math.cos(a), sa = Math.sin(a);
            const rm = cr + pl * 0.42;
            const b = pos.length / 3;
            pos.push(ca * cr, 0.0009, sa * cr);                     // base
            pos.push(ca * rm - sa * pw, ph, sa * rm + ca * pw);      // side
            pos.push(ca * (cr + pl), 0.0006, sa * (cr + pl));        // tip
            pos.push(ca * rm + sa * pw, ph, sa * rm - ca * pw);      // side
            for (let k = 0; k < 4; k++) col.push(petal.r, petal.g, petal.b);
            idx.push(b, b + 1, b + 2, b, b + 2, b + 3);
        }

        const g = new THREE.BufferGeometry();
        g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
        g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
        g.setIndex(idx);
        g.computeVertexNormals();
        return g;
    }

    function initThree(canvas) {
        renderer = new THREE.WebGLRenderer({ canvas: canvas, antialias: true });
        renderer.setPixelRatio(window.devicePixelRatio);
        renderer.setClearColor(new THREE.Color(0x000000), 1);

        scene = new THREE.Scene();
        camera = new THREE.PerspectiveCamera(42, 1, 0.1, 100);
        camera.position.set(0, 0, 2.9);   // dollied in: the planet fills the frame

        scene.add(new THREE.AmbientLight(0xffffff, 0.18));

        // The host star sits in the background, upper left, and is the key
        // light. Anything visible in frame necessarily lies behind the planet's
        // plane, so the star backlights it; a gentle fill from the camera side
        // keeps the daisies on the near face legible.
        STAR_POS = new THREE.Vector3(-4.5, 2.4, -9);
        starLight = new THREE.DirectionalLight(0xfff0d0, 1.0);
        starLight.position.copy(STAR_POS);
        scene.add(starLight);

        const fill = new THREE.DirectionalLight(0xdfe8ff, 0.45);
        fill.position.set(2, 1, 6);
        scene.add(fill);

        const cv = document.createElement('canvas');
        cv.width = cv.height = 128;
        const c2 = cv.getContext('2d');
        const g = c2.createRadialGradient(64, 64, 0, 64, 64, 64);
        g.addColorStop(0.00, 'rgba(255,253,245,1)');
        g.addColorStop(0.12, 'rgba(255,246,214,0.95)');
        g.addColorStop(0.30, 'rgba(255,214,130,0.55)');
        g.addColorStop(0.60, 'rgba(255,178,80,0.16)');
        g.addColorStop(1.00, 'rgba(255,160,60,0)');
        c2.fillStyle = g;
        c2.fillRect(0, 0, 128, 128);
        starGlow = new THREE.Sprite(new THREE.SpriteMaterial({
            map: new THREE.CanvasTexture(cv), transparent: true,
            depthWrite: false, blending: THREE.AdditiveBlending
        }));
        starGlow.position.copy(STAR_POS);
        scene.add(starGlow);

        // Offset right and down, so the planet sits low in the frame and leaves
        // the upper left to the star.
        planet = new THREE.Group();
        planet.position.set(0.52, -0.38, 0);
        scene.add(planet);

        const ground = new THREE.Mesh(
            new THREE.SphereGeometry(1, 64, 48),
            new THREE.MeshPhongMaterial({ color: 0x9a8b73, shininess: 4 })
        );
        planet.add(ground);

        // Two instanced meshes of flower heads. Instance i of the black mesh
        // sits at site perm[i]; instance i of the white mesh at perm[N-1-i].
        // Matrices are written once; each frame we only set .count.
        const flowerMat = () => new THREE.MeshPhongMaterial({
            vertexColors: true, shininess: 22, side: THREE.DoubleSide
        });
        blackMesh = new THREE.InstancedMesh(makeFlower(0x15150f), flowerMat(), N_SITES);
        whiteMesh = new THREE.InstancedMesh(makeFlower(0xf6f4ec), flowerMat(), N_SITES);

        // Each flower is stood up along its site's outward normal, spun by a
        // random angle about that normal, and given a little size variation.
        const m = new THREE.Matrix4();
        const up = new THREE.Vector3(0, 1, 0);
        const qA = new THREE.Quaternion(), qS = new THREE.Quaternion();
        const sc = new THREE.Vector3();

        function place(mesh, i, siteIdx) {
            const n = sites[siteIdx];
            qA.setFromUnitVectors(up, n);                             // align +Y to normal
            qS.setFromAxisAngle(up, Math.random() * Math.PI * 2);     // spin about it
            qA.multiply(qS);
            const s = 0.8 + Math.random() * 0.5;
            sc.set(s, s, s);
            m.compose(n.clone().multiplyScalar(1.004), qA, sc);
            mesh.setMatrixAt(i, m);
        }

        for (let i = 0; i < N_SITES; i++) {
            place(blackMesh, i, perm[i]);
            place(whiteMesh, i, perm[N_SITES - 1 - i]);
        }
        blackMesh.instanceMatrix.needsUpdate = true;
        whiteMesh.instanceMatrix.needsUpdate = true;
        blackMesh.count = 0;
        whiteMesh.count = 0;
        planet.add(blackMesh);
        planet.add(whiteMesh);
    }

    function resize(canvas) {
        const w = canvas.clientWidth, h = canvas.clientHeight;
        if (!w || !h) return;
        renderer.setSize(w, h, false);
        camera.aspect = w / h;
        camera.updateProjectionMatrix();
    }

    /* ── Phase plot: T against L, with the biome-free curve underneath ──── */
    const L_LO = 0.3, L_HI = 2.4, T_LO = 0.6, T_HI = 1.35;

    function drawPhase(cv) {
        const ctx = cv.getContext('2d');
        const W = cv.width, H = cv.height;
        // Axis labels now live in HTML outside the canvas, so the plot only
        // needs a hairline margin to keep the marker from clipping.
        const px = L => (L - L_LO) / (L_HI - L_LO) * (W - 10) + 5;
        const py = T => H - 6 - (T - T_LO) / (T_HI - T_LO) * (H - 12);

        // Left transparent: the overlay is translucent, so the planet shows
        // through the plot.
        ctx.clearRect(0, 0, W, H);

        // The tolerable band. The growth window w = exp(-alpha (T-1)^4) falls to
        // a fraction q of its peak at |T - 1| = (-ln q)^(1/4) / alpha^(1/4), and
        // alpha = 8 (T_opt/dT)^4, so every contour's half-width is proportional
        // to the bandwidth dT. Stacking q = 0.05 ... 0.95 at low alpha builds the
        // same nested strips drawn in Fig. 1: widen dT and the band visibly
        // thickens.
        const a4 = Math.pow(8 * Math.pow(P.T_opt / P.dT, 4), 0.25);
        ctx.fillStyle = 'rgba(70, 200, 70, 0.028)';
        for (let q = 0.05; q < 0.96; q += 0.05) {
            const hw = Math.pow(-Math.log(q), 0.25) / a4;
            const yT = py(1 + hw), yB = py(1 - hw);
            ctx.fillRect(px(L_LO), yT, px(L_HI) - px(L_LO), yB - yT);
        }

        // T = T_opt
        ctx.strokeStyle = '#555'; ctx.setLineDash([3, 3]); ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(px(L_LO), py(1)); ctx.lineTo(px(L_HI), py(1)); ctx.stroke();
        ctx.setLineDash([]);

        // Bare planet: T = L^(1/4)
        ctx.strokeStyle = '#888'; ctx.lineWidth = 1.4;
        ctx.beginPath();
        for (let i = 0; i <= 100; i++) {
            const L = L_LO + (L_HI - L_LO) * i / 100;
            const T = Math.pow(L, 0.25);
            i ? ctx.lineTo(px(L), py(T)) : ctx.moveTo(px(L), py(T));
        }
        ctx.stroke();

        // Trail, then the live point.
        ctx.strokeStyle = 'rgba(0,150,0,0.55)'; ctx.lineWidth = 1.2;
        ctx.beginPath();
        for (let i = 0; i < trail.length; i++) {
            const p = trail[i];
            i ? ctx.lineTo(px(p[0]), py(p[1])) : ctx.moveTo(px(p[0]), py(p[1]));
        }
        ctx.stroke();

        ctx.fillStyle = '#2fbf2f';
        ctx.beginPath(); ctx.arc(px(X[3]), py(X[2]), 3.5, 0, 6.283); ctx.fill();
        // Axis labels are HTML (#edw-ylab, #edw-xlab): canvas text cannot be
        // typeset by MathJax.
    }

    /* ── Loop ───────────────────────────────────────────────────────────── */
    let canvas, phase;

    function frame() {
        if (running) {
            const th = theta();
            for (let i = 0; i < SUBSTEPS; i++) X = srk1(X, th);

            trail.push([X[3], X[2]]);
            if (trail.length > 220) trail.shift();

            blackMesh.count = Math.min(N_SITES, Math.round(X[0] * N_SITES));
            whiteMesh.count = Math.min(N_SITES - blackMesh.count,
                                       Math.round(X[1] * N_SITES));

            // The star's apparent size and brightness track the instantaneous
            // luminosity, so it swells as you push the slider and flickers with
            // the Ornstein-Uhlenbeck noise.
            starLight.intensity = 0.22 + 0.80 * X[3];
            const s = 0.6 + 0.85 * X[3];
            starGlow.scale.set(s, s, 1);

            updateReadouts();
            drawPhase(phase);
        }
        if (!dragging) spinY += 0.0018;
        planet.rotation.y = spinY;
        planet.rotation.x = spinX;
        renderer.render(scene, camera);
        frameId = requestAnimationFrame(frame);
    }

    function set(id, v) {
        const el = document.getElementById(id);
        if (el) el.textContent = v;
    }

    function updateReadouts() {
        set('edw-fb', X[0].toFixed(3));
        set('edw-fw', X[1].toFixed(3));
        set('edw-T', X[2].toFixed(3));
        set('edw-L', X[3].toFixed(3));
        set('edw-A', albedo().toFixed(3));
        set('edw-V', ((X[0] + X[1]) / P.f).toFixed(3));
    }

    /* ── UI ─────────────────────────────────────────────────────────────── */
    const UI = `
      <div id="edw-head">
        <span class="edw-title">Exo-Daisy World</span>
        <button class="edw-hbtn" id="edw-reset">Reset</button>
        <button class="edw-hbtn" id="edw-pause">Pause</button>
      </div>
      <div id="edw-stage">
        <canvas id="edw-canvas"></canvas>

        <div id="edw-readouts" class="edw-ov">
          <div class="edw-readout"><span class="k">\\(f_B\\)</span><span class="v" id="edw-fb">—</span></div>
          <div class="edw-readout"><span class="k">\\(T/T_{opt}\\)</span><span class="v" id="edw-T">—</span></div>
          <div class="edw-readout"><span class="k">\\(f_W\\)</span><span class="v" id="edw-fw">—</span></div>
          <div class="edw-readout"><span class="k">\\(L/L_{opt}\\)</span><span class="v" id="edw-L">—</span></div>
          <div class="edw-readout"><span class="k">\\(A\\)</span><span class="v" id="edw-A">—</span></div>
          <div class="edw-readout"><span class="k">\\(V\\)</span><span class="v" id="edw-V">—</span></div>
        </div>

        <div id="edw-plotwrap" class="edw-ov">
          <div id="edw-plotrow">
            <div id="edw-ylab"><span>\\(T/T_{opt}\\)</span></div>
            <canvas id="edw-phase" width="176" height="104"></canvas>
          </div>
          <div id="edw-xlab">\\(L/L_{opt}\\)</div>
        </div>

        <div id="edw-sliders" class="edw-ov">
          <div class="edw-row">
            <label for="edw-lam">\\(\\langle L\\rangle / L_{opt}\\)
              <span class="val" id="edw-lam-val">1.00</span></label>
            <input type="range" id="edw-lam" min="-0.7" max="1.4" step="0.005" value="0">
          </div>
          <div class="edw-row">
            <label for="edw-bw">\\(\\Delta T / T_{opt}\\)
              <span class="val" id="edw-bw-val">0.100</span></label>
            <input type="range" id="edw-bw" min="0.01" max="0.3" step="0.002" value="0.1">
          </div>
        </div>
      </div>
    `;

    /* The docs load MathJax (sphinx.ext.mathjax), but it has already swept the
     * page by the time we inject this markup, so typeset our labels by hand.
     * MathJax may still be loading, hence the startup promise. */
    function typesetMath(el) {
        const M = window.MathJax;
        if (!M) return;
        if (M.startup && M.startup.promise) {
            M.startup.promise.then(() => M.typesetPromise([el])).catch(() => {});
        } else if (M.typesetPromise) {
            M.typesetPromise([el]).catch(() => {});
        }
    }

    function boot() {
        const root = document.getElementById('edw-app');
        if (!root || typeof THREE === 'undefined') return;
        root.innerHTML = UI;
        typesetMath(root);

        canvas = document.getElementById('edw-canvas');
        phase = document.getElementById('edw-phase');

        buildSites();
        seed();
        initThree(canvas);
        resize(canvas);
        window.addEventListener('resize', () => resize(canvas));

        // Drag to spin.
        canvas.addEventListener('pointerdown', e => {
            dragging = true; lastX = e.clientX; lastY = e.clientY;
            canvas.setPointerCapture(e.pointerId);
        });
        canvas.addEventListener('pointermove', e => {
            if (!dragging) return;
            spinY += (e.clientX - lastX) * 0.006;
            spinX += (e.clientY - lastY) * 0.006;
            spinX = Math.max(-1.3, Math.min(1.3, spinX));
            lastX = e.clientX; lastY = e.clientY;
        });
        canvas.addEventListener('pointerup', e => {
            dragging = false;
            canvas.releasePointerCapture(e.pointerId);
        });

        document.getElementById('edw-lam').addEventListener('input', function () {
            P.lam = parseFloat(this.value);
            document.getElementById('edw-lam-val').textContent = (1 + P.lam).toFixed(2);
        });
        document.getElementById('edw-bw').addEventListener('input', function () {
            P.dT = parseFloat(this.value) * P.T_opt;
            document.getElementById('edw-bw-val').textContent = parseFloat(this.value).toFixed(3);
        });
        document.getElementById('edw-reset').addEventListener('click', seed);
        document.getElementById('edw-pause').addEventListener('click', function () {
            running = !running;
            this.textContent = running ? 'Pause' : 'Resume';
        });

        frameId = requestAnimationFrame(frame);
    }

    /* Three.js is loaded from the CDN by the page; wait for it if need be. */
    function whenReady() {
        if (typeof THREE !== 'undefined') { boot(); return; }
        setTimeout(whenReady, 60);
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', whenReady);
    } else {
        whenReady();
    }
})();
