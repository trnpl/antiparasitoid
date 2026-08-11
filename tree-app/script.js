// =============================================================================
// Drosophila Phylogeny Explorer -- static site logic
//
// App built using Claude Sonnet 5, last modified by RLT 260811
// =============================================================================

const TREE_VARIANTS = {
  plain:        "layers/tree_plain.svg",
  asr:          "layers/tree_asr_branches.svg",
  cafe_ppo1:    "layers/cafe_ppo1.svg",
  cafe_ppo234:  "layers/cafe_ppo234.svg"
};

const LAYER_ORDER = [
  { key: "asr_pies",      files: ["layers/asr_pies.svg"],                                   controlledBy: "show_pies" },
  { key: "states",        files: ["layers/ring_states.svg"],                                 controlledBy: "ring:states" },
  { key: "ppo1",          files: ["layers/ring_ppo1.svg", "layers/ring_ppo1_qc.svg"],        controlledBy: "ring:ppo1" },
  { key: "ppo2",          files: ["layers/ring_ppo2.svg", "layers/ring_ppo2_qc.svg"],        controlledBy: "ring:ppo2" },
  { key: "ppo3",          files: ["layers/ring_ppo3.svg", "layers/ring_ppo3_qc.svg"],        controlledBy: "ring:ppo3" },
  { key: "ppo4",          files: ["layers/ring_ppo4.svg", "layers/ring_ppo4_qc.svg"],        controlledBy: "ring:ppo4" },
  { key: "cdtb",          files: ["layers/ring_cdtb.svg", "layers/ring_cdtb_qc.svg"],        controlledBy: "ring:cdtb" },
  { key: "clade_labels",  files: ["layers/clade_labels.svg"],                                 controlledBy: "show_clade_labels" },
  { key: "bootstrap",     files: ["layers/bootstrap.svg"],                                    controlledBy: "show_bootstrap" },
  { key: "scale_bar",     files: ["layers/scale_bar.svg"],                                    controlledBy: "show_scale_bar" },
  { key: "tip_labels",    files: ["layers/tip_labels.svg"],                                   controlledBy: "show_tip_labels" }
];

const treeStack = document.getElementById("tree-stack");
const imageEls = {}; // key -> [<img>, ...]

function buildImageStack() {
  // The .tree-stack container has no explicit height, and every layer below
  // is position:absolute (removed from normal flow) -- so without this,
  // the container collapses to 0px tall, and every hover cell's top/height
  // (given as %) resolves against that zero, making them all unhoverable.
  // This sizer stays in NORMAL flow (not absolute) so its own intrinsic
  // aspect ratio gives the container real height; every other layer is
  // still absolutely positioned on top of it as before.
  const sizer = document.createElement("img");
  sizer.className = "sizer";
  sizer.src = TREE_VARIANTS.plain; // any layer works, they all share one aspect ratio
  sizer.alt = "";
  treeStack.appendChild(sizer);

  // Base tree variants first (bottom of stack) -- only one visible at a time.
  Object.keys(TREE_VARIANTS).forEach(key => {
    imageEls[key] = [makeImg(TREE_VARIANTS[key])];
  });
  Object.keys(TREE_VARIANTS).forEach(key => {
    imageEls[key][0].style.display = (key === "plain") ? "block" : "none";
  });

  LAYER_ORDER.forEach(layer => {
    imageEls[layer.key] = layer.files.map(f => makeImg(f));
  });

  // asr_pies has no dedicated checkbox anymore (it follows the ASR radio
  // button) so it needs its initial visibility set explicitly here,
  // matching "Plain" being the default selected tree variant.
  setLayerVisible("asr_pies", false);
}

function makeImg(src) {
  const img = document.createElement("img");
  img.className = "layer";
  img.src = src;
  img.alt = "";
  treeStack.appendChild(img);
  return img;
}

function setLayerVisible(key, visible) {
  (imageEls[key] || []).forEach(img => { img.style.display = visible ? "block" : "none"; });
}

// ---------------------------------------------------------------------------
// Wire up controls
// ---------------------------------------------------------------------------

function initControls() {
  // Base tree variant (radio)
  document.querySelectorAll('input[name="tree_variant"]').forEach(radio => {
    radio.addEventListener("change", () => {
      Object.keys(TREE_VARIANTS).forEach(key => {
        setLayerVisible(key, radio.value === key && radio.checked);
      });
      // Pies always follow ASR branch coloring specifically, not any other variant.
      setLayerVisible("asr_pies", radio.value === "asr" && radio.checked);
      renderLegend();
    });
  });

  // Simple checkbox -> layer toggles
  const simpleToggles = [
    ["show_tip_labels", "tip_labels"],
    ["show_clade_labels", "clade_labels"],
    ["show_bootstrap", "bootstrap"],
    ["show_scale_bar", "scale_bar"]
  ];
  simpleToggles.forEach(([inputId, layerKey]) => {
    const el = document.getElementById(inputId);
    setLayerVisible(layerKey, el.checked);
    el.addEventListener("change", () => setLayerVisible(layerKey, el.checked));
  });

  // Tile ring toggles
  document.querySelectorAll(".ring-toggle").forEach(cb => {
    setLayerVisible(cb.value, cb.checked);
    cb.addEventListener("change", () => {
      setLayerVisible(cb.value, cb.checked);
      updateActiveCountBadge();
      renderLegend();
    });
  });

  // Hover grid on/off
  const hoverToggle = document.getElementById("show_hover");
  hoverToggle.addEventListener("change", () => {
    const grid = document.getElementById("hover-grid");
    if (grid) grid.style.display = hoverToggle.checked ? "block" : "none";
  });

  updateActiveCountBadge();
}

function updateActiveCountBadge() {
  const rings = document.querySelectorAll(".ring-toggle");
  const activeCount = Array.from(rings).filter(r => r.checked).length;
  document.getElementById("active-count-badge").textContent =
    `${activeCount} of ${rings.length} characters shown`;
}

// ---------------------------------------------------------------------------
// Hover tooltip grid -- built once from tip_coordinates.json, positioned
// with simple CSS percentages. Tippy is attached immediately after each
// div is created.
// ---------------------------------------------------------------------------

function buildHoverGrid(coords) {
  const grid = document.createElement("div");
  grid.id = "hover-grid";
  grid.style.position = "absolute";
  grid.style.top = "0";
  grid.style.left = "0";
  grid.style.width = "100%";
  grid.style.height = "100%";

  const rowH = coords.row_height_pct;

  // Map ring key -> x-range, and ring key -> which tooltip field to use.
  const ringXRanges = {};
  coords.rings.forEach(r => {
    ringXRanges[r.ring] = { start: r.x_start_pct, width: r.x_end_pct - r.x_start_pct };
  });
  const ringTooltipField = {
    states: "tt_states", ppo1: "tt_ppo1", ppo2: "tt_ppo2",
    ppo3: "tt_ppo3", ppo4: "tt_ppo4", cdtb: "tt_cdtb"
  };

  // These rings hide their box entirely when the value is 0/NA (see
  // hide_na_zero in build_tree_layers.R) -- suppress the hover target
  // there too, since there's nothing drawn to hover over. ppo1/ppo2 are
  // deliberately NOT in this list; their boxes (and hover) show at 0.
  const ringValueField = { states: "val_states", ppo3: "val_ppo3", ppo4: "val_ppo4", cdtb: "val_cdtb" };
  const suppressZeroRings = new Set(Object.keys(ringValueField));

  const labelStart = coords.tip_label_region.x_start_pct;
  const labelWidth = coords.tip_label_region.x_end_pct - coords.tip_label_region.x_start_pct;

  coords.tips.forEach(tip => {
    const top = tip.y_pct - rowH / 2;

    // One cell over the tip-label area -> full summary.
    addHoverCell(grid, labelStart, top, labelWidth, rowH, tip.tt_full);

    // One cell PER RING -> that ring's specific detail (e.g. PPO1 copy number).
    Object.keys(ringXRanges).forEach(ringKey => {
      if (suppressZeroRings.has(ringKey)) {
        const v = tip[ringValueField[ringKey]];
        if (v === 0 || v === null || v === undefined) return; // no box drawn there -- skip the hover target too
      }
      const xr = ringXRanges[ringKey];
      const field = ringTooltipField[ringKey];
      addHoverCell(grid, xr.start, top, xr.width, rowH, tip[field]);
    });
  });

  treeStack.appendChild(grid);
}

function addHoverCell(grid, xStartPct, yStartPct, widthPct, heightPct, tooltipHtml) {
  const cell = document.createElement("div");
  cell.className = "hover-cell";
  cell.style.left = xStartPct + "%";
  cell.style.top = yStartPct + "%";
  cell.style.width = widthPct + "%";
  cell.style.height = heightPct + "%";
  grid.appendChild(cell);

  tippy(cell, {
    content: tooltipHtml,
    allowHTML: true,
    placement: "right"
  });
}

// ---------------------------------------------------------------------------
// Legend / key panel -- reflects only whichever rings/ASR mode are
// currently toggled on. Built from coords.legend (categorical palettes +
// continuous value ranges), fetched once alongside the hover-grid data.
// ---------------------------------------------------------------------------

let legendData = null;

const RING_LABELS = {
  states: "Defense mechanism",
  ppo1: "PPO1 copy number",
  ppo2: "PPO2 copy number",
  ppo3: "PPO3 copy number",
  ppo4: "PPO4 copy number",
  cdtb: "cdtB copy number"
};

function renderLegend() {
  const container = document.getElementById("legend-content");
  container.innerHTML = "";
  if (!legendData) return;
  
  // QC legend is shared across PPO1-4/cdtb (not per-ring) -- show it once
  // if ANY of those rings is currently active, since that's when QC
  // markers could actually be visible on the tree.
  const qcBearingRings = ["ppo1", "ppo2", "ppo3", "ppo4", "cdtb"];
  const anyQcActive = qcBearingRings.some(key => {
    const cb = document.querySelector(`.ring-toggle[value="${key}"]`);
    return cb && cb.checked;
  });
  if (anyQcActive && legendData.qc) {
    container.appendChild(buildLegendBlock("QC", legendData.qc, "qc"));
  }

  // Tile rings: one legend block per currently-checked ring toggle.
  document.querySelectorAll(".ring-toggle").forEach(cb => {
    if (!cb.checked) return;
    const key = cb.value;
    const entry = legendData[key];
    if (!entry) return;
    container.appendChild(buildLegendBlock(RING_LABELS[key] || key, entry, key));
  });

  // ASR branch coloring: only relevant when that tree variant is selected.
  const selectedVariant = document.querySelector('input[name="tree_variant"]:checked').value;
  const variantLegendKey = {
    asr: "asr",
    cafe_ppo1: "cafe_ppo1",
    cafe_ppo234: "cafe_ppo234"
  }[selectedVariant];
  const variantLegendTitle = {
    asr: "ASR branch state",
    cafe_ppo1: "PPO copy number (CAFE5)",
    cafe_ppo234: "PPO copy number (CAFE5)"
  }[selectedVariant];
  // Icon key is separate from the data key: both CAFE5 variants share one
  // set of icon files (icons/cafe_ppo_copy_*.png / icons/cafe_ppo_copy.png)
  // even though their underlying legend data is tracked separately.
  const variantIconKey = {
    asr: "asr",
    cafe_ppo1: "cafe_ppo_copy",
    cafe_ppo234: "cafe_ppo_copy"
  }[selectedVariant];
  if (variantLegendKey && legendData[variantLegendKey]) {
    container.appendChild(buildLegendBlock(variantLegendTitle, legendData[variantLegendKey], variantIconKey));
  }

  if (!container.children.length) {
    const empty = document.createElement("p");
    empty.className = "legend-empty";
    empty.textContent = "No colored layers active.";
    container.appendChild(empty);
  }
}

// Icons override color swatches by naming convention: drop a file at
// icons/{ringKey}_{value}.png and it's used automatically -- no code
// changes needed per icon. Falls back to the color swatch if the image
// 404s, so partial icon sets (some values have custom art, others don't)
// work without extra configuration.
// Two icon conventions, checked independently:
//   1. icons/{ringKey}_{value}.png -- per-value, overrides one swatch.
//   2. icons/{ringKey}.png         -- one static image for the WHOLE
//      group, overrides every swatch at once (e.g. a single reference
//      diagram instead of uploading one icon per QC code).
// Both fall back gracefully to the normal color swatches/gradient if the
// corresponding image isn't present.
function buildLegendBlock(title, entry, ringKey) {
  const block = document.createElement("div");
  block.className = "legend-block";

  const heading = document.createElement("h3");
  heading.textContent = title;
  block.appendChild(heading);

  // Everything below goes into `content` so a whole-group static image
  // can hide it as one unit, rather than needing to know its internal shape.
  const content = document.createElement("div");

  if (entry.type === "categorical") {
    const list = document.createElement("div");
    list.className = "legend-swatch-list";
    entry.items.forEach(item => {
      const row = document.createElement("div");
      row.className = "legend-swatch-row";

      const swatch = document.createElement("span");
      swatch.className = "legend-swatch";
      swatch.style.background = item.color;
      row.appendChild(swatch);

      if (ringKey) {
        const icon = document.createElement("img");
        icon.className = "legend-swatch-icon";
        icon.src = `icons/${ringKey}_${item.value}.png`;
        icon.alt = "";
        icon.onerror = () => { icon.remove(); };              // no custom icon -- keep the color swatch
        icon.onload = () => { swatch.style.display = "none"; }; // custom icon loaded -- hide the swatch it replaces
        row.appendChild(icon);
      }

      const label = document.createElement("span");
      label.className = "legend-swatch-label";
      label.textContent = item.label || item.value;

      row.appendChild(label);
      list.appendChild(row);
    });
    content.appendChild(list);
  } else if (entry.type === "continuous") {
    const bar = document.createElement("div");
    bar.className = "legend-gradient-bar";
    bar.style.background = `linear-gradient(to right, ${entry.low_color}, ${entry.high_color})`;

    const labels = document.createElement("div");
    labels.className = "legend-gradient-labels";
    const minEl = document.createElement("span");
    minEl.textContent = entry.min != null ? entry.min : "\u2013";
    const maxEl = document.createElement("span");
    maxEl.textContent = entry.max != null ? entry.max : "\u2013";
    labels.appendChild(minEl);
    labels.appendChild(maxEl);

    content.appendChild(bar);
    content.appendChild(labels);
  }

  block.appendChild(content);

  if (ringKey) {
    const staticImg = document.createElement("img");
    staticImg.className = "legend-static-image";
    staticImg.src = `icons/${ringKey}.png`;
    staticImg.alt = "";
    staticImg.onerror = () => { staticImg.remove(); };            // no whole-group image -- keep normal swatches
    staticImg.onload = () => { content.style.display = "none"; }; // whole-group image found -- replace swatches entirely
    block.appendChild(staticImg);
  }

  return block;
}

async function init() {
  buildImageStack();
  initControls();

  let coords;
  try {
    const res = await fetch("tip_coordinates.json");
    if (!res.ok) throw new Error(`HTTP ${res.status} fetching tip_coordinates.json`);
    coords = await res.json();
  } catch (err) {
    console.error("Could not fetch/parse tip_coordinates.json:", err);
    console.warn(
      "If you're opening index.html directly (file://), fetch() will be " +
      "blocked by the browser's CORS policy. Serve this folder with a local " +
      "server instead, e.g.: python3 -m http.server, then open " +
      "http://localhost:8000"
    );
    return;
  }

  // Schema check -- catches a stale tip_coordinates.json left over from an
  // older version of build_tree_layers.R before it causes a confusing
  // silent failure deep inside buildHoverGrid.
  const requiredTipFields = ["tt_full", "tt_states", "tt_ppo1", "tt_ppo2", "tt_ppo3", "tt_ppo4", "tt_cdtb"];
  const sampleTip = coords.tips && coords.tips[0];
  const missingTopLevel = !coords.tip_label_region || !coords.rings || !sampleTip;
  const missingTipFields = sampleTip && requiredTipFields.filter(f => !(f in sampleTip));

  if (missingTopLevel || (missingTipFields && missingTipFields.length)) {
    console.error(
      "tip_coordinates.json is missing expected fields " +
      "(tip_label_region / tt_states / tt_ppo1 etc). This usually means " +
      "it's a stale file from an older version of build_tree_layers.R -- " +
      "re-run build_tree_layers.R to regenerate it, then reload this page.",
      { missingTopLevel, missingTipFields }
    );
    return;
  }

  try {
    buildHoverGrid(coords);
  } catch (err) {
    console.error("Hover grid failed to build:", err);
  }

  legendData = coords.legend || null;
  renderLegend();
}

init();
