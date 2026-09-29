"use strict";

const PERIOD_LABELS = { J1: "Mid-century", J2: "Late century" };
const SSP_LABELS = {
  "ssp126soc-adapt": "SSP1-2.6",
  "ssp370soc-adapt": "SSP3-7.0",
  "ssp585soc-adapt": "SSP5-8.5"
};
const MODEL_LABELS = {
  combined: "IMAGE-MAgPIE consensus",
  image: "IMAGE",
  magpie: "MAgPIE"
};
const WORLD_EXTENT = [-17243959.06, -8392927.60, 17243959.06, 8392927.60];

function firstElement(ids, fallback = null) {
  for (const id of ids) {
    const element = document.getElementById(id);
    if (element) return element;
  }
  return fallback;
}

const selects = Array.from(document.querySelectorAll("select"));
const statusBox = firstElement(["status", "viewer-status", "status-box"]);
const scenarioSelect = firstElement(["scenario", "scenario-select", "ssp-select"], selects[0]);
const periodSelect = firstElement(["period", "period-select", "jump-select"], selects[1]);
const modelSelect = firstElement(["model", "model-select", "lum-select"], selects[2]);
const layerSelect = firstElement(["layer", "layer-select", "variable-select"], selects[3]);
const opacityInput = firstElement(["opacity", "opacity-slider"], document.querySelector('input[type="range"]'));
const searchInput = firstElement(["search-input", "searchInput", "country-search"], document.querySelector('input[placeholder*="Mexico" i]'));
const searchButton = firstElement(["search-button", "searchButton", "find-button"], Array.from(document.querySelectorAll("button")).find(b => b.textContent.trim().toLowerCase() === "find"));
const searchResults = firstElement(["search-results", "searchResults", "country-results"]);
const globalButton = firstElement(["global-view", "globalView", "global-view-button"], Array.from(document.querySelectorAll("button")).find(b => b.textContent.trim().toLowerCase() === "global view"));
const legend = firstElement(["legend", "map-legend", "legend-container"]);
const modelHelp = firstElement(["model-help", "modelHelp", "model-description"]);

function setStatus(message, type = "ok") {
  if (!statusBox) return;
  statusBox.textContent = message;
  statusBox.className = `status ${type}`;
}

proj4.defs("EPSG:8857", "+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs +type=crs");
ol.proj.proj4.register(proj4);
const equalEarth = ol.proj.get("EPSG:8857");
equalEarth.setExtent(WORLD_EXTENT);
equalEarth.setWorldExtent([-180, -90, 180, 90]);

const view = new ol.View({
  projection: equalEarth,
  center: [0, 0],
  resolution: 52000,
  maxResolution: 150000,
  minResolution: 180,
  extent: WORLD_EXTENT,
  constrainOnlyCenter: true,
  showFullExtent: true,
  smoothExtentConstraint: false
});

const oceanLayer = new ol.layer.Vector({
  source: new ol.source.Vector({
    features: [
      new ol.Feature(new ol.geom.Polygon([[
        [WORLD_EXTENT[0], WORLD_EXTENT[1]], [WORLD_EXTENT[2], WORLD_EXTENT[1]],
        [WORLD_EXTENT[2], WORLD_EXTENT[3]], [WORLD_EXTENT[0], WORLD_EXTENT[3]],
        [WORLD_EXTENT[0], WORLD_EXTENT[1]]
      ]]))
    ]
  }),
  style: new ol.style.Style({ fill: new ol.style.Fill({ color: "#EAF1F4" }) }),
  zIndex: 0
});

const countrySource = new ol.source.Vector();
const countriesLayer = new ol.layer.Vector({
  source: countrySource,
  style: new ol.style.Style({
    fill: new ol.style.Fill({ color: "rgba(255,255,255,0)" }),
    stroke: new ol.style.Stroke({ color: "#394842", width: 1.05 })
  }),
  zIndex: 20
});

const map = new ol.Map({ target: "map", layers: [oceanLayer, countriesLayer], view });
let catalog = null;
let rasterLayer = null;
let searchLayer = null;
let boundariesReady = false;

function globalView() {
  if (searchLayer) {
    map.removeLayer(searchLayer);
    searchLayer = null;
  }
  map.updateSize();
  view.fit(WORLD_EXTENT, {
    size: map.getSize(),
    padding: [25, 25, 25, 25],
    duration: 350,
    nearest: false
  });
}

async function loadCountryBoundaries() {
  const url = `data_boundaries/ne_110m_admin_0_countries_8857.geojson?v=${Date.now()}`;
  try {
    const response = await fetch(url, { cache: "no-store" });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const json = await response.json();
    const features = new ol.format.GeoJSON().readFeatures(json, {
      dataProjection: "EPSG:8857",
      featureProjection: "EPSG:8857"
    });
    if (!features.length) throw new Error("GeoJSON contains no features");
    countrySource.clear(true);
    countrySource.addFeatures(features);
    boundariesReady = true;
    countriesLayer.changed();
    if (searchResults) searchResults.textContent = "Country search ready.";
    console.log(`Loaded ${features.length} country boundaries`);
  } catch (error) {
    boundariesReady = false;
    console.error("Country boundary loading failed:", error);
    if (searchResults) searchResults.textContent = `Country boundaries failed: ${error.message}`;
  }
}

function unique(values) { return [...new Set(values)]; }
function addOptions(select, values, labels = {}) {
  select.innerHTML = "";
  values.forEach(value => {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = labels[value] || value;
    select.appendChild(option);
  });
}
function matchingLayers() {
  return catalog.layers.filter(item => item.ssp === scenarioSelect.value && item.period === periodSelect.value && item.model === modelSelect.value);
}
function selectedMetadata() { return matchingLayers()[Number(layerSelect.value || 0)]; }

function updateHelp() {
  if (!modelHelp) return;
  modelHelp.textContent = modelSelect.value === "combined"
    ? "Consensus-adjusted IMAGE-MAgPIE score. An ineligible model contributes zero."
    : `Normalised within-action percentile calculated from ${MODEL_LABELS[modelSelect.value]} outputs only.`;
}
function updateLegend(item) {
  if (!legend || !item) return;
  if (item.kind === "categorical") {
    legend.innerHTML = `<strong>${item.title}</strong>` + Object.entries(catalog.classes).map(([code, label]) => `<div><span style="display:inline-block;width:16px;height:12px;margin-right:6px;background:${catalog.colors[code]};border:1px solid #555"></span>${label}</div>`).join("");
  } else {
    const text = modelSelect.value === "combined" ? "Consensus-adjusted score" : `${MODEL_LABELS[modelSelect.value]} normalized within-action percentile`;
    legend.innerHTML = `<strong>${item.title}</strong><div>${text}</div><div style="display:flex;justify-content:space-between"><span>${item.min}</span><span>${item.max}</span></div><div style="height:14px;background:linear-gradient(90deg,#440154,#31688e,#35b779,#fde725)"></div>`;
  }
}
function loadSelectedLayer() {
  const item = selectedMetadata();
  if (!item) return;
  if (rasterLayer) map.removeLayer(rasterLayer);
  const source = new ol.source.ImageStatic({
    url: `${item.file}?v=2`,
    imageExtent: item.imageExtent,
    projection: "EPSG:8857",
    interpolate: item.kind !== "categorical",
    crossOrigin: "anonymous"
  });
  rasterLayer = new ol.layer.Image({ source, opacity: Number(opacityInput.value || 1), zIndex: 2 });
  map.addLayer(rasterLayer);
  countriesLayer.setZIndex(20);
  source.on("imageloaderror", () => setStatus(`Layer failed to load: ${item.file}`, "error"));
  source.on("imageloadend", () => setStatus(`${SSP_LABELS[item.ssp]} | ${PERIOD_LABELS[item.period]} | ${MODEL_LABELS[item.model]} | ${item.title}`));
  updateLegend(item);
}
function refreshLayerOptions() {
  const layers = matchingLayers();
  layerSelect.innerHTML = "";
  layers.forEach((item, index) => {
    const option = document.createElement("option");
    option.value = String(index);
    option.textContent = item.title;
    layerSelect.appendChild(option);
  });
  updateHelp();
  loadSelectedLayer();
}

function normalizeText(value) {
  return value.normalize("NFD").replace(/[\u0300-\u036f]/g, "").toLowerCase().trim();
}
function countryNames(feature) {
  const p = feature.getProperties();
  return [p.ADMIN, p.NAME, p.NAME_EN, p.SOVEREIGNT, p.BRK_NAME, p.FORMAL_EN].filter(Boolean).map(String);
}
function fitCountry(feature) {
  const extent = feature.getGeometry().getExtent();
  if (!extent.every(Number.isFinite) || extent[0] >= extent[2] || extent[1] >= extent[3]) return;
  map.updateSize();
  view.fit(extent, { size: map.getSize(), padding: [55, 55, 55, 55], maxZoom: 6.5, duration: 400 });
  if (searchLayer) map.removeLayer(searchLayer);
  const highlighted = feature.clone();
  highlighted.setStyle(new ol.style.Style({ fill: new ol.style.Fill({ color: "rgba(213,94,0,0.08)" }), stroke: new ol.style.Stroke({ color: "#D55E00", width: 2.4 }) }));
  searchLayer = new ol.layer.Vector({ source: new ol.source.Vector({ features: [highlighted] }), zIndex: 30 });
  map.addLayer(searchLayer);
}
function searchCountry() {
  const query = normalizeText(searchInput.value || "");
  if (!query) return;
  if (!boundariesReady) {
    if (searchResults) searchResults.textContent = "Country boundaries are not ready. Reload the page if this persists.";
    return;
  }
  const matches = countrySource.getFeatures().map(feature => {
    const names = countryNames(feature);
    const normalized = names.map(normalizeText);
    let score = 99;
    if (normalized.some(name => name === query)) score = 0;
    else if (normalized.some(name => name.startsWith(query))) score = 1;
    else if (normalized.some(name => name.includes(query))) score = 2;
    return { feature, label: names[0] || "Country", score };
  }).filter(item => item.score < 99).sort((a, b) => a.score - b.score || a.label.localeCompare(b.label));
  if (!matches.length) {
    if (searchResults) searchResults.textContent = "Country not found.";
    return;
  }
  fitCountry(matches[0].feature);
  if (searchResults) {
    searchResults.innerHTML = "";
    matches.slice(0, 5).forEach((item, index) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "search-result";
      button.textContent = `${index === 0 ? "Best match: " : ""}${item.label}`;
      button.addEventListener("click", () => { fitCountry(item.feature); searchResults.innerHTML = ""; });
      searchResults.appendChild(button);
    });
  }
}

async function initialise() {
  try {
    const [catalogResponse] = await Promise.all([
      fetch(`catalog_png.json?v=${Date.now()}`, { cache: "no-store" }),
      loadCountryBoundaries()
    ]);
    if (!catalogResponse.ok) throw new Error(`catalog_png.json returned ${catalogResponse.status}`);
    catalog = await catalogResponse.json();
    addOptions(scenarioSelect, unique(catalog.layers.map(item => item.ssp)), SSP_LABELS);
    addOptions(periodSelect, unique(catalog.layers.map(item => item.period)), PERIOD_LABELS);
    addOptions(modelSelect, unique(catalog.layers.map(item => item.model)), MODEL_LABELS);
    scenarioSelect.addEventListener("change", refreshLayerOptions);
    periodSelect.addEventListener("change", refreshLayerOptions);
    modelSelect.addEventListener("change", refreshLayerOptions);
    layerSelect.addEventListener("change", loadSelectedLayer);
    opacityInput.addEventListener("input", () => { if (rasterLayer) rasterLayer.setOpacity(Number(opacityInput.value)); });
    searchButton.addEventListener("click", searchCountry);
    searchInput.addEventListener("keydown", event => { if (event.key === "Enter") searchCountry(); });
    globalButton.addEventListener("click", globalView);
    refreshLayerOptions();
    setTimeout(globalView, 100);
  } catch (error) {
    console.error(error);
    setStatus(`Viewer initialisation failed: ${error.message}`, "error");
  }
}

initialise();
