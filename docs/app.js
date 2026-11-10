"use strict";

const WORLD_EXTENT = [-17250000, -8400000, 17250000, 8400000];
const ENSEMBLE_SIZE = 180;

const GROUP_LABELS = {
  priority: "Prioritised strategic actions",
  targets: "Strategic-action ensemble agreement",
  calibrated_target_scores: "Calibrated strategic-action scores",
  priority_diagnostics: "Priority diagnostics",
  connectivity: "Connectivity-class agreement",
  target2: "Restoration need and recovery"
};

const PALETTE_GRADIENTS = {
  agreement: "linear-gradient(90deg,#FFFFFF 0%,#DEEBF7 20%,#9ECAE1 40%,#6BAED6 60%,#3182BD 78%,#08519C 92%,#08306B 100%)",
  plasma: "linear-gradient(90deg,#FFFFFF 0%,#F4EFF8 7%,#0D0887 12%,#7E03A8 36%,#CC4778 59%,#F89540 81%,#F0F921 100%)",
  magma: "linear-gradient(90deg,#FFFFFF 0%,#F6EEF5 7%,#000004 12%,#3B0F70 36%,#8C2981 59%,#DE4968 81%,#FCFDBF 100%)",
  ylgn: "linear-gradient(90deg,#FFFFFF 0%,#FFFFE5 8%,#D9F0A3 30%,#78C679 55%,#238443 79%,#004529 100%)"
};

proj4.defs("EPSG:8857", "+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs +type=crs");
ol.proj.proj4.register(proj4);
const equalEarth = ol.proj.get("EPSG:8857");
equalEarth.setExtent(WORLD_EXTENT);
equalEarth.setWorldExtent([-180, -90, 180, 90]);

const $ = id => document.getElementById(id);
const sspSelect = $("ssp"), radiusSelect = $("radius"), groupSelect = $("group"), layerSelect = $("layer");
const legend = $("legend"), statusBox = $("status"), opacityInput = $("opacity");
const searchInput = $("place-search"), searchButton = $("search-button"), searchResults = $("search-results"), globalButton = $("global-view");
let catalog = null, rasterLayer = null, searchLayer = null;

function setStatus(message, error = false) { statusBox.textContent = message; statusBox.className = error ? "status error" : "status"; }
function unique(values) { return [...new Set(values)]; }
function setOptions(element, values, labelFunction) { const previous = element.value; element.innerHTML = ""; values.forEach(value => { const option = document.createElement("option"); option.value = String(value); option.textContent = labelFunction(value); element.appendChild(option); }); if (values.map(String).includes(previous)) element.value = previous; }

const oceanFeature = new ol.Feature(new ol.geom.Polygon([[[WORLD_EXTENT[0], WORLD_EXTENT[1]], [WORLD_EXTENT[2], WORLD_EXTENT[1]], [WORLD_EXTENT[2], WORLD_EXTENT[3]], [WORLD_EXTENT[0], WORLD_EXTENT[3]], [WORLD_EXTENT[0], WORLD_EXTENT[1]]]]));
const oceanLayer = new ol.layer.Vector({ source: new ol.source.Vector({ features: [oceanFeature] }), style: new ol.style.Style({ fill: new ol.style.Fill({ color: "#EAF1F4" }) }), zIndex: 0 });
const countrySource = new ol.source.Vector({ url: "data_boundaries/ne_110m_admin_0_countries_8857.geojson", format: new ol.format.GeoJSON({ dataProjection: "EPSG:8857", featureProjection: "EPSG:8857" }) });
const countriesLayer = new ol.layer.Vector({ source: countrySource, style: new ol.style.Style({ fill: new ol.style.Fill({ color: "rgba(255,255,255,0)" }), stroke: new ol.style.Stroke({ color: "rgba(45,55,55,0.68)", width: 0.85 }) }), zIndex: 20 });
const view = new ol.View({ projection: equalEarth, center: [0, 0], resolution: 52000, maxResolution: 150000, minResolution: 180, extent: WORLD_EXTENT, constrainOnlyCenter: true, showFullExtent: true });
const map = new ol.Map({ target: "map", layers: [oceanLayer, countriesLayer], view });

function globalView() { if (searchLayer) { map.removeLayer(searchLayer); searchLayer = null; } map.updateSize(); view.fit(WORLD_EXTENT, { size: map.getSize(), padding: [24, 24, 24, 24], duration: 250 }); }
function layersForContext() { return catalog.layers.filter(item => item.ssp === sspSelect.value && String(item.radius) === radiusSelect.value); }
function filteredLayers() { return layersForContext().filter(item => item.group === groupSelect.value); }
function selectedLayer() { return filteredLayers()[Number(layerSelect.value || 0)]; }
function isAgreementLayer(item) { return item.metric === "ensemble_agreement" || item.group === "connectivity" || item.file.includes("selection_percent") || item.file.includes("agreement_percent"); }
function actionAgreementGradient(item) { if (item.action === undefined || item.metric !== "ensemble_agreement") return null; const color = catalog.colors[String(item.action)]; return color ? `linear-gradient(90deg,#FFFFFF 0%,${color} 100%)` : null; }

function layerDescription(item) {
  if (item.kind === "categorical") {
    return (
      "Strategic action selected most frequently across the available " +
      "ensemble members in each grid cell."
    );
  }

  if (item.metric === "modal_agreement") {
    return (
      "Percentage of available ensemble members supporting the modal " +
      "strategic action in each grid cell. Higher values indicate greater " +
      "agreement about the most frequently selected action."
    );
  }

  if (item.metric === "ensemble_agreement") {
    return (
      "Percentage of available ensemble members selecting this specific " +
      "strategic action in each grid cell, whether or not it is the modal action."
    );
  }

  if (item.metric === "mean_calibrated_score") {
    return (
      "Mean calibrated score across available ensemble members. For each " +
      "strategic action, positive scores are expressed as empirical percentiles " +
      "relative to all positive scores for that action across SSPs, connectivity " +
      "radii and ensemble members. Higher values indicate relatively high " +
      "scores for that action."
    );
  }


  if (item.metric === "winning_margin") {
    return (
      "Mean difference between the highest and second-highest calibrated " +
      "scores across available ensemble members. Higher values indicate a " +
      "larger difference between the selected action and the next-ranked action."
    );
  }


  if (item.metric === "connectivity_agreement") {
    return (
      "Percentage of available connectivity ensemble members assigning " +
      "the grid cell to this connectivity class."
    );
  }

  if (item.metric === "recovery_score") {
    return (
      "Mean projected secondary-land recovery score across available ensemble " +
      "members. Values range from 0 to 1, with higher values indicating greater " +
      "projected recovery of secondary land."
    );
  }

  if (item.metric === "recovery_alignment") {
    return (
      "Indicates the spatial correspondence between Target 2 restoration " +
      "priority and projected secondary-land recovery. Higher values occur " +
      "where high restoration priority coincides with high projected recovery."
    );
  }

  if (item.metric === "restoration_gap") {
    return (
      "Indicates where Target 2 restoration priority exceeds projected " +
      "secondary-land recovery. Higher values occur where restoration priority " +
      "is high but projected recovery is low."
    );
  }

  return "Continuous ensemble summary.";
}

function legendLabels(item) { const maximum = Number(item.max ?? 100); if (isAgreementLayer(item)) return ["0%", "25%", "50%", "75%", "100% "]; if (maximum <= 1) return ["0", "0.25", "0.50", "0.75", "1.00"]; return ["0", (maximum * .25).toFixed(0), (maximum * .5).toFixed(0), (maximum * .75).toFixed(0), maximum.toFixed(0)]; }

function updateLegend(item) {
  if (!item) { legend.innerHTML = ""; return; }
  let html = `<h3>${item.title}</h3><p class="legend-description">${layerDescription(item)}</p>`;
  if (item.kind === "categorical") {
    Object.entries(catalog.classes).forEach(([code, label]) => { const color = catalog.colors[code]; if (color) html += `<div class="legend-row"><span class="swatch" style="background:${color}"></span><span> ${label}</span></div>`; });
  } else {
    const agreement = isAgreementLayer(item); const actionGradient = actionAgreementGradient(item); const palette = agreement ? "agreement" : String(item.palette || "plasma").toLowerCase(); const gradient = actionGradient || PALETTE_GRADIENTS[palette] || PALETTE_GRADIENTS.plasma; const labels = legendLabels(item); const zeroText = agreement ? "0% = no ensemble agreement" : "0 = no positive evidence";
    html += `<div class="continuous-legend"><div class="gradient-zero"><span class="zero-chip"></span><span>${zeroText}</span></div><div class="gradient-bar" style="background:${gradient}"></div><div class="gradient-labels">${labels.map(value => `<span>${value}</span>`).join("")}</div></div>`;
  }
  legend.innerHTML = html;
}

function showLayer() { const item = selectedLayer(); if (rasterLayer) { map.removeLayer(rasterLayer); rasterLayer = null; } if (!item) { setStatus("No layer available for this selection.", true); updateLegend(null); return; } const extent = item.imageExtent.map(Number); const source = new ol.source.ImageStatic({ url: `${item.file}?v=${Date.now()}`, projection: equalEarth, imageExtent: extent, interpolate: false }); source.on("imageloadstart", () => setStatus(`Loading ${item.title}...`)); source.on("imageloadend", () => setStatus(`${catalog.sspLabels[item.ssp]} | ${catalog.radiusLabels[String(item.radius)]} | ${item.title}`)); source.on("imageloaderror", () => setStatus(`Layer failed to load: ${item.file}`, true)); rasterLayer = new ol.layer.Image({ source, opacity: Number(opacityInput.value || 1), zIndex: 2 }); map.addLayer(rasterLayer); countriesLayer.setZIndex(20); map.updateSize(); view.fit(extent, { size: map.getSize(), padding: [20, 20, 20, 20], duration: 0 }); updateLegend(item); }
function refreshLayerMenu() { const available = filteredLayers(); layerSelect.innerHTML = ""; available.forEach((item, index) => { const option = document.createElement("option"); option.value = String(index); option.textContent = item.title; layerSelect.appendChild(option); }); showLayer(); }
function refreshGroups() { const groups = unique(layersForContext().map(item => item.group)); setOptions(groupSelect, groups, value => GROUP_LABELS[value] || value); refreshLayerMenu(); }
function normalizeText(value) { return String(value || "").normalize("NFD").replace(/[\u0300-\u036f]/g, "").toLowerCase().trim(); }
function countryNames(feature) { const p = feature.getProperties(); return [p.ADMIN, p.NAME, p.NAME_EN, p.SOVEREIGNT, p.BRK_NAME, p.FORMAL_EN].filter(Boolean).map(String); }
function fitCountry(feature) { const geometry = feature.getGeometry(); if (!geometry) return; view.fit(geometry.getExtent(), { size: map.getSize(), padding: [55, 55, 55, 55], maxZoom: 6.5, duration: 350 }); if (searchLayer) map.removeLayer(searchLayer); const clone = feature.clone(); clone.setStyle(new ol.style.Style({ fill: new ol.style.Fill({ color: "rgba(213,94,0,0.08)" }), stroke: new ol.style.Stroke({ color: "#D55E00", width: 2.4 }) })); searchLayer = new ol.layer.Vector({ source: new ol.source.Vector({ features: [clone] }), zIndex: 30 }); map.addLayer(searchLayer); }
function searchCountry() { const query = normalizeText(searchInput.value); if (!query) { searchResults.textContent = "Enter a country name."; return; } const matches = countrySource.getFeatures().map(feature => { const names = countryNames(feature), normalised = names.map(normalizeText); let score = 99; if (normalised.some(name => name === query)) score = 0; else if (normalised.some(name => name.startsWith(query))) score = 1; else if (normalised.some(name => name.includes(query))) score = 2; return { feature, label: names[0] || "Country", score }; }).filter(item => item.score < 99).sort((a, b) => a.score - b.score || a.label.localeCompare(b.label)); if (!matches.length) { searchResults.textContent = "Country not found."; return; } fitCountry(matches[0].feature); searchResults.innerHTML = ""; matches.slice(0, 5).forEach((match, index) => { const button = document.createElement("button"); button.type = "button"; button.className = "search-result"; button.textContent = `${index === 0 ? "Best match: " : ""}${match.label}`; button.addEventListener("click", () => { fitCountry(match.feature); searchResults.innerHTML = ""; }); searchResults.appendChild(button); }); }
async function initialise() { try { setStatus("Loading catalogue..."); const response = await fetch(`catalog_png.json?v=${Date.now()}`, { cache: "no-store" }); if (!response.ok) throw new Error(`catalog_png.json returned HTTP ${response.status}`); catalog = await response.json(); if (!catalog.layers?.length) throw new Error("Catalogue contains no layers"); setOptions(sspSelect, Object.keys(catalog.sspLabels), value => catalog.sspLabels[value]); setOptions(radiusSelect, Object.keys(catalog.radiusLabels), value => catalog.radiusLabels[value]); if (Object.keys(catalog.sspLabels).includes("ssp126")) sspSelect.value = "ssp126"; if (Object.keys(catalog.radiusLabels).includes("6")) radiusSelect.value = "6"; refreshGroups(); sspSelect.addEventListener("change", refreshGroups); radiusSelect.addEventListener("change", refreshGroups); groupSelect.addEventListener("change", refreshLayerMenu); layerSelect.addEventListener("change", showLayer); opacityInput.addEventListener("input", () => { if (rasterLayer) rasterLayer.setOpacity(Number(opacityInput.value)); }); searchButton.addEventListener("click", searchCountry); searchInput.addEventListener("keydown", event => { if (event.key === "Enter") searchCountry(); }); globalButton.addEventListener("click", globalView); setTimeout(() => { map.updateSize(); globalView(); showLayer(); }, 100); } catch (error) { console.error(error); setStatus(`Viewer initialisation failed: ${error.message}`, true); } }
initialise();
