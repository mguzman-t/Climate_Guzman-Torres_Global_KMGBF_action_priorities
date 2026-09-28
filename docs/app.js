"use strict";

const statusBox = document.getElementById("status");
const sspSelect = document.getElementById("ssp");
const periodSelect = document.getElementById("period");
const layerSelect = document.getElementById("layer");
const opacityInput = document.getElementById("opacity");
const legend = document.getElementById("legend");

const PERIOD_LABELS = {
  J1: "2041–2070",
  J2: "2071–2100"
};

const SCENARIO_LABELS = {
  "ssp126soc-adapt": "SSP1-2.6",
  "ssp370soc-adapt": "SSP3-7.0",
  "ssp585soc-adapt": "SSP5-8.5"
};

const palettes = {
  viridis: [
    [0.00, "#440154"],
    [0.25, "#3b528b"],
    [0.50, "#21918c"],
    [0.75, "#5ec962"],
    [1.00, "#fde725"]
  ],
  magma: [
    [0.00, "#000004"],
    [0.25, "#51127c"],
    [0.50, "#b73779"],
    [0.75, "#fc8961"],
    [1.00, "#fcfdbf"]
  ],
  ylgn: [
    [0.00, "#ffffe5"],
    [0.50, "#78c679"],
    [1.00, "#005a32"]
  ]
};

let catalog = null;
let currentLayer = null;

function setStatus(message, state = "") {
  statusBox.textContent = message;
  statusBox.className = `status ${state}`.trim();
}

function assertLibraries() {
  const missing = [];

  if (typeof L === "undefined") missing.push("Leaflet");
  if (typeof parseGeoraster === "undefined") missing.push("GeoRaster");
  if (typeof GeoRasterLayer === "undefined") missing.push("GeoRasterLayer");

  if (missing.length) {
    throw new Error(`Missing map libraries: ${missing.join(", ")}`);
  }
}

assertLibraries();

const map = L.map("map", {
  worldCopyJump: true,
  minZoom: 2
}).setView([15, 0], 2);

L.tileLayer("https://tile.openstreetmap.org/{z}/{x}/{y}.png", {
  attribution: "© OpenStreetMap contributors",
  maxZoom: 19
}).addTo(map);

function addCountrySearch() {
  if (!L.Control || !L.Control.geocoder) {
    console.warn("Country search plugin did not load.");
    return;
  }

  L.Control.geocoder({
    defaultMarkGeocode: false,
    placeholder: "Search country or place",
    errorMessage: "Location not found",
    collapsed: false,
    position: "topright",
    geocoder: L.Control.Geocoder.nominatim({
      geocodingQueryParams: {
        addressdetails: 1
      }
    })
  })
    .on("markgeocode", event => {
      const result = event.geocode;

      if (result.bbox) {
        map.fitBounds(result.bbox, {
          padding: [20, 20],
          maxZoom: 6
        });
      } else if (result.center) {
        map.setView(result.center, 5);
      }
    })
    .addTo(map);
}

function addGlobalViewButton() {
  const GlobalViewControl = L.Control.extend({
    options: { position: "topright" },

    onAdd() {
      const container = L.DomUtil.create(
        "div",
        "leaflet-bar leaflet-control"
      );

      const button = L.DomUtil.create(
        "a",
        "global-view-button",
        container
      );

      button.href = "#";
      button.title = "Return to global view";
      button.setAttribute("aria-label", "Return to global view");
      button.textContent = "🌍";

      L.DomEvent.disableClickPropagation(container);
      L.DomEvent.on(button, "click", event => {
        L.DomEvent.preventDefault(event);
        map.setView([15, 0], 2);
      });

      return container;
    }
  });

  map.addControl(new GlobalViewControl());
}

addCountrySearch();
addGlobalViewButton();

function interpolateColor(stops, value) {
  let lower = stops[0];
  let upper = stops[stops.length - 1];

  for (let index = 1; index < stops.length; index += 1) {
    if (value <= stops[index][0]) {
      lower = stops[index - 1];
      upper = stops[index];
      break;
    }
  }

  const fraction = (value - lower[0]) / (upper[0] - lower[0] || 1);
  const parseHex = color => parseInt(color.slice(1), 16);
  const lowerRgb = parseHex(lower[1]);
  const upperRgb = parseHex(upper[1]);

  return "#" + [16, 8, 0].map(shift => {
    const start = (lowerRgb >> shift) & 255;
    const end = (upperRgb >> shift) & 255;
    return Math.round(start * (1 - fraction) + end * fraction)
      .toString(16)
      .padStart(2, "0");
  }).join("");
}

function scenarioLabel(value) {
  return SCENARIO_LABELS[value] || value;
}

function periodLabel(value) {
  return PERIOD_LABELS[value] || value;
}

function populateSelectors() {
  const scenarios = [...new Set(catalog.layers.map(item => item.ssp))];
  const periods = [...new Set(catalog.layers.map(item => item.period))];

  sspSelect.innerHTML = scenarios
    .map(value => `<option value="${value}">${scenarioLabel(value)}</option>`)
    .join("");

  periodSelect.innerHTML = periods
    .map(value => `<option value="${value}">${periodLabel(value)}</option>`)
    .join("");

  updateLayerOptions();
}

function selectedLayers() {
  return catalog.layers.filter(item =>
    item.ssp === sspSelect.value && item.period === periodSelect.value
  );
}

function updateLayerOptions() {
  const layers = selectedLayers();

  layerSelect.innerHTML = layers
    .map((item, index) => (
      `<option value="${index}">${item.title}</option>`
    ))
    .join("");

  if (!layers.length) {
    legend.innerHTML = "<b>No layers match this scenario and period.</b>";
    setStatus("No layers available for this selection.", "error");

    if (currentLayer) {
      map.removeLayer(currentLayer);
      currentLayer = null;
    }

    return;
  }

  renderSelectedLayer();
}

async function loadGeoRaster(url) {
  const response = await fetch(url);

  if (!response.ok) {
    throw new Error(
      `Raster request failed: ${response.status} ${response.statusText} for ${url}`
    );
  }

  const arrayBuffer = await response.arrayBuffer();
  return parseGeoraster(arrayBuffer);
}

async function renderSelectedLayer() {
  const layers = selectedLayers();
  const metadata = layers[Number(layerSelect.value)];

  if (!metadata) return;

  if (currentLayer) {
    map.removeLayer(currentLayer);
    currentLayer = null;
  }

  setStatus(`Loading ${metadata.title}…`);

  try {
    const georaster = await loadGeoRaster(metadata.file);

    currentLayer = new GeoRasterLayer({
      georaster,
      opacity: Number(opacityInput.value),
      resolution: 256,
      resampleMethod: metadata.kind === "categorical" ? "nearest" : "bilinear",
      pixelValuesToColorFn: values => valueToColor(metadata, values[0])
    });

    currentLayer.addTo(map);
    legend.innerHTML = legendHtml(metadata);

    setStatus(
      `${scenarioLabel(metadata.ssp)} | ${periodLabel(metadata.period)} | ${metadata.title}`,
      "ready"
    );
  } catch (error) {
    console.error("Could not render raster layer", metadata, error);
    setStatus("Layer failed to load. See the message below.", "error");
    legend.innerHTML = `
      <b>Layer failed to load.</b>
      <p>${error.message}</p>
      <p>Confirm that <code>${metadata.file}</code> exists under <code>docs/</code>.</p>
    `;
  }
}

function valueToColor(metadata, value) {
  if (
    value === undefined ||
    value === null ||
    Number.isNaN(value) ||
    value === 255 ||
    value === -9999
  ) {
    return null;
  }

  if (metadata.kind === "categorical") {
    return catalog.colors[String(Math.round(value))] || null;
  }

  const minimum = Number(metadata.min);
  const maximum = Number(metadata.max);
  const scaled = Math.max(
    0,
    Math.min(1, (value - minimum) / (maximum - minimum))
  );

  return interpolateColor(
    palettes[metadata.palette] || palettes.viridis,
    scaled
  );
}

function legendHtml(metadata) {
  if (metadata.kind === "categorical") {
    return "<h3>Priority classes</h3>" + Object.entries(catalog.classes)
      .map(([code, label]) => `
        <div class="legend-row">
          <span class="swatch" style="background:${catalog.colors[code]}"></span>
          <span>${code}: ${label}</span>
        </div>
      `)
      .join("");
  }

  const palette = palettes[metadata.palette] || palettes.viridis;

  return `
    <h3>${metadata.title}</h3>
    <div>${metadata.min}<span style="float:right">${metadata.max}</span></div>
    <div style="height:14px;background:linear-gradient(90deg,${palette.map(item => item[1]).join(",")})"></div>
  `;
}

sspSelect.addEventListener("change", updateLayerOptions);
periodSelect.addEventListener("change", updateLayerOptions);
layerSelect.addEventListener("change", renderSelectedLayer);
opacityInput.addEventListener("input", () => {
  if (currentLayer) {
    currentLayer.setOpacity(Number(opacityInput.value));
  }
});

fetch("catalog.json")
  .then(response => {
    if (!response.ok) {
      throw new Error(
        `catalog.json request failed: ${response.status} ${response.statusText}`
      );
    }

    return response.json();
  })
  .then(data => {
    catalog = data;

    if (!Array.isArray(catalog.layers) || catalog.layers.length === 0) {
      throw new Error(
        "catalog.json contains no layers. Run scripts/05_build_web_repository.py."
      );
    }

    populateSelectors();
  })
  .catch(error => {
    console.error("Could not initialise viewer", error);
    setStatus("Viewer initialisation failed.", "error");
    legend.innerHTML = `<b>Viewer initialisation failed.</b><p>${error.message}</p>`;
  });
