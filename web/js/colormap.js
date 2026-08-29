/*
 * Colour maps.
 *
 * The stability plots read a signed quantity, log |R(z)|, whose zero is the
 * only value that matters. A diverging map centred on that zero puts the
 * stability boundary exactly where the colours turn, so the picture carries the
 * information even before the contour line is drawn on top.
 */

/** ColorBrewer RdYlBu, from red through pale yellow to blue. */
const RD_YL_BU = [
  [165, 0, 38],
  [215, 48, 39],
  [244, 109, 67],
  [253, 174, 97],
  [254, 224, 144],
  [255, 255, 191],
  [224, 243, 248],
  [171, 217, 233],
  [116, 173, 209],
  [69, 117, 180],
  [49, 54, 149],
];

/** Sequential blue map for magnitudes that have no meaningful zero. */
const BLUES = [
  [12, 20, 32],
  [20, 44, 74],
  [26, 72, 112],
  [33, 102, 148],
  [56, 135, 178],
  [96, 168, 200],
  [148, 199, 219],
  [200, 227, 238],
];

function interpolate(stops, t) {
  const clamped = Math.min(1, Math.max(0, t));
  const scaled = clamped * (stops.length - 1);
  const index = Math.min(stops.length - 2, Math.floor(scaled));
  const frac = scaled - index;
  const a = stops[index];
  const b = stops[index + 1];
  return [
    Math.round(a[0] + (b[0] - a[0]) * frac),
    Math.round(a[1] + (b[1] - a[1]) * frac),
    Math.round(a[2] + (b[2] - a[2]) * frac),
  ];
}

/**
 * Diverging map for log magnitudes: blue where the mode decays, red where it
 * grows. `t` runs from zero at the low end to one at the high end.
 */
export function diverging(t) {
  return interpolate([...RD_YL_BU].reverse(), t);
}

export function sequential(t) {
  return interpolate(BLUES, t);
}

export function rgb(color) {
  return `rgb(${color[0]} ${color[1]} ${color[2]})`;
}

/**
 * Quantize a value into one of `levels` bands and return its colour.
 *
 * Banding rather than a smooth ramp is deliberate: discrete steps let the eye
 * read a value off the colour bar, which a continuous gradient does not.
 */
export function band(value, low, high, levels, map = diverging) {
  if (!Number.isFinite(value)) return map(1);
  const span = high - low;
  if (span <= 0) return map(0.5);
  const normalized = (value - low) / span;
  const index = Math.min(levels - 1, Math.max(0, Math.floor(normalized * levels)));
  return map((index + 0.5) / levels);
}

/**
 * Level boundaries for a diverging map that must place a band edge exactly at
 * zero, so the stability boundary coincides with a colour change.
 */
export function symmetricLevels(low, high, count) {
  const reach = Math.max(Math.abs(low), Math.abs(high));
  const step = (2 * reach) / count;
  const levels = [];
  for (let i = 0; i <= count; i += 1) {
    levels.push(-reach + i * step);
  }
  return levels;
}

/** Colours for the series in a line plot, taken from the token set. */
export function seriesColor(index) {
  const style = getComputedStyle(document.documentElement);
  const name = `--color-series-${(index % 6) + 1}`;
  return style.getPropertyValue(name).trim() || "#6fb2e0";
}
