/*
 * A small canvas plotting layer.
 *
 * Written rather than pulled in, because the one plot that matters here is a
 * banded contour field with an explicit unit contour drawn over it, and a
 * general charting library gets in the way of that more than it helps. What is
 * needed is narrow: linear and logarithmic axes, a filled contour field, a
 * marching squares contour line, a colour bar, and line series.
 */

import { band, diverging, rgb, seriesColor } from "./colormap.js";

const MARGIN = { top: 26, right: 14, bottom: 42, left: 54 };
const COLORBAR_WIDTH = 12;
const COLORBAR_GAP = 8;
const COLORBAR_LABELS = 46;

function token(name, fallback) {
  const value = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return value || fallback;
}

/** A canvas with axes, scales and a plotting area. */
export class Figure {
  constructor(canvas, options = {}) {
    this.canvas = canvas;
    this.context = canvas.getContext("2d");
    this.title = options.title ?? "";
    this.xLabel = options.xLabel ?? "";
    this.yLabel = options.yLabel ?? "";
    this.xScale = options.xScale ?? "linear";
    this.yScale = options.yScale ?? "linear";
    this.colorbar = options.colorbar ?? null;
    this.domain = { x0: 0, x1: 1, y0: 0, y1: 1 };
    this.resize();
  }

  resize() {
    const ratio = window.devicePixelRatio || 1;
    const rect = this.canvas.getBoundingClientRect();
    const width = Math.max(160, Math.round(rect.width));
    const height = Math.max(120, Math.round(rect.height));
    this.canvas.width = Math.round(width * ratio);
    this.canvas.height = Math.round(height * ratio);
    this.context.setTransform(ratio, 0, 0, ratio, 0, 0);
    this.width = width;
    this.height = height;

    const right = this.colorbar ? MARGIN.right + COLORBAR_WIDTH + COLORBAR_GAP + COLORBAR_LABELS : MARGIN.right;
    this.plot = {
      left: MARGIN.left,
      top: MARGIN.top,
      right: width - right,
      bottom: height - MARGIN.bottom,
    };
    this.plot.width = this.plot.right - this.plot.left;
    this.plot.height = this.plot.bottom - this.plot.top;
  }

  setDomain(x0, x1, y0, y1) {
    this.domain = { x0, x1, y0, y1 };
  }

  toPixelX(x) {
    const { x0, x1 } = this.domain;
    const value = this.xScale === "log" ? Math.log10(x) : x;
    const low = this.xScale === "log" ? Math.log10(x0) : x0;
    const high = this.xScale === "log" ? Math.log10(x1) : x1;
    return this.plot.left + ((value - low) / (high - low)) * this.plot.width;
  }

  toPixelY(y) {
    const { y0, y1 } = this.domain;
    const value = this.yScale === "log" ? Math.log10(y) : y;
    const low = this.yScale === "log" ? Math.log10(y0) : y0;
    const high = this.yScale === "log" ? Math.log10(y1) : y1;
    return this.plot.bottom - ((value - low) / (high - low)) * this.plot.height;
  }

  clear() {
    const context = this.context;
    context.clearRect(0, 0, this.width, this.height);
    context.fillStyle = token("--color-surface", "#1e2023");
    context.fillRect(0, 0, this.width, this.height);
  }

  /** Frame, ticks, labels and the dotted axes through the origin. */
  drawAxes(options = {}) {
    const context = this.context;
    const border = token("--color-border-strong", "#474c53");
    const text = token("--color-text-muted", "#9aa1a9");
    const dim = token("--color-text-dim", "#6b7279");
    const font = token("--font-sans", "sans-serif");

    context.save();
    context.strokeStyle = border;
    context.lineWidth = 1;
    context.strokeRect(this.plot.left + 0.5, this.plot.top + 0.5, this.plot.width, this.plot.height);

    context.fillStyle = text;
    context.font = `10px ${font}`;
    context.textAlign = "center";
    context.textBaseline = "top";

    const xTicks = options.xTicks ?? this.ticks("x");
    for (const value of xTicks) {
      const px = this.toPixelX(value);
      if (px < this.plot.left - 1 || px > this.plot.right + 1) continue;
      context.beginPath();
      context.moveTo(px, this.plot.bottom);
      context.lineTo(px, this.plot.bottom + 4);
      context.stroke();
      context.fillText(formatTick(value, this.xScale), px, this.plot.bottom + 7);
    }

    context.textAlign = "right";
    context.textBaseline = "middle";
    const yTicks = options.yTicks ?? this.ticks("y");
    for (const value of yTicks) {
      const py = this.toPixelY(value);
      if (py < this.plot.top - 1 || py > this.plot.bottom + 1) continue;
      context.beginPath();
      context.moveTo(this.plot.left - 4, py);
      context.lineTo(this.plot.left, py);
      context.stroke();
      context.fillText(formatTick(value, this.yScale), this.plot.left - 7, py);
    }

    // Dotted lines through the origin, the reference every stability plot needs.
    if (options.origin !== false && this.xScale === "linear" && this.yScale === "linear") {
      context.save();
      context.strokeStyle = dim;
      context.setLineDash([1, 3]);
      const { x0, x1, y0, y1 } = this.domain;
      if (x0 < 0 && x1 > 0) {
        const px = this.toPixelX(0);
        context.beginPath();
        context.moveTo(px, this.plot.top);
        context.lineTo(px, this.plot.bottom);
        context.stroke();
      }
      if (y0 < 0 && y1 > 0) {
        const py = this.toPixelY(0);
        context.beginPath();
        context.moveTo(this.plot.left, py);
        context.lineTo(this.plot.right, py);
        context.stroke();
      }
      context.restore();
    }

    if (this.title) {
      context.fillStyle = token("--color-text", "#e6e8ea");
      context.font = `600 11px ${font}`;
      context.textAlign = "center";
      context.textBaseline = "alphabetic";
      context.fillText(this.title, (this.plot.left + this.plot.right) / 2, MARGIN.top - 10);
    }
    if (this.xLabel) {
      context.fillStyle = text;
      context.font = `10px ${font}`;
      context.textAlign = "center";
      context.textBaseline = "alphabetic";
      context.fillText(this.xLabel, (this.plot.left + this.plot.right) / 2, this.height - 6);
    }
    if (this.yLabel) {
      context.save();
      context.fillStyle = text;
      context.font = `10px ${font}`;
      context.translate(11, (this.plot.top + this.plot.bottom) / 2);
      context.rotate(-Math.PI / 2);
      context.textAlign = "center";
      context.textBaseline = "alphabetic";
      context.fillText(this.yLabel, 0, 0);
      context.restore();
    }
    context.restore();
  }

  ticks(axis) {
    const scale = axis === "x" ? this.xScale : this.yScale;
    const low = axis === "x" ? this.domain.x0 : this.domain.y0;
    const high = axis === "x" ? this.domain.x1 : this.domain.y1;
    return scale === "log" ? decadeTicks(low, high) : niceTicks(low, high, 6);
  }

  /**
   * Filled contour bands from a row major grid, rows running bottom to top.
   *
   * The grid is resampled bilinearly to pixel resolution before banding, so the
   * band edges land where the data says and not on the grid raster.
   */
  drawField(data, gridWidth, gridHeight, options = {}) {
    const levels = options.levels ?? 14;
    const low = options.low ?? -3;
    const high = options.high ?? 3;
    const map = options.map ?? diverging;

    const context = this.context;
    const width = Math.round(this.plot.width);
    const height = Math.round(this.plot.height);
    const image = context.createImageData(width, height);

    for (let py = 0; py < height; py += 1) {
      // Canvas rows run downward, the grid runs upward.
      const gy = ((height - 1 - py) / Math.max(1, height - 1)) * (gridHeight - 1);
      const y0 = Math.min(gridHeight - 1, Math.floor(gy));
      const y1 = Math.min(gridHeight - 1, y0 + 1);
      const fy = gy - y0;
      for (let px = 0; px < width; px += 1) {
        const gx = (px / Math.max(1, width - 1)) * (gridWidth - 1);
        const x0 = Math.min(gridWidth - 1, Math.floor(gx));
        const x1 = Math.min(gridWidth - 1, x0 + 1);
        const fx = gx - x0;

        const v00 = data[y0 * gridWidth + x0];
        const v10 = data[y0 * gridWidth + x1];
        const v01 = data[y1 * gridWidth + x0];
        const v11 = data[y1 * gridWidth + x1];
        const value =
          v00 * (1 - fx) * (1 - fy) + v10 * fx * (1 - fy) + v01 * (1 - fx) * fy + v11 * fx * fy;

        const color = band(value, low, high, levels, map);
        const offset = (py * width + px) * 4;
        image.data[offset] = color[0];
        image.data[offset + 1] = color[1];
        image.data[offset + 2] = color[2];
        image.data[offset + 3] = 255;
      }
    }
    context.putImageData(image, this.plot.left, this.plot.top);
  }

  /**
   * Marching squares contour at one level, in grid coordinates mapped through
   * the domain. This is the line that says where stability actually ends.
   */
  drawContour(data, gridWidth, gridHeight, level, options = {}) {
    const context = this.context;
    context.save();
    context.beginPath();
    context.strokeStyle = options.color ?? "#000000";
    context.lineWidth = options.lineWidth ?? 1.4;
    context.lineJoin = "round";

    const { x0, x1, y0, y1 } = this.domain;
    const toX = (gx) => this.toPixelX(x0 + ((x1 - x0) * gx) / (gridWidth - 1));
    const toY = (gy) => this.toPixelY(y0 + ((y1 - y0) * gy) / (gridHeight - 1));

    const at = (ix, iy) => data[iy * gridWidth + ix];
    const cross = (a, b) => (level - a) / (b - a);

    for (let iy = 0; iy < gridHeight - 1; iy += 1) {
      for (let ix = 0; ix < gridWidth - 1; ix += 1) {
        const v0 = at(ix, iy);
        const v1 = at(ix + 1, iy);
        const v2 = at(ix + 1, iy + 1);
        const v3 = at(ix, iy + 1);
        if (![v0, v1, v2, v3].every(Number.isFinite)) continue;

        let index = 0;
        if (v0 > level) index |= 1;
        if (v1 > level) index |= 2;
        if (v2 > level) index |= 4;
        if (v3 > level) index |= 8;
        if (index === 0 || index === 15) continue;

        const bottom = () => [ix + cross(v0, v1), iy];
        const right = () => [ix + 1, iy + cross(v1, v2)];
        const top = () => [ix + cross(v3, v2), iy + 1];
        const left = () => [ix, iy + cross(v0, v3)];

        const segments = [];
        switch (index) {
          case 1:
          case 14:
            segments.push([left(), bottom()]);
            break;
          case 2:
          case 13:
            segments.push([bottom(), right()]);
            break;
          case 3:
          case 12:
            segments.push([left(), right()]);
            break;
          case 4:
          case 11:
            segments.push([right(), top()]);
            break;
          case 5:
            segments.push([left(), top()], [bottom(), right()]);
            break;
          case 6:
          case 9:
            segments.push([bottom(), top()]);
            break;
          case 7:
          case 8:
            segments.push([left(), top()]);
            break;
          case 10:
            segments.push([left(), bottom()], [right(), top()]);
            break;
          default:
            break;
        }
        for (const [a, b] of segments) {
          context.moveTo(toX(a[0]), toY(a[1]));
          context.lineTo(toX(b[0]), toY(b[1]));
        }
      }
    }
    context.stroke();
    context.restore();
  }

  /** The colour bar for a field, with its own ticks. */
  drawColorbar(low, high, options = {}) {
    if (!this.colorbar) return;
    const context = this.context;
    const levels = options.levels ?? 14;
    const map = options.map ?? diverging;
    const left = this.plot.right + COLORBAR_GAP;
    const top = this.plot.top;
    const height = this.plot.height;

    for (let i = 0; i < levels; i += 1) {
      const y = top + height - ((i + 1) * height) / levels;
      const color = map((i + 0.5) / levels);
      context.fillStyle = rgb(color);
      context.fillRect(left, y, COLORBAR_WIDTH, Math.ceil(height / levels) + 0.5);
    }
    context.strokeStyle = token("--color-border-strong", "#474c53");
    context.lineWidth = 1;
    context.strokeRect(left + 0.5, top + 0.5, COLORBAR_WIDTH, height);

    context.fillStyle = token("--color-text-muted", "#9aa1a9");
    context.font = `10px ${token("--font-sans", "sans-serif")}`;
    context.textAlign = "left";
    context.textBaseline = "middle";
    const ticks = niceTicks(low, high, 5);
    for (const value of ticks) {
      const y = top + height - ((value - low) / (high - low)) * height;
      if (y < top - 1 || y > top + height + 1) continue;
      context.beginPath();
      context.moveTo(left + COLORBAR_WIDTH, y);
      context.lineTo(left + COLORBAR_WIDTH + 3, y);
      context.stroke();
      context.fillText(formatTick(value, "linear"), left + COLORBAR_WIDTH + 6, y);
    }

    if (this.colorbar.label) {
      context.save();
      context.translate(left + COLORBAR_WIDTH + COLORBAR_LABELS - 4, top + height / 2);
      context.rotate(-Math.PI / 2);
      context.textAlign = "center";
      context.textBaseline = "alphabetic";
      context.fillText(this.colorbar.label, 0, 0);
      context.restore();
    }
  }

  /** Line series, each `{ x: [], y: [], label, color, marker }`. */
  drawSeries(series, options = {}) {
    const context = this.context;
    context.save();
    context.beginPath();
    context.rect(this.plot.left, this.plot.top, this.plot.width, this.plot.height);
    context.clip();

    series.forEach((entry, index) => {
      const color = entry.color ?? seriesColor(index);
      context.strokeStyle = color;
      context.fillStyle = color;
      context.lineWidth = entry.lineWidth ?? 1.6;
      context.beginPath();
      let started = false;
      for (let i = 0; i < entry.x.length; i += 1) {
        const x = entry.x[i];
        const y = entry.y[i];
        if (!Number.isFinite(x) || !Number.isFinite(y)) {
          started = false;
          continue;
        }
        if (this.xScale === "log" && x <= 0) continue;
        if (this.yScale === "log" && y <= 0) continue;
        const px = this.toPixelX(x);
        const py = this.toPixelY(y);
        if (started) context.lineTo(px, py);
        else {
          context.moveTo(px, py);
          started = true;
        }
      }
      context.stroke();

      if (entry.marker !== false) {
        for (let i = 0; i < entry.x.length; i += 1) {
          const x = entry.x[i];
          const y = entry.y[i];
          if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
          if (this.xScale === "log" && x <= 0) continue;
          if (this.yScale === "log" && y <= 0) continue;
          context.beginPath();
          context.arc(this.toPixelX(x), this.toPixelY(y), 2.4, 0, Math.PI * 2);
          context.fill();
        }
      }
    });
    context.restore();

    if (options.legend !== false && series.some((s) => s.label)) {
      this.drawLegend(series);
    }
  }

  drawLegend(series) {
    const context = this.context;
    const font = token("--font-sans", "sans-serif");
    context.save();
    context.font = `10px ${font}`;
    context.textAlign = "left";
    context.textBaseline = "middle";
    const entries = series.filter((s) => s.label);
    const widths = entries.map((s) => context.measureText(s.label).width);
    const boxWidth = Math.max(...widths, 0) + 26;
    const boxHeight = entries.length * 14 + 8;
    const x = this.plot.right - boxWidth - 8;
    const y = this.plot.top + 8;

    context.fillStyle = token("--color-bg", "#17181a");
    context.globalAlpha = 0.82;
    context.fillRect(x, y, boxWidth, boxHeight);
    context.globalAlpha = 1;
    context.strokeStyle = token("--color-border", "#33373c");
    context.strokeRect(x + 0.5, y + 0.5, boxWidth, boxHeight);

    entries.forEach((entry, index) => {
      const color = entry.color ?? seriesColor(index);
      const cy = y + 11 + index * 14;
      context.strokeStyle = color;
      context.lineWidth = 2;
      context.beginPath();
      context.moveTo(x + 7, cy);
      context.lineTo(x + 19, cy);
      context.stroke();
      context.fillStyle = token("--color-text-muted", "#9aa1a9");
      context.fillText(entry.label, x + 23, cy);
    });
    context.restore();
  }

  /** A guide line of given slope, for reading convergence order off a plot. */
  drawSlopeGuide(slope, anchor, options = {}) {
    if (this.xScale !== "log" || this.yScale !== "log") return;
    const context = this.context;
    const { x0, x1 } = this.domain;
    const y = (x) => anchor.y * Math.pow(x / anchor.x, slope);
    context.save();
    context.strokeStyle = options.color ?? token("--color-text-dim", "#6b7279");
    context.setLineDash([4, 4]);
    context.lineWidth = 1;
    context.beginPath();
    context.moveTo(this.toPixelX(x0), this.toPixelY(y(x0)));
    context.lineTo(this.toPixelX(x1), this.toPixelY(y(x1)));
    context.stroke();
    context.restore();
  }
}

function formatTick(value, scale) {
  if (scale === "log") {
    const exponent = Math.round(Math.log10(value));
    return `1e${exponent}`;
  }
  if (value === 0) return "0";
  const magnitude = Math.abs(value);
  if (magnitude >= 1e4 || magnitude < 1e-3) return value.toExponential(0);
  const digits = magnitude >= 100 ? 0 : magnitude >= 10 ? 1 : magnitude >= 1 ? 1 : 2;
  return value.toFixed(digits);
}

function niceTicks(low, high, target) {
  if (!(high > low)) return [low];
  const raw = (high - low) / target;
  const magnitude = Math.pow(10, Math.floor(Math.log10(raw)));
  const normalized = raw / magnitude;
  const step = (normalized >= 5 ? 5 : normalized >= 2 ? 2 : 1) * magnitude;
  const first = Math.ceil(low / step) * step;
  const ticks = [];
  for (let value = first; value <= high + step * 1e-9; value += step) {
    ticks.push(Math.abs(value) < step * 1e-9 ? 0 : value);
  }
  return ticks;
}

function decadeTicks(low, high) {
  if (!(low > 0) || !(high > low)) return [];
  const first = Math.ceil(Math.log10(low));
  const last = Math.floor(Math.log10(high));
  const ticks = [];
  for (let exponent = first; exponent <= last; exponent += 1) {
    ticks.push(Math.pow(10, exponent));
  }
  return ticks;
}

/** Extent of a set of values, padded, for setting a plot domain. */
export function extent(values, { pad = 0.05, log = false } = {}) {
  const finite = values.filter((v) => Number.isFinite(v) && (!log || v > 0));
  if (finite.length === 0) return log ? [1e-16, 1] : [0, 1];
  let low = Math.min(...finite);
  let high = Math.max(...finite);
  if (log) {
    const l = Math.log10(low);
    const h = Math.log10(high);
    const span = Math.max(h - l, 0.5);
    return [Math.pow(10, l - span * pad), Math.pow(10, h + span * pad)];
  }
  const span = Math.max(high - low, 1e-12);
  low -= span * pad;
  high += span * pad;
  return [low, high];
}
