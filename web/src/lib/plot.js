/*
 * Plotly, bundled down to the four trace types this interface uses, and themed
 * once so every figure in the app looks like it came from the same place.
 *
 * The stability plots read a signed quantity, log |R(z)|, whose zero is the
 * only value that matters. The colour scale is therefore diverging, banded
 * rather than smooth so a value can be read off the bar, and symmetric so the
 * turn in the colours lands exactly on the stability boundary.
 */

import Plotly from 'plotly.js/lib/core';
import heatmap from 'plotly.js/lib/heatmap';
import contour from 'plotly.js/lib/contour';
import scattergl from 'plotly.js/lib/scattergl';
import bar from 'plotly.js/lib/bar';

Plotly.register([heatmap, contour, scattergl, bar]);

export default Plotly;

export const COLORS = {
	charcoal: '#0f0f0f',
	surface: '#1a1a18',
	cream: '#f0efe9',
	accent: '#969591',
	line: 'rgba(240, 239, 233, 0.10)',
	zero: 'rgba(240, 239, 233, 0.22)'
};

/** Series colours, taken from the same palette as the rest of the site. */
export const SERIES = [
	'#00d9c0',
	'#f5a623',
	'#3399d6',
	'#818cf8',
	'#d9513c',
	'#d4d3cd',
	'#00b3a0',
	'#f7b84d'
];

export function seriesColor(index) {
	return SERIES[index % SERIES.length];
}

/** ColorBrewer RdYlBu, blue at the low end where a mode decays. */
const RD_YL_BU = [
	'#313695',
	'#4575b4',
	'#74add1',
	'#abd9e9',
	'#e0f3f8',
	'#ffffbf',
	'#fee090',
	'#fdae61',
	'#f46d43',
	'#d73027',
	'#a50026'
];

function mix(a, b, t) {
	const parse = (hex) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));
	const [r1, g1, b1] = parse(a);
	const [r2, g2, b2] = parse(b);
	const channel = (x, y) => Math.round(x + (y - x) * t);
	return `rgb(${channel(r1, r2)},${channel(g1, g2)},${channel(b1, b2)})`;
}

function sample(stops, t) {
	const scaled = Math.min(1, Math.max(0, t)) * (stops.length - 1);
	const index = Math.min(stops.length - 2, Math.floor(scaled));
	return mix(stops[index], stops[index + 1], scaled - index);
}

/**
 * A banded colour scale: each band is one flat colour, so the plot reads as
 * discrete levels the way a contour plot does rather than as a gradient.
 */
export function bandedScale(levels = 12, stops = RD_YL_BU) {
	const scale = [];
	for (let i = 0; i < levels; i += 1) {
		const color = sample(stops, (i + 0.5) / levels);
		scale.push([i / levels, color], [(i + 1) / levels, color]);
	}
	return scale;
}

const AXIS = {
	color: COLORS.accent,
	gridcolor: COLORS.line,
	linecolor: COLORS.line,
	zerolinecolor: COLORS.zero,
	tickfont: { family: 'JetBrains Mono, monospace', size: 9 },
	titlefont: { family: 'JetBrains Mono, monospace', size: 10 },
	ticks: 'outside',
	ticklen: 3,
	// Powers of ten rather than SI prefixes: an error axis reads in decades.
	exponentformat: 'power',
	tickcolor: COLORS.line,
	showline: true,
	mirror: true,
	automargin: false
};

/**
 * Layout shared by every figure. `compact` strips it down to the data alone,
 * which is what a card wants.
 *
 * A figure carries no title of its own: the heading above it is a heading like
 * any other on the page, set in the same type as the rest of the interface
 * rather than in whatever the plotting library draws.
 */
/**
 * `equalScale` ties one unit on the vertical axis to one unit on the horizontal
 * one. Without it a window computed for one aspect ratio and drawn in a
 * container of another stretches the picture, and a stability region is a shape
 * whose proportions carry the meaning. With it the frame shows whatever extra
 * plane it has room for instead.
 */
export function layout({ compact = false, equalScale = false, x = {}, y = {} } = {}) {
	return {
		paper_bgcolor: COLORS.surface,
		plot_bgcolor: COLORS.charcoal,
		font: { family: 'JetBrains Mono, monospace', size: 10, color: COLORS.accent },
		margin: compact ? { l: 1, r: 1, t: 1, b: 1 } : { l: 62, r: 18, t: 10, b: 44 },
		xaxis: {
			...AXIS,
			showticklabels: !compact,
			showgrid: !compact,
			ticks: compact ? '' : 'outside',
			...x
		},
		yaxis: {
			...AXIS,
			showticklabels: !compact,
			showgrid: !compact,
			ticks: compact ? '' : 'outside',
			...(equalScale ? { scaleanchor: 'x', scaleratio: 1 } : {}),
			...y
		},
		showlegend: false,
		hovermode: compact ? false : 'closest',
		legend: {
			bgcolor: 'rgba(15,15,15,0.85)',
			bordercolor: COLORS.line,
			borderwidth: 1,
			font: { family: 'JetBrains Mono, monospace', size: 9, color: COLORS.cream }
		}
	};
}

export function config({ compact = false } = {}) {
	return {
		displayModeBar: false,
		responsive: true,
		staticPlot: compact,
		scrollZoom: !compact,
		// Suppresses the "double click to zoom back out" toast Plotly shows on
		// the first interaction with a figure.
		showTips: false
	};
}

/**
 * The unit circle, as a reference rather than as data: drawn exactly like the
 * axes, because it is the same kind of mark. Both axes of these figures carry
 * one scale, so it is a circle on screen and not an ellipse.
 */
export function unitCircle({ color = '#000000' } = {}) {
	return {
		type: 'circle',
		xref: 'x',
		yref: 'y',
		layer: 'above',
		x0: -1,
		y0: -1,
		x1: 1,
		y1: 1,
		line: { color, width: 1, dash: 'dot' }
	};
}

/**
 * The two axes through the origin, drawn over the data.
 *
 * Plotly puts its own zero lines under the traces, where a filled contour or a
 * heatmap hides them completely. These are shapes on the upper layer instead,
 * which is the only way to keep the reference visible on the plots that need it
 * most.
 */
export function originLines({ x = true, y = true, compact = false, color = '#000000' } = {}) {
	const line = { color, width: 1, dash: 'dot' };
	const shapes = [];
	if (y) {
		shapes.push({
			type: 'line',
			layer: 'above',
			xref: 'paper',
			x0: 0,
			x1: 1,
			yref: 'y',
			y0: 0,
			y1: 0,
			opacity: compact ? 0.5 : 0.65,
			line
		});
	}
	if (x) {
		shapes.push({
			type: 'line',
			layer: 'above',
			yref: 'paper',
			y0: 0,
			y1: 1,
			xref: 'x',
			x0: 0,
			x1: 0,
			opacity: compact ? 0.5 : 0.65,
			line
		});
	}
	return shapes;
}

/** Axis range on a log scale, in the exponents Plotly wants. */
export function logRange(values, pad = 0.06) {
	const finite = values.filter((v) => Number.isFinite(v) && v > 0);
	if (finite.length === 0) return [-16, 0];
	const low = Math.log10(Math.min(...finite));
	const high = Math.log10(Math.max(...finite));
	const span = Math.max(high - low, 0.5);
	return [low - span * pad, high + span * pad];
}

export function linearRange(values, pad = 0.05) {
	const finite = values.filter((v) => Number.isFinite(v));
	if (finite.length === 0) return [0, 1];
	const low = Math.min(...finite);
	const high = Math.max(...finite);
	const span = Math.max(high - low, 1e-12);
	return [low - span * pad, high + span * pad];
}
