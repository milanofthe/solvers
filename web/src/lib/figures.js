/*
 * What a card draws.
 *
 * A card in the library is a small window onto one method, and which window it
 * is should be the reader's choice: the same grid answers "which of these is
 * stable here", "which of these is cheap", and "which of these actually
 * converges at the order it claims", depending on the mode.
 *
 * Each mode declares what it needs from the engine and turns the answer into a
 * figure. Nothing here computes; the engine does, on its own schedule.
 */

import {
	bandedScale,
	COLORS,
	config,
	layout,
	logRange,
	linearRange,
	originLines,
	seriesColor
} from './plot.js';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** The flat row major grid the core returns, as the nested array Plotly wants. */
function reshape(data, width, height) {
	const rows = new Array(height);
	for (let r = 0; r < height; r += 1) {
		rows[r] = Array.from(data.subarray(r * width, (r + 1) * width));
	}
	return rows;
}

function axisValues(low, high, count) {
	const values = new Array(count);
	for (let i = 0; i < count; i += 1) values[i] = low + ((high - low) * i) / (count - 1);
	return values;
}

/**
 * The panel is wider than it is tall, and the window is cut to that ratio so a
 * region fills the frame without a shape being stretched out of true. Where the
 * window sits is not guessed here: the core probes the method and reports back
 * a box that actually contains its region.
 */
const PANEL_ASPECT = 1.7;

/** The order star lives around the origin and reaches about as far as the order. */
function orderStarWindow(method) {
	const reach = Math.max(3.5, (method.order ?? 4) * 1.1);
	const half = reach / PANEL_ASPECT;
	return { re: [-reach, reach], im: [-half, half] };
}

const GRID = { runge_kutta: 140, linear_multistep: 72 };

export function empty(message) {
	return {
		data: [],
		layout: {
			...layout({ compact: true }),
			annotations: [
				{
					text: message,
					showarrow: false,
					font: { family: 'JetBrains Mono, monospace', size: 9, color: '#969591' },
					xref: 'paper',
					yref: 'paper',
					x: 0.5,
					y: 0.5
				}
			]
		},
		config: config({ compact: true })
	};
}

// ---------------------------------------------------------------------------
// Method modes
// ---------------------------------------------------------------------------

const stability = {
	key: 'stability',
	label: 'stability',
	note: 'Banded log |R(z)| on the complex plane with the unit contour on top. The region it encloses is where the method does not amplify.',
	request(engine, method, priority) {
		const resolution = GRID[method.class] ?? 96;
		return engine.request(
			'stabilityGridAuto',
			{ id: method.id, aspect: PANEL_ASPECT, width: resolution, height: resolution },
			{ key: `stability:${method.id}`, priority }
		);
	},
	figure(result, method, { compact = true } = {}) {
		const { data, width, height, re, im } = result;
		// The colour range follows the data rather than a fixed span. A method
		// whose region is a half plane varies by a fraction of a decade across
		// the whole window, and a fixed range would render it as one flat wash.
		const sample = [];
		for (let i = 0; i < data.length; i += Math.max(1, Math.floor(data.length / 2000))) {
			if (Number.isFinite(data[i])) sample.push(Math.abs(data[i]));
		}
		sample.sort((a, b) => a - b);
		const quantile = sample[Math.floor(sample.length * 0.97)] ?? 1;
		const reach = Math.min(Math.max(quantile, 0.3), 4);
		const z = reshape(data, width, height);
		const x = axisValues(re[0], re[1], width);
		const y = axisValues(im[0], im[1], height);
		const levels = compact ? 12 : 16;

		return {
			data: [
				{
					type: 'heatmap',
					z,
					x,
					y,
					zmin: -reach,
					zmax: reach,
					colorscale: bandedScale(levels),
					showscale: !compact,
					zsmooth: 'best',
					hoverinfo: 'skip',
					colorbar: compact
						? undefined
						: {
								title: { text: '|R| log', side: 'right', font: { size: 9 } },
								thickness: 10,
								len: 0.9,
								outlinewidth: 1,
								outlinecolor: 'rgba(240,239,233,0.1)',
								tickfont: { size: 9 }
							}
				},
				{
					type: 'contour',
					z,
					x,
					y,
					contours: { start: 0, end: 0, size: 1, coloring: 'none' },
					line: { color: '#000000', width: compact ? 1 : 1.4 },
					showscale: false,
					hoverinfo: 'skip'
				}
			],
			layout: {
				...layout({
					compact,
					x: { title: compact ? undefined : { text: 'Re(z)' }, range: re },
					y: { title: compact ? undefined : { text: 'Im(z)' }, range: im }
				}),
				shapes: originLines({ compact })
			},
			config: config({ compact })
		};
	}
};

const structure = {
	key: 'structure',
	label: 'structure',
	note: 'The coefficient matrix itself. A strictly lower triangle is explicit, a filled diagonal needs one solve per stage, a full matrix couples all of them.',
	request(engine, method, priority) {
		return engine.request(
			'methodDetail',
			{ id: method.id },
			{ key: `detail:${method.id}`, priority }
		);
	},
	figure(detail, method, { compact = true } = {}) {
		const coefficients = detail.coefficients;
		let rows;
		let labels;
		if (coefficients.kind === 'runge_kutta') {
			rows = coefficients.a.map((row) => row.map((entry) => entry.value));
			rows.push(coefficients.b.map((entry) => entry.value));
			labels = coefficients.a.map((_, i) => `stage ${i + 1}`).concat(['b']);
			if (coefficients.bEmbedded) {
				rows.push(coefficients.bEmbedded.map((entry) => entry.value));
				labels.push('b hat');
			}
		} else {
			rows = [coefficients.alpha ?? [], coefficients.beta ?? []];
			labels = ['alpha', 'beta'];
		}

		const columns = Math.max(...rows.map((row) => row.length), 1);
		const padded = rows.map((row) => {
			const filled = new Array(columns).fill(0);
			row.forEach((value, i) => {
				filled[i] = value ?? 0;
			});
			return filled;
		});
		// A signed logarithm, so a coefficient three decades smaller than its
		// neighbours is still visible instead of reading as zero.
		const shaped = padded.map((row) =>
			row.map((v) => (v === 0 ? 0 : Math.sign(v) * Math.log10(1 + Math.abs(v) * 20)))
		);
		const reach = Math.max(...shaped.flat().map(Math.abs), 1);

		return {
			data: [
				{
					type: 'heatmap',
					// Plotly draws the first row at the bottom.
					z: shaped.slice().reverse(),
					y: labels.slice().reverse(),
					zmin: -reach,
					zmax: reach,
					colorscale: bandedScale(9, ['#d9513c', '#8a3a2c', '#1a1a18', '#0a6e63', '#00d9c0']),
					showscale: false,
					xgap: 1,
					ygap: 1,
					hoverinfo: compact ? 'skip' : 'y+z'
				}
			],
			layout: layout({
				compact,
				x: { showticklabels: false, showgrid: false, ticks: '' },
				y: { showgrid: false, ticks: '' }
			}),
			config: config({ compact })
		};
	}
};

const damping = {
	key: 'damping',
	label: 'damping',
	note: 'How far a decaying mode is damped, along the negative real axis. Staying under one is A-stability; falling away at the right is L-stability.',
	request(engine, method, priority) {
		return engine.request(
			'methodSummary',
			{ id: method.id },
			{ key: `summary:${method.id}`, priority }
		);
	},
	figure(summary, method, { compact = true } = {}) {
		if (!summary.damping?.length) return empty('no stability function');
		const x = summary.dampingAxis.map((v) => Math.abs(v));
		const y = summary.damping.map((v) => Math.max(v, 1e-8));
		return {
			data: [
				{
					type: 'scatter',
					mode: 'lines',
					x,
					y,
					line: { color: seriesColor(0), width: 1.6 },
					hoverinfo: compact ? 'skip' : 'x+y'
				}
			],
			layout: {
				...layout({
					compact,
					x: { type: 'log', range: logRange(x, 0.02), title: compact ? undefined : { text: '-z' } },
					y: {
						type: 'log',
						range: [-6, Math.min(2, Math.log10(Math.max(...y)) + 0.2)],
						title: compact ? undefined : { text: '|R(z)|' }
					}
				}),
				shapes: [
					{
						type: 'line',
						xref: 'paper',
						x0: 0,
						x1: 1,
						yref: 'y',
						y0: 1,
						y1: 1,
						line: { color: '#969591', width: 1, dash: 'dot' }
					}
				]
			},
			config: config({ compact })
		};
	}
};

const order = {
	key: 'order',
	label: 'order',
	note: 'Fixed step error against step size. The dotted line is the slope the coefficients promise; the measurement lying on it is the method converging as advertised.',
	request(engine, method, priority, context) {
		const p = method.order;
		const coarse = p >= 7 ? 0.5 : p >= 5 ? 0.4 : 0.3;
		const ratio = p >= 5 ? 0.7 : 0.6;
		return engine.request(
			'convergence',
			{ id: method.id, problem: context.problem, coarse, ratio, count: 7 },
			{ key: `convergence:${method.id}:${context.problem}`, priority }
		);
	},
	figure(study, method, { compact = true } = {}) {
		const points = study.points.filter((p) => Number.isFinite(p.error) && p.error > 1e-15);
		if (points.length < 2) return empty('no usable points');
		const x = points.map((p) => p.h);
		const y = points.map((p) => p.error);

		// The guides are anchored on the finest step, where the method is
		// asymptotic, so the line for its own order should lie along the data all
		// the way back to the coarse end. The neighbouring orders fan out from the
		// same point, which turns "what slope is this" into a reading rather than
		// an estimate.
		const anchor = points.reduce((best, p) => (p.h < best.h ? p : best), points[0]);
		const span = [Math.min(...x), Math.max(...x)];
		const guideOrders = compact
			? [method.order]
			: [method.order - 2, method.order - 1, method.order, method.order + 1].filter((k) => k >= 1);

		const guides = guideOrders.map((k) => {
			const own = k === method.order;
			return {
				type: 'scatter',
				mode: 'lines',
				x: span,
				y: span.map((h) => anchor.error * (h / anchor.h) ** k),
				line: {
					color: own ? '#f0efe9' : '#969591',
					width: own ? 1.2 : 0.8,
					dash: own ? 'dash' : 'dot'
				},
				opacity: own ? 0.8 : 0.45,
				hoverinfo: 'skip'
			};
		});

		const annotations = compact
			? []
			: guideOrders.map((k) => ({
					x: Math.log10(span[1]),
					y: Math.log10(anchor.error * (span[1] / anchor.h) ** k),
					text: `h<sup>${k}</sup>`,
					showarrow: false,
					xanchor: 'left',
					yanchor: 'middle',
					xshift: 4,
					font: {
						family: 'JetBrains Mono, monospace',
						size: 9,
						color: k === method.order ? '#f0efe9' : '#969591'
					}
				}));

		// Below the round off floor the measurement is reporting the arithmetic,
		// not the method.
		const floor = 1e-14;
		const shapes =
			compact || Math.min(...y) > floor * 30
				? []
				: [
						{
							type: 'line',
							xref: 'paper',
							x0: 0,
							x1: 1,
							yref: 'y',
							y0: floor,
							y1: floor,
							line: { color: '#d9513c', width: 1, dash: 'dot' }
						}
					];

		// The guide ends have to fit, or the labels sit on the frame. They are
		// allowed to push the view up by a decade and a half and no further, so a
		// steep guide cannot shrink the data it is there to be compared against.
		const guideEnds = compact ? [] : guides.flatMap((guide) => guide.y);
		const ceiling = Math.max(...y) * 10 ** 1.5;
		const yRange = logRange([...y, ...guideEnds.filter((v) => v <= ceiling)], 0.16);

		// Room on the right for the exponent labels.
		const xRange = logRange(x, compact ? 0.06 : 0.1);
		if (!compact) xRange[1] += 0.12;

		return {
			data: [
				...guides,
				{
					type: 'scatter',
					mode: compact ? 'lines' : 'lines+markers',
					x,
					y,
					line: { color: seriesColor(0), width: 1.8 },
					marker: { size: 5 },
					hovertemplate: compact ? undefined : 'h %{x:.3e}  error %{y:.3e}<extra></extra>',
					hoverinfo: compact ? 'skip' : undefined
				}
			],
			layout: {
				...layout({
					compact,
					x: {
						type: 'log',
						range: xRange,
						title: compact ? undefined : { text: 'step size h' }
					},
					y: {
						type: 'log',
						range: yRange,
						title: compact ? undefined : { text: 'relative error' }
					}
				}),
				annotations,
				shapes
			},
			config: config({ compact })
		};
	}
};

const cost = {
	key: 'cost',
	label: 'cost',
	note: 'Achieved accuracy against the work it took, one adaptive run per tolerance. Accuracy improves to the right, so the lower curve is the cheaper method.',
	request(engine, method, priority, context) {
		return engine.request(
			'workPrecision',
			{ id: method.id, problem: context.problem, from: -3, to: -10 },
			{ key: `cost:${method.id}:${context.problem}`, priority }
		);
	},
	figure(result, method, { compact = true } = {}) {
		const points = result.points.filter(
			(p) => p.succeeded && p.error > 0 && Number.isFinite(p.error)
		);
		if (points.length < 2) return empty('did not complete');
		const x = points.map((p) => p.error);
		const y = points.map((p) => p.rhs_evals);
		const range = logRange(x);
		return {
			data: [
				{
					type: 'scatter',
					mode: compact ? 'lines' : 'lines+markers',
					x,
					y,
					line: { color: seriesColor(1), width: 1.6 },
					marker: { size: 5 },
					hoverinfo: compact ? 'skip' : 'x+y'
				}
			],
			layout: layout({
				compact,
				x: {
					type: 'log',
					range: [range[1], range[0]],
					title: compact ? undefined : { text: 'achieved error' }
				},
				y: {
					type: 'log',
					range: logRange(y, 0.12),
					title: compact ? undefined : { text: 'rhs evaluations' }
				}
			}),
			config: config({ compact })
		};
	}
};

const orderStar = {
	key: 'orderStar',
	label: 'order star',
	note: 'log |R(z) exp(-z)|. Where the stability region says whether a mode decays, the star says whether the method tracks the exact solution: the sectors meeting at the origin count one more than the order, and a method is A-acceptable exactly when no sector reaches into the left half plane.',
	request(engine, method, priority) {
		if (method.class !== 'runge_kutta') return Promise.resolve(null);
		const view = orderStarWindow(method);
		const resolution = 160;
		return engine.request(
			'orderStar',
			{ id: method.id, re: view.re, im: view.im, width: resolution, height: resolution },
			{ key: `orderStar:${method.id}`, priority }
		);
	},
	figure(result, method, { compact = true } = {}) {
		if (!result) return empty('Runge-Kutta only');
		const { data, width, height, re, im } = result;
		const z = reshape(data, width, height);
		const x = axisValues(re[0], re[1], width);
		const y = axisValues(im[0], im[1], height);
		// The star is a two colour object: inside or outside the unit modulus.
		// A narrow symmetric range keeps the boundary crisp instead of washing
		// it out against values that run to thirty decades.
		const reach = 1.2;
		return {
			data: [
				{
					type: 'heatmap',
					z,
					x,
					y,
					zmin: -reach,
					zmax: reach,
					colorscale: bandedScale(2, ['#17313a', '#17313a', '#5a3417', '#5a3417']),
					showscale: false,
					zsmooth: 'best',
					hoverinfo: 'skip'
				},
				{
					type: 'contour',
					z,
					x,
					y,
					contours: { start: 0, end: 0, size: 1, coloring: 'none' },
					line: { color: '#f5a623', width: compact ? 1 : 1.4 },
					showscale: false,
					hoverinfo: 'skip'
				}
			],
			layout: {
				...layout({
					compact,
					x: { title: compact ? undefined : { text: 'Re(z)' }, range: re },
					y: { title: compact ? undefined : { text: 'Im(z)' }, range: im }
				}),
				shapes: originLines({ compact })
			},
			config: config({ compact })
		};
	}
};

/**
 * The leading error coefficients, one bar per rooted tree.
 *
 * Two methods of the same order are not equally accurate. What separates them
 * is how large the residuals of the first unsatisfied order conditions are and
 * which elementary differentials they sit on, which is exactly what this shows.
 */
export const errorCoefficients = {
	key: 'errorCoefficients',
	label: 'error coefficients',
	request(engine, method, priority) {
		if (method.class !== 'runge_kutta') return Promise.resolve(null);
		return engine.request(
			'errorCoefficients',
			{ id: method.id },
			{ key: `errors:${method.id}`, priority }
		);
	},
	figure(result, method, { compact = false } = {}) {
		if (!result) return empty('Runge-Kutta only');
		const all = result.conditions.map((condition) => ({
			...condition,
			magnitude: Math.abs(condition.residual)
		}));
		const largest = all.reduce((acc, c) => Math.max(acc, c.magnitude), 0);
		// Anything eight decades below the largest is zero to rounding and only
		// stretches the axis.
		const floor = Math.max(largest * 1e-8, 1e-18);
		const conditions = all
			.filter((c) => c.magnitude >= floor)
			.sort((a, b) => a.magnitude - b.magnitude)
			.slice(-12);
		if (conditions.length === 0) return empty('no conditions at this order');
		return {
			data: [
				{
					type: 'bar',
					orientation: 'h',
					x: conditions.map((c) => Math.max(c.magnitude, floor)),
					y: conditions.map((c) => c.tree),
					marker: { color: seriesColor(1) },
					hovertemplate: '%{y}  %{x:.3e}<extra></extra>'
				}
			],
			layout: layout({
				compact,
				x: {
					type: 'log',
					range: logRange(
						conditions.map((c) => Math.max(c.magnitude, floor)),
						0.08
					),
					title: compact ? undefined : { text: `residual at order ${result.atOrder}` }
				},
				y: { type: 'category', automargin: true, ticks: '' }
			}),
			config: config({ compact })
		};
	}
};

export const METHOD_MODES = [stability, orderStar, structure, damping, order, cost];
export const MODES_NEEDING_A_PROBLEM = new Set(['order', 'cost']);

// ---------------------------------------------------------------------------
// Problem modes
// ---------------------------------------------------------------------------

function profileRequest(engine, problem, priority) {
	return engine.request(
		'problemProfile',
		{ id: problem.id, samples: 600 },
		{ key: `profile:${problem.id}`, priority }
	);
}

const solution = {
	key: 'solution',
	label: 'solution',
	note: 'The components over the interval.',
	request: profileRequest,
	figure(profile, problem, { compact = true } = {}) {
		const components = Math.min(profile.dim, 6);
		const data = [];
		for (let c = 0; c < components; c += 1) {
			data.push({
				type: 'scatter',
				mode: 'lines',
				x: profile.t,
				y: profile.y.map((row) => row[c]),
				name: `y${c + 1}`,
				line: { color: seriesColor(c), width: 1.3 },
				hoverinfo: compact ? 'skip' : 'x+y+name'
			});
		}
		return {
			data,
			layout: {
				...layout({
					compact,
					x: { range: profile.tSpan, title: compact ? undefined : { text: 't' } },
					y: {
						range: linearRange(data.flatMap((s) => s.y)),
						title: compact ? undefined : { text: 'y' }
					}
				}),
				shapes: originLines({ x: false, compact, color: COLORS.accent })
			},
			config: config({ compact })
		};
	}
};

const phase = {
	key: 'phase',
	label: 'phase',
	note: 'The first two components against each other.',
	request: profileRequest,
	figure(profile, problem, { compact = true } = {}) {
		if (profile.dim < 2) return empty('scalar problem');
		const x = profile.y.map((row) => row[0]);
		const y = profile.y.map((row) => row[1]);
		return {
			data: [
				{
					type: 'scatter',
					mode: 'lines',
					x,
					y,
					line: { color: seriesColor(3), width: 1.3 },
					hoverinfo: compact ? 'skip' : 'x+y'
				}
			],
			layout: {
				...layout({
					compact,
					x: { range: linearRange(x), title: compact ? undefined : { text: 'y1' } },
					y: { range: linearRange(y), title: compact ? undefined : { text: 'y2' } }
				}),
				shapes: originLines({ compact, color: COLORS.accent })
			},
			config: config({ compact })
		};
	}
};

const stiffness = {
	key: 'stiffness',
	label: 'stiffness',
	note: 'The spectral abscissa of the Jacobian along the solution. Where it spikes is where an explicit method has to take tiny steps.',
	request: profileRequest,
	figure(profile, problem, { compact = true } = {}) {
		const y = profile.decayRate.map((v) => Math.max(v, 1e-6));
		return {
			data: [
				{
					type: 'scatter',
					mode: 'lines',
					x: profile.t,
					y,
					line: { color: seriesColor(1), width: 1.3 },
					hoverinfo: compact ? 'skip' : 'x+y'
				}
			],
			layout: layout({
				compact,
				x: { range: profile.tSpan, title: compact ? undefined : { text: 't' } },
				y: {
					type: 'log',
					range: logRange(y, 0.1),
					title: compact ? undefined : { text: 'max(-Re lambda)' }
				}
			}),
			config: config({ compact })
		};
	}
};

const steps = {
	key: 'steps',
	label: 'step size',
	note: 'The step size a high accuracy reference run needed, which is the same information seen from the solver side.',
	request: profileRequest,
	figure(profile, problem, { compact = true } = {}) {
		const times = [];
		const sizes = [];
		let t = profile.tSpan[0];
		for (const h of profile.steps) {
			times.push(t);
			sizes.push(Math.abs(h));
			t += h;
		}
		if (sizes.length < 2) return empty('no steps recorded');
		return {
			data: [
				{
					type: 'scatter',
					mode: 'lines',
					x: times,
					y: sizes,
					line: { color: seriesColor(0), width: 1.2 },
					hoverinfo: compact ? 'skip' : 'x+y'
				}
			],
			layout: layout({
				compact,
				x: { range: profile.tSpan, title: compact ? undefined : { text: 't' } },
				y: {
					type: 'log',
					range: logRange(sizes, 0.1),
					title: compact ? undefined : { text: 'h' }
				}
			}),
			config: config({ compact })
		};
	}
};

export const PROBLEM_MODES = [solution, phase, stiffness, steps];
