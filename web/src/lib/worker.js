/*
 * The compute worker.
 *
 * Hosts one instance of the solver core and answers one job at a time. Every
 * number the interface shows is produced here, which keeps the main thread free
 * to scroll a grid of fifty four cards while stability regions and work
 * precision runs are still being computed behind it.
 *
 * JSON is parsed here rather than on the main thread, so what crosses back is
 * already a structure. A grid crosses back as a transferable buffer and is not
 * copied at all.
 */

import init, * as wasm from './wasm/solvers_wasm.js';
import wasmUrl from './wasm/solvers_wasm_bg.wasm?url';

const ready = init({ module_or_path: wasmUrl });

/**
 * A window grown so that a panel of the given shape is covered by it.
 *
 * The frame is the measured box, but both axes carry one scale, so a panel
 * whose shape differs from the box widens one of the two until they agree.
 * Growing to cover every shape a panel could have wastes most of the grid on a
 * region far from square: the stability region of Adams-Bashforth 5 is twice as
 * tall as it is wide, and in a landscape panel four fifths of the samples would
 * land outside the frame. So the caller says what shape it actually is, with a
 * tenth of slack either way so that nudging a window is not a resample.
 */
function covering({ re, im }, widest = 2.2, tallest = 0.7) {
	const width = re[1] - re[0];
	const height = im[1] - im[0];
	const middleX = (re[0] + re[1]) / 2;
	const middleY = (im[0] + im[1]) / 2;
	// Two percent over what is needed, so a window that matches the panel
	// exactly still has data at its very edge.
	const reach = (1.02 * Math.max(width, height * widest)) / 2;
	const rise = (1.02 * Math.max(height, width / tallest)) / 2;
	return {
		re: [middleX - reach, middleX + reach],
		im: [middleY - rise, middleY + rise]
	};
}

const HANDLERS = {
	methodCatalog: () => JSON.parse(wasm.method_catalog()),
	methodCount: () => wasm.method_count(),
	methodCatalogSlice: ({ offset, count }) =>
		JSON.parse(wasm.method_catalog_slice(offset, count)),
	problemCatalog: () => JSON.parse(wasm.problem_catalog()),
	optionsCatalog: () => JSON.parse(wasm.options_catalog()),
	methodDetail: ({ id }) => JSON.parse(wasm.method_detail(id)),
	methodSummary: ({ id }) => JSON.parse(wasm.method_summary(id)),
	stabilityFunction: ({ id }) => JSON.parse(wasm.stability_function(id)),
	errorCoefficients: ({ id }) => JSON.parse(wasm.error_coefficients(id)),
	orderStar: ({ id, re, im, width, height }) => {
		const data = wasm.order_star_grid(id, re[0], re[1], im[0], im[1], width, height);
		return { data, width, height, re, im };
	},
	problemProfile: ({ id, samples }) => JSON.parse(wasm.problem_profile(id, samples)),
	convergence: ({ id, problem, coarse, ratio, count }) =>
		JSON.parse(wasm.convergence_study(id, problem, coarse, ratio, count)),
	workPrecision: ({ id, problem, from, to }) =>
		JSON.parse(wasm.work_precision(id, problem, from, to)),
	trajectory: ({ id, problem, rtol, atol, samples }) =>
		JSON.parse(wasm.trajectory(id, problem, rtol, atol, samples)),
	stepHistory: ({ id, problem, rtol, atol, controller }) =>
		JSON.parse(wasm.step_history(id, problem, rtol, atol, controller)),
	stabilityGrid: ({ id, re, im, width, height }) => {
		const data = wasm.stability_grid(id, re[0], re[1], im[0], im[1], width, height);
		return { data, width, height, re, im };
	},
	// The same figure over a window the caller names, for a reader who has
	// panned or zoomed somewhere the measured window does not reach.
	stabilityWindow: ({ id, re, im, view, width, height }) => ({
		data: wasm.stability_grid(id, re[0], re[1], im[0], im[1], width, height),
		width,
		height,
		re,
		im,
		view,
		boundary: wasm.region_boundary(id, re[0], re[1], im[0], im[1])
	}),
	stabilityGridAuto: ({ id, width, height, aspect = 1.6 }) => {
		const view = JSON.parse(wasm.stability_window(id));
		// Sample wider than the frame. Both axes of the picture carry the same
		// scale, so a panel of any other shape than the measured box widens one
		// of the two; sampling only the box would leave that strip empty and
		// read as a hole in the figure rather than as more of the plane.
		const sample = covering(view, aspect * 1.1, aspect / 1.1);
		const data = wasm.stability_grid(
			id,
			sample.re[0],
			sample.re[1],
			sample.im[0],
			sample.im[1],
			width,
			height
		);
		// The boundary comes with the grid but is not read off it. The colours
		// are a smooth field that a coarse grid renders perfectly well; the
		// boundary is a curve, and a curve wants its samples on itself. Two
		// resolutions for two different objects, one round trip.
		return {
			data,
			width,
			height,
			re: sample.re,
			im: sample.im,
			view,
			boundary: wasm.region_boundary(
				id,
				sample.re[0],
				sample.re[1],
				sample.im[0],
				sample.im[1]
			)
		};
	}
};

self.onmessage = async (event) => {
	const { id, kind, args } = event.data;
	try {
		await ready;
		const handler = HANDLERS[kind];
		if (!handler) throw new Error(`unknown job: ${kind}`);
		const result = handler(args ?? {});
		const transfer = Object.values(result ?? {})
			.filter((v) => v instanceof Float64Array)
			.map((v) => v.buffer);
		self.postMessage({ id, ok: true, result }, transfer);
	} catch (error) {
		self.postMessage({ id, ok: false, error: String(error?.message ?? error) });
	}
};
