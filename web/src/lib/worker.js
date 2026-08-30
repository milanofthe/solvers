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

/** A window grown about its own centre, for sampling past what is framed. */
function grown({ re, im }, factor) {
	const stretch = ([low, high]) => {
		const middle = (low + high) / 2;
		const half = ((high - low) / 2) * factor;
		return [middle - half, middle + half];
	};
	return { re: stretch(re), im: stretch(im) };
}

const HANDLERS = {
	methodCatalog: () => JSON.parse(wasm.method_catalog()),
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
	stabilityGridAuto: ({ id, width, height, locus = 0 }) => {
		const view = JSON.parse(wasm.stability_window(id));
		// Sample wider than the frame. Both axes of the picture carry the same
		// scale, so a panel of any other shape than the measured box widens one
		// of the two; sampling only the box would leave that strip empty and
		// read as a hole in the figure rather than as more of the plane.
		const sample = grown(view, 1.7);
		const data = wasm.stability_grid(
			id,
			sample.re[0],
			sample.re[1],
			sample.im[0],
			sample.im[1],
			width,
			height
		);
		// The boundary of a multistep region is available in closed form, so it
		// is fetched with the grid rather than traced out of it. A region with
		// no interior has nothing for a contour to find, and that is not a rare
		// case: it is what the whole weakly stable half of the library looks
		// like.
		return {
			data,
			width,
			height,
			re: sample.re,
			im: sample.im,
			view,
			locus: locus ? wasm.region_boundary(id, locus) : null
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
