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
	boundaryLocus: ({ id, samples }) => wasm.boundary_locus(id, samples),
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
	}
};

self.onmessage = async (event) => {
	const { id, kind, args } = event.data;
	try {
		await ready;
		const handler = HANDLERS[kind];
		if (!handler) throw new Error(`unknown job: ${kind}`);
		const result = handler(args ?? {});
		const transfer = result?.data instanceof Float64Array ? [result.data.buffer] : [];
		self.postMessage({ id, ok: true, result }, transfer);
	} catch (error) {
		self.postMessage({ id, ok: false, error: String(error?.message ?? error) });
	}
};
