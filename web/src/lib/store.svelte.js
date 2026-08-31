/*
 * Shared state.
 *
 * One engine for the whole session, one copy of the catalogues, and the
 * handful of choices that persist as the reader moves between the library, the
 * problems and the comparison.
 */

import { Engine, PRIORITY } from './engine.js';

let instance = null;

export function engine() {
	if (!instance) instance = new Engine();
	return instance;
}

/// The catalogue is fetched in slices of this many methods, which is a
/// compromise: small enough that the first tiles appear quickly, large enough
/// that the per job overhead stays negligible against the analysis itself.
const CATALOG_SLICE = 16;

export const catalog = $state({
	methods: [],
	problems: [],
	options: { controllers: [], nonlinearSolvers: [] },
	/// How many methods the library holds, known before any of them arrive.
	total: 0,
	/// How many have arrived so far.
	loaded: 0,
	ready: false,
	error: ''
});

export const ui = $state({
	// Library
	mode: 'stability',
	problem: 'nonlinear_decay',
	query: '',
	tags: [],
	sort: 'family',
	selecting: false,
	selection: [],
	// Problems
	problemMode: 'solution',
	// Comparison
	compareView: 'cost',
	compareProblem: 'lotka_volterra',
	compareController: 'pi4020',
	presetKey: 'explicit-pairs'
});

export const progress = $state({ outstanding: 0 });

export async function loadCatalogues() {
	if (catalog.ready) return;
	const jobs = engine();
	jobs.onProgress((outstanding) => {
		progress.outstanding = outstanding;
	});
	try {
		const [problems, options, total] = await Promise.all([
			jobs.request('problemCatalog', {}, { key: 'catalog:problems', priority: PRIORITY.immediate }),
			jobs.request('optionsCatalog', {}, { key: 'catalog:options', priority: PRIORITY.immediate }),
			jobs.request('methodCount', {}, { key: 'catalog:count', priority: PRIORITY.immediate })
		]);
		catalog.problems = problems;
		catalog.options = options;
		catalog.total = total;

		// Deriving the properties of the whole library is seconds of work, and
		// asking for it in one job means an empty page for all of them. In
		// slices the pool runs several at once and the grid fills as they land,
		// which is both faster and something to look at while it happens.
		const slices = [];
		for (let offset = 0; offset < total; offset += CATALOG_SLICE) {
			const count = Math.min(CATALOG_SLICE, total - offset);
			slices.push(
				jobs
					.request(
						'methodCatalogSlice',
						{ offset, count },
						{ key: `catalog:methods:${offset}`, priority: PRIORITY.immediate }
					)
					.then((entries) => {
						// Slices can land out of order; the grid sorts itself,
						// but the underlying list stays in library order so that
						// a rerun of the same filters gives the same answer.
						catalog.methods = [...catalog.methods, ...entries].sort(
							(left, right) => left.id.localeCompare(right.id)
						);
						catalog.loaded = catalog.methods.length;
						if (!catalog.ready) catalog.ready = true;
					})
			);
		}
		await Promise.all(slices);
		catalog.ready = true;
	} catch (error) {
		catalog.error = String(error);
	}
}

export function toggleSelection(id) {
	const index = ui.selection.indexOf(id);
	if (index >= 0) ui.selection.splice(index, 1);
	else ui.selection.push(id);
}

export function toggleTag(name) {
	const index = ui.tags.indexOf(name);
	if (index >= 0) ui.tags.splice(index, 1);
	else ui.tags.push(name);
}
