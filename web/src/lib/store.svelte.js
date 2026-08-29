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

export const catalog = $state({
	methods: [],
	problems: [],
	options: { controllers: [], nonlinearSolvers: [] },
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
		const [methods, problems, options] = await Promise.all([
			jobs.request('methodCatalog', {}, { key: 'catalog:methods', priority: PRIORITY.immediate }),
			jobs.request('problemCatalog', {}, { key: 'catalog:problems', priority: PRIORITY.immediate }),
			jobs.request('optionsCatalog', {}, { key: 'catalog:options', priority: PRIORITY.immediate })
		]);
		catalog.methods = methods;
		catalog.problems = problems;
		catalog.options = options;
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
