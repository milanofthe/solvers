<script>
	/*
	 * Comparison.
	 *
	 * Two ways in: pick cards out of the library, or take a preset. The presets
	 * are the comparisons worth making, the ones that answer a question somebody
	 * actually has: which explicit pair to use on a smooth problem, whether the
	 * extra stages of a higher order method pay for themselves, what BDF costs
	 * against ESDIRK on something stiff.
	 */
	import { getContext } from 'svelte';
	import Plot from '$lib/components/Plot.svelte';
	import { catalog, engine, ui } from '$lib/store.svelte.js';
	import { isCancelled, PRIORITY } from '$lib/engine.js';
	import { config, layout, logRange, linearRange, seriesColor } from '$lib/plot.js';
	import { num, stabilityLabel } from '$lib/format.js';

	const rail = getContext('rail');

	const PRESETS = [
		{
			key: 'explicit-pairs',
			label: 'explicit pairs',
			problem: 'lotka_volterra',
			methods: ['rkbs32', 'rkf45', 'rkck54', 'rkdp54', 'rkv65', 'rkdp87'],
			note: 'The embedded explicit pairs against each other on a smooth problem, which is the choice for anything non stiff.'
		},
		{
			key: 'order-pays',
			label: 'does order pay',
			problem: 'arenstorf',
			methods: ['rkbs32', 'rkdp54', 'rkv65', 'rkdp87'],
			note: 'The same family at four orders. On a demanding orbit the high order methods reach an accuracy the low order ones cannot buy at any price.'
		},
		{
			key: 'stiff-families',
			label: 'stiff families',
			problem: 'van_der_pol_stiff',
			methods: ['esdirk32', 'esdirk43', 'esdirk54', 'bdf3', 'bdf4', 'radau_iia_5'],
			note: 'One step against multistep against fully implicit, on the standard stiff benchmark.'
		},
		{
			key: 'bdf-order',
			label: 'BDF by order',
			problem: 'robertson',
			methods: ['bdf2', 'bdf3', 'bdf4', 'bdf5', 'bdf6'],
			note: 'Higher order BDF is cheaper per digit but its stability wedge closes. On Robertson that trade is visible directly.'
		},
		{
			key: 'adams',
			label: 'Adams both ways',
			problem: 'brusselator',
			methods: ['adams_bashforth_3', 'adams_bashforth_5', 'adams_moulton_3', 'adams_moulton_5'],
			note: 'One evaluation per step against one solve per step, at matched orders.'
		},
		{
			key: 'collocation',
			label: 'collocation',
			problem: 'kaps',
			methods: [
				'radau_iia_3',
				'radau_iia_5',
				'gauss_legendre_4',
				'gauss_legendre_6',
				'lobatto_iiic_4'
			],
			note: 'Gauss, Radau and Lobatto on a problem built to expose order reduction.'
		},
		{
			key: 'order-reduction',
			label: 'order reduction',
			problem: 'prothero_robinson',
			methods: ['esdirk43', 'sdirk3', 'radau_iia_5', 'trapezoidal', 'bdf4'],
			note: 'Prothero-Robinson punishes a low stage order. What each method converges at here is often not its classical order.'
		}
	];

	const VIEWS = [
		{ key: 'cost', label: 'work precision' },
		{ key: 'order', label: 'convergence' },
		{ key: 'steps', label: 'step size' },
		{ key: 'table', label: 'summary' }
	];

	let figure = $state.raw(null);
	let status = $state('');

	const chosen = $derived(
		ui.selection.map((id) => catalog.methods.find((m) => m.id === id)).filter(Boolean)
	);
	const preset = $derived(PRESETS.find((p) => p.key === ui.presetKey));

	function applyPreset(entry) {
		engine().drop(() => true);
		ui.presetKey = entry.key;
		ui.selection = entry.methods.filter((id) => catalog.methods.some((m) => m.id === id));
		ui.compareProblem = entry.problem;
	}

	if (ui.selection.length === 0 && preset) applyPreset(preset);

	function axisTitle(text) {
		return { text, font: { family: 'JetBrains Mono, monospace', size: 10 } };
	}

	async function build() {
		if (chosen.length === 0) {
			figure = null;
			status = 'Nothing selected. Take a preset, or choose cards in the library.';
			return;
		}
		if (ui.compareView === 'table') {
			figure = null;
			status = '';
			return;
		}
		status = 'computing';
		const jobs = engine();
		const view = ui.compareView;
		const problem = ui.compareProblem;

		const results = await Promise.all(
			chosen.map((method) => {
				if (view === 'cost') {
					return jobs
						.request(
							'workPrecision',
							{ id: method.id, problem, from: -3, to: -11 },
							{ key: `cost:${method.id}:${problem}`, priority: PRIORITY.immediate }
						)
						.then((data) => ({ method, data }))
						.catch((error) => (isCancelled(error) ? null : { method, error }));
				}
				if (view === 'order') {
					const p = method.order;
					return jobs
						.request(
							'convergence',
							{
								id: method.id,
								problem,
								coarse: p >= 7 ? 0.5 : p >= 5 ? 0.4 : 0.3,
								ratio: p >= 5 ? 0.7 : 0.6,
								count: 7
							},
							{ key: `convergence:${method.id}:${problem}`, priority: PRIORITY.immediate }
						)
						.then((data) => ({ method, data }))
						.catch((error) => (isCancelled(error) ? null : { method, error }));
				}
				return jobs
					.request(
						'stepHistory',
						{
							id: method.id,
							problem,
							rtol: 1e-6,
							atol: 1e-8,
							controller: ui.compareController
						},
						{
							key: `steps:${method.id}:${problem}:${ui.compareController}`,
							priority: PRIORITY.immediate
						}
					)
					.then((data) => ({ method, data }))
					.catch((error) => (isCancelled(error) ? null : { method, error }));
			})
		);

		// The reader may have moved on while this was running.
		if (view !== ui.compareView || problem !== ui.compareProblem) return;

		const traces = [];
		results.forEach((entry, index) => {
			if (!entry?.data) return;
			const color = seriesColor(index);
			if (view === 'cost') {
				const points = entry.data.points.filter(
					(p) => p.succeeded && p.error > 0 && Number.isFinite(p.error)
				);
				if (points.length < 2) return;
				traces.push({
					type: 'scatter',
					mode: 'lines+markers',
					name: entry.method.id,
					x: points.map((p) => p.error),
					y: points.map((p) => p.rhs_evals),
					line: { color, width: 1.6 },
					marker: { size: 5 }
				});
			} else if (view === 'order') {
				const points = entry.data.points.filter((p) => Number.isFinite(p.error) && p.error > 1e-15);
				if (points.length < 2) return;
				traces.push({
					type: 'scatter',
					mode: 'lines+markers',
					name: entry.method.id,
					x: points.map((p) => p.h),
					y: points.map((p) => p.error),
					line: { color, width: 1.6 },
					marker: { size: 5 }
				});
			} else {
				if (!entry.data.h?.length) return;
				traces.push({
					type: 'scattergl',
					mode: 'lines',
					name: entry.method.id,
					x: entry.data.t,
					y: entry.data.h.map(Math.abs),
					line: { color, width: 1.2 }
				});
			}
		});

		if (traces.length === 0) {
			figure = null;
			status = 'No method produced a usable run on this problem.';
			return;
		}

		const allX = traces.flatMap((t) => t.x);
		const allY = traces.flatMap((t) => t.y);
		let base;
		if (view === 'cost') {
			const range = logRange(allX);
			base = layout({
				title: `work precision on ${problem}`,
				x: { type: 'log', range: [range[1], range[0]], title: axisTitle('achieved relative error') },
				y: { type: 'log', range: logRange(allY, 0.12), title: axisTitle('rhs evaluations') }
			});
			status =
				'Accuracy improves to the right, so the lower curve is the cheaper method at the same accuracy.';
		} else if (view === 'order') {
			base = layout({
				title: `convergence on ${problem}`,
				x: { type: 'log', range: logRange(allX), title: axisTitle('step size h') },
				y: { type: 'log', range: logRange(allY, 0.12), title: axisTitle('relative error') }
			});
			const reduced = results
				.filter(
					(entry) =>
						entry?.data &&
						Number.isFinite(entry.data.local_order) &&
						entry.method.order - entry.data.local_order > 0.5
				)
				.map((entry) => `${entry.method.id} at ${num(entry.data.local_order, 3)} instead of ${entry.method.order}`);
			status = reduced.length
				? `Order reduction here: ${reduced.join(', ')}.`
				: 'Every selected method converges at the order its coefficients promise.';
		} else {
			base = layout({
				title: `step size on ${problem}, controller ${ui.compareController}`,
				x: { range: linearRange(allX, 0.01), title: axisTitle('t') },
				y: { type: 'log', range: logRange(allY, 0.1), title: axisTitle('step size h') }
			});
			status = 'Where the step size collapses is where the solution does something the method has to resolve.';
		}

		figure = { data: traces, layout: { ...base, showlegend: true }, config: config({}) };
	}

	$effect(() => {
		// Rebuild whenever any of these change.
		ui.compareView;
		ui.compareProblem;
		ui.compareController;
		ui.selection.length;
		build();
	});

	$effect(() => {
		rail.panel = controls;
		return () => {
			rail.panel = null;
		};
	});
</script>

{#snippet controls()}
	<div>
		<p class="label mb-1">show</p>
		{#each VIEWS as view}
			<button
				type="button"
				class="facet {ui.compareView === view.key ? 'facet-on' : ''}"
				onclick={() => (ui.compareView = view.key)}
			>
				<span>{view.label}</span>
			</button>
		{/each}
	</div>

	<div>
		<p class="label mb-1">presets</p>
		{#each PRESETS as entry}
			<button
				type="button"
				class="facet {ui.presetKey === entry.key ? 'facet-on' : ''}"
				onclick={() => applyPreset(entry)}
			>
				<span class="truncate">{entry.label}</span>
			</button>
		{/each}
	</div>

	<div>
		<p class="label mb-1">problem</p>
		{#each catalog.problems as problem}
			<button
				type="button"
				class="facet {ui.compareProblem === problem.id ? 'facet-on' : ''}"
				onclick={() => {
					ui.presetKey = '';
					ui.compareProblem = problem.id;
				}}
			>
				<span class="truncate">{problem.id}</span>
			</button>
		{/each}
	</div>

	{#if ui.compareView === 'steps'}
		<div>
			<p class="label mb-1">controller</p>
			{#each catalog.options.controllers as name}
				<button
					type="button"
					class="facet {ui.compareController === name ? 'facet-on' : ''}"
					onclick={() => (ui.compareController = name)}
				>
					<span>{name}</span>
				</button>
			{/each}
		</div>
	{/if}

	<div>
		<p class="label mb-1">selected</p>
		{#each chosen as method}
			<button
				type="button"
				class="facet facet-on"
				title="remove"
				onclick={() => {
					ui.presetKey = '';
					ui.selection = ui.selection.filter((id) => id !== method.id);
				}}
			>
				<span class="truncate">{method.id}</span>
				<span class="text-accent">&times;</span>
			</button>
		{/each}
		<a href="/" class="action mt-2 inline-block">[ choose in the library ]</a>
	</div>
{/snippet}

<header class="mb-6 border-b border-cream/10 pb-4">
	<h1 class="font-display text-lg text-cream">Comparison</h1>
	<p class="mt-1 max-w-[68ch] text-xs text-cream/60">
		{preset ? preset.note : 'Whatever is selected in the library, on the problem chosen in the rail.'}
	</p>
</header>

{#if ui.compareView === 'table'}
	<div class="overflow-x-auto border border-cream/10">
		<table class="w-full border-collapse font-mono text-xs">
			<thead>
				<tr class="border-b border-cream/15 text-accent">
					<th class="px-3 py-2 text-left font-normal uppercase tracking-[0.15em]">method</th>
					<th class="px-3 py-2 text-right font-normal uppercase tracking-[0.15em]">order</th>
					<th class="px-3 py-2 text-right font-normal uppercase tracking-[0.15em]">size</th>
					<th class="px-3 py-2 text-right font-normal uppercase tracking-[0.15em]">cost</th>
					<th class="px-3 py-2 text-right font-normal uppercase tracking-[0.15em]">stage order</th>
					<th class="px-3 py-2 text-left font-normal uppercase tracking-[0.15em]">stability</th>
					<th class="px-3 py-2 text-left font-normal uppercase tracking-[0.15em]">kind</th>
				</tr>
			</thead>
			<tbody>
				{#each chosen as method}
					<tr class="border-b border-cream/10">
						<td class="px-3 py-1"
							><a href="/methods/{method.id}" class="text-cream hover:text-amber">{method.name}</a
							></td
						>
						<td class="px-3 py-1 text-right text-cream">{method.order}</td>
						<td class="px-3 py-1 text-right text-cream">{method.size}</td>
						<td class="px-3 py-1 text-right text-cream">{method.stageCost}</td>
						<td class="px-3 py-1 text-right text-cream">{method.stageOrder ?? '–'}</td>
						<td class="px-3 py-1 text-cream/70">{stabilityLabel(method)}</td>
						<td class="px-3 py-1 text-cream/70">{method.implicit ? 'implicit' : 'explicit'}</td>
					</tr>
				{/each}
			</tbody>
		</table>
	</div>
{:else}
	<div class="border border-cream/10 bg-charcoal-light">
		<Plot {figure} class="h-[32rem] w-full" />
	</div>
{/if}

<p class="mt-3 max-w-[80ch] text-xs text-cream/50">{status}</p>
