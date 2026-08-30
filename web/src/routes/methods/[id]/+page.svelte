<script>
	/*
	 * One method in full.
	 *
	 * Two columns. The figures take the width, because they are what there is to
	 * look at; the numbers run alongside in a rail that stays put while the
	 * figures scroll, because they are what you keep looking back at. Below the
	 * width where a third column fits, the rail becomes a block under the title
	 * instead of disappearing.
	 *
	 * The figures are grouped by the question they answer rather than by size:
	 * where the method is stable, how accurate it is, and what it is made of.
	 * Within a group every figure is the same height, which is the only thing
	 * that makes a grid of plots read as a grid.
	 */
	import { getContext } from 'svelte';
	import { page } from '$app/state';
	import Panel from '$lib/components/Panel.svelte';
	import { catalog, engine, ui } from '$lib/store.svelte.js';
	import { PRIORITY } from '$lib/engine.js';
	import { METHOD_MODES, errorCoefficients } from '$lib/figures.js';
	import { limit, num } from '$lib/format.js';

	const rail = getContext('rail');

	const id = $derived(page.params.id);
	const method = $derived(catalog.methods.find((m) => m.id === id));

	let detail = $state.raw(null);
	let stabilityFunction = $state.raw(null);

	const MODES = Object.fromEntries(METHOD_MODES.map((mode) => [mode.key, mode]));

	/** A loader for one figure, bound to the method currently on screen. */
	function panel(mode, current, problem) {
		return (box, view) =>
			(view && mode.resample
				? mode.resample(engine(), current, view, box, PRIORITY.immediate)
				: mode.request(engine(), current, PRIORITY.immediate, { problem, box })
			).then((result) => mode.figure(result, current, { compact: false }));
	}

	$effect(() => {
		const current = method;
		if (!current) return;
		let stale = false;
		detail = null;
		stabilityFunction = null;

		engine()
			.request(
				'methodDetail',
				{ id: current.id },
				{ key: `detail:${current.id}`, priority: PRIORITY.immediate }
			)
			.then((result) => {
				if (!stale) detail = result;
			})
			.catch(() => {});

		if (current.class !== 'linear_multistep') {
			engine()
				.request(
					'stabilityFunction',
					{ id: current.id },
					{ key: `function:${current.id}`, priority: PRIORITY.prefetch }
				)
				.then((result) => {
					if (!stale) stabilityFunction = result;
				})
				.catch(() => {});
		}

		return () => {
			stale = true;
		};
	});

	$effect(() => {
		rail.panel = controls;
		return () => {
			rail.panel = null;
		};
	});

	/**
	 * The numbers, in three groups: what the method is, how it behaves in the
	 * complex plane, and how the two were established.
	 */
	const facts = $derived.by(() => {
		if (!detail) return [];
		const report = detail.report;
		const size = report.class === 'linear_multistep' ? 'steps' : 'stages';
		return [
			[
				'method',
				[
					['order', report.computed_order],
					['embedded', report.computed_embedded_order ?? '–'],
					[size, report.size],
					['stage order', report.stage_order ?? '–'],
					['cost per step', report.stage_cost]
				]
			],
			[
				'stability',
				[
					['A-stable', report.a_stable ? 'yes' : 'no'],
					['L-stable', report.l_stable ? 'yes' : 'no'],
					[
						'stiffly accurate',
						report.stiffly_accurate == null ? '–' : report.stiffly_accurate ? 'yes' : 'no'
					],
					['R(inf)', damping(report.damping_at_infinity)],
					['A(alpha)', report.alpha_angle == null ? '–' : `${num(report.alpha_angle, 4)}°`],
					['real limit', boundless(report.real_stability_limit)],
					['imaginary limit', boundless(report.imaginary_stability_limit)]
				]
			],
			['checked', [['order conditions', report.exact_arithmetic ? 'exactly' : 'numerically']]]
		];
	});

	/** A limit with no finite value is worth a word rather than a glyph. */
	function boundless(value) {
		if (value == null) return '–';
		return value === 'unbounded' ? 'unbounded' : num(value);
	}

	function damping(value) {
		if (value == null) return '–';
		return value === 'unbounded' ? 'unbounded' : num(value);
	}

	const structureLine = $derived.by(() => {
		if (!detail) return '';
		const c = detail.coefficients;
		if (c.kind === 'runge_kutta') {
			return [
				c.structure?.replace('_', ' '),
				c.singlyDiagonal ? 'singly diagonal' : null,
				c.explicitFirstStage ? 'explicit first stage' : null,
				c.fsal ? 'FSAL' : null
			]
				.filter(Boolean)
				.join(' · ');
		}
		if (c.kind === 'rosenbrock') {
			return [
				'linearly implicit',
				c.singlyDiagonal ? 'singly diagonal' : null,
				c.stifflyAccurate ? 'stiffly accurate' : null
			]
				.filter(Boolean)
				.join(' · ');
		}
		return [`${c.steps} steps`, c.startup ? `starts with ${c.startup}` : 'self starting'].join(
			' · '
		);
	});

	function polynomial(coefficients) {
		const terms = coefficients
			.map((value, power) => ({ value, power }))
			.filter((term) => Math.abs(term.value) > 1e-14)
			.map(({ value, power }) => {
				const unit = Math.abs(value - 1) < 1e-12;
				const factor = unit && power > 0 ? '' : num(value, 6);
				const variable = power === 0 ? '' : power === 1 ? 'z' : `z^${power}`;
				return `${factor}${factor && variable ? ' ' : ''}${variable}`;
			});
		return terms.join('  +  ') || '0';
	}

	// A multistep method has no rational stability function, so two of the six
	// figures have nothing to draw. They are left out rather than shown empty.
	const hasRationalFunction = $derived(method && method.class !== 'linear_multistep');

	/** One coefficient, exactly as the file holds it. */
	function entry(value) {
		if (value.value === 0) return '0';
		return value.exact ? value.text : num(value.value, 8);
	}
</script>

{#snippet controls()}
	<div class="rail-group">
		<p class="label mb-1">run on</p>
		{#each catalog.problems as problem}
			<button
				type="button"
				class="facet {ui.problem === problem.id ? 'facet-on' : ''}"
				onclick={() => (ui.problem = problem.id)}
			>
				<span class="truncate">{problem.id}</span>
			</button>
		{/each}
	</div>
	<div class="rail-group">
		<a href="/" class="rail-link">back to the library</a>
	</div>
{/snippet}

{#snippet numbers()}
	{#each facts as [group, rows]}
		<div class="border-t border-cream/10 pt-3">
			<p class="label mb-2">{group}</p>
			{#each rows as [name, value]}
				<div class="flex items-baseline justify-between gap-3 py-[2px] font-mono text-xs">
					<span class="truncate text-cream/50">{name}</span>
					<span class="shrink-0 text-cream">{value}</span>
				</div>
			{/each}
		</div>
	{/each}
{/snippet}

{#snippet section(title, body)}
	<section class="mt-9 border-t border-cream/10 pt-4">
		<p class="label mb-4">{title}</p>
		{@render body()}
	</section>
{/snippet}

{#if !method}
	<p class="label">unknown method: {id}</p>
{:else}
	<div class="flex gap-12">
		<div class="min-w-0 flex-1">
			<header>
				<div class="flex items-baseline justify-between gap-6">
					<h1 class="font-display text-lg text-cream">{method.name}</h1>
					<span class="shrink-0 font-mono text-xs text-accent">{method.id}</span>
				</div>
				{#if structureLine}
					<p class="label mt-2">{structureLine}</p>
				{/if}
				{#each detail?.report.discrepancies ?? [] as issue}
					<p class="mt-2 font-mono text-xs text-[#d9513c]">{issue}</p>
				{/each}
			</header>

			<!-- Below the width where the rail fits beside the figures, the same
			     groups run across instead of down. -->
			<div class="mt-8 grid items-start gap-x-10 gap-y-6 sm:grid-cols-2 lg:grid-cols-3 xl:hidden">
				{@render numbers()}
			</div>

			{#key method.id}
				{#snippet stability()}
					<div class="grid gap-x-10 gap-y-8 lg:grid-cols-2">
						<Panel
							label="stability region"
							height="h-[20rem]"
							load={panel(MODES.stability, method, ui.problem)}
						/>
						{#if hasRationalFunction}
							<Panel
								label="order star"
								height="h-[20rem]"
								load={panel(MODES.orderStar, method, ui.problem)}
							/>
						{/if}
						<Panel
							label="damping"
							height="h-[20rem]"
							load={panel(MODES.damping, method, ui.problem)}
						/>
						{#if hasRationalFunction}
							<Panel
								label="error coefficients"
								height="h-[20rem]"
								load={panel(errorCoefficients, method, ui.problem)}
							/>
						{/if}
					</div>
				{/snippet}
				{@render section('stability', stability)}

				{#snippet accuracy()}
					<div class="grid gap-x-10 gap-y-8 lg:grid-cols-2">
						<Panel
							label="convergence on {ui.problem}"
							height="h-[20rem]"
							load={panel(MODES.order, method, ui.problem)}
						/>
						<Panel
							label="work precision on {ui.problem}"
							height="h-[20rem]"
							load={panel(MODES.cost, method, ui.problem)}
						/>
					</div>
				{/snippet}
				{@render section('accuracy', accuracy)}
			{/key}

			{#if detail}
				{#snippet coefficients()}
					<div class="overflow-x-auto">
						{#if detail.coefficients.kind === 'runge_kutta'}
							<table class="border-collapse font-mono text-xs">
								<tbody>
									{#each detail.coefficients.a as row, i}
										<tr>
											<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent">
												{entry(detail.coefficients.c[i])}
											</td>
											{#each row as value}
												<td
													class="py-[2px] pl-4 text-right {value.value === 0
														? 'text-cream/20'
														: 'text-cream'}">{entry(value)}</td
												>
											{/each}
										</tr>
									{/each}
									<tr class="border-t border-cream/20">
										<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent">b</td>
										{#each detail.coefficients.b as value}
											<td
												class="py-[2px] pl-4 text-right {value.value === 0
													? 'text-cream/20'
													: 'text-cream'}">{entry(value)}</td
											>
										{/each}
									</tr>
									{#if detail.coefficients.bEmbedded}
										<tr>
											<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent"
												>b hat</td
											>
											{#each detail.coefficients.bEmbedded as value}
												<td
													class="py-[2px] pl-4 text-right {value.value === 0
														? 'text-cream/20'
														: 'text-cream'}">{entry(value)}</td
												>
											{/each}
										</tr>
									{/if}
								</tbody>
							</table>
						{:else if detail.coefficients.kind === 'rosenbrock'}
							<table class="border-collapse font-mono text-xs">
								<tbody>
									{#each detail.coefficients.alpha as row, i}
										<tr>
											<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent"
												>{i === 0 ? 'alpha' : ''}</td
											>
											{#each row as value}
												<td
													class="py-[2px] pl-4 text-right {value.value === 0
														? 'text-cream/20'
														: 'text-cream'}">{entry(value)}</td
												>
											{/each}
										</tr>
									{/each}
									{#each detail.coefficients.gamma as row, i}
										<tr class={i === 0 ? 'border-t border-cream/20' : ''}>
											<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent"
												>{i === 0 ? 'gamma' : ''}</td
											>
											{#each row as value}
												<td
													class="py-[2px] pl-4 text-right {value.value === 0
														? 'text-cream/20'
														: 'text-cream'}">{entry(value)}</td
												>
											{/each}
										</tr>
									{/each}
									<tr class="border-t border-cream/20">
										<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent">b</td>
										{#each detail.coefficients.b as value}
											<td
												class="py-[2px] pl-4 text-right {value.value === 0
													? 'text-cream/20'
													: 'text-cream'}">{entry(value)}</td
											>
										{/each}
									</tr>
									{#if detail.coefficients.bEmbedded}
										<tr>
											<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent"
												>b hat</td
											>
											{#each detail.coefficients.bEmbedded as value}
												<td
													class="py-[2px] pl-4 text-right {value.value === 0
														? 'text-cream/20'
														: 'text-cream'}">{entry(value)}</td
												>
											{/each}
										</tr>
									{/if}
								</tbody>
							</table>
						{:else}
							<table class="border-collapse font-mono text-xs">
								<tbody>
									<tr>
										<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent"
											>alpha</td
										>
										{#each detail.coefficients.alpha ?? [] as value}
											<td class="py-[2px] pl-4 text-right text-cream">{num(value, 6)}</td>
										{/each}
									</tr>
									<tr>
										<td class="border-r border-cream/20 py-[2px] pr-3 text-right text-accent">beta</td
										>
										{#each detail.coefficients.beta ?? [] as value}
											<td class="py-[2px] pl-4 text-right text-cream">{num(value, 6)}</td>
										{/each}
									</tr>
								</tbody>
							</table>
						{/if}
					</div>

					{#if stabilityFunction}
						<div class="mt-8 font-mono text-xs">
							<p class="text-cream/50">
								N(z) <span class="text-cream">= {polynomial(stabilityFunction.numerator)}</span>
							</p>
							<p class="mt-1 text-cream/50">
								D(z) <span class="text-cream">= {polynomial(stabilityFunction.denominator)}</span>
							</p>
						</div>
					{/if}

					<div class="mt-8">
						<Panel
							label="coefficient structure"
							height="h-[18rem]"
							load={panel(MODES.structure, method, ui.problem)}
						/>
					</div>
				{/snippet}
				{@render section(
					detail.coefficients.kind === 'linear_multistep'
						? 'coefficients on a uniform grid'
						: 'coefficients',
					coefficients
				)}

				{#if detail.references.length}
					{#snippet references()}
						<ul class="flex flex-col gap-4">
							{#each detail.references as reference}
								<li>
									{#if reference.link}
										<a
											href={reference.link}
											target="_blank"
											rel="noreferrer"
											class="text-xs text-cream underline decoration-cream/20 underline-offset-4 hover:decoration-cream"
											>{reference.title ?? reference.link}</a
										>
									{:else}
										<span class="text-xs text-cream">{reference.title ?? 'untitled'}</span>
									{/if}
									<p class="mt-1 font-mono text-xs text-accent">
										{[reference.authors, reference.year, reference.source]
											.filter(Boolean)
											.join(' · ')}
									</p>
								</li>
							{/each}
						</ul>
					{/snippet}
					{@render section('published in', references)}
				{/if}
			{/if}
		</div>

		<aside class="sticky top-6 hidden h-fit w-52 shrink-0 space-y-5 xl:block">
			{@render numbers()}
		</aside>
	</div>
{/if}
