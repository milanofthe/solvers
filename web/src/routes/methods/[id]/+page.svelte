<script>
	/*
	 * One method in full.
	 *
	 * Laid out as the answers to four questions in order: what is it, how does it
	 * behave in the complex plane, how large is the error it leaves behind, and
	 * what is it made of. The coefficients come last on purpose; they are the
	 * reference, not the summary.
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
		return () =>
			mode
				.request(engine(), current, PRIORITY.immediate, { problem })
				.then((result) => mode.figure(result, current, { compact: false }));
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

		if (current.class === 'runge_kutta') {
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

	const properties = $derived.by(() => {
		if (!detail) return [];
		const report = detail.report;
		return [
			['order', report.computed_order],
			['embedded', report.computed_embedded_order ?? '–'],
			[report.class === 'runge_kutta' ? 'stages' : 'steps', report.size],
			['stage order', report.stage_order ?? '–'],
			['cost per step', report.stage_cost],
			['A-stable', report.a_stable ? 'yes' : 'no'],
			['L-stable', report.l_stable ? 'yes' : 'no'],
			[
				'stiffly accurate',
				report.stiffly_accurate == null ? '–' : report.stiffly_accurate ? 'yes' : 'no'
			],
			['R(inf)', report.damping_at_infinity == null ? '–' : limit(report.damping_at_infinity)],
			['A(alpha)', report.alpha_angle == null ? '–' : `${num(report.alpha_angle, 4)} deg`],
			['real limit', limit(report.real_stability_limit)],
			['imaginary limit', limit(report.imaginary_stability_limit)],
			['order check', report.exact_arithmetic ? 'exact' : 'numeric']
		];
	});

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

{#if !method}
	<p class="label">unknown method: {id}</p>
{:else}
	<header class="mb-6">
		<div class="flex items-baseline justify-between gap-4">
			<h1 class="font-display text-lg text-cream">{method.name}</h1>
			<span class="font-mono text-xs text-accent">{method.id}</span>
		</div>
		{#if structureLine}
			<p class="label mt-3">{structureLine}</p>
		{/if}
		{#each detail?.report.discrepancies ?? [] as issue}
			<p class="mt-1 font-mono text-xs text-[#d9513c]">{issue}</p>
		{/each}
	</header>

	{#if detail}
		<section class="mb-8">
			<p class="label mb-2">properties</p>
			<dl
				class="grid grid-cols-2 gap-x-6 font-mono text-xs sm:grid-cols-3 lg:grid-cols-4 xl:grid-cols-5"
			>
				{#each properties as [name, value]}
					<div class="flex items-baseline justify-between gap-2 border-b border-cream/10 py-[3px]">
						<dt class="truncate text-cream/50">{name}</dt>
						<dd class="text-cream">{value}</dd>
					</div>
				{/each}
			</dl>
		</section>
	{/if}

	{#key method.id}
		<div class="mb-8 grid gap-6 xl:grid-cols-2">
			<Panel
				label="stability region"
				height="h-[26rem]"
				load={panel(MODES.stability, method, ui.problem)}
				note="Banded log |R(z)| with the unit contour on top. The region it encloses is where the method does not amplify."
			/>
			<Panel
				label="order star"
				height="h-[26rem]"
				load={panel(MODES.orderStar, method, ui.problem)}
				note="log |R(z) exp(-z)|, with the amber line where the two agree in modulus. The sectors meeting at the origin count one more than the order."
			/>
		</div>

		<div class="mb-8 grid gap-6 lg:grid-cols-3">
			<Panel
				label="damping"
				height="h-64"
				load={panel(MODES.damping, method, ui.problem)}
				note="A decaying mode along the negative real axis."
			/>
			<Panel
				label="error coefficients"
				height="h-64"
				load={panel(errorCoefficients, method, ui.problem)}
				note="The residual of each order condition the method does not satisfy."
			/>
			<Panel
				label="convergence on {ui.problem}"
				height="h-64"
				load={panel(MODES.order, method, ui.problem)}
				note="Measured against the slope the coefficients promise."
			/>
		</div>

		<div class="mb-8 grid gap-6 lg:grid-cols-2">
			<Panel
				label="coefficient structure"
				height="h-64"
				load={panel(MODES.structure, method, ui.problem)}
				note="Teal positive, red negative, shaded by magnitude."
			/>
			<Panel
				label="work precision on {ui.problem}"
				height="h-64"
				load={panel(MODES.cost, method, ui.problem)}
				note="Accuracy improves to the right."
			/>
		</div>
	{/key}

	{#if detail}
		<section class="mb-8">
			<p class="label mb-2">
				{detail.coefficients.kind === 'runge_kutta'
					? 'butcher tableau'
					: 'coefficients on a uniform grid'}
			</p>
			<div class="overflow-x-auto border border-cream/10 bg-charcoal-light p-3">
				{#if detail.coefficients.kind === 'runge_kutta'}
					<table class="border-collapse font-mono text-xs">
						<tbody>
							{#each detail.coefficients.a as row, i}
								<tr>
									<td class="border-r border-cream/20 px-3 text-right text-accent">
										{detail.coefficients.c[i].exact
											? detail.coefficients.c[i].text
											: num(detail.coefficients.c[i].value, 8)}
									</td>
									{#each row as entry}
										<td class="px-3 text-right {entry.value === 0 ? 'text-cream/20' : 'text-cream'}">
											{entry.value === 0 ? '0' : entry.exact ? entry.text : num(entry.value, 8)}
										</td>
									{/each}
								</tr>
							{/each}
							<tr class="border-t border-cream/20">
								<td class="border-r border-cream/20 px-3 text-right text-accent">b</td>
								{#each detail.coefficients.b as entry}
									<td class="px-3 text-right {entry.value === 0 ? 'text-cream/20' : 'text-cream'}">
										{entry.value === 0 ? '0' : entry.exact ? entry.text : num(entry.value, 8)}
									</td>
								{/each}
							</tr>
							{#if detail.coefficients.bEmbedded}
								<tr>
									<td class="border-r border-cream/20 px-3 text-right text-accent">b hat</td>
									{#each detail.coefficients.bEmbedded as entry}
										<td class="px-3 text-right {entry.value === 0 ? 'text-cream/20' : 'text-cream'}">
											{entry.value === 0 ? '0' : entry.exact ? entry.text : num(entry.value, 8)}
										</td>
									{/each}
								</tr>
							{/if}
						</tbody>
					</table>
				{:else}
					<table class="border-collapse font-mono text-xs">
						<tbody>
							<tr>
								<td class="border-r border-cream/20 px-3 text-right text-accent">alpha</td>
								{#each detail.coefficients.alpha ?? [] as value}
									<td class="px-3 text-right text-cream">{num(value, 6)}</td>
								{/each}
							</tr>
							<tr>
								<td class="border-r border-cream/20 px-3 text-right text-accent">beta</td>
								{#each detail.coefficients.beta ?? [] as value}
									<td class="px-3 text-right text-cream">{num(value, 6)}</td>
								{/each}
							</tr>
						</tbody>
					</table>
				{/if}
			</div>
		</section>

		{#if stabilityFunction}
			<section class="mb-8">
				<p class="label mb-2">stability function</p>
				<p class="font-mono text-xs text-cream/80">
					N(z) = {polynomial(stabilityFunction.numerator)}
				</p>
				<p class="font-mono text-xs text-cream/80">
					D(z) = {polynomial(stabilityFunction.denominator)}
				</p>
			</section>
		{/if}

		{#if detail.references.length}
			<section class="mb-8">
				<p class="label mb-2">references</p>
				<ul class="flex flex-col gap-3">
					{#each detail.references as reference}
						<li class="border-l border-cream/20 pl-3">
							{#if reference.link}
								<a
									href={reference.link}
									target="_blank"
									rel="noreferrer"
									class="text-xs text-cream hover:text-amber">{reference.title ?? reference.link}</a
								>
							{:else}
								<span class="text-xs text-cream">{reference.title ?? 'untitled'}</span>
							{/if}
							<p class="font-mono text-xs text-accent">
								{[reference.authors, reference.year, reference.source].filter(Boolean).join(' · ')}
							</p>
						</li>
					{/each}
				</ul>
			</section>
		{/if}
	{/if}
{/if}
