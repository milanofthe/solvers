<script>
	/*
	 * One method in full: the coefficients it is made of, what those imply, its
	 * stability region at a size worth looking at, and where it was published.
	 */
	import { getContext } from 'svelte';
	import { page } from '$app/state';
	import Plot from '$lib/components/Plot.svelte';
	import { catalog, engine, ui } from '$lib/store.svelte.js';
	import { isCancelled, PRIORITY } from '$lib/engine.js';
	import { METHOD_MODES } from '$lib/figures.js';
	import { limit, num } from '$lib/format.js';

	const rail = getContext('rail');

	const id = $derived(page.params.id);
	const method = $derived(catalog.methods.find((m) => m.id === id));

	let detail = $state.raw(null);
	let figure = $state.raw(null);
	let stabilityFunction = $state.raw(null);
	let panelMode = $state('stability');

	$effect(() => {
		const current = method;
		if (!current) return;
		let stale = false;
		detail = null;
		figure = null;
		stabilityFunction = null;

		engine()
			.request('methodDetail', { id: current.id }, { key: `detail:${current.id}`, priority: PRIORITY.immediate })
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

	// The large panel, in whichever mode is chosen in the rail.
	$effect(() => {
		const current = method;
		const mode = METHOD_MODES.find((m) => m.key === panelMode) ?? METHOD_MODES[0];
		const problem = ui.problem;
		if (!current) return;
		let stale = false;
		mode
			.request(engine(), current, PRIORITY.immediate, { problem })
			.then((result) => {
				if (!stale) figure = mode.figure(result, current, { compact: false });
			})
			.catch((error) => {
				if (!stale && !isCancelled(error)) figure = null;
			});
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
		<p class="label mb-1">panel</p>
		{#each METHOD_MODES as entry}
			<button
				type="button"
				class="facet {panelMode === entry.key ? 'facet-on' : ''}"
				onclick={() => (panelMode = entry.key)}
			>
				<span>{entry.label}</span>
			</button>
		{/each}
	</div>
	<a href="/" class="action">[ back to the library ]</a>
{/snippet}

{#if !method}
	<p class="label">unknown method: {id}</p>
{:else}
	<header class="mb-6">
		<div class="flex items-baseline justify-between gap-4">
			<h1 class="font-display text-lg text-cream">{method.name}</h1>
			<span class="font-mono text-xs text-accent">{method.id}</span>
		</div>
		{#if method.description}
			<p class="mt-2 max-w-[72ch] text-xs text-cream/60">{method.description}</p>
		{/if}
		{#if detail}
			<p class="label mt-3">
				{#if detail.coefficients.kind === 'runge_kutta'}
					{detail.coefficients.structure?.replace('_', ' ')}
					{#if detail.coefficients.singlyDiagonal}&middot; singly diagonal{/if}
					{#if detail.coefficients.explicitFirstStage}&middot; explicit first stage{/if}
					{#if detail.coefficients.fsal}&middot; FSAL{/if}
				{:else}
					{detail.coefficients.steps} steps &middot;
					{detail.coefficients.startup
						? `starts with ${detail.coefficients.startup}`
						: 'self starting'}
				{/if}
			</p>
			{#each detail.report.discrepancies as issue}
				<p class="mt-1 font-mono text-xs text-[#d9513c]">{issue}</p>
			{/each}
		{/if}
	</header>

	<div class="grid gap-6 lg:grid-cols-[1fr_minmax(20rem,26rem)]">
		<div class="border border-cream/10 bg-charcoal-light">
			<Plot {figure} class="h-[30rem] w-full" />
		</div>

		<div>
			{#if detail}
				<p class="label mb-2">properties</p>
				<dl class="grid grid-cols-2 gap-x-4 font-mono text-xs">
					{#each properties as [name, value]}
						<div class="flex items-baseline justify-between border-b border-cream/10 py-[2px]">
							<dt class="text-cream/50">{name}</dt>
							<dd class="text-cream">{value}</dd>
						</div>
					{/each}
				</dl>
			{/if}
		</div>
	</div>

	{#if detail}
		<section class="mt-8">
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
			<p class="mt-2 max-w-[72ch] text-xs text-cream/50">
				{detail.coefficients.kind === 'runge_kutta'
					? 'An exact fraction is shown as the file writes it. A decimal means the coefficient is only known as a double, and the order conditions were then checked numerically rather than as identities.'
					: 'A multistep file stores which coefficients are free, not their values. These are the values on a uniform grid; on a varying step size they are solved for again at every step.'}
			</p>
		</section>

		{#if stabilityFunction}
			<section class="mt-8">
				<p class="label mb-2">stability function</p>
				<p class="font-mono text-xs text-cream/80">N(z) = {polynomial(stabilityFunction.numerator)}</p>
				<p class="font-mono text-xs text-cream/80">
					D(z) = {polynomial(stabilityFunction.denominator)}
				</p>
				<p class="mt-2 max-w-[72ch] text-xs text-cream/50">
					R agrees with exp to order {stabilityFunction.orderOfConsistency}. R(inf) = {limit(
						stabilityFunction.atInfinity
					)}.
					{#if stabilityFunction.poles.length}
						Poles at {stabilityFunction.poles
							.map((p) => `${num(p.re, 4)}${p.im >= 0 ? '+' : ''}${num(p.im, 4)}i`)
							.join(', ')}.
					{:else}
						No poles, the method is explicit.
					{/if}
				</p>
			</section>
		{/if}

		{#if detail.references.length}
			<section class="mt-8">
				<p class="label mb-2">references</p>
				<ul class="flex flex-col gap-3">
					{#each detail.references as reference}
						<li class="border-l border-cream/20 pl-3">
							{#if reference.link}
								<a href={reference.link} target="_blank" rel="noreferrer" class="text-xs text-cream hover:text-amber"
									>{reference.title ?? reference.link}</a
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
