<script>
	/*
	 * One method, drawn in whichever mode the library is currently showing.
	 *
	 * The card does not compute anything itself. It asks the pipeline when it
	 * comes into view, and asks again when the mode changes; anything it asked
	 * for and no longer wants is dropped by the pipeline rather than cancelled
	 * here.
	 */
	import Plot from './Plot.svelte';
	import { viewport } from '$lib/viewport.js';
	import { engine, ui, toggleSelection } from '$lib/store.svelte.js';
	import { isCancelled, PRIORITY } from '$lib/engine.js';
	import { METHOD_MODES, MODES_NEEDING_A_PROBLEM, empty } from '$lib/figures.js';
	import { limit, num, stabilityLabel } from '$lib/format.js';

	let { method } = $props();

	let visible = $state(false);
	let figure = $state.raw(null);

	const selected = $derived(ui.selection.includes(method.id));
	const mode = $derived(METHOD_MODES.find((m) => m.key === ui.mode) ?? METHOD_MODES[0]);
	const context = $derived({ problem: MODES_NEEDING_A_PROBLEM.has(ui.mode) ? ui.problem : null });

	$effect(() => {
		const active = mode;
		const problem = context.problem;
		if (!visible) return;
		let stale = false;
		figure = empty('…');
		active
			.request(engine(), method, PRIORITY.visible, { problem })
			.then((result) => {
				if (stale) return;
				figure = active.figure(result, method, { compact: true });
			})
			.catch((error) => {
				if (stale || isCancelled(error)) return;
				figure = empty('unavailable');
			});
		return () => {
			stale = true;
		};
	});

	// The facts line already carries the numbers and the stability, so the card
	// shows only the two tags that say something neither of those does.
	const INTERESTING = [
		'FSAL',
		'symplectic',
		'symmetric',
		'strong stability preserving',
		'stiffly accurate',
		'singly diagonal',
		'fully implicit',
		'explicit first stage',
		'multistep',
		'explicit'
	];

	const shownTags = $derived.by(() => {
		const names = new Set((method.tags ?? []).map((tag) => tag.name));
		return INTERESTING.filter((name) => names.has(name)).slice(0, 2);
	});
</script>

<article
	use:viewport={(value) => (visible = value)}
	class="flex flex-col border bg-charcoal-light transition-colors {selected
		? 'border-amber'
		: 'border-cream/10 hover:border-cream/25'}"
>
	<!-- The card is the method, so it links to the method. While the library is
	     in selection mode the same click picks it for a comparison instead; the
	     href stays put, which keeps opening it in a new tab working either way. -->
	<a
		href="/methods/{method.id}"
		class="block"
		data-selected={selected}
		onclick={(event) => {
			if (!ui.selecting) return;
			event.preventDefault();
			toggleSelection(method.id);
		}}
	>
		<header class="flex items-baseline justify-between gap-2 px-3 pt-2">
			<h3 class="truncate font-display text-sm text-cream" title={method.name}>{method.name}</h3>
			<span class="shrink-0 font-mono text-xs {selected ? 'text-amber' : 'text-accent'}"
				>{method.id}</span
			>
		</header>

		<Plot {figure} class="mt-2 h-40 w-full" />

		<dl class="flex flex-wrap gap-x-3 px-3 pt-2 font-mono text-xs text-cream/50">
			<div><dt class="inline">order</dt> <dd class="inline text-cream">{method.order}</dd></div>
			<div>
				<dt class="inline">{method.class === 'runge_kutta' ? 'stages' : 'steps'}</dt>
				<dd class="inline text-cream">{method.size}</dd>
			</div>
			<div><dt class="inline">cost</dt> <dd class="inline text-cream">{method.stageCost}</dd></div>
			{#if method.dampingAtInfinity != null}
				<div>
					<dt class="inline">R&#8734;</dt>
					<dd class="inline text-cream">{limit(method.dampingAtInfinity)}</dd>
				</div>
			{/if}
		</dl>

		<p class="px-3 pb-3 pt-1 font-mono text-xs text-accent">
			{stabilityLabel(method)}{#each shownTags as name}<span class="text-cream/25"> / </span>{name}{/each}
		</p>

		{#if method.discrepancies.length}
			<p class="px-3 pb-3 font-mono text-xs text-[#d9513c]">{method.discrepancies[0]}</p>
		{/if}
	</a>
</article>
