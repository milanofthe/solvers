<script>
	/*
	 * One reference problem. Selecting it makes it the problem the method cards
	 * and the comparison run on, which is the only reason the problems need a
	 * grid of their own.
	 */
	import Plot from './Plot.svelte';
	import { viewport } from '$lib/viewport.js';
	import { engine, ui } from '$lib/store.svelte.js';
	import { isCancelled, PRIORITY } from '$lib/engine.js';
	import { PROBLEM_MODES, empty } from '$lib/figures.js';
	import { num } from '$lib/format.js';

	let { problem } = $props();

	let visible = $state(false);
	let figure = $state.raw(null);
	let arrived = $state(false);

	const selected = $derived(ui.problem === problem.id);
	const mode = $derived(PROBLEM_MODES.find((m) => m.key === ui.problemMode) ?? PROBLEM_MODES[0]);

	$effect(() => {
		const active = mode;
		if (!visible) return;
		let stale = false;
		figure = empty('…');
		active
			.request(engine(), problem, PRIORITY.visible)
			.then((profile) => {
				if (stale) return;
				figure = active.figure(profile, problem, { compact: true });
				arrived = true;
			})
			.catch((error) => {
				if (stale || isCancelled(error)) return;
				figure = empty('unavailable');
				arrived = true;
			});
		return () => {
			stale = true;
		};
	});
</script>

<article
	use:viewport={(value) => (visible = value)}
	class="tile flex aspect-square flex-col overflow-hidden border bg-charcoal-light {arrived
		? 'tile-in'
		: ''} {selected
		? 'border-amber'
		: 'border-cream/10 hover:border-cream/25'}"
>
	<header class="flex items-baseline justify-between gap-2 px-3 pt-2">
		<h3 class="truncate font-display text-sm text-cream" title={problem.name}>{problem.name}</h3>
		<span class="shrink-0 font-mono text-xs {selected ? 'text-amber' : 'text-accent'}"
			>{problem.id}</span
		>
	</header>

	<button
		type="button"
		class="mt-2 flex min-h-0 flex-1 cursor-pointer flex-col text-left"
		onclick={() => (ui.problem = problem.id)}
		aria-pressed={selected}
		aria-label="use {problem.name}"
	>
		<Plot {figure} class="min-h-0 w-full flex-1" />
	</button>

	<p class="flex flex-wrap gap-x-3 px-3 pt-2 font-mono text-xs text-cream/50">
		<span>dim <span class="text-cream">{problem.dim}</span></span>
		<span>t <span class="text-cream">{num(problem.tStart, 3)}…{num(problem.tEnd, 3)}</span></span>
		<span>{problem.stiff ? 'stiff' : 'non stiff'}</span>
		{#if problem.hasExact}<span>closed form</span>{/if}
	</p>

	<div class="px-3 pb-3 pt-2">
		<p class="line-clamp-2 text-xs text-cream/50">{problem.description}</p>
	</div>
</article>
