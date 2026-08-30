<script>
	/*
	 * A titled figure. The heading is a heading like every other one on the
	 * page, set in the interface's own type rather than drawn by the plotting
	 * library, and the panel takes care of asking the pipeline and of the state
	 * while it waits.
	 */
	import Plot from './Plot.svelte';
	import { isCancelled } from '$lib/engine.js';
	import { empty } from '$lib/figures.js';

	let { label, load, height = 'h-72', note = '' } = $props();

	let figure = $state.raw(null);
	// Outside the reactive graph on purpose: reading `figure` in the effect that
	// writes it is a loop, and all this decides is whether a resize blanks the
	// panel or redraws it in place.
	let drawn = false;

	// The figure is sampled to fit its container, so the container is measured
	// before it is asked for and again when it changes. Snapping to a step keeps
	// a dragged window from asking for a new grid on every frame.
	let box = $state.raw(null);

	function measure(next) {
		const snapped = {
			width: Math.round(next.width / 24) * 24,
			height: Math.round(next.height / 24) * 24
		};
		if (box && box.width === snapped.width && box.height === snapped.height) return;
		box = snapped;
	}

	$effect(() => {
		const request = load;
		const measured = box;
		if (!measured) return;
		let stale = false;
		if (!drawn) figure = empty('');
		Promise.resolve()
			.then(() => request(measured))
			.then((result) => {
				if (stale) return;
				figure = result;
				drawn = true;
			})
			.catch((error) => {
				if (!stale && !isCancelled(error)) figure = empty('unavailable');
			});
		return () => {
			stale = true;
		};
	});
</script>

<section>
	<p class="label mb-2">{label}</p>
	<div class="border border-cream/10 bg-charcoal-light">
		<Plot {figure} class="{height} w-full" onsize={measure} />
	</div>
	{#if note}
		<p class="mt-2 max-w-[70ch] text-xs text-cream/40">{note}</p>
	{/if}
</section>
