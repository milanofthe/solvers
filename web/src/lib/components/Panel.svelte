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

	$effect(() => {
		const request = load;
		let stale = false;
		figure = empty('');
		Promise.resolve()
			.then(() => request())
			.then((result) => {
				if (!stale) figure = result;
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
		<Plot {figure} class="{height} w-full" />
	</div>
	{#if note}
		<p class="mt-2 max-w-[70ch] text-xs text-cream/40">{note}</p>
	{/if}
</section>
