<script>
	/*
	 * A Plotly figure that redraws in place when the figure it is given changes.
	 * `Plotly.react` diffs against what is already there, so switching a card
	 * from one mode to another reuses the same canvas instead of tearing it down.
	 */
	import { onDestroy } from 'svelte';
	import Plotly from '$lib/plot.js';

	let { figure = null, class: className = '' } = $props();

	let node;
	let drawn = false;

	$effect(() => {
		// Depend on the figure's identity alone. Plotly writes its own
		// bookkeeping into the objects it is handed, so anything reactive passed
		// in here would report those writes back as a change and loop; figures
		// are therefore always raw state, replaced rather than edited.
		const current = figure;
		if (!node || !current) return;
		Plotly.react(node, current.data, current.layout, current.config);
		drawn = true;
	});

	onDestroy(() => {
		if (node && drawn) Plotly.purge(node);
	});
</script>

<div bind:this={node} class={className}></div>
