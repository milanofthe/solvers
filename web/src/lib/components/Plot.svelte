<script>
	/*
	 * A Plotly figure that redraws in place when the figure it is given changes.
	 * `Plotly.react` diffs against what is already there, so switching a card
	 * from one mode to another reuses the same canvas instead of tearing it down.
	 */
	import { onDestroy } from 'svelte';
	import Plotly from '$lib/plot.js';

	let { figure = null, class: className = '', onsize = null } = $props();

	let node;
	let drawn = false;

	// What the figure will occupy, reported to whoever asked for it. A figure
	// of the complex plane is sampled to fit its container, so the container has
	// to say how large it is before the sampling and again whenever it changes.
	$effect(() => {
		if (!node || !onsize) return;
		const observer = new ResizeObserver(([entry]) => {
			const { width, height } = entry.contentRect;
			if (width > 0 && height > 0) onsize({ width, height });
		});
		observer.observe(node);
		return () => observer.disconnect();
	});

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
