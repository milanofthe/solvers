<script>
	/*
	 * A Plotly figure that redraws in place when the figure it is given changes.
	 * `Plotly.react` diffs against what is already there, so switching a card
	 * from one mode to another reuses the same canvas instead of tearing it down.
	 */
	import { onDestroy } from 'svelte';
	import Plotly from '$lib/plot.js';

	let { figure = null, class: className = '', onsize = null, onview = null } = $props();

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

	// Panning and zooming change what the figure has to cover, and the data it
	// was given covers what it used to cover. The new window is reported so the
	// owner can ask for the data that reaches it. Settling first, because a
	// scroll wheel emits one relayout per notch.
	let timer;
	let listening = false;

	function listen() {
		if (listening || !onview || typeof node?.on !== 'function') return;
		listening = true;
		node.on('plotly_relayout', () => {
			clearTimeout(timer);
			timer = setTimeout(() => {
				const layout = node?._fullLayout;
				if (!layout) return;
				onview({ re: [...layout.xaxis.range], im: [...layout.yaxis.range] });
			}, 220);
		});
	}

	$effect(() => {
		// Depend on the figure's identity alone. Plotly writes its own
		// bookkeeping into the objects it is handed, so anything reactive passed
		// in here would report those writes back as a change and loop; figures
		// are therefore always raw state, replaced rather than edited.
		const current = figure;
		if (!node || !current) return;
		// A plot only gains its event emitter once it has been drawn once, so
		// the listener is attached on the way out rather than on mount.
		Plotly.react(node, current.data, current.layout, current.config).then(listen);
		drawn = true;
	});

	onDestroy(() => {
		clearTimeout(timer);
		if (node && drawn) Plotly.purge(node);
	});
</script>

<div bind:this={node} class={className}></div>
