<script>
	/*
	 * One piece of mathematics, set by KaTeX.
	 *
	 * The library and its stylesheet are imported here rather than globally, so
	 * they are part of whichever route uses them and cost the rest of the site
	 * nothing.
	 */
	import katex from 'katex';
	import 'katex/dist/katex.min.css';

	let { formula, display = false } = $props();

	const rendered = $derived(
		katex.renderToString(formula, {
			displayMode: display,
			throwOnError: false,
			output: 'html'
		})
	);
</script>

<span class={display ? 'block overflow-x-auto py-1 text-cream' : 'text-cream'}>
	<!-- eslint-disable-next-line svelte/no-at-html-tags -->
	{@html rendered}
</span>
