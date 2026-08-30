<script>
	/*
	 * The frame: a quiet rail on the left carrying the navigation and whatever
	 * controls the current page needs, and the content to the right of it. No
	 * tabs, no chrome around the controls; a control is a line of monospace that
	 * turns from grey to cream when it is on.
	 */
	import '../app.css';
	import { setContext } from 'svelte';
	import { page } from '$app/state';
	import { catalog, loadCatalogues, progress } from '$lib/store.svelte.js';

	let { children } = $props();

	let settled = $state(false);

	// Pages hand their own controls to the rail as a snippet, and a page with
	// numbers worth keeping in view hands those to a second rail on the other
	// side. Both are siblings of the content rather than part of it, so their
	// rules run the height of the window like the navigation does.
	const rail = $state({ panel: null, right: null });
	setContext('rail', rail);

	loadCatalogues();

	const links = [
		{ href: '/', label: 'methods' },
		{ href: '/problems', label: 'problems' },
		{ href: '/compare', label: 'compare' }
	];

	// After the first paint, so the frame fades in rather than appearing. The
	// loading line is inside it and fades in too, which is the point: there is
	// always something on screen, it just is not abrupt.
	$effect(() => {
		settled = true;
	});

	function active(href) {
		if (href === '/') return page.url.pathname === '/' || page.url.pathname.startsWith('/methods');
		return page.url.pathname.startsWith(href);
	}
</script>

<div class="flex min-h-screen w-full">
	<aside
		class="sticky top-0 hidden h-screen w-72 shrink-0 flex-col overflow-y-auto border-r border-cream/10 px-5 py-6 md:flex"
	>
		<a href="/" class="font-display text-lg leading-none text-cream">solvers</a>
		<p class="label mt-2">numerical integration</p>

		<nav class="rail-group flex flex-col">
			{#each links as link}
				<a
					href={link.href}
					class="py-[1px] font-mono text-xs transition-colors {active(link.href)
						? 'text-cream'
						: 'text-cream/50 hover:text-cream'}">{link.label}</a
				>
			{/each}
		</nav>

		{#if rail.panel}
			<div class="flex flex-col">
				{@render rail.panel()}
			</div>
		{/if}

		<div class="rail-group mt-auto">
			<p class="font-mono text-xs text-accent">
				{catalog.ready ? `${catalog.methods.length} methods` : 'loading'}
				{#if progress.outstanding > 0}
					<span class="text-amber">&middot; {progress.outstanding} computing</span>
				{/if}
			</p>
			<a
				href="https://github.com/milanofthe/solvers"
				class="rail-link mt-2"
				target="_blank"
				rel="noreferrer">source</a
			>
			<a href="/api/methods.json" class="rail-link" download>download as json</a>
		</div>
	</aside>

	<main class="page min-w-0 flex-1 px-6 py-6 {settled ? 'page-in' : ''}">
		{#if catalog.error}
			<p class="font-mono text-xs text-[#d9513c]">{catalog.error}</p>
		{:else if !catalog.ready}
			<p class="label">starting the compute pipeline</p>
		{:else}
			{@render children()}
		{/if}
	</main>

	{#if rail.right}
		<aside
			class="sticky top-0 hidden h-screen w-60 shrink-0 overflow-y-auto border-l border-cream/10 px-5 py-6 xl:block"
		>
			{@render rail.right()}
		</aside>
	{/if}
</div>
