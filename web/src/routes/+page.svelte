<script>
	/*
	 * The library.
	 *
	 * Fifty four methods is too many to read as a table and too few to need a
	 * search engine. What it needs is faceting: narrow by the properties that
	 * decide a choice, and see the consequence of that choice drawn on every
	 * card at once.
	 */
	import { flip } from 'svelte/animate';
	import { cubicOut } from 'svelte/easing';
	import { getContext } from 'svelte';
	import { page } from '$app/state';
	import MethodCard from '$lib/components/MethodCard.svelte';
	import { catalog, engine, ui, toggleTag } from '$lib/store.svelte.js';
	import { METHOD_MODES, MODES_NEEDING_A_PROBLEM } from '$lib/figures.js';

	const rail = getContext('rail');

	const SORTS = [
		{ key: 'family', label: 'family', of: (m) => `${m.family}/${m.id}` },
		{ key: 'name', label: 'name', of: (m) => m.name.toLowerCase() },
		{ key: 'order', label: 'order', of: (m) => -m.order },
		// Undated entries sort to the end rather than to 1900.
		{ key: 'year', label: 'year', of: (m) => m.year ?? 9999 },
		{ key: 'size', label: 'size', of: (m) => m.size },
		{ key: 'cost', label: 'cost', of: (m) => m.stageCost },
		// Undated and unmeasured entries sort to the end rather than to zero.
		{ key: 'error', label: 'error constant', of: (m) => m.errorConstant ?? Infinity },
		{
			key: 'ssp',
			label: 'SSP coefficient',
			of: (m) => -(m.sspCoefficient === 'unbounded' ? 1e9 : (m.sspCoefficient ?? 0))
		},
		{
			key: 'stability',
			label: 'stability',
			of: (m) => -(m.lStable ? 200 : m.aStable ? 100 : (m.alphaAngle ?? 0))
		}
	];

	const GROUP_ORDER = [
		'structure',
		'stability',
		'order',
		'control',
		'use',
		'geometry',
		'provenance',
		'other'
	];

	const groups = $derived.by(() => {
		const counts = new Map();
		for (const method of catalog.methods) {
			for (const tag of method.tags ?? []) {
				if (!counts.has(tag.group)) counts.set(tag.group, new Map());
				const group = counts.get(tag.group);
				group.set(tag.name, (group.get(tag.name) ?? 0) + 1);
			}
		}
		return [...counts.entries()]
			.sort((a, b) => GROUP_ORDER.indexOf(a[0]) - GROUP_ORDER.indexOf(b[0]))
			.map(([group, names]) => ({
				group,
				names: [...names.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
			}));
	});

	const groupOf = $derived.by(() => {
		const index = new Map();
		for (const { group, names } of groups) {
			for (const [name] of names) index.set(name, group);
		}
		return index;
	});

	// A filter can arrive in the address, which is what lets the reference page
	// link a definition to the methods that satisfy it. Applied once: after that
	// the rail owns the selection, and a back navigation should not undo it.
	let arrived = false;
	$effect(() => {
		if (arrived) return;
		arrived = true;
		const wanted = page.url.searchParams.get('tag');
		const within = page.url.searchParams.get('family');
		if (wanted) ui.tags = [wanted];
		if (within) ui.query = within;
	});

	const shown = $derived.by(() => {
		const query = ui.query.trim().toLowerCase();
		const wanted = new Map();
		for (const name of ui.tags) {
			const group = groupOf.get(name) ?? 'other';
			if (!wanted.has(group)) wanted.set(group, []);
			wanted.get(group).push(name);
		}

		const rows = catalog.methods.filter((method) => {
			const names = new Set((method.tags ?? []).map((t) => t.name));
			// Any of the chosen tags inside a group, every group together.
			for (const [, group] of wanted) {
				if (!group.some((name) => names.has(name))) return false;
			}
			if (!query) return true;
			return [method.id, method.name, method.family, ...names]
				.join(' ')
				.toLowerCase()
				.includes(query);
		});

		const sort = SORTS.find((s) => s.key === ui.sort) ?? SORTS[0];
		return [...rows].sort((a, b) => {
			const x = sort.of(a);
			const y = sort.of(b);
			if (x === y) return a.id.localeCompare(b.id);
			return x > y ? 1 : -1;
		});
	});

	// Switching mode makes everything still queued for the old one pointless.
	function setMode(key) {
		if (ui.mode === key) return;
		engine().drop(() => true);
		ui.mode = key;
	}

	$effect(() => {
		rail.panel = controls;
		return () => {
			rail.panel = null;
		};
	});
</script>

{#snippet controls()}
	<div class="rail-group">
		<p class="label mb-1">select</p>
		<button
			type="button"
			class="facet {ui.selecting ? 'facet-on' : ''}"
			onclick={() => (ui.selecting = !ui.selecting)}
		>
			<span>{ui.selecting ? 'selecting' : 'off'}</span>
			<span class="text-accent">{ui.selection.length || ''}</span>
		</button>
		{#if ui.selection.length > 0}
			<button type="button" class="facet" onclick={() => (ui.selection = [])}>
				<span>clear</span>
			</button>
			<a href="/compare" class="facet facet-on"><span>compare</span></a>
		{/if}
	</div>

	<div class="rail-group">
		<p class="label mb-1">show</p>
		{#each METHOD_MODES as entry}
			<button type="button" class="facet {ui.mode === entry.key ? 'facet-on' : ''}" onclick={() => setMode(entry.key)}>
				<span>{entry.label}</span>
			</button>
		{/each}
	</div>

	{#if MODES_NEEDING_A_PROBLEM.has(ui.mode)}
		<div class="rail-group">
			<p class="label mb-1">run on</p>
			{#each catalog.problems as problem}
				<button
					type="button"
					class="facet {ui.problem === problem.id ? 'facet-on' : ''}"
					onclick={() => {
						engine().drop(() => true);
						ui.problem = problem.id;
					}}
				>
					<span class="truncate">{problem.id}</span>
				</button>
			{/each}
		</div>
	{/if}

	<div class="rail-group">
		<p class="label mb-1">sort</p>
		{#each SORTS as sort}
			<button type="button" class="facet {ui.sort === sort.key ? 'facet-on' : ''}" onclick={() => (ui.sort = sort.key)}>
				<span>{sort.label}</span>
			</button>
		{/each}
	</div>

	<div class="rail-group">
		<p class="label mb-1">search</p>
		<input
			type="search"
			bind:value={ui.query}
			placeholder="name, family, tag"
			class="w-full border border-cream/10 bg-charcoal px-2 py-1 font-mono text-xs text-cream placeholder:text-cream/30 focus:border-cream/30 focus:outline-none"
		/>
	</div>

	{#each groups as { group, names }}
		<div class="rail-group">
			<p class="label mb-1">{group}</p>
			{#each names as [name, count]}
				<button
					type="button"
					class="facet {ui.tags.includes(name) ? 'facet-on' : ''}"
					onclick={() => toggleTag(name)}
				>
					<span class="truncate">{name}</span>
					<span class="text-accent">{count}</span>
				</button>
			{/each}
		</div>
	{/each}

	{#if ui.tags.length > 0}
		<button type="button" class="rail-link" onclick={() => (ui.tags = [])}>clear filters</button>
	{/if}
{/snippet}

<header class="mb-4">
	<p class="label">
		{shown.length} of {catalog.total || catalog.methods.length}
		<!-- The library arrives in slices, so while they are still landing the
		     count says how far along it is rather than pretending the grid is
		     already the whole thing. -->
		{#if catalog.total && catalog.loaded < catalog.total}
			<span class="text-cream/40">&middot; deriving {catalog.loaded}/{catalog.total}</span>
		{/if}
		{#if ui.selection.length > 0}
			<span class="text-amber">&middot; {ui.selection.length} selected</span>
		{/if}
	</p>
</header>

<!-- Keyed by id, so re-sorting moves the tiles that are already there rather
     than rebuilding them, and `flip` takes them to their new places instead of
     letting them jump. The wrapper exists because an animation goes on an
     element and a card is a component. -->
<div class="grid grid-cols-[repeat(auto-fill,minmax(17rem,1fr))] gap-4">
	{#each shown as method (method.id)}
		<div animate:flip={{ duration: 320, easing: cubicOut }}>
			<MethodCard {method} />
		</div>
	{/each}
</div>
