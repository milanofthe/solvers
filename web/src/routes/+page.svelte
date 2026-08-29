<script>
	/*
	 * The library.
	 *
	 * Fifty four methods is too many to read as a table and too few to need a
	 * search engine. What it needs is faceting: narrow by the properties that
	 * decide a choice, and see the consequence of that choice drawn on every
	 * card at once.
	 */
	import { getContext } from 'svelte';
	import MethodCard from '$lib/components/MethodCard.svelte';
	import { catalog, engine, ui, toggleTag } from '$lib/store.svelte.js';
	import { METHOD_MODES, MODES_NEEDING_A_PROBLEM } from '$lib/figures.js';

	const rail = getContext('rail');

	const SORTS = [
		{ key: 'family', label: 'family', of: (m) => `${m.family}/${m.id}` },
		{ key: 'name', label: 'name', of: (m) => m.name.toLowerCase() },
		{ key: 'order', label: 'order', of: (m) => -m.order },
		{ key: 'size', label: 'size', of: (m) => m.size },
		{ key: 'cost', label: 'cost', of: (m) => m.stageCost },
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
			return [method.id, method.name, method.family, method.description ?? '', ...names]
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
		<button type="button" class="action text-left" onclick={() => (ui.tags = [])}>[ clear filters ]</button>
	{/if}
{/snippet}

<header class="mb-4">
	<p class="label">
		{shown.length} of {catalog.methods.length}
		{#if ui.selection.length > 0}
			<span class="text-amber">&middot; {ui.selection.length} selected</span>
		{/if}
	</p>
</header>

<div class="grid grid-cols-[repeat(auto-fill,minmax(17rem,1fr))] gap-4">
	{#each shown as method (method.id)}
		<MethodCard {method} />
	{/each}
</div>

{#if ui.selection.length > 0}
	<div
		class="sticky bottom-4 mt-6 flex items-center gap-4 border border-cream/15 bg-charcoal-light px-3 py-2"
	>
		<span class="font-mono text-xs text-cream">{ui.selection.length} selected</span>
		<span class="truncate font-mono text-xs text-accent">{ui.selection.join('  ')}</span>
		<button type="button" class="action ml-auto" onclick={() => (ui.selection = [])}>[ clear ]</button>
		<a href="/compare" class="action-on">[ compare ]</a>
	</div>
{/if}
