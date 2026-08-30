<script>
	/*
	 * The reference problems.
	 *
	 * A method is only fast or accurate relative to something, so the problems
	 * get the same treatment as the methods. Picking one makes it the problem the
	 * method cards and the comparison run on.
	 */
	import { flip } from 'svelte/animate';
	import { cubicOut } from 'svelte/easing';
	import { getContext } from 'svelte';
	import ProblemCard from '$lib/components/ProblemCard.svelte';
	import { catalog, engine, ui } from '$lib/store.svelte.js';
	import { PROBLEM_MODES } from '$lib/figures.js';

	const rail = getContext('rail');
	const stiff = $derived(catalog.problems.filter((p) => p.stiff));
	const nonStiff = $derived(catalog.problems.filter((p) => !p.stiff));

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
		{#each PROBLEM_MODES as entry}
			<button
				type="button"
				class="facet {ui.problemMode === entry.key ? 'facet-on' : ''}"
				onclick={() => {
					engine().drop(() => true);
					ui.problemMode = entry.key;
				}}
			>
				<span>{entry.label}</span>
			</button>
		{/each}
	</div>

	<div class="rail-group">
		<p class="label mb-1">selected</p>
		<p class="font-mono text-xs text-cream">{ui.problem}</p>
		<a href="/compare" class="rail-link mt-2">compare on it</a>
	</div>
{/snippet}

<section class="mb-8">
	<p class="label mb-2">non stiff</p>
	<div class="grid grid-cols-[repeat(auto-fill,minmax(17rem,1fr))] gap-4">
		{#each nonStiff as problem (problem.id)}
			<div animate:flip={{ duration: 320, easing: cubicOut }}>
				<ProblemCard {problem} />
			</div>
		{/each}
	</div>
</section>

<section>
	<p class="label mb-2">stiff</p>
	<div class="grid grid-cols-[repeat(auto-fill,minmax(17rem,1fr))] gap-4">
		{#each stiff as problem (problem.id)}
			<div animate:flip={{ duration: 320, easing: cubicOut }}>
				<ProblemCard {problem} />
			</div>
		{/each}
	</div>
</section>
