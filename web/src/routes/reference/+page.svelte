<script>
	/*
	 * What the words mean.
	 *
	 * The library can be filtered by a dozen properties and says nowhere what any
	 * of them is. This page is the other half: the defining form of each family
	 * and the definition of each property, with nothing added. Every count and
	 * every figure on it comes from the same library and the same pipeline as the
	 * rest of the site, so it cannot drift away from what is actually there.
	 */
	import { getContext } from 'svelte';
	import Math from '$lib/components/Math.svelte';
	import Panel from '$lib/components/Panel.svelte';
	import { catalog, engine, ui } from '$lib/store.svelte.js';
	import { PRIORITY } from '$lib/engine.js';
	import { METHOD_MODES } from '$lib/figures.js';

	const rail = getContext('rail');

	$effect(() => {
		rail.panel = controls;
		return () => {
			rail.panel = null;
		};
	});

	const STABILITY = METHOD_MODES.find((mode) => mode.key === 'stability');

	/** A stability region, drawn from the library rather than pasted in. */
	function region(id) {
		return (box, view) => {
			const method = catalog.methods.find((m) => m.id === id);
			if (!method) return Promise.resolve(null);
			return (
				view && STABILITY.resample
					? STABILITY.resample(engine(), method, view, box, PRIORITY.visible)
					: STABILITY.request(engine(), method, PRIORITY.visible, { box })
			).then((result) => STABILITY.figure(result, method, { compact: true }));
		};
	}

	const counts = $derived.by(() => {
		const byFamily = new Map();
		const byTag = new Map();
		for (const method of catalog.methods) {
			byFamily.set(method.family, (byFamily.get(method.family) ?? 0) + 1);
			for (const tag of method.tags ?? []) {
				byTag.set(tag.name, (byTag.get(tag.name) ?? 0) + 1);
			}
		}
		return { byFamily, byTag };
	});

	const family = (...names) => names.reduce((sum, n) => sum + (counts.byFamily.get(n) ?? 0), 0);
	const tagged = (name) => counts.byTag.get(name) ?? 0;

	/**
	 * The families, each with the sparsity of its coefficient matrix.
	 *
	 * The glyph is not decoration. Explicit, diagonally implicit and fully
	 * implicit are not three degrees of anything: they are three answers to which
	 * entries of `A` may be nonzero, and that is the whole distinction.
	 */
	const GROUPS = [
		{
			heading: 'explicit',
			families: [
				{ name: 'Runge-Kutta', fill: 'strict', count: () => family('erk'), href: '/?tag=explicit' },
				{
					name: 'strong stability preserving',
					fill: 'strict',
					count: () => family('ssprk'),
					href: '/?tag=strong+stability+preserving'
				},
				{
					name: 'Chebyshev',
					fill: 'strict',
					count: () => family('rkc'),
					href: '/?family=rkc'
				}
			]
		},
		{
			heading: 'implicit',
			families: [
				{
					name: 'diagonally implicit',
					fill: 'lower',
					count: () => family('dirk'),
					href: '/?tag=singly+diagonal'
				},
				{
					name: 'explicit first stage',
					fill: 'esdirk',
					count: () => family('esdirk'),
					href: '/?tag=explicit+first+stage'
				},
				{
					name: 'fully implicit',
					fill: 'full',
					count: () => family('irk'),
					href: '/?tag=fully+implicit'
				}
			]
		},
		{
			heading: 'linearly implicit',
			families: [
				{
					name: 'Rosenbrock-Wanner',
					fill: 'lower',
					count: () => family('rosenbrock'),
					href: '/?tag=linearly+implicit'
				}
			]
		},
		{
			heading: 'multistep',
			families: [
				{
					name: 'Adams-Bashforth',
					fill: 'past-explicit',
					count: () => family('adams_bashforth'),
					href: '/?family=adams_bashforth'
				},
				{
					name: 'Adams-Moulton',
					fill: 'past',
					count: () => family('adams_moulton'),
					href: '/?family=adams_moulton'
				},
				{ name: 'backward differentiation', fill: 'past', count: () => family('bdf'), href: '/?family=bdf' },
				{
					name: 'Nystrom and Milne-Simpson',
					fill: 'past-explicit',
					count: () => family('nystrom', 'milne_simpson'),
					href: '/?family=nystrom'
				}
			]
		}
	];

	const SIDE = 5;

	/** Which cells of the glyph are filled, for one structure. */
	function cells(fill) {
		const out = [];
		for (let i = 0; i < SIDE; i += 1) {
			for (let j = 0; j < SIDE; j += 1) {
				const keep =
					fill === 'strict'
						? j < i
						: fill === 'lower'
							? j <= i
							: fill === 'esdirk'
								? j <= i && i > 0
								: fill === 'full'
									? true
									: fill === 'past'
										? i === SIDE - 1
										: /* past-explicit */ i === SIDE - 1 && j < SIDE - 1;
				if (keep) out.push([i, j]);
			}
		}
		return out;
	}
</script>

{#snippet controls()}
	<div class="rail-group">
		<a href="/" class="rail-link">back to the library</a>
	</div>
{/snippet}

{#snippet section(title, body)}
	<section class="-mx-6 mt-9 border-t border-cream/10 px-6 pt-4">
		<p class="label mb-4">{title}</p>
		{@render body()}
	</section>
{/snippet}

{#snippet term(name, formula, sentence, tag)}
	<div class="border-t border-cream/10 py-4 first:border-t-0 first:pt-0">
		<div class="flex flex-wrap items-baseline justify-between gap-x-6 gap-y-1">
			<h3 class="font-mono text-xs uppercase tracking-[0.2em] text-cream">{name}</h3>
			{#if tag}
				<a
					href="/?tag={encodeURIComponent(tag)}"
					class="font-mono text-xs text-accent hover:text-cream"
					>{tagged(tag)}
					{tag}</a
				>
			{/if}
		</div>
		{#if formula}
			<div class="mt-2 text-sm"><Math {formula} display /></div>
		{/if}
		<p class="mt-2 max-w-[78ch] text-xs text-cream/60">{sentence}</p>
	</div>
{/snippet}

<header>
	<div class="flex items-baseline justify-between gap-6">
		<h1 class="font-display text-lg text-cream">reference</h1>
		<span class="shrink-0 font-mono text-xs text-accent">{catalog.methods.length} methods</span>
	</div>
	<p class="label mt-2">the forms and the properties, as they are defined</p>
</header>

{#snippet families()}
	<div class="grid gap-x-10 gap-y-8 sm:grid-cols-2 xl:grid-cols-4">
		{#each GROUPS as group}
			<div>
				<p class="label mb-3">{group.heading}</p>
				<div class="flex flex-col gap-4">
					{#each group.families as entry}
						<a href={entry.href} class="group flex items-center gap-3">
							<svg viewBox="0 0 {SIDE} {SIDE}" class="h-9 w-9 shrink-0" aria-hidden="true">
								<rect width={SIDE} height={SIDE} fill="#141413" />
								{#each cells(entry.fill) as [i, j]}
									<rect
										x={j + 0.1}
										y={i + 0.1}
										width="0.8"
										height="0.8"
										class="fill-accent group-hover:fill-teal"
									/>
								{/each}
							</svg>
							<span class="min-w-0">
								<span class="block truncate text-xs text-cream group-hover:text-teal"
									>{entry.name}</span
								>
								<span class="block font-mono text-xs text-cream/40">{entry.count()}</span>
							</span>
						</a>
					{/each}
				</div>
			</div>
		{/each}
	</div>
	<p class="mt-6 max-w-[78ch] text-xs text-cream/60">
		The square is the coefficient matrix and the filled cells are the entries a family allows to be
		nonzero. Explicit, diagonally implicit and fully implicit are not three degrees of the same
		thing; they are three answers to that one question, and everything else about a family follows
		from it. A multistep family has no such matrix, so the row shown is what it reaches back into.
	</p>
{/snippet}
{@render section('families', families)}

{#snippet forms()}
	<div class="grid gap-x-12 gap-y-8 lg:grid-cols-3">
		<div>
			<p class="label mb-3">Runge-Kutta</p>
			<Math
				formula={'y_{n+1} = y_n + h\\sum_{i=1}^{s} b_i k_i'}
				display
			/>
			<Math
				formula={'k_i = f\\!\\left(t_n + c_i h,\\; y_n + h\\sum_{j=1}^{s} a_{ij} k_j\\right)'}
				display
			/>
			<p class="mt-2 text-xs text-cream/60">
				The whole method is the array of numbers, written as
				<Math formula={'\\begin{array}{c|c} c & A \\\\ \\hline & b^{\\top}\\end{array}'} />.
			</p>
		</div>
		<div>
			<p class="label mb-3">linear multistep</p>
			<Math
				formula={'\\sum_{j=0}^{k} \\alpha_j\\, y_{n-j} \\;=\\; h \\sum_{j=0}^{k} \\beta_j\\, f_{n-j}'}
				display
			/>
			<p class="mt-2 text-xs text-cream/60">
				No stages: the past is the memory. A file here states which coefficients are free rather
				than their values, and they are solved for at every step from the step sizes actually
				taken, so a family is defined on a nonuniform grid rather than only on an even one.
			</p>
		</div>
		<div>
			<p class="label mb-3">Rosenbrock-Wanner</p>
			<Math
				formula={'\\begin{aligned}(I - h\\gamma_{ii} J)\\,K_i &= h f\\!\\left(t_n + c_i h,\\, y_n + \\sum_j \\alpha_{ij} K_j\\right) \\\\ &\\quad + hJ\\sum_j \\gamma_{ij} K_j\\end{aligned}'}
				display
			/>
			<p class="mt-2 text-xs text-cream/60">
				<Math formula={'J = f\'(y_n)'} /> is part of the formula, not a device for solving one. There
				is nothing to iterate: a step is
				<Math formula={'s'} /> linear solves decided in advance.
			</p>
		</div>
	</div>
{/snippet}
{@render section('what a method is', forms)}

{#snippet order()}
	<div>
		{@render term(
			'order',
			'\\Phi(t) = \\sum_j b_j \\psi_j(t) = \\frac{1}{\\gamma(t)} \\quad \\text{for every rooted tree } |t| \\le p',
			'One condition per rooted tree. Where a file states its coefficients as fractions these are checked as identities rather than against a tolerance, and the report says which of the two happened.',
			null
		)}
		{@render term(
			'simplifying assumptions',
			'B(p):\\; \\sum_i b_i c_i^{k-1} = \\tfrac{1}{k} \\qquad C(q):\\; \\sum_j a_{ij} c_j^{k-1} = \\tfrac{c_i^k}{k} \\qquad D(r):\\; \\sum_i b_i c_i^{k-1} a_{ij} = \\tfrac{b_j (1 - c_j^k)}{k}',
			'Butcher: a method satisfying all three has order at least min(p, q + r + 1, 2q + 2), with no reference to trees. This is the only affordable route above order ten, where the trees number in the hundreds of thousands.',
			null
		)}
		{@render term(
			'stage order',
			'q = \\max\\{k : C(k) \\text{ holds}\\}',
			'How well the internal stages themselves approximate the solution. On a stiff problem a method converges at its stage order rather than at its order until the step is small enough, which is the order reduction that makes a high order method behave like a low order one.',
			null
		)}
	</div>
{/snippet}
{@render section('how order is established', order)}

{#snippet linear()}
	<div class="grid gap-10 lg:grid-cols-[1fr_20rem]">
		<div>
			{@render term(
				'stability function',
				'R(z) = 1 + z\\, b^{\\top} (I - zA)^{-1} \\mathbb{1} = \\frac{\\det(I - zA + z\\,\\mathbb{1}b^{\\top})}{\\det(I - zA)}',
				'What one step does to the scalar problem whose solution decays or turns at a rate lambda, with z the step size times that rate. The stability region is where the modulus is at most one, and the black curve on every figure here is its boundary.',
				null
			)}
			{@render term(
				'A-stable',
				'|R(z)| \\le 1 \\quad \\text{for all } \\operatorname{Re} z \\le 0',
				'No step size is too large for a decaying mode. Necessary for a stiff problem and not sufficient: it says nothing about how fast the stiff modes are damped.',
				'A-stable'
			)}
			{@render term(
				'L-stable',
				'\\text{A-stable and } \\lim_{z \\to -\\infty} R(z) = 0',
				'The infinitely stiff modes are killed outright rather than merely not amplified. Without it a mode at the far end of the spectrum survives every step at almost full size.',
				'L-stable'
			)}
			{@render term(
				'A(alpha)-stable',
				'\\{ z : |\\arg(-z)| \\le \\alpha \\} \\subseteq \\{ z : |R(z)| \\le 1 \\}',
				'A wedge about the negative real axis instead of the whole half plane. The backward differentiation formulae lose this angle as their order rises, from 86 degrees at order three to 18 at order six, which is why the family stops there.',
				'A(alpha)-stable'
			)}
			{@render term(
				'stiffly accurate',
				'a_{sj} = b_j \\quad \\text{for all } j',
				'The step ends on the solution point, so the last stage is the answer. It is what gives an implicit method R(inf) = 0 and what makes it usable on a differential algebraic problem.',
				'stiffly accurate'
			)}
		</div>
		<div class="flex flex-col gap-6">
			<Panel label="A-stable and L-stable: Radau IIA 5" height="h-40" load={region('radau_iia_5')} />
			<Panel label="A-stable, not L-stable: Gauss-Legendre 6" height="h-40" load={region('gauss_legendre_6')} />
			<Panel label="conditionally stable: classical Runge-Kutta" height="h-40" load={region('rk4')} />
			<Panel label="A(alpha)-stable: BDF 5" height="h-40" load={region('bdf5')} />
		</div>
	</div>
{/snippet}
{@render section('linear stability', linear)}

{#snippet nonlinear()}
	<div>
		{@render term(
			'algebraically stable',
			'b \\ge 0 \\quad \\text{and} \\quad M = \\operatorname{diag}(b) A + A^{\\top} \\operatorname{diag}(b) - b b^{\\top} \\succeq 0',
			'Two solutions of a contractive problem never move apart under the method, however large the step. Linear stability says nothing about this: it is a statement about the problem being nonlinear.',
			'algebraically stable'
		)}
		{@render term(
			'strong stability preserving',
			'\\mathcal{C} = \\max\\{ r : (I + rA)^{-1} \\text{ exists},\\; rA(I + rA)^{-1} \\ge 0,\\; r b^{\\top}(I + rA)^{-1} \\ge 0,\\; (I + rA)^{-1}\\mathbb{1} \\ge 0 \\}',
			'The radius of absolute monotonicity: how many forward Euler steps of size h/C the step is a convex combination of. A method with C > 0 keeps whatever bound forward Euler keeps, at a step C times as large.',
			'strong stability preserving'
		)}
	</div>
{/snippet}
{@render section('nonlinear stability', nonlinear)}

{#snippet oscillation()}
	<div>
		{@render term(
			'dispersion and dissipation order',
			'W(y) = R(iy)\\,e^{-iy} = 1 + \\sum_{m \\ge 1} w_m y^m',
			'On the imaginary axis the exact solution turns without growing. The first term of W whose imaginary part does not vanish gives the order to which the method gets the phase right; the first term of |W| that departs from one gives the order for the amplitude.',
			null
		)}
		{@render term(
			'non dissipative',
			'|R(iy)| = 1 \\quad \\text{for all real } y',
			'No amplitude error at any order. The Gauss methods have this exactly, which is what makes them worth their cost on a problem that has to run for many periods.',
			'non dissipative'
		)}
	</div>
{/snippet}
{@render section('oscillation', oscillation)}

{#snippet control()}
	<div>
		{@render term(
			'embedded pair',
			'\\hat{y}_{n+1} = y_n + h \\sum_i \\hat{b}_i k_i, \\qquad \\text{err} = \\left\\| \\sum_i (b_i - \\hat{b}_i) k_i \\right\\|',
			'A second set of weights on the same stages, one order lower, whose difference from the first estimates the error at no extra evaluations. A method without one is stepped twice at half the size instead, which costs about three times as much.',
			'embedded pair'
		)}
		{@render term(
			'FSAL',
			'a_{sj} = b_j \\quad \\text{and} \\quad c_s = 1',
			'First same as last: an explicit method whose final stage is the first stage of the next step, so a stage costs nothing on every step that is accepted.',
			'FSAL'
		)}
	</div>
{/snippet}
{@render section('error control', control)}

{#snippet algebraic()}
	<div>
		<p class="max-w-[78ch] text-xs text-cream/60">
			A differential algebraic problem <Math formula={'M \\dot{y} = f(t, y)'} /> with
			<Math formula={'M'} /> singular is not covered by this library, which integrates
			<Math formula={'\\dot{y} = f(t, y)'} /> only. Two of the properties above are nonetheless the
			ones that decide a method's behaviour there, and are worth reading with that in mind.
		</p>
		<div class="mt-4">
			{@render term(
				'index',
				'\\text{the number of differentiations of } f \\text{ needed to recover } \\dot{y}',
				'Index one is a constrained system whose constraint can be eliminated; index two and above cannot, and the order a method attains falls with the index rather than staying at its classical value.',
				null
			)}
			{@render term(
				'what carries over',
				null,
				'A method that is stiffly accurate satisfies the constraint exactly at the end of each step, which is why the Radau IIA, Lobatto IIIC and RODAS families are the ones used there. Stage order bounds the order actually attained on the algebraic part, which is why a method of high order and low stage order gains nothing.',
				'stiffly accurate'
			)}
		</div>
	</div>
{/snippet}
{@render section('differential algebraic problems', algebraic)}
