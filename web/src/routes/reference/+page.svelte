<script>
	/*
	 * What the words mean.
	 *
	 * The library can be filtered by a dozen properties and says nowhere what any
	 * of them is. This page is the other half: what is being solved, the defining
	 * form of each family, and the definition of each property. Every count and
	 * every figure comes from the same library and the same pipeline as the rest
	 * of the site, so it cannot drift away from what is actually there.
	 *
	 * A term is one row of three columns: the name, the definition, and what the
	 * definition says. Three columns rather than one because a formula centred in
	 * a wide page with a paragraph capped beside it wastes the width twice over.
	 * Every rule runs the width of the content, so the page reads as a stack of
	 * rows rather than as a column floating in a frame.
	 */
	import { getContext } from 'svelte';
	import Math from '$lib/components/Math.svelte';
	import Panel from '$lib/components/Panel.svelte';
	import { catalog, engine } from '$lib/store.svelte.js';
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
	 * The families, each with the sparsity of its coefficient matrix. The glyph
	 * is the definition: which entries of `A` a family allows to be nonzero.
	 */
	const GROUPS = [
		{
			heading: 'explicit',
			families: [
				{
					name: 'Runge-Kutta',
					fill: 'strict',
					count: () => family('erk'),
					href: '/?tag=explicit',
					what: 'Nothing on or above the diagonal, so each stage is just a formula in the ones before it. One right hand side per stage, nothing to solve. The catch is stability: your step size is capped by the stability region no matter how much accuracy you actually need.'
				},
				{
					name: 'strong stability preserving',
					fill: 'strict',
					count: () => family('ssprk'),
					href: '/?tag=strong+stability+preserving',
					what: 'Explicit methods whose step is a convex combination of forward Euler steps. Whatever holds for Euler (positivity, a maximum principle, total variation) still holds after a full step, up to some step size limit. How large that limit is is the number these are optimised for.'
				},
				{
					name: 'Chebyshev',
					fill: 'strict',
					count: () => family('rkc'),
					href: '/?family=rkc',
					what: 'Explicit, but with a lot of stages spent stretching the stability region down the negative real axis instead of raising the order. Twenty stages get you past -700, which is what makes an explicit method workable on a problem whose stiffness is real rather than oscillatory.'
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
					href: '/?tag=singly+diagonal',
					what: 'Diagonal filled, nothing above it, so the stages solve one after the other. Each one is a system the size of the problem, not a multiple of it. If the diagonal is a single repeated value you factorize once and reuse it for every stage of every step.'
				},
				{
					name: 'explicit first stage',
					fill: 'esdirk',
					count: () => family('esdirk'),
					href: '/?tag=explicit+first+stage',
					what: 'The same, but the first stage stays explicit. That one free stage is what lets the rest share a diagonal and still land on the solution point, and it is why these are the usual pick for a stiff nonlinear problem today.'
				},
				{
					name: 'fully implicit',
					fill: 'full',
					count: () => family('irk'),
					href: '/?tag=fully+implicit',
					what: 'Nothing is zero, so all the stages solve together as one system of size s times n. You get the most order a stage count can carry (Gauss gets you twice the stages) and the best stability available, and you pay for both in the size of that solve.'
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
					href: '/?tag=linearly+implicit',
					what: 'One linear solve per stage against the Jacobian, and nothing to iterate. On a linear problem it is exact. On a nonlinear one the Jacobian shows up in the order conditions, so the tableau carries a second matrix and the method only does what it says with a Jacobian that is current.'
				}
			]
		},
		{
			heading: 'additive',
			families: [
				{
					name: 'implicit-explicit',
					fill: 'additive',
					count: () => family('imex'),
					href: '/?family=imex',
					what: 'Two tableaux sharing one set of stages: the explicit one takes the cheap part of the right hand side, the implicit one takes the stiff part. Its order is not the order of either half. Every way the two interleave has to come out right too, which is a much longer list of conditions.'
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
					href: '/?family=adams_bashforth',
					what: 'Keeps the last k right hand sides and integrates the polynomial through them. Order k for a single evaluation per step, which no one step method gets anywhere near. In exchange the stability region shrinks as the order goes up, and you need something else to start the history.'
				},
				{
					name: 'Adams-Moulton',
					fill: 'past',
					count: () => family('adams_moulton'),
					href: '/?family=adams_moulton',
					what: 'Same construction, but the new point is in the polynomial too. One order more out of the same number of past points, and a solve to pay for it.'
				},
				{
					name: 'backward differentiation',
					fill: 'past',
					count: () => family('bdf'),
					href: '/?family=bdf',
					what: 'Differentiates the polynomial through the last k solution values instead of integrating the one through the derivatives. A-stable at orders one and two, the wedge closes as the order goes up, and past six it is not zero stable at all.'
				},
				{
					name: 'Nystrom and Milne-Simpson',
					fill: 'past-explicit',
					count: () => family('nystrom', 'milne_simpson'),
					href: '/?family=nystrom',
					what: 'These step across two intervals instead of one, which puts a second root of the characteristic polynomial on the unit circle. That root never decays, so they are only weakly stable and cannot be run on a decaying problem at all.'
				}
			]
		}
	];

	const SIDE = 5;

	/**
	 * Which cells of the glyph are filled, for one structure, and which of them
	 * belong to the second matrix of a family that has two.
	 */
	function cells(fill) {
		const out = [];
		for (let i = 0; i < SIDE; i += 1) {
			for (let j = 0; j < SIDE; j += 1) {
				if (fill === 'additive') {
					// The explicit half below the diagonal, the implicit half on
					// it: the pair is what the two together allow.
					if (j < i) out.push([i, j, false]);
					else if (j === i && i > 0) out.push([i, j, true]);
					continue;
				}
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
				if (keep) out.push([i, j, false]);
			}
		}
		return out;
	}

	const SOURCES = [
		{
			authors: 'Hairer, E., Noersett, S. P., & Wanner, G.',
			title: 'Solving Ordinary Differential Equations I: Nonstiff Problems',
			where: 'Springer, 2nd edition, 1993',
			doi: '10.1007/978-3-540-78862-1'
		},
		{
			authors: 'Hairer, E., & Wanner, G.',
			title: 'Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems',
			where: 'Springer, 2nd edition, 1996',
			doi: '10.1007/978-3-642-05221-7'
		},
		{
			authors: 'Butcher, J. C.',
			title: 'Numerical Methods for Ordinary Differential Equations',
			where: 'Wiley, 3rd edition, 2016',
			doi: '10.1002/9781119121534'
		},
		{
			authors: 'Dahlquist, G.',
			title: 'A special stability problem for linear multistep methods',
			where: 'BIT Numerical Mathematics 3, 1963',
			doi: '10.1007/BF01963532'
		},
		{
			authors: 'Burrage, K., & Butcher, J. C.',
			title: 'Stability criteria for implicit Runge-Kutta methods',
			where: 'SIAM Journal on Numerical Analysis 16(1), 1979',
			doi: '10.1137/0716037'
		},
		{
			authors: 'Gottlieb, S., Ketcheson, D. I., & Shu, C.-W.',
			title: 'Strong Stability Preserving Runge-Kutta and Multistep Time Discretizations',
			where: 'World Scientific, 2011',
			doi: '10.1142/7498'
		},
		{
			authors: 'van der Houwen, P. J., & Sommeijer, B. P.',
			title:
				'Explicit Runge-Kutta (-Nystroem) methods with reduced phase errors for computing oscillating solutions',
			where: 'SIAM Journal on Numerical Analysis 24(3), 1987',
			doi: '10.1137/0724041'
		},
		{
			authors: 'Soederlind, G.',
			title: 'Digital filters in adaptive time-stepping',
			where: 'ACM Transactions on Mathematical Software 29(1), 2003',
			doi: '10.1145/641876.641877'
		},
		{
			authors: 'Ehle, B. L.',
			title:
				'On Pade approximations to the exponential function and A-stable methods for the numerical solution of initial value problems',
			where: 'Research Report CSRR 2010, University of Waterloo, 1969'
		}
	];
</script>

{#snippet controls()}
	<div class="rail-group">
		<a href="/" class="rail-link">back to the library</a>
	</div>
{/snippet}

<!-- A section rule runs the width of the content, out to both panels. -->
{#snippet section(title, body)}
	<section class="-mx-6 mt-10 border-t border-cream/10 px-6 pt-5">
		<p class="label mb-5 text-cream">{title}</p>
		{@render body()}
	</section>
{/snippet}

{#snippet term(name, formula, sentence, tag)}
	<div
		class="-mx-6 grid gap-x-8 gap-y-2 border-t border-cream/10 px-6 py-4 first:border-t-0 first:pt-0 lg:grid-cols-[12rem_minmax(0,1fr)_minmax(0,32rem)] lg:items-baseline"
	>
		<div>
			<h3 class="font-mono text-xs uppercase tracking-[0.2em] text-cream">{name}</h3>
			{#if tag}
				<a
					href="/?tag={encodeURIComponent(tag)}"
					class="mt-1 block font-mono text-xs text-accent hover:text-cream">{tagged(tag)} {tag}</a
				>
			{/if}
		</div>
		{#if formula}
			<div class="min-w-0 overflow-x-auto"><Math {formula} display /></div>
			<p class="text-xs text-cream/60">{sentence}</p>
		{:else}
			<!-- Nothing to set, so the sentence takes the width the formula would
			     have had rather than leaving a column empty. -->
			<p class="text-xs text-cream/60 lg:col-span-2">{sentence}</p>
		{/if}
	</div>
{/snippet}

<header>
	<div class="flex items-baseline justify-between gap-6">
		<h1 class="font-display text-lg text-cream">reference</h1>
		<span class="shrink-0 font-mono text-xs text-accent">{catalog.methods.length} methods</span>
	</div>
	<p class="label mt-2">what you are solving, and what the properties of a solver actually mean</p>
</header>

{#snippet problem()}
	<div>
		{@render term(
			'initial value problem',
			'\\dot{y} = f(t, y), \\qquad y(t_0) = y_0, \\qquad y(t) \\in \\mathbb{R}^n',
			'A state and a rule for how fast it changes. Everything in this library is a way of moving that state forward by one step, and solving the problem is doing that over and over.',
			null
		)}
		{@render term(
			'step and error',
			'y_{n+1} \\approx y(t_n + h_n), \\qquad \\|y_{n+1} - y(t_n + h_n)\\| = O(h^{p+1})',
			'The step size gets picked again every step. An error of order p + 1 in a single step piles up to order p over a fixed interval, and p is the number the method is named after.',
			null
		)}
		{@render term(
			'stiffness',
			'\\dot{y} = J y, \\qquad \\frac{\\max_i |\\operatorname{Re}\\lambda_i(J)|}{\\min_i |\\operatorname{Re}\\lambda_i(J)|} \\gg 1',
			'Rates in the problem sitting orders of magnitude apart. The fast ones die out almost immediately but still cap the step an explicit method can take, so your step size ends up set by stability and not by the accuracy you wanted. That is the whole reason the implicit families exist.',
			'stiff problems'
		)}
		{@render term(
			'differential algebraic',
			'M \\dot{y} = f(t, y), \\qquad M \\text{ singular}',
			'Part of the system is a constraint instead of a derivative. The index counts how often you have to differentiate f to get that derivative back, and from index two on, methods drop below their classical order. This library does not integrate these, it solves the explicit form.',
			null
		)}
		{@render term(
			'what a solver is judged on',
			null,
			'Three things, and they trade against each other: the order it converges at, where in the complex plane it stays stable, and what a single step costs you. The rest of this page is the vocabulary for those three.',
			null
		)}
	</div>
{/snippet}
{@render section('the problem', problem)}

{#snippet families()}
	<p class="mb-6 max-w-[52rem] text-xs text-cream/60">
		The square next to each name is the coefficient matrix, and the filled cells are the entries
		that family lets be nonzero. That pattern is the family: it decides what a step costs, whether
		you have to solve anything and how big that solve is, and how much stability is on the table at
		all. Multistep families have no such matrix, so what you see is the past they read instead. The
		pair has two matrices, and the lighter cells are the ones only its implicit half fills.
	</p>
	{#each GROUPS as group}
		<div class="-mx-6 border-t border-cream/10 px-6 pt-5 first:border-t-0 first:pt-0">
			<p class="label mb-4">{group.heading}</p>
			{#each group.families as entry}
				<div
					class="mb-5 grid gap-x-8 gap-y-2 lg:grid-cols-[14rem_minmax(0,40rem)] lg:items-baseline"
				>
					<a href={entry.href} class="group flex items-center gap-3">
						<svg viewBox="0 0 {SIDE} {SIDE}" class="h-9 w-9 shrink-0" aria-hidden="true">
							<rect width={SIDE} height={SIDE} fill="#141413" />
							{#each cells(entry.fill) as [i, j, second]}
								<rect
									x={j + 0.1}
									y={i + 0.1}
									width="0.8"
									height="0.8"
									class="{second
										? 'fill-accent/40'
										: 'fill-accent'} group-hover:fill-teal"
								/>
							{/each}
						</svg>
						<span class="min-w-0">
							<span class="block text-xs text-cream group-hover:text-teal">{entry.name}</span>
							<span class="block font-mono text-xs text-cream/40">{entry.count()}</span>
						</span>
					</a>
					<p class="text-xs text-cream/60">{entry.what}</p>
				</div>
			{/each}
		</div>
	{/each}
{/snippet}
{@render section('families', families)}

{#snippet forms()}
	<!-- Two to a row rather than four: the Rosenbrock and additive forms are
	     wide enough that a quarter of the page cuts them off. -->
	<div class="grid gap-x-12 gap-y-10 lg:grid-cols-2">
		<div>
			<p class="label mb-3">Runge-Kutta</p>
			<Math formula={'y_{n+1} = y_n + h\\sum_{i=1}^{s} b_i k_i'} display />
			<Math
				formula={'k_i = f\\!\\left(t_n + c_i h,\\; y_n + h\\sum_{j=1}^{s} a_{ij} k_j\\right)'}
				display
			/>
			<p class="mt-3 text-xs text-cream/60">
				Written as the array
				<Math formula={'\\begin{array}{c|c} c & A \\\\ \\hline & b^{\\top}\\end{array}'} />, which is
				the whole method.
			</p>
		</div>
		<div>
			<p class="label mb-3">linear multistep</p>
			<Math
				formula={'\\sum_{j=0}^{k} \\alpha_j\\, y_{n-j} \\;=\\; h \\sum_{j=0}^{k} \\beta_j\\, f_{n-j}'}
				display
			/>
			<p class="mt-3 text-xs text-cream/60">
				No stages at all. A file here says which coefficients are free instead of giving their
				values, and they get solved for at every step from the step sizes actually taken, so a
				family is defined on a nonuniform grid.
			</p>
		</div>
		<div>
			<p class="label mb-3">Rosenbrock-Wanner</p>
			<Math
				formula={'\\begin{aligned}(I - h\\gamma_{ii} J)\\,K_i &= h f\\!\\left(t_n + c_i h,\\, y_n + \\sum_j \\alpha_{ij} K_j\\right) \\\\ &\\quad + hJ\\sum_j \\gamma_{ij} K_j\\end{aligned}'}
				display
			/>
			<p class="mt-3 text-xs text-cream/60">
				<Math formula={'J = f\'(y_n)'} /> is part of the formula itself. A step is
				<Math formula={'s'} /> linear solves decided in advance, with nothing to iterate.
			</p>
		</div>
		<div>
			<p class="label mb-3">additive Runge-Kutta</p>
			<Math
				formula={'y_{n+1} = y_n + h\\sum_{i=1}^{s}\\left( b^{E}_i f_E(Y_i) + b^{I}_i f_I(Y_i) \\right)'}
				display
			/>
			<Math
				formula={'Y_i = y_n + h\\sum_{j=1}^{s}\\left( a^{E}_{ij} f_E(Y_j) + a^{I}_{ij} f_I(Y_j) \\right)'}
				display
			/>
			<p class="mt-3 text-xs text-cream/60">
				Two arrays for one right hand side split in two, on shared stages. The abscissae are not
				shared though: each half evaluates at its own <Math formula={'c_i'} />, and several
				published pairs really do differ there.
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
			'One condition per rooted tree. Fractions get checked as identities, decimals against a tolerance, and the report tells you which of the two happened.',
			null
		)}
		{@render term(
			'simplifying assumptions',
			'\\begin{aligned} B(p) &: \\sum_i b_i c_i^{k-1} = \\tfrac{1}{k} \\\\ C(q) &: \\sum_j a_{ij} c_j^{k-1} = \\tfrac{c_i^k}{k} \\\\ D(r) &: \\sum_i b_i c_i^{k-1} a_{ij} = \\tfrac{b_j (1 - c_j^k)}{k} \\end{aligned}',
			'Butcher\'s theorem: the three together give you order at least min(p, q + r + 1, 2q + 2). Past order ten there are hundreds of thousands of rooted trees, so this is the only route left.',
			null
		)}
		{@render term(
			'stage order',
			'q = \\max\\{k : C(k) \\text{ holds}\\}',
			'The order the internal stages hit on their own. On a stiff problem you see the stage order, not the real one, until the step gets small enough.',
			null
		)}
	</div>
{/snippet}
{@render section('how order is established', order)}

{#snippet linear()}
	<div class="mb-8 grid gap-x-8 gap-y-6 sm:grid-cols-2 xl:grid-cols-4">
		<Panel label="Radau IIA 5 · A-stable, L-stable" height="h-44" load={region('radau_iia_5')} />
		<Panel
			label="Gauss-Legendre 6 · A-stable only"
			height="h-44"
			load={region('gauss_legendre_6')}
		/>
		<Panel label="BDF 5 · A(alpha)-stable" height="h-44" load={region('bdf5')} />
		<Panel label="RK4 · conditionally stable" height="h-44" load={region('rk4')} />
	</div>
	<div>
		{@render term(
			'stability function',
			'R(z) = 1 + z\\, b^{\\top} (I - zA)^{-1} \\mathbb{1} = \\frac{\\det(I - zA + z\\,\\mathbb{1}b^{\\top})}{\\det(I - zA)}',
			'One step applied to a single scalar mode, with z the step size times the rate. The stability region is where the modulus stays at or below one, and the black curve on every figure here is that boundary.',
			null
		)}
		{@render term(
			'A-stable',
			'|R(z)| \\le 1 \\quad \\text{for all } \\operatorname{Re} z \\le 0',
			'No step size is too big for a decaying mode. It says nothing about how fast a stiff mode actually gets damped.',
			'A-stable'
		)}
		{@render term(
			'L-stable',
			'\\text{A-stable and } \\lim_{z \\to -\\infty} R(z) = 0',
			'The infinitely stiff modes get damped to nothing in a single step, not just kept from growing.',
			'L-stable'
		)}
		{@render term(
			'A(alpha)-stable',
			'\\{ z : |\\arg(-z)| \\le \\alpha \\} \\subseteq \\{ z : |R(z)| \\le 1 \\}',
			'A wedge around the negative real axis instead of the whole half plane. For BDF the angle drops from 86 degrees at order three to 18 at order six.',
			'A(alpha)-stable'
		)}
		{@render term(
			'stiffly accurate',
			'a_{sj} = b_j \\quad \\text{for all } j',
			'The step ends on the solution point, so the last stage is already the answer. You get R at infinity equal to zero out of it, and the constraint of a differential algebraic problem holds exactly at the end of every step.',
			'stiffly accurate'
		)}
	</div>
{/snippet}
{@render section('linear stability', linear)}

{#snippet nonlinear()}
	<div>
		{@render term(
			'algebraically stable',
			'b \\ge 0 \\quad \\text{and} \\quad M = \\operatorname{diag}(b) A + A^{\\top} \\operatorname{diag}(b) - b b^{\\top} \\succeq 0',
			'Two solutions of a contractive problem never drift apart under the method, at any step size at all. It implies B-stability, and linear stability tells you nothing about whether you have it.',
			'algebraically stable'
		)}
		{@render term(
			'strong stability preserving',
			'\\begin{aligned} \\mathcal{C} = \\max\\{\\, r :\\;& (I + rA)^{-1} \\text{ exists}, \\\\ & rA(I + rA)^{-1} \\ge 0, \\\\ & r\\, b^{\\top}(I + rA)^{-1} \\ge 0, \\\\ & (I + rA)^{-1}\\mathbb{1} \\ge 0 \\,\\} \\end{aligned}',
			'The radius of absolute monotonicity. The step is a convex combination of forward Euler steps of size h over C, so anything Euler respects, you keep at a step C times as large.',
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
			'The first term whose imaginary part does not vanish gives you the order of the phase error, and the first one where the modulus leaves one gives the order of the amplitude error.',
			null
		)}
		{@render term(
			'non dissipative',
			'|R(iy)| = 1 \\quad \\text{for all real } y',
			'No amplitude error at all, at any order. The Gauss methods have this exactly.',
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
			'A second set of weights on the same stages, one order lower. The difference between them estimates the error and costs you no extra evaluations. Without one you step twice at half the size instead, for roughly three times the work.',
			'embedded pair'
		)}
		{@render term(
			'FSAL',
			'a_{sj} = b_j \\quad \\text{and} \\quad c_s = 1',
			'First same as last. The final stage of an explicit method is the first stage of the next one, so an accepted step costs you one evaluation less.',
			'FSAL'
		)}
	</div>
{/snippet}
{@render section('error control', control)}

{#snippet sources()}
	<ul class="grid gap-x-12 gap-y-5 lg:grid-cols-2">
		{#each SOURCES as source}
			<li>
				{#if source.doi}
					<a
						href="https://doi.org/{source.doi}"
						target="_blank"
						rel="noreferrer"
						class="text-xs text-cream underline decoration-cream/20 underline-offset-4 hover:decoration-cream"
						>{source.title}</a
					>
				{:else}
					<span class="text-xs text-cream">{source.title}</span>
				{/if}
				<p class="mt-1 font-mono text-xs text-accent">{source.authors} · {source.where}</p>
			</li>
		{/each}
	</ul>
{/snippet}
{@render section('sources', sources)}
