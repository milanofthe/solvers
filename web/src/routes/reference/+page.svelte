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
					what: 'Nothing on or above the diagonal, so every stage is a formula in the ones before it and a step costs one right hand side per stage. There is no solve and no A-stability: the step size is bounded by the stability region, whatever the accuracy asks for.'
				},
				{
					name: 'strong stability preserving',
					fill: 'strict',
					count: () => family('ssprk'),
					href: '/?tag=strong+stability+preserving',
					what: 'Explicit methods whose step is a convex combination of forward Euler steps. Whatever bound holds for Euler, positivity or a maximum principle or total variation, then survives the whole step up to a step size limit, and how large that limit is is the number these are optimised for.'
				},
				{
					name: 'Chebyshev',
					fill: 'strict',
					count: () => family('rkc'),
					href: '/?family=rkc',
					what: 'Explicit, with many stages spent stretching the stability region along the negative real axis instead of raising the order. Twenty stages reach past -700, which is what makes an explicit method usable on a problem whose stiffness is real and not oscillatory.'
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
					what: 'The diagonal is filled and nothing above it, so the stages solve one after another, each of them a system the size of the problem rather than a multiple of it. Where the diagonal is one repeated value, a single factorization serves every stage of every step.'
				},
				{
					name: 'explicit first stage',
					fill: 'esdirk',
					count: () => family('esdirk'),
					href: '/?tag=explicit+first+stage',
					what: 'The same, with the first stage left explicit. That one free stage is what lets the rest share a diagonal and still end the step on the solution point, which is why these are the usual choice for a stiff nonlinear problem today.'
				},
				{
					name: 'fully implicit',
					fill: 'full',
					count: () => family('irk'),
					href: '/?tag=fully+implicit',
					what: 'Nothing is zero, so the stages are one coupled system of size s times n. It buys the most order a stage count can carry, twice the stages for Gauss, and the best stability there is, and it pays for both in the size of the solve.'
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
					what: 'One linear solve per stage against the Jacobian and no iteration anywhere. On a linear problem that is exact; on a nonlinear one the Jacobian enters the order conditions, so the tableau carries a second matrix and the method is only what it claims to be with a Jacobian that is current.'
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
					what: 'Two tableaux on one set of stages: the explicit one integrates the part of the right hand side that is cheap, the implicit one the part that is stiff. Its order is not the order of either half. Every way the two interleave has to be right as well, which is a longer list of conditions than either half alone.'
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
					what: 'Keeps the last k right hand sides and integrates the polynomial through them. Order k for one evaluation a step, which no one step method comes close to, and in exchange a stability region that shrinks as the order rises and a history that has to be started by something else.'
				},
				{
					name: 'Adams-Moulton',
					fill: 'past',
					count: () => family('adams_moulton'),
					href: '/?family=adams_moulton',
					what: 'The same construction with the new point included in the polynomial. One order more from the same number of past points, and a solve to pay for it.'
				},
				{
					name: 'backward differentiation',
					fill: 'past',
					count: () => family('bdf'),
					href: '/?family=bdf',
					what: 'Differentiates the polynomial through the last k solution values instead of integrating the one through the derivatives. A-stable at orders one and two, a stability wedge that closes as the order rises, and not zero stable at all past six.'
				},
				{
					name: 'Nystrom and Milne-Simpson',
					fill: 'past-explicit',
					count: () => family('nystrom', 'milne_simpson'),
					href: '/?family=nystrom',
					what: 'Step across two intervals rather than one, which puts a second root of the characteristic polynomial on the unit circle. That root never decays, so these are only weakly stable and cannot be run on a decaying problem at all.'
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
	<p class="label mt-2">what is being solved, and what the properties of a solver are</p>
</header>

{#snippet problem()}
	<div>
		{@render term(
			'initial value problem',
			'\\dot{y} = f(t, y), \\qquad y(t_0) = y_0, \\qquad y(t) \\in \\mathbb{R}^n',
			'A state and a rule for how fast it changes. Everything in this library is a rule for advancing that state by one step; solving the problem is applying it repeatedly.',
			null
		)}
		{@render term(
			'step and error',
			'y_{n+1} \\approx y(t_n + h_n), \\qquad \\|y_{n+1} - y(t_n + h_n)\\| = O(h^{p+1})',
			'The step size is chosen anew each step. An error of order p + 1 in one step accumulates to order p over a fixed interval, and p is what a method is called.',
			null
		)}
		{@render term(
			'stiffness',
			'\\dot{y} = J y, \\qquad \\frac{\\max_i |\\operatorname{Re}\\lambda_i(J)|}{\\min_i |\\operatorname{Re}\\lambda_i(J)|} \\gg 1',
			'Rates in the problem separated by orders of magnitude. The fast ones die out but still bound the step an explicit method may take, so the step is set by stability rather than by accuracy. That is what the implicit families exist for.',
			'stiff problems'
		)}
		{@render term(
			'differential algebraic',
			'M \\dot{y} = f(t, y), \\qquad M \\text{ singular}',
			'Part of the system is a constraint rather than a derivative. The index is the number of differentiations of f needed to recover the derivative; from index two the order a method attains falls below its classical one. Not integrated by this library, which solves the explicit form.',
			null
		)}
		{@render term(
			'what a solver is judged on',
			null,
			'Three things, and they trade against each other: the order it converges at, the region of the complex plane in which it is stable, and what one step costs. The rest of this page is the vocabulary for those three.',
			null
		)}
	</div>
{/snippet}
{@render section('the problem', problem)}

{#snippet families()}
	<p class="mb-6 max-w-[52rem] text-xs text-cream/60">
		The square beside each name is its coefficient matrix, and the filled cells are the entries the
		family allows to be nonzero. That pattern is the family: it decides what a step costs, whether
		there is a solve and how large it is, and how much stability is available at all. A multistep
		family has no such matrix, so the row shown is the past it reads instead. The pair has two
		matrices, and the lighter cells are what only its implicit half fills.
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
				No stages. A file here states which coefficients are free rather than their values; they are
				solved for at every step from the step sizes actually taken, so a family is defined on a
				nonuniform grid.
			</p>
		</div>
		<div>
			<p class="label mb-3">Rosenbrock-Wanner</p>
			<Math
				formula={'\\begin{aligned}(I - h\\gamma_{ii} J)\\,K_i &= h f\\!\\left(t_n + c_i h,\\, y_n + \\sum_j \\alpha_{ij} K_j\\right) \\\\ &\\quad + hJ\\sum_j \\gamma_{ij} K_j\\end{aligned}'}
				display
			/>
			<p class="mt-3 text-xs text-cream/60">
				<Math formula={'J = f\'(y_n)'} /> is part of the formula. A step is
				<Math formula={'s'} /> linear solves fixed in advance, with nothing to iterate.
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
				shared: each half evaluates at its own <Math formula={'c_i'} />, and several published
				pairs differ there.
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
			'One condition per rooted tree. Coefficients stated as fractions are checked as identities, decimals against a tolerance, and the report says which.',
			null
		)}
		{@render term(
			'simplifying assumptions',
			'\\begin{aligned} B(p) &: \\sum_i b_i c_i^{k-1} = \\tfrac{1}{k} \\\\ C(q) &: \\sum_j a_{ij} c_j^{k-1} = \\tfrac{c_i^k}{k} \\\\ D(r) &: \\sum_i b_i c_i^{k-1} a_{ij} = \\tfrac{b_j (1 - c_j^k)}{k} \\end{aligned}',
			'All three together give order at least min(p, q + r + 1, 2q + 2), by Butcher. Above order ten the rooted trees number in the hundreds of thousands and this is the route that remains.',
			null
		)}
		{@render term(
			'stage order',
			'q = \\max\\{k : C(k) \\text{ holds}\\}',
			'The order the internal stages themselves attain. On a stiff problem a method converges at its stage order until the step is small enough, not at its order.',
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
			'One step applied to a scalar mode, with z the step size times the rate. The stability region is where its modulus is at most one; the black curve on every figure here is that boundary.',
			null
		)}
		{@render term(
			'A-stable',
			'|R(z)| \\le 1 \\quad \\text{for all } \\operatorname{Re} z \\le 0',
			'No step size is too large for a decaying mode. Says nothing about how fast a stiff mode is damped.',
			'A-stable'
		)}
		{@render term(
			'L-stable',
			'\\text{A-stable and } \\lim_{z \\to -\\infty} R(z) = 0',
			'The infinitely stiff modes are damped to nothing in one step rather than merely not amplified.',
			'L-stable'
		)}
		{@render term(
			'A(alpha)-stable',
			'\\{ z : |\\arg(-z)| \\le \\alpha \\} \\subseteq \\{ z : |R(z)| \\le 1 \\}',
			'A wedge about the negative real axis instead of the half plane. For the backward differentiation formulae the angle falls from 86 degrees at order three to 18 at order six.',
			'A(alpha)-stable'
		)}
		{@render term(
			'stiffly accurate',
			'a_{sj} = b_j \\quad \\text{for all } j',
			'The step ends on the solution point, so the last stage is the answer. Gives R at infinity equal to zero, and satisfies the constraint of a differential algebraic problem exactly at the end of each step.',
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
			'Two solutions of a contractive problem do not move apart under the method, at any step size. Implies B-stability. Linear stability does not decide it.',
			'algebraically stable'
		)}
		{@render term(
			'strong stability preserving',
			'\\begin{aligned} \\mathcal{C} = \\max\\{\\, r :\\;& (I + rA)^{-1} \\text{ exists}, \\\\ & rA(I + rA)^{-1} \\ge 0, \\\\ & r\\, b^{\\top}(I + rA)^{-1} \\ge 0, \\\\ & (I + rA)^{-1}\\mathbb{1} \\ge 0 \\,\\} \\end{aligned}',
			'The radius of absolute monotonicity: the step is a convex combination of forward Euler steps of size h over C. Any bound forward Euler keeps is kept at a step C times as large.',
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
			'The first term whose imaginary part does not vanish gives the order of the phase error; the first at which the modulus departs from one gives the order of the amplitude error.',
			null
		)}
		{@render term(
			'non dissipative',
			'|R(iy)| = 1 \\quad \\text{for all real } y',
			'No amplitude error at any order. Exact for the Gauss methods.',
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
			'A second set of weights on the same stages, one order lower. Their difference estimates the error at no further evaluations. Without one a method is stepped twice at half the size instead, at about three times the cost.',
			'embedded pair'
		)}
		{@render term(
			'FSAL',
			'a_{sj} = b_j \\quad \\text{and} \\quad c_s = 1',
			'First same as last. The final stage of an explicit method is the first of the next step, so an accepted step costs one evaluation fewer.',
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
