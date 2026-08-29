/*
 * Convergence and cost.
 *
 * Both views run the solvers in the page. Convergence is fixed step against a
 * known or tightly integrated solution, and reports the slope that comes out
 * next to the order the coefficients promise. Cost is the adaptive counterpart:
 * one run per tolerance, plotting the accuracy actually achieved against the
 * work it took, which is the question that decides which method to use.
 */

import { el, num, replace } from "../dom.js";
import { extent, Figure } from "../plot.js";
import { seriesColor } from "../colormap.js";

const DEFAULT_METHODS = ["rkdp54", "rkbs32", "rk4", "esdirk43", "bdf4"];

function methodPicker(methods, selected, onChange) {
  const list = el(
    "select",
    {
      multiple: true,
      size: "10",
      onchange: (event) => {
        const chosen = [...event.target.selectedOptions].map((option) => option.value);
        onChange(chosen);
      },
    },
    methods.map((method) =>
      el(
        "option",
        { value: method.id, selected: selected.includes(method.id) },
        `${method.name}  (order ${method.order})`,
      ),
    ),
  );
  return list;
}

function problemPicker(problems, selected, onChange) {
  return el(
    "select",
    { onchange: (event) => onChange(event.target.value) },
    problems.map((problem) =>
      el(
        "option",
        { value: problem.id, selected: problem.id === selected },
        `${problem.name}${problem.stiff ? " (stiff)" : ""}`,
      ),
    ),
  );
}

export function convergenceView(context) {
  const { wasm, methods, problems, state } = context;
  state.convergenceMethods ??= DEFAULT_METHODS.filter((id) => methods.some((m) => m.id === id));
  state.convergenceProblem ??= "nonlinear_decay";

  const canvas = el("canvas", { class: "figure__canvas figure__canvas--tall" });
  const table = el("tbody");
  const status = el("p", { class: "card__meta" });

  function run() {
    const chosen = state.convergenceMethods;
    if (chosen.length === 0) {
      status.textContent = "Select at least one method.";
      return;
    }
    status.textContent = "Running…";

    const studies = [];
    for (const id of chosen) {
      const method = methods.find((m) => m.id === id);
      if (!method) continue;
      // A high order method needs a coarser ladder to stay above round off.
      const coarse = method.order >= 7 ? 0.5 : method.order >= 5 ? 0.4 : 0.3;
      const ratio = method.order >= 5 ? 0.7 : 0.6;
      try {
        const study = JSON.parse(
          wasm.convergence_study(id, state.convergenceProblem, coarse, ratio, 7),
        );
        studies.push({ method, study });
      } catch (error) {
        studies.push({ method, error: String(error) });
      }
    }

    const figure = new Figure(canvas, {
      title: `convergence on ${state.convergenceProblem}`,
      xLabel: "step size h",
      yLabel: "relative error at the end point",
      xScale: "log",
      yScale: "log",
    });
    const allH = studies.flatMap((s) => s.study?.points.map((p) => p.h) ?? []);
    const allE = studies.flatMap(
      (s) => s.study?.points.map((p) => p.error).filter((e) => e > 1e-16 && Number.isFinite(e)) ?? [],
    );
    const [hLow, hHigh] = extent(allH, { log: true, pad: 0.08 });
    const [eLow, eHigh] = extent(allE, { log: true, pad: 0.12 });
    figure.setDomain(hLow, hHigh, eLow, eHigh);
    figure.clear();
    figure.drawAxes({ origin: false });
    figure.drawSeries(
      studies
        .filter((s) => s.study)
        .map((entry, index) => ({
          x: entry.study.points.map((p) => p.h),
          y: entry.study.points.map((p) => p.error),
          label: entry.method.id,
          color: seriesColor(index),
        })),
    );

    replace(
      table,
      studies.map((entry) =>
        el(
          "tr",
          {},
          el("td", { class: "table__name" }, entry.method.name),
          el("td", { class: "table__numeric" }, String(entry.method.order)),
          el("td", { class: "table__numeric" }, entry.study ? num(entry.study.estimated_order, 3) : "–"),
          el("td", { class: "table__numeric" }, entry.study ? num(entry.study.local_order, 3) : "–"),
          el(
            "td",
            {},
            entry.error
              ? el("span", { class: "badge badge--error" }, "failed")
              : Math.abs((entry.study?.local_order ?? 0) - entry.method.order) < 0.5
                ? el("span", { class: "badge badge--ok" }, "matches")
                : el("span", { class: "badge badge--warn" }, "differs"),
          ),
        ),
      ),
    );

    const reference = studies.find((s) => s.study)?.study;
    status.textContent = reference?.exact_reference
      ? "Errors are measured against the closed form solution."
      : "Errors are measured against a reference run at tolerance 1e-13.";
  }

  const picker = methodPicker(methods, state.convergenceMethods, (chosen) => {
    state.convergenceMethods = chosen;
    run();
  });

  const container = el(
    "section",
    { class: "section" },
    el("h2", { class: "section__title" }, "Convergence order"),
    el(
      "p",
      { class: "section__note" },
      "Fixed step runs over a ladder of step sizes. The slope of the fit is the order the method " +
        "actually converges at, which is not always the order its coefficients promise: a stiff " +
        "problem can reduce it, and a linear one can flatter it.",
    ),
    el(
      "div",
      { class: "toolbar" },
      el("label", { class: "field" }, el("span", { class: "field__label" }, "problem"),
        problemPicker(problems, state.convergenceProblem, (id) => {
          state.convergenceProblem = id;
          run();
        })),
      el("label", { class: "field" }, el("span", { class: "field__label" }, "methods"), picker),
    ),
    el("div", { class: "figure" }, canvas),
    el(
      "div",
      { class: "card", style: "margin-top:1rem" },
      el(
        "div",
        { class: "scroll-x" },
        el(
          "table",
          { class: "table" },
          el(
            "thead",
            {},
            el(
              "tr",
              {},
              el("th", {}, "method"),
              el("th", {}, "stated order"),
              el("th", {}, "fitted"),
              el("th", {}, "asymptotic"),
              el("th", {}, "agreement"),
            ),
          ),
          table,
        ),
      ),
      el("div", { class: "card__body" }, status),
    ),
  );

  requestAnimationFrame(run);
  return container;
}

export function costView(context) {
  const { wasm, methods, problems, state } = context;
  state.costMethods ??= DEFAULT_METHODS.filter((id) => methods.some((m) => m.id === id));
  state.costProblem ??= "lotka_volterra";
  state.costMetric ??= "rhs_evals";

  const canvas = el("canvas", { class: "figure__canvas figure__canvas--tall" });
  const status = el("p", { class: "card__meta" });

  const METRICS = [
    { key: "rhs_evals", label: "right hand side evaluations" },
    { key: "steps", label: "attempted steps" },
    { key: "lu_decompositions", label: "factorizations" },
    { key: "nonlinear_iterations", label: "nonlinear iterations" },
  ];

  function run() {
    const chosen = state.costMethods.filter((id) => methods.some((m) => m.id === id && m.adaptive));
    if (chosen.length === 0) {
      status.textContent = "Select at least one adaptive method.";
      return;
    }
    status.textContent = "Running…";

    const results = [];
    for (const id of chosen) {
      try {
        results.push({
          id,
          data: JSON.parse(wasm.work_precision(id, state.costProblem, -3, -11)),
        });
      } catch (error) {
        results.push({ id, error: String(error) });
      }
    }

    const series = results
      .filter((entry) => entry.data)
      .map((entry, index) => {
        const points = entry.data.points.filter(
          (p) => p.succeeded && Number.isFinite(p.error) && p.error > 0,
        );
        return {
          x: points.map((p) => p.error),
          y: points.map((p) => p[state.costMetric]),
          label: entry.id,
          color: seriesColor(index),
        };
      })
      .filter((s) => s.x.length > 1);

    const figure = new Figure(canvas, {
      title: `work precision on ${state.costProblem}`,
      xLabel: "achieved relative error",
      yLabel: METRICS.find((m) => m.key === state.costMetric).label,
      xScale: "log",
      yScale: "log",
    });
    const [xLow, xHigh] = extent(series.flatMap((s) => s.x), { log: true, pad: 0.1 });
    const [yLow, yHigh] = extent(series.flatMap((s) => s.y), { log: true, pad: 0.15 });
    figure.setDomain(xHigh, xLow, yLow, yHigh);
    figure.clear();
    figure.drawAxes({ origin: false });
    figure.drawSeries(series);

    status.textContent =
      "Accuracy improves to the right of the plot, so a method sitting lower is doing the same " +
      "job for less work. Each point is one adaptive run at one tolerance.";
  }

  const picker = methodPicker(
    methods.filter((m) => m.adaptive),
    state.costMethods,
    (chosen) => {
      state.costMethods = chosen;
      run();
    },
  );

  const container = el(
    "section",
    { class: "section" },
    el("h2", { class: "section__title" }, "Work precision"),
    el(
      "p",
      { class: "section__note" },
      "One adaptive run per tolerance, from 1e-3 down to 1e-11. The axis is reversed so that " +
        "accuracy increases to the right and the better method is the lower curve.",
    ),
    el(
      "div",
      { class: "toolbar" },
      el("label", { class: "field" }, el("span", { class: "field__label" }, "problem"),
        problemPicker(problems, state.costProblem, (id) => {
          state.costProblem = id;
          run();
        })),
      el(
        "label",
        { class: "field" },
        el("span", { class: "field__label" }, "cost measure"),
        el(
          "select",
          {
            onchange: (event) => {
              state.costMetric = event.target.value;
              run();
            },
          },
          METRICS.map((metric) =>
            el("option", { value: metric.key, selected: metric.key === state.costMetric }, metric.label),
          ),
        ),
      ),
      el("label", { class: "field" }, el("span", { class: "field__label" }, "methods"), picker),
    ),
    el("div", { class: "figure" }, canvas),
    el("div", { class: "card", style: "margin-top:1rem" }, el("div", { class: "card__body" }, status)),
  );

  requestAnimationFrame(run);
  return container;
}
