/*
 * Step size control.
 *
 * The controller is the part of an adaptive solver that is easiest to get
 * wrong and hardest to see. Plotting the step size against time makes the
 * difference between the presets legible: a deadbeat controller chatters on a
 * problem with sharp transitions, a filtered one rides through them, and a
 * predictive one recovers faster after the switch.
 */

import { el, num, replace } from "../dom.js";
import { extent, Figure } from "../plot.js";
import { seriesColor } from "../colormap.js";

const DEFAULT_CONTROLLERS = ["i", "pi4020", "h211b", "gustafsson"];

export function controllerView(context) {
  const { wasm, methods, problems, state } = context;
  const options = JSON.parse(wasm.options_catalog());

  state.controllerMethod ??= "esdirk43";
  state.controllerProblem ??= "van_der_pol_stiff";
  state.controllerTolerance ??= 1e-6;
  state.controllers ??= DEFAULT_CONTROLLERS.filter((name) => options.controllers.includes(name));

  const stepCanvas = el("canvas", { class: "figure__canvas figure__canvas--tall" });
  const solutionCanvas = el("canvas", { class: "figure__canvas" });
  const table = el("tbody");
  const status = el("p", { class: "card__meta" });

  function run() {
    const method = methods.find((m) => m.id === state.controllerMethod);
    if (!method || state.controllers.length === 0) {
      status.textContent = "Select a method and at least one controller.";
      return;
    }
    status.textContent = "Running…";

    const runs = [];
    for (const name of state.controllers) {
      try {
        runs.push({
          name,
          data: JSON.parse(
            wasm.step_history(
              method.id,
              state.controllerProblem,
              state.controllerTolerance,
              state.controllerTolerance * 1e-2,
              name,
            ),
          ),
        });
      } catch (error) {
        runs.push({ name, error: String(error) });
      }
    }

    const usable = runs.filter((entry) => entry.data && entry.data.h.length > 1);
    const stepFigure = new Figure(stepCanvas, {
      title: `step size on ${state.controllerProblem}`,
      xLabel: "t",
      yLabel: "step size h",
      yScale: "log",
    });
    const [tLow, tHigh] = extent(usable.flatMap((entry) => entry.data.t), { pad: 0.02 });
    const [hLow, hHigh] = extent(
      usable.flatMap((entry) => entry.data.h.map(Math.abs)),
      { log: true, pad: 0.1 },
    );
    stepFigure.setDomain(tLow, tHigh, hLow, hHigh);
    stepFigure.clear();
    stepFigure.drawAxes({ origin: false });
    stepFigure.drawSeries(
      usable.map((entry, index) => ({
        x: entry.data.t,
        y: entry.data.h.map(Math.abs),
        label: entry.name,
        color: seriesColor(index),
        marker: false,
        lineWidth: 1.2,
      })),
    );

    // The solution itself, so the step size plot can be read against it.
    try {
      const solution = JSON.parse(
        wasm.trajectory(
          method.id,
          state.controllerProblem,
          state.controllerTolerance,
          state.controllerTolerance * 1e-2,
          600,
        ),
      );
      const components = solution.y[0]?.length ?? 0;
      const solutionFigure = new Figure(solutionCanvas, {
        title: `solution of ${state.controllerProblem}`,
        xLabel: "t",
        yLabel: "y",
      });
      const series = [];
      for (let c = 0; c < Math.min(components, 4); c += 1) {
        series.push({
          x: solution.t,
          y: solution.y.map((row) => row[c]),
          label: `y${c + 1}`,
          color: seriesColor(c),
          marker: false,
        });
      }
      const [yLow, yHigh] = extent(series.flatMap((s) => s.y), { pad: 0.06 });
      solutionFigure.setDomain(tLow, tHigh, yLow, yHigh);
      solutionFigure.clear();
      solutionFigure.drawAxes({ origin: false });
      solutionFigure.drawSeries(series);
    } catch (error) {
      /* the solution panel is a nicety, the step size plot is the point */
    }

    replace(
      table,
      runs.map((entry) =>
        el(
          "tr",
          {},
          el("td", { class: "table__name" }, entry.name),
          el("td", { class: "table__numeric" }, entry.data ? String(entry.data.stats.accepted) : "–"),
          el("td", { class: "table__numeric" }, entry.data ? String(entry.data.stats.rejected) : "–"),
          el(
            "td",
            { class: "table__numeric" },
            entry.data && entry.data.stats.steps
              ? `${num((100 * entry.data.stats.rejected) / entry.data.stats.steps, 3)}%`
              : "–",
          ),
          el("td", { class: "table__numeric" }, entry.data ? String(entry.data.stats.rhs_evals) : "–"),
          el("td", { class: "table__numeric" }, entry.data ? String(entry.data.stats.lu_decompositions) : "–"),
          el(
            "td",
            {},
            entry.error
              ? el("span", { class: "badge badge--error" }, "failed")
              : entry.data.status === "success"
                ? el("span", { class: "badge badge--ok" }, "success")
                : el("span", { class: "badge badge--warn" }, entry.data.status),
          ),
        ),
      ),
    );

    status.textContent =
      "The rejection rate is the direct read on how well a controller suits the problem: a good " +
      "one keeps it low without giving up step size.";
  }

  const container = el(
    "section",
    { class: "section" },
    el("h2", { class: "section__title" }, "Step size control"),
    el(
      "p",
      { class: "section__note" },
      "Every controller here is the same digital filter with different gains, so the comparison " +
        "isolates the control law from everything else in the solver.",
    ),
    el(
      "div",
      { class: "toolbar" },
      el(
        "label",
        { class: "field" },
        el("span", { class: "field__label" }, "method"),
        el(
          "select",
          {
            onchange: (event) => {
              state.controllerMethod = event.target.value;
              run();
            },
          },
          methods
            .filter((m) => m.adaptive)
            .map((m) =>
              el("option", { value: m.id, selected: m.id === state.controllerMethod }, m.name),
            ),
        ),
      ),
      el(
        "label",
        { class: "field" },
        el("span", { class: "field__label" }, "problem"),
        el(
          "select",
          {
            onchange: (event) => {
              state.controllerProblem = event.target.value;
              run();
            },
          },
          problems.map((p) =>
            el(
              "option",
              { value: p.id, selected: p.id === state.controllerProblem },
              `${p.name}${p.stiff ? " (stiff)" : ""}`,
            ),
          ),
        ),
      ),
      el(
        "label",
        { class: "field" },
        el("span", { class: "field__label" }, "tolerance"),
        el(
          "select",
          {
            onchange: (event) => {
              state.controllerTolerance = Number(event.target.value);
              run();
            },
          },
          [1e-4, 1e-6, 1e-8, 1e-10].map((value) =>
            el(
              "option",
              { value: String(value), selected: value === state.controllerTolerance },
              value.toExponential(0),
            ),
          ),
        ),
      ),
      el(
        "label",
        { class: "field" },
        el("span", { class: "field__label" }, "controllers"),
        el(
          "select",
          {
            multiple: true,
            size: "6",
            onchange: (event) => {
              state.controllers = [...event.target.selectedOptions].map((o) => o.value);
              run();
            },
          },
          options.controllers.map((name) =>
            el("option", { value: name, selected: state.controllers.includes(name) }, name),
          ),
        ),
      ),
    ),
    el("div", { class: "figure" }, stepCanvas),
    el("div", { class: "figure", style: "margin-top:1rem" }, solutionCanvas),
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
              el("th", {}, "controller"),
              el("th", {}, "accepted"),
              el("th", {}, "rejected"),
              el("th", {}, "rejection rate"),
              el("th", {}, "rhs evaluations"),
              el("th", {}, "factorizations"),
              el("th", {}, "status"),
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
