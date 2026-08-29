/*
 * Stability regions.
 *
 * A gallery of `|R(z)|` on the complex plane, one panel per method: filled log
 * magnitude bands with the unit contour drawn on top, which is the boundary of
 * the stability region. For a multistep method the same picture is the largest
 * root modulus of `rho(w) - z sigma(w)`, which plays the same role.
 */

import { el, replace } from "../dom.js";
import { Figure } from "../plot.js";

const RESOLUTION = { runge_kutta: 220, linear_multistep: 96 };

/** Default window on the complex plane, widened for the larger regions. */
function windowFor(method) {
  if (method.class === "linear_multistep") {
    const reach = method.aStable ? 10 : 6;
    return { re: [-reach, reach * 0.4], im: [-reach * 0.7, reach * 0.7] };
  }
  if (method.implicit) {
    return { re: [-12, 6], im: [-9, 9] };
  }
  const reach = Math.max(4, (method.size ?? 4) * 1.1);
  return { re: [-reach * 1.6, reach * 0.5], im: [-reach, reach] };
}

/** Draw one panel into a canvas. */
export function drawStabilityPanel(canvas, wasm, method, options = {}) {
  const view = options.window ?? windowFor(method);
  const resolution = options.resolution ?? RESOLUTION[method.class] ?? 160;

  const figure = new Figure(canvas, {
    title: options.title ?? `stability region ${method.id}`,
    xLabel: "Re(z)",
    yLabel: "Im(z)",
    colorbar: { label: "|R| (log)" },
  });
  figure.setDomain(view.re[0], view.re[1], view.im[0], view.im[1]);

  let data;
  try {
    data = wasm.stability_grid(
      method.id,
      view.re[0],
      view.re[1],
      view.im[0],
      view.im[1],
      resolution,
      resolution,
    );
  } catch (error) {
    figure.clear();
    figure.drawAxes();
    return;
  }

  // Clip the range so a single pole does not flatten the whole colour scale.
  let low = 0;
  let high = 0;
  for (const value of data) {
    if (!Number.isFinite(value)) continue;
    if (value < low) low = value;
    if (value > high) high = value;
  }
  const reach = Math.min(Math.max(Math.max(-low, high), 1), 5);

  figure.clear();
  figure.drawField(data, resolution, resolution, {
    low: -reach,
    high: reach,
    levels: options.levels ?? 14,
  });
  // The unit contour: log |R| = 0.
  figure.drawContour(data, resolution, resolution, 0, { color: "#000000", lineWidth: 1.3 });
  figure.drawAxes();
  figure.drawColorbar(-reach, reach, { levels: options.levels ?? 14 });
  return figure;
}

export function stabilityView(context) {
  const { wasm, methods, state, navigate } = context;
  const families = [...new Set(methods.map((m) => m.family))].sort();

  const container = el("div");
  const gallery = el("div", { class: "grid grid--wide" });

  const familySelect = el(
    "select",
    {
      onchange: (event) => {
        state.stabilityFamily = event.target.value;
        render();
      },
    },
    el("option", { value: "" }, "all families"),
    families.map((family) =>
      el("option", { value: family, selected: state.stabilityFamily === family }, family),
    ),
  );

  const levelsInput = el("input", {
    type: "range",
    min: "6",
    max: "24",
    step: "2",
    value: String(state.stabilityLevels ?? 14),
    oninput: (event) => {
      state.stabilityLevels = Number(event.target.value);
      render();
    },
  });

  const toolbar = el(
    "div",
    { class: "toolbar" },
    el("label", { class: "field" }, el("span", { class: "field__label" }, "family"), familySelect),
    el("label", { class: "field" }, el("span", { class: "field__label" }, "contour levels"), levelsInput),
  );

  function render() {
    const selected = state.stabilityFamily
      ? methods.filter((m) => m.family === state.stabilityFamily)
      : methods;
    replace(
      gallery,
      selected.map((method) => {
        const canvas = el("canvas", { class: "figure__canvas figure__canvas--tall" });
        const card = el(
          "div",
          { class: "figure card--link", onclick: () => navigate(`#/method/${method.id}`) },
          canvas,
          el(
            "p",
            { class: "figure__caption" },
            `${method.name} · order ${method.order} · ${method.size} ${
              method.class === "runge_kutta" ? "stages" : "steps"
            }${method.aStable ? " · A-stable" : ""}${method.lStable ? " · L-stable" : ""}`,
          ),
        );
        // The canvas needs its layout size before the figure can size itself.
        requestAnimationFrame(() =>
          drawStabilityPanel(canvas, wasm, method, { levels: state.stabilityLevels ?? 14 }),
        );
        return card;
      }),
    );
  }

  render();
  container.append(
    el(
      "section",
      { class: "section" },
      el("h2", { class: "section__title" }, "Stability of the integration methods"),
      el(
        "p",
        { class: "section__note" },
        "Filled bands show log |R(z)| for a Runge-Kutta method, and the largest root modulus of " +
          "rho(w) - z sigma(w) for a multistep method. The black line is the unit contour, so the " +
          "region it encloses is where the method does not amplify. Every panel is computed in the " +
          "page from the tableau.",
      ),
      toolbar,
      gallery,
    ),
  );
  return container;
}
