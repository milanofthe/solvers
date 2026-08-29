/*
 * Entry point.
 *
 * Loads the WebAssembly build of the solver core, pulls the method and problem
 * catalogues out of it, and routes between the views. Nothing is fetched from a
 * server: the method files are compiled into the module, and every number on
 * screen is produced in the page.
 */

import init, * as wasm from "../pkg/solvers_wasm.js";
import { clear, el } from "./dom.js";
import { libraryView } from "./views/library.js";
import { stabilityView } from "./views/stability.js";
import { methodView } from "./views/method.js";
import { convergenceView, costView } from "./views/analysis.js";
import { controllerView } from "./views/controller.js";

const ROUTES = [
  { path: "#/library", label: "Library", view: (context) => libraryView(context) },
  { path: "#/stability", label: "Stability", view: (context) => stabilityView(context) },
  { path: "#/convergence", label: "Convergence", view: (context) => convergenceView(context) },
  { path: "#/cost", label: "Cost", view: (context) => costView(context) },
  { path: "#/control", label: "Step control", view: (context) => controllerView(context) },
];

const state = {
  sortKey: "family",
  sortAscending: true,
};

async function main() {
  const root = document.getElementById("view");
  const nav = document.getElementById("nav");

  await init();

  const methods = JSON.parse(wasm.method_catalog());
  const problems = JSON.parse(wasm.problem_catalog());
  const context = {
    wasm,
    methods,
    problems,
    state,
    navigate: (hash) => {
      window.location.hash = hash;
    },
  };

  for (const route of ROUTES) {
    nav.append(
      el(
        "button",
        {
          class: "nav__item",
          type: "button",
          dataset: { path: route.path },
          onclick: () => context.navigate(route.path),
        },
        route.label,
      ),
    );
  }

  function render() {
    const hash = window.location.hash || ROUTES[0].path;
    for (const button of nav.children) {
      const active = hash === button.dataset.path || (hash.startsWith("#/method/") && button.dataset.path === "#/library");
      if (active) button.setAttribute("aria-current", "page");
      else button.removeAttribute("aria-current");
    }

    clear(root);
    if (hash.startsWith("#/method/")) {
      const id = decodeURIComponent(hash.slice("#/method/".length));
      root.append(
        el(
          "p",
          { class: "card__meta", style: "margin-bottom:1rem" },
          el("a", { href: "#/library" }, "← back to the library"),
        ),
        methodView(context, id),
      );
      return;
    }
    const route = ROUTES.find((r) => r.path === hash) ?? ROUTES[0];
    root.append(route.view(context));
  }

  window.addEventListener("hashchange", render);
  // Redraw the canvases when the layout changes.
  let resizeTimer = null;
  window.addEventListener("resize", () => {
    window.clearTimeout(resizeTimer);
    resizeTimer = window.setTimeout(render, 200);
  });

  render();
}

main().catch((error) => {
  const root = document.getElementById("view");
  clear(root);
  root.append(el("p", { class: "status" }, `Failed to start: ${error}`));
  console.error(error);
});
