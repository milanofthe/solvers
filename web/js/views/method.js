/*
 * One method in full: the coefficients it is made of, what those coefficients
 * imply, its stability region, and where it was published.
 */

import { el, limit, num } from "../dom.js";
import { drawStabilityPanel } from "./stability.js";

function coefficientCell(entry) {
  if (!entry) return el("td", { class: "tableau__zero" }, "0");
  const zero = entry.value === 0;
  const text = entry.exact ? entry.text : num(entry.value, 8);
  return el(
    "td",
    { class: zero ? "tableau__zero" : "", title: entry.exact ? "exact" : "floating point" },
    zero ? "0" : text,
  );
}

function butcherTableau(coefficients) {
  const stages = coefficients.stages;
  const rows = [];
  for (let i = 0; i < stages; i += 1) {
    rows.push(
      el(
        "tr",
        {},
        el("td", { class: "tableau__c" }, coefficients.c[i].exact ? coefficients.c[i].text : num(coefficients.c[i].value, 8)),
        coefficients.a[i].map((entry) => coefficientCell(entry)),
      ),
    );
  }
  rows.push(
    el(
      "tr",
      { class: "tableau__rule" },
      el("td", { class: "tableau__c tableau__label" }, "b"),
      coefficients.b.map((entry) => coefficientCell(entry)),
    ),
  );
  if (coefficients.bEmbedded) {
    rows.push(
      el(
        "tr",
        {},
        el("td", { class: "tableau__c tableau__label" }, "b̂"),
        coefficients.bEmbedded.map((entry) => coefficientCell(entry)),
      ),
    );
  }
  return el("table", { class: "tableau" }, el("tbody", {}, rows));
}

function multistepCoefficients(coefficients) {
  const rows = [];
  const header = ["j"];
  for (let j = 0; j <= coefficients.steps; j += 1) header.push(String(j));
  rows.push(el("tr", {}, header.map((label, index) =>
    el("td", { class: index === 0 ? "tableau__c tableau__label" : "tableau__label" }, label),
  )));
  const line = (label, values) =>
    el(
      "tr",
      {},
      el("td", { class: "tableau__c tableau__label" }, label),
      (values ?? []).map((value) => el("td", {}, num(value, 6))),
    );
  rows.push(line("α", coefficients.alpha));
  rows.push(line("β", coefficients.beta));
  return el("table", { class: "tableau" }, el("tbody", {}, rows));
}

function properties(report) {
  const entries = [
    ["order", report.computed_order],
    ["embedded order", report.computed_embedded_order ?? "–"],
    [report.class === "runge_kutta" ? "stages" : "steps", report.size],
    ["stage order", report.stage_order ?? "–"],
    ["cost per step", report.stage_cost],
    ["A-stable", report.a_stable ? "yes" : "no"],
    ["L-stable", report.l_stable ? "yes" : "no"],
    ["stiffly accurate", report.stiffly_accurate == null ? "–" : report.stiffly_accurate ? "yes" : "no"],
    ["R(∞)", report.damping_at_infinity == null ? "–" : limit(report.damping_at_infinity)],
    ["A(α) angle", report.alpha_angle == null ? "–" : `${num(report.alpha_angle, 4)}°`],
    ["real limit", limit(report.real_stability_limit)],
    ["imaginary limit", limit(report.imaginary_stability_limit)],
    ["order check", report.exact_arithmetic ? "exact" : "numeric"],
  ];
  return el(
    "dl",
    { class: "properties" },
    entries.map(([label, value]) => el("div", {}, el("dt", {}, label), el("dd", {}, String(value)))),
  );
}

function stabilityFunction(wasm, id) {
  let data;
  try {
    data = JSON.parse(wasm.stability_function(id));
  } catch (error) {
    return null;
  }
  const polynomial = (coefficients) =>
    coefficients
      .map((value, power) => ({ value, power }))
      .filter((term) => Math.abs(term.value) > 1e-14)
      .map((term) => `${num(term.value, 6)}${term.power ? ` z^${term.power}` : ""}`)
      .join(" + ") || "0";

  return el(
    "div",
    { class: "card" },
    el(
      "div",
      { class: "card__body" },
      el("h3", { class: "card__title" }, "Stability function"),
      el("p", { class: "mono" }, `N(z) = ${polynomial(data.numerator)}`),
      el("p", { class: "mono" }, `D(z) = ${polynomial(data.denominator)}`),
      el(
        "p",
        { class: "card__meta" },
        `R agrees with exp to order ${data.orderOfConsistency}. R(∞) = ${limit(data.atInfinity)}. ` +
          (data.poles.length
            ? `Poles at ${data.poles.map((p) => `${num(p.re, 4)}${p.im >= 0 ? "+" : ""}${num(p.im, 4)}i`).join(", ")}.`
            : "No poles, the method is explicit."),
      ),
    ),
  );
}

export function methodView(context, id) {
  const { wasm } = context;
  let detail;
  try {
    detail = JSON.parse(wasm.method_detail(id));
  } catch (error) {
    return el("p", { class: "status" }, `Unknown method: ${id}`);
  }

  const summary = context.methods.find((m) => m.id === id);
  const canvas = el("canvas", { class: "figure__canvas figure__canvas--tall" });
  requestAnimationFrame(() => {
    if (summary) drawStabilityPanel(canvas, wasm, summary, { title: `stability region ${id}` });
  });

  const coefficients = detail.coefficients;
  const structure =
    coefficients.kind === "runge_kutta"
      ? [
          coefficients.structure?.replace("_", " "),
          coefficients.singlyDiagonal ? "singly diagonal" : null,
          coefficients.explicitFirstStage ? "explicit first stage" : null,
          coefficients.fsal ? "FSAL" : null,
        ].filter(Boolean)
      : [
          `${coefficients.steps} steps`,
          coefficients.startup ? `starts with ${coefficients.startup}` : "self starting",
        ];

  return el(
    "article",
    {},
    el(
      "section",
      { class: "section" },
      el("h2", { class: "section__title" }, detail.name),
      detail.description ? el("p", { class: "section__note" }, detail.description) : null,
      el("div", { class: "badges" }, structure.map((item) => el("span", { class: "badge" }, item))),
      detail.report.discrepancies.length
        ? el(
            "div",
            { class: "badges" },
            detail.report.discrepancies.map((issue) => el("span", { class: "badge badge--error" }, issue)),
          )
        : null,
    ),
    el("section", { class: "section" }, properties(detail.report)),
    el(
      "section",
      { class: "section" },
      el("div", { class: "grid grid--wide" }, [
        el(
          "div",
          { class: "card" },
          el(
            "div",
            { class: "card__body" },
            el("h3", { class: "card__title" }, coefficients.kind === "runge_kutta" ? "Butcher tableau" : "Coefficients on a uniform grid"),
            el(
              "div",
              { class: "scroll-x" },
              coefficients.kind === "runge_kutta"
                ? butcherTableau(coefficients)
                : multistepCoefficients(coefficients),
            ),
            el(
              "p",
              { class: "card__meta" },
              coefficients.kind === "runge_kutta"
                ? "Exact fractions are shown as written in the method file; a decimal means the coefficient is only known as a double."
                : "A multistep file stores which coefficients are free, not their values. These are the values that come out on a uniform grid; on a varying step size they are solved for again at every step.",
            ),
          ),
        ),
        el("div", { class: "figure" }, canvas, el("p", { class: "figure__caption" }, "Stability region")),
      ]),
    ),
    coefficients.kind === "runge_kutta"
      ? el("section", { class: "section" }, stabilityFunction(wasm, id))
      : null,
    detail.references.length
      ? el(
          "section",
          { class: "section" },
          el("h3", { class: "section__title" }, "References"),
          el(
            "ul",
            { class: "references" },
            detail.references.map((reference) =>
              el(
                "li",
                {},
                reference.link
                  ? el("a", { href: reference.link, target: "_blank", rel: "noreferrer" }, reference.title ?? reference.link)
                  : el("span", {}, reference.title ?? "untitled"),
                el(
                  "div",
                  { class: "muted" },
                  [reference.authors, reference.year, reference.source].filter(Boolean).join(" · "),
                ),
              ),
            ),
          ),
        )
      : null,
  );
}
