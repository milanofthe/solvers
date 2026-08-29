/*
 * The method library.
 *
 * Every column here is derived from the coefficients rather than read from the
 * file, so the table doubles as a check on the library: a method whose file
 * claims something the analysis disagrees with is flagged rather than believed.
 */

import { el, limit, num, replace } from "../dom.js";

const COLUMNS = [
  { key: "name", label: "method", sort: (m) => m.name.toLowerCase() },
  { key: "family", label: "family", sort: (m) => m.family },
  { key: "order", label: "order", numeric: true, sort: (m) => m.order },
  {
    key: "embeddedOrder",
    label: "embedded",
    numeric: true,
    sort: (m) => m.embeddedOrder ?? -1,
    render: (m) => (m.embeddedOrder == null ? "–" : String(m.embeddedOrder)),
  },
  {
    key: "size",
    label: "stages",
    numeric: true,
    sort: (m) => m.size,
    render: (m) => `${m.size}`,
  },
  {
    key: "stageOrder",
    label: "stage order",
    numeric: true,
    sort: (m) => m.stageOrder ?? -1,
    render: (m) => (m.stageOrder == null ? "–" : String(m.stageOrder)),
  },
  {
    key: "stageCost",
    label: "cost/step",
    numeric: true,
    sort: (m) => m.stageCost,
    render: (m) => `${m.stageCost}`,
  },
  {
    key: "stability",
    label: "stability",
    sort: (m) => (m.lStable ? 2 : m.aStable ? 1 : 0),
    render: (m) => {
      if (m.lStable) return "L-stable";
      if (m.aStable) return "A-stable";
      if (m.alphaAngle != null && m.alphaAngle > 0) return `A(${num(m.alphaAngle, 4)}°)`;
      return "conditional";
    },
  },
  {
    key: "dampingAtInfinity",
    label: "R(∞)",
    numeric: true,
    sort: (m) =>
      typeof m.dampingAtInfinity === "number" ? Math.abs(m.dampingAtInfinity) : Number.POSITIVE_INFINITY,
    render: (m) => (m.dampingAtInfinity == null ? "–" : limit(m.dampingAtInfinity)),
  },
  {
    key: "doi",
    label: "source",
    sort: (m) => m.doi ?? "",
    render: (m) =>
      m.doi
        ? el("a", { href: `https://doi.org/${m.doi}`, target: "_blank", rel: "noreferrer" }, "DOI")
        : "–",
  },
];

export function libraryView(context) {
  const { methods, state, navigate } = context;

  const body = el("tbody");
  const summary = el("p", { class: "card__meta" });

  const search = el("input", {
    type: "search",
    placeholder: "name, family or property",
    value: state.query ?? "",
    oninput: (event) => {
      state.query = event.target.value;
      render();
    },
  });

  const filters = [
    { key: "implicitOnly", label: "implicit only" },
    { key: "adaptiveOnly", label: "adaptive only" },
    { key: "stiffOnly", label: "A-stable only" },
  ];

  const checklist = el(
    "ul",
    { class: "checklist" },
    filters.map((filter) =>
      el(
        "li",
        {},
        el(
          "label",
          {},
          el("input", {
            type: "checkbox",
            checked: Boolean(state[filter.key]),
            onchange: (event) => {
              state[filter.key] = event.target.checked;
              render();
            },
          }),
          filter.label,
        ),
      ),
    ),
  );

  const head = el(
    "tr",
    {},
    COLUMNS.map((column) =>
      el(
        "th",
        {
          "aria-sort": state.sortKey === column.key ? (state.sortAscending ? "ascending" : "descending") : "none",
          onclick: () => {
            if (state.sortKey === column.key) state.sortAscending = !state.sortAscending;
            else {
              state.sortKey = column.key;
              state.sortAscending = true;
            }
            render();
          },
        },
        column.label,
      ),
    ),
  );

  function render() {
    const query = (state.query ?? "").trim().toLowerCase();
    let rows = methods.filter((method) => {
      if (state.implicitOnly && !method.implicit) return false;
      if (state.adaptiveOnly && !method.adaptive) return false;
      if (state.stiffOnly && !method.aStable) return false;
      if (!query) return true;
      return [method.id, method.name, method.family, method.class, method.description ?? ""]
        .join(" ")
        .toLowerCase()
        .includes(query);
    });

    const column = COLUMNS.find((c) => c.key === state.sortKey) ?? COLUMNS[1];
    const direction = state.sortAscending ? 1 : -1;
    rows = [...rows].sort((a, b) => {
      const x = column.sort(a);
      const y = column.sort(b);
      if (x === y) return a.id.localeCompare(b.id);
      return x > y ? direction : -direction;
    });

    replace(
      body,
      rows.map((method) =>
        el(
          "tr",
          { onclick: () => navigate(`#/method/${method.id}`), style: "cursor:pointer" },
          COLUMNS.map((c) =>
            el(
              "td",
              { class: c.numeric ? "table__numeric" : c.key === "name" ? "table__name" : "" },
              c.render ? c.render(method) : String(method[c.key] ?? "–"),
            ),
          ),
        ),
      ),
    );

    const flagged = rows.filter((m) => m.discrepancies.length > 0).length;
    summary.textContent =
      `${rows.length} of ${methods.length} methods` +
      (flagged ? ` · ${flagged} disagree with their file` : " · all agree with their analysis");
  }

  render();

  return el(
    "section",
    { class: "section" },
    el("h2", { class: "section__title" }, "Method library"),
    el(
      "p",
      { class: "section__note" },
      "Order, stage order and stability are computed from the coefficients, not taken from the " +
        "file. The cost column counts right hand side evaluations per step, with FSAL reuse " +
        "already deducted.",
    ),
    el(
      "div",
      { class: "toolbar" },
      el("label", { class: "field field--grow" }, el("span", { class: "field__label" }, "search"), search),
      el("div", { class: "field" }, el("span", { class: "field__label" }, "filter"), checklist),
    ),
    el("div", { class: "scroll-x card" }, el("table", { class: "table" }, el("thead", {}, head), body)),
    el("div", { style: "margin-top:0.5rem" }, summary),
  );
}
