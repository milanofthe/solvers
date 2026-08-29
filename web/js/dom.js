/* Minimal DOM helpers, so the views read as structure rather than as string
 * concatenation and nothing has to reach for innerHTML. */

export function el(tag, props = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(props ?? {})) {
    if (value === null || value === undefined || value === false) continue;
    if (key === "class") node.className = value;
    else if (key === "text") node.textContent = value;
    else if (key === "dataset") Object.assign(node.dataset, value);
    else if (key.startsWith("on") && typeof value === "function") {
      node.addEventListener(key.slice(2).toLowerCase(), value);
    } else node.setAttribute(key, value === true ? "" : String(value));
  }
  for (const child of children.flat()) {
    if (child === null || child === undefined || child === false) continue;
    node.append(child instanceof Node ? child : document.createTextNode(String(child)));
  }
  return node;
}

export function clear(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
}

export function replace(node, ...children) {
  clear(node);
  for (const child of children.flat()) {
    if (child) node.append(child);
  }
}

/** Fixed significant digits, with a dash for anything not a number. */
export function num(value, digits = 3) {
  if (value === null || value === undefined || !Number.isFinite(value)) return "–";
  if (value === 0) return "0";
  const magnitude = Math.abs(value);
  if (magnitude >= 1e5 || magnitude < 1e-3) return value.toExponential(digits - 1);
  return Number(value.toPrecision(digits)).toString();
}

/** Renders `unbounded` for the limits that legitimately have no finite value. */
export function limit(value) {
  if (value === "unbounded") return "∞";
  return num(value);
}

export function boolBadge(value, trueLabel, falseLabel) {
  if (value === null || value === undefined) return null;
  return el("span", { class: value ? "badge badge--ok" : "badge" }, value ? trueLabel : falseLabel);
}
