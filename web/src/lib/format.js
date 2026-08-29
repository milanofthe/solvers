/** Number formatting shared by the tables and the card facts. */
export function num(value, digits = 3) {
	if (value === null || value === undefined || !Number.isFinite(value)) return '\u2013';
	if (value === 0) return '0';
	const magnitude = Math.abs(value);
	if (magnitude >= 1e5 || magnitude < 1e-3) return value.toExponential(digits - 1);
	return String(Number(value.toPrecision(digits)));
}

/** Renders the limits that legitimately have no finite value. */
export function limit(value) {
	if (value === 'unbounded') return '\u221e';
	return num(value);
}

export function stabilityLabel(method) {
	if (method.lStable) return 'L-stable';
	if (method.aStable) return 'A-stable';
	if (method.alphaAngle > 0) return `A(${num(method.alphaAngle, 3)})`;
	return 'conditional';
}
