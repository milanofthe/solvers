/**
 * Tells a card whether it is on screen, so the pipeline can compute what the
 * reader is actually looking at before anything else.
 */
export function viewport(node, onChange) {
	const observer = new IntersectionObserver(
		(entries) => onChange(entries[0].isIntersecting),
		{ rootMargin: '400px 0px' }
	);
	observer.observe(node);
	return {
		destroy() {
			observer.disconnect();
		}
	};
}
