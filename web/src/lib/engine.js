/*
 * The job pipeline.
 *
 * A grid of cards asks for far more work than it can display at once: fifty
 * four stability regions, and in the order or cost modes fifty four solver runs
 * of a hundred steps each. Issuing those as they are asked for would stall the
 * page for seconds and compute most of them for nothing, since the reader only
 * ever looks at a screenful.
 *
 * So requests do not run, they queue. The queue is ordered by priority, cards
 * raise their priority when they scroll into view, anything still queued for a
 * mode nobody is looking at any more is dropped, and identical requests share
 * one result. A small pool of workers takes jobs off the front, which keeps the
 * interactive thread free and uses a few cores without trying to use all of
 * them.
 */

export const PRIORITY = {
	immediate: 0,
	visible: 1,
	prefetch: 2
};

export class CancelledError extends Error {
	constructor(key) {
		super(`cancelled: ${key}`);
		this.name = 'CancelledError';
		this.cancelled = true;
	}
}

/** True for the rejection a dropped job produces, which callers ignore. */
export function isCancelled(error) {
	return Boolean(error?.cancelled);
}

export class Engine {
	constructor({ size } = {}) {
		const cores = globalThis.navigator?.hardwareConcurrency ?? 4;
		// Leave the machine something to run the interface with.
		this.size = size ?? Math.max(1, Math.min(4, cores - 2));
		this.workers = [];
		this.idle = [];
		this.queue = [];
		this.inFlight = new Map();
		this.cache = new Map();
		this.pending = new Map();
		this.nextId = 1;
		this.listeners = new Set();
	}

	get outstanding() {
		return this.queue.length + this.inFlight.size;
	}

	onProgress(listener) {
		this.listeners.add(listener);
		return () => this.listeners.delete(listener);
	}

	announce() {
		for (const listener of this.listeners) listener(this.outstanding);
	}

	spawn() {
		const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
		worker.onmessage = (event) => this.finish(worker, event.data);
		worker.onerror = (event) => {
			const job = this.inFlight.get(worker);
			if (job) {
				this.inFlight.delete(worker);
				job.reject(new Error(event.message ?? 'worker failed'));
			}
			this.idle.push(worker);
			this.pump();
		};
		this.workers.push(worker);
		return worker;
	}

	/**
	 * Ask for a result. Two requests with the same `key` share one computation
	 * and the answer is remembered.
	 */
	request(kind, args, { key, priority = PRIORITY.visible } = {}) {
		const identity = key ?? `${kind}:${JSON.stringify(args ?? {})}`;
		if (this.cache.has(identity)) return Promise.resolve(this.cache.get(identity));

		const existing = this.pending.get(identity);
		if (existing) {
			if (priority < existing.priority) {
				existing.priority = priority;
				this.sort();
			}
			return existing.promise;
		}

		const job = { id: this.nextId++, kind, args, key: identity, priority };
		job.promise = new Promise((resolve, reject) => {
			job.resolve = (value) => {
				this.cache.set(identity, value);
				this.pending.delete(identity);
				resolve(value);
			};
			job.reject = (error) => {
				this.pending.delete(identity);
				reject(error);
			};
		});
		this.pending.set(identity, job);
		this.queue.push(job);
		this.sort();
		this.pump();
		this.announce();
		return job.promise;
	}

	/** Raise the priority of anything still queued under a matching key. */
	raise(match, priority) {
		let changed = false;
		for (const job of this.queue) {
			if (match(job) && priority < job.priority) {
				job.priority = priority;
				changed = true;
			}
		}
		if (changed) this.sort();
	}

	/**
	 * Drop queued jobs nobody is waiting for any more. Jobs already running are
	 * left alone; stopping one would mean tearing down a worker.
	 */
	drop(match) {
		const kept = [];
		for (const job of this.queue) {
			if (match(job)) {
				this.pending.delete(job.key);
				job.reject(new CancelledError(job.key));
			} else {
				kept.push(job);
			}
		}
		this.queue = kept;
		this.announce();
	}

	sort() {
		this.queue.sort((a, b) => a.priority - b.priority || a.id - b.id);
	}

	pump() {
		while (this.queue.length > 0) {
			let worker = this.idle.pop();
			if (!worker && this.workers.length < this.size) worker = this.spawn();
			if (!worker) return;
			const job = this.queue.shift();
			this.inFlight.set(worker, job);
			worker.postMessage({ id: job.id, kind: job.kind, args: job.args });
		}
	}

	finish(worker, message) {
		const job = this.inFlight.get(worker);
		this.inFlight.delete(worker);
		this.idle.push(worker);
		if (job) {
			if (message.ok) job.resolve(message.result);
			else job.reject(new Error(message.error));
		}
		this.pump();
		this.announce();
	}
}
