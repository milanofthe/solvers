import adapter from '@sveltejs/adapter-static';
import { vitePreprocess } from '@sveltejs/vite-plugin-svelte';

/** @type {import('@sveltejs/kit').Config} */
const config = {
	preprocess: vitePreprocess(),
	kit: {
		adapter: adapter({
			pages: 'build',
			assets: 'build',
			fallback: 'index.html',
			precompress: false,
			strict: false
		}),
		paths: { base: '' },
		// Nothing is server rendered, so nothing links to anything at build time
		// and the crawler finds nothing. The fixed routes are named instead; the
		// method pages are deliberately absent and fall back to the shell.
		prerender: { entries: ['/', '/problems', '/compare', '/reference'] }
	}
};

export default config;
