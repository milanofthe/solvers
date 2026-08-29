import { sveltekit } from '@sveltejs/kit/vite';

export default {
	plugins: [sveltekit()],
	worker: {
		format: 'es'
	},
	server: {
		port: 5180
	}
};
