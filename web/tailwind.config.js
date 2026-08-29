/**
 * The palette and type scale are taken from milanrother.com so the two read as
 * one body of work. The corners are the one departure: everything here is
 * square.
 */
export default {
	content: ['./src/**/*.{html,js,svelte}'],
	theme: {
		// Three sizes, and no more.
		fontSize: {
			xs: ['0.6875rem', { lineHeight: '1.45' }],
			sm: ['0.8125rem', { lineHeight: '1.55' }],
			lg: ['1.125rem', { lineHeight: '1.3' }]
		},
		borderRadius: {
			none: '0',
			DEFAULT: '0'
		},
		extend: {
			colors: {
				charcoal: {
					DEFAULT: '#0f0f0f',
					warm: '#121210',
					light: '#1a1a18'
				},
				cream: {
					DEFAULT: '#f0efe9',
					dark: '#d4d3cd',
					light: '#f7f6f2'
				},
				accent: {
					DEFAULT: '#969591',
					light: '#b3b2ae'
				},
				amber: {
					DEFAULT: '#f5a623',
					dark: '#d4900f',
					light: '#f7b84d'
				},
				teal: {
					DEFAULT: '#00d9c0',
					dark: '#00b3a0'
				},
				pathsim: {
					DEFAULT: '#0070c0',
					light: '#3399d6'
				}
			},
			fontFamily: {
				display: ['Space Grotesk', 'system-ui', 'sans-serif'],
				body: ['Inter', 'system-ui', 'sans-serif'],
				mono: ['JetBrains Mono', 'Fira Code', 'monospace']
			}
		}
	},
	plugins: []
};
