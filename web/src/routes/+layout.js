// The whole interface runs on WebAssembly in the browser, so there is nothing
// to render on a server. Prerendering still earns its keep: with no server
// rendering it writes out one empty shell per route, which is what turns a
// direct hit on `/reference` into a page rather than into the 404 the host
// serves for a path it has no file for. The method pages keep the fallback,
// because their addresses are not known until the library is read.
export const ssr = false;
export const prerender = true;
