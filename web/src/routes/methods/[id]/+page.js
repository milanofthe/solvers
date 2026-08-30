// One page per method, and the ids are not known until the library is read out
// of the WebAssembly module. So this route is the one that keeps the fallback:
// the host answers an address it has no file for with the shell, and the router
// takes it from there.
export const prerender = false;
