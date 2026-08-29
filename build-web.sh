#!/usr/bin/env sh
# Build the browser bundle and start the interface.
#
# The method files are compiled into the WebAssembly module, so the page needs
# no server side; the production build is static and any file server will do.
set -e
root=$(cd "$(dirname "$0")" && pwd)

wasm-pack build "$root/crates/solvers-wasm"     --target web     --out-dir "$root/web/src/lib/wasm"     --release     --no-typescript

cd "$root/web"
[ -d node_modules ] || npm install
npm run dev
