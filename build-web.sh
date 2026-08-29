#!/usr/bin/env sh
# Build the browser bundle and serve the interface.
#
# The method files are compiled into the wasm module, so the page needs no
# server side; any static file server will do.
set -e
root=$(cd "$(dirname "$0")" && pwd)

wasm-pack build "$root/crates/solvers-wasm" \
    --target web \
    --out-dir "$root/web/pkg" \
    --release \
    --no-typescript

echo
echo "serving http://127.0.0.1:8099"
python -m http.server 8099 --directory "$root/web"
