#!/bin/sh
# run_tests.sh -- regression tests for collharness itself.
# Usage: ./run_tests.sh   (run from the harness root, after `make`)
set -e

BIN=./collharness
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

fail() { echo "FAIL: $1" >&2; exit 1; }

echo "== edge-case-only exact values =="
$BIN gen --out "$TMP/edge.collcase" --n 29 --seed 7 > "$TMP/edge.log"
grep -q "found_count=50 candidate_count=24 dedup_count=8 value_match_count=24" \
	"$TMP/edge.log" || fail "edge-case stats do not match hand-verified values"

echo "== round trip: strip expected output, recompute via ref, cmp =="
for n in 500 5000; do
	$BIN gen --out "$TMP/g$n.collcase" --n "$n" --seed "$n" >/dev/null
	# a tiny helper program isn't available standalone here, so use
	# `ref` directly on the freshly generated file (it recomputes and
	# overwrites the expected section unconditionally) and diff the
	# reported numbers against the original `gen` output instead of
	# a binary strip -- equivalent coverage, no extra binary needed.
	$BIN ref --in "$TMP/g$n.collcase" --out "$TMP/r$n.collcase" > "$TMP/r$n.log"
	$BIN cmp "$TMP/g$n.collcase" "$TMP/r$n.collcase" | grep -q "^MATCH" \
		|| fail "round trip mismatch for n=$n"
done

echo "== root_bytes=8, bucket_hash=1, skew (both modes) run without error =="
$BIN gen --out "$TMP/r8.collcase" --n 2000 --root-bytes 8 --key-bits 40 --seed 1 >/dev/null
$BIN gen --out "$TMP/h1.collcase" --n 2000 --bucket-hash 1 --seed 2 >/dev/null
$BIN gen --out "$TMP/sk0.collcase" --n 1000 --skew 200 --bucket-hash 0 --seed 3 >/dev/null
$BIN gen --out "$TMP/sk1.collcase" --n 1000 --skew 200 --bucket-hash 1 --seed 4 >/dev/null
for f in r8 h1 sk0 sk1; do
	$BIN info "$TMP/$f.collcase" >/dev/null || fail "info failed for $f"
done

echo "== suite generation (tiny + med, no --large) =="
mkdir -p "$TMP/suite"
$BIN suite --out-dir "$TMP/suite" > "$TMP/suite.log"
test -f "$TMP/suite/tiny_basic4.collcase" || fail "suite did not produce tiny_basic4"
grep -q "skipping large_default4" "$TMP/suite.log" || fail "large case should be skipped by default"

echo "ALL TESTS PASSED"
