#!/bin/bash
# run the test_dlp_composites() harness (yafu "special()"), summary lines only
cd /mnt/c/Users/bbuhr/projects/factoring/yafu-3.1.8
BIN=${1:-./yafu}
$BIN "special()" 2>&1 | grep -E "testing|got|average time"
