#!/usr/bin/env bash
# Build the paper as PDF and HTML from generation_codes.md. Run from paper/.
set -euo pipefail
cd "$(dirname "$0")"
(cd .. && python3 paper/make_figures.py)
pandoc generation_codes.md -o generation_codes.pdf --pdf-engine=xelatex \
  --default-image-extension=pdf -H header.tex \
  -V geometry:margin=1in -V fontsize=11pt -V colorlinks=true -V linkcolor=teal -V urlcolor=teal \
  -V fontfamily=fontspec -V mainfont="Latin Modern Roman" -V mathfont="Latin Modern Math"
pandoc generation_codes.md -o generation_codes.html -s --mathml --embed-resources \
  --default-image-extension=svg --css paper.css \
  --metadata pagetitle="Beyond One Generation" \
  -V include-before='<p class="download"><a href="generation_codes.pdf">PDF version</a></p>'
echo "built generation_codes.pdf and generation_codes.html"
