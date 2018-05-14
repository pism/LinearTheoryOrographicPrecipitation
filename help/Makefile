MATHJAX_URL=https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML

help.html: help.rst
	rst2html5 --math-output="MathJax ${MATHJAX_URL}" $^ $@
