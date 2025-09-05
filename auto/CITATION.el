;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "CITATION"
 (lambda ()
   (LaTeX-add-bibitems
    "ePowerSim.jl"))
 '(or :bibtex :latex))

