(TeX-add-style-hook
 "readme"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-environments-local "minted")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "inputenc"
    "fontenc"
    "graphicx"
    "grffile"
    "longtable"
    "wrapfig"
    "rotating"
    "ulem"
    "amsmath"
    "textcomp"
    "amssymb"
    "capt-of"
    "hyperref"
    "minted")
   (LaTeX-add-labels
    "sec:org5b38cc3"
    "sec:org6d95a01"
    "sec:org46af9f0"
    "sec:org0e31fc4"
    "sec:org028dfb7"
    "sec:org04d2a7d"
    "sec:org5e7249f"
    "sec:orga952248"
    "sec:orgad8991a"
    "sec:org715b463"
    "sec:org6e7e368"
    "sec:org664a940"
    "sec:orga40ca8c"
    "sec:orgf1a23b1"
    "sec:orgcc15a15"
    "sec:orgf899ab8"
    "sec:org1f5f3da"
    "sec:orgaa61298"
    "sec:org59c08b9"
    "sec:org4d5ae3b"
    "sec:org17ab98d"
    "sec:org5ec0355"
    "sec:org15010ea"
    "sec:org18b491a"
    "sec:orgec1d30d"
    "sec:org83e5060"
    "sec:org04b268f"
    "sec:org8da3aed"
    "sec:org6eb6ca4"
    "sec:orgb0a6848"
    "sec:orgebc566e"
    "sec:orgeed60b3"
    "sec:orgf3d08d5"
    "sec:org0669f4f"
    "sec:orgbf523ec"
    "sec:org4ab39c8"
    "sec:org77138c2"
    "sec:org82f6eeb"
    "sec:orgc22eeb7"
    "sec:org9971c98"
    "sec:orgb0feddd"
    "sec:org3302c69"
    "sec:org1f125b7"
    "sec:orgff454f0"
    "sec:org8a43dbf"
    "sec:orgbb99309"
    "sec:org014bb6d"
    "sec:orgf4426f7"
    "sec:org9cdbc28"
    "sec:orgbc969e9"
    "sec:org98ef15c"
    "sec:org36ec673"
    "sec:org7b8765f"
    "sec:org2ace33d"
    "sec:org198c7d8"
    "sec:orgc1f2cdc"
    "sec:org2db8d92"
    "sec:orgbbc7a5a"
    "sec:orgbc35629"
    "sec:orgf599c6c"
    "sec:org311825d"))
 :latex)

