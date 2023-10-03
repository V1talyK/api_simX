using Documenter
using ITPM_SimX  # your package name here

 makedocs(
     sitename = "ITPM_SimX",  # your package name here
     format = Documenter.HTML(prettyurls = false),  # optional
     pages = [
         "Introduction" => "index.md"
     ]
 )

 # Documenter can also automatically deploy documentation to gh-pages.
 # See "Hosting Documentation" and deploydocs() in the Documenter manual
 # for more information.
 deploydocs(
     repo = "github.com/V1talyK/ITPM_SimX.jl.git",
 )
